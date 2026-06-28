/*
  Copyright (c) 2026 Ronald de Man

  This file is distributed under the terms of the GNU GPL, version 2.
*/

#include <inttypes.h>
#include <stdint.h>
#include <stdlib.h>

#include "defs.h"
#include "generatep.h"
#include "index.h"
#include "kslicep.h"
#include "movegen.h"
#include "probe.h"
#include "tb8gen.h"
#include "threads.h"
#include "types.h"
#include "util.h"

const char *wdl_name[7] = {
  "loss", "bloss", "draw", "cwin", "win", "nobloss", "noloss"
};

const char side[2][6] = { "white", "black" };
const char clr_L[4][3] = { "31", "32", "33", "34" };
const char clr_W[4][3] = { "94", "93", "92", "91" };

uint64_t capt_cnt[2][5], sub_cnt[2][5];
uint64_t psub_cnt[5], pcapt_cnt[5], pawn_cnt[5];
int max_iteration;
static uint8_t *merged_table, *merged_table2;

INLINE uint64_t sum16(uint64_t *num)
{
  uint64_t cnt = 0;
  for (int r = 0; r < 16; r++)
    cnt += num[r];
  return cnt;
}

INLINE void mark_king_uncaptures(int stm, int k, Bitboard occ,
    struct IdxState *is)
{
  uint8_t ksq[2] = { is->sq[0], is->sq[1] };

  // The stm king uncaptures, so we need to place a piece where the king was.
  is->bb[k + 1] |= bit(ksq[stm]);
  Bitboard occ0 = is->bb[0] ^ bit(ksq[stm]);

  Bitboard b = king_attacks(ksq[stm]) & ~king_attacks(ksq[stm ^ 1]) & ~occ;
  while (b) {
    ksq[stm] = pop_lsb(&b);
    is->bb[0] = occ0 | bit(ksq[stm]);
    uint64_t idx = rank_bb(is->bb, &ri);
    uint8_t *p = kslice_get_address(ksq);
    kslice_bit_set_atomic(p, idx);
  }

  // Restore the set occupancy bitboards.
  is->bb[0] = is->occ[0];
  is->bb[k + 1] ^= bit(is->sq[stm]);
}

// Uncapture a piece in set k by a piece in set j.
INLINE void mark_uncaptures(int stm, const int pc, int k,
    uint8_t *restrict const p, Bitboard occ, struct IdxState *is)
{
  int j = g_piece_set[stm][pc];
  if (j < 0) return;

  int l = min(j, k);
  uint64_t idx0 = 0;
  for (int i = 0; i < l; i++)
    idx0 = idx0 * ri.factor[i] + is->sub[i];

  Bitboard bb = is->bb[j + 1];
  while (bb) {
    Bitboard to_bb = bb & -bb;
    is->bb[j + 1] ^= to_bb;
    is->bb[k + 1] ^= to_bb;
    bb ^= to_bb;
    Bitboard b = piece_moves(pc, lsb(to_bb), occ);
    while (b) {
      Bitboard from_bb = b & -b;
      is->bb[j + 1] ^= from_bb;
      uint64_t idx = rank_bb_from(is->bb, idx0, l, is->occ[l], &ri);
      kslice_bit_set_atomic(p, idx);
      is->bb[j + 1] ^= from_bb;
      b ^= from_bb;
    }
    is->bb[j + 1] ^= to_bb;
    is->bb[k + 1] ^= to_bb;
  }
}

INLINE void mark_king_unmoves(int stm, Bitboard occ, struct IdxState *is)
{
  uint8_t ksq[2] = { is->sq[0], is->sq[1] };
  Bitboard occ0 = is->bb[0] ^ bit(ksq[stm]);

  Bitboard b = king_attacks(ksq[stm]) & ~king_attacks(ksq[stm ^ 1]) & ~occ;
  while (b) {
    ksq[stm] = pop_lsb(&b);
    is->bb[0] = occ0 | bit(ksq[stm]);
    uint64_t idx = rank_bb(is->bb, &ri);
    uint8_t *p = kslice_get_address(ksq);
    kslice_bit_set_atomic(p, idx);
  }

  is->bb[0] = is->occ[0];
}

INLINE void mark_unmoves(int stm, const int pc, uint8_t *restrict const p,
    Bitboard occ, struct IdxState *is)
{
  int k = g_piece_set[stm][pc];
  if (k < 0) return;

  uint64_t idx0 = 0;
  for (int i = 0; i < k; i++)
    idx0 = idx0 * ri.factor[i] + is->sub[i];

  Bitboard bb = is->bb[k + 1];
  while (bb) {
    Bitboard from_bb = bb & -bb;
    is->bb[k + 1] ^= from_bb;
    bb ^= from_bb;
    Bitboard b = piece_moves(pc, lsb(from_bb), occ);
    while (b) {
      Bitboard to_bb = b & -b;
      is->bb[k + 1] ^= to_bb;
      uint64_t idx = rank_bb_from(is->bb, idx0, k, is->occ[k], &ri);
      kslice_bit_set_atomic(p, idx);
      is->bb[k + 1] ^= to_bb;
      b ^= to_bb;
    }
    is->bb[k + 1] ^= from_bb;
  }
}

// Uncapture stm^1 king or pawn by an stm piece being added to set k.
INLINE void mark_uncapture_king_or_pawn(int stm, const int pc, int capsq,
    uint8_t *restrict const p, Bitboard occ, struct IdxState *is)
{
  int k = g_piece_set[stm][pc];
  if (k < 0) return;

  uint64_t idx0 = 0;
  for (int i = 0; i < k; i++)
    idx0 = idx0 * ri.factor[i] + is->sub[i];

  Bitboard b = piece_moves(pc, capsq, occ);
  while (b) {
    Bitboard from_bb = b & -b;
    is->bb[k + 1] ^= from_bb;
    uint64_t idx = rank_bb_from(is->bb, idx0, k, is->occ[k], &ri);
    kslice_bit_set_atomic(p, idx);
    is->bb[k + 1] ^= from_bb;
    b ^= from_bb;
  }
}

INLINE bool check_king_moves(int stm, Bitboard occ, struct IdxState *is)
{
  uint8_t ksq[2] = { is->sq[0], is->sq[1] };
  Bitboard occ0 = is->bb[0] ^ bit(ksq[stm]);

  Bitboard b = king_attacks(ksq[stm]) & ~king_attacks(ksq[stm ^ 1]) & ~occ;
  while (b) {
    ksq[stm] = pop_lsb(&b);
    is->bb[0] = occ0 | bit(ksq[stm]);
    uint64_t idx = rank_bb(is->bb, &ri);
    uint8_t *p = kslice_get_address(ksq);
    if (!kslice_bit_test(p, idx)) {
      is->bb[0] = is->occ[0];
      return false;
    }
  }

  is->bb[0] = is->occ[0];
  return true;
}

INLINE bool check_moves(int stm, const int pc, uint8_t *restrict const p,
    Bitboard occ, struct IdxState *is)
{
  int k = g_piece_set[stm][pc];
  if (k < 0) return true;

  uint64_t idx0 = 0;
  for (int i = 0; i < k; i++)
    idx0 = idx0 * ri.factor[i] + is->sub[i];

  Bitboard bb = is->bb[k + 1];
  while (bb) {
    Bitboard from_bb = bb & -bb;
    is->bb[k + 1] ^= from_bb;
    bb ^= from_bb;

    Bitboard b = piece_moves(pc, lsb(from_bb), occ);
    while (b) {
      Bitboard to_bb = b & -b;
      is->bb[k + 1] ^= to_bb;
      uint64_t idx = rank_bb_from(is->bb, idx0, k, is->occ[k], &ri);
      is->bb[k + 1] ^= to_bb;
      if (!kslice_bit_test(p, idx)) {
        is->bb[k + 1] ^= from_bb;
        return false;
      }
      b ^= to_bb;
    }
    is->bb[k + 1] ^= from_bb;
  }

  return true;
}

static int work_set, work_r, work_lower, work_upper;

static struct Work work_g_dynamic, work_g_static, work_capt_dynamic[MAX_SETS];

static constexpr uint64_t GENERATE_MIN_DYNAMIC_CHUNK = 1ULL << 18;
static constexpr int GENERATE_DYNAMIC_FACTOR = 4;

void init_generation_work(void)
{
  work_init(&work_g_dynamic, kslice_size, 0x1ff, WORK_DYNAMIC,
      GENERATE_DYNAMIC_FACTOR, GENERATE_MIN_DYNAMIC_CHUNK);
  work_init(&work_g_static, kslice_size, 0x1ff, WORK_STATIC, 1,
      GENERATE_MIN_DYNAMIC_CHUNK);
  for (int k = 0; k < ri.numsets; k++)
    work_init(&work_capt_dynamic[k], capt_ri[k].sizes[0], 0x1ff, WORK_DYNAMIC,
        GENERATE_DYNAMIC_FACTOR, GENERATE_MIN_DYNAMIC_CHUNK);
}

static void calc_sub_worker(struct ThreadData *thread)
{
  struct IdxState is;
  int k = work_set;
  int r = work_r;
  int stm = g_slice.stm;

  Position pos = g_pos;
  pos.sq[0] = g_slice.sq[0];
  pos.sq[1] = g_slice.sq[1];
  pos.sq[2] = g_slice.sq[2];
  pos.stm = stm;
  int n = --pos.num;
  int m = ri.last[k];
  pos.pt[m] = pos.pt[n];
  uint8_t *restrict p[5];
  for (int i = 0; i < 5; i++)
    p[i] = k16slice_sub_buf[i] + sub_offset[k] + r * kslice_sub_alloc_size[k];

  Bitboard occ = idx_state_init(&is, thread->begin, g_slice.sq, &capt_ri[k],
      false);
  for (uint64_t idx = thread->begin, end = thread->end; idx < end;
      idx++, occ = idx_state_inc(&is, &capt_ri[k]))
  {
    if (!idx_state_legal(&is, stm, occ))
      continue;
    pos.occ = occ;
    idx_state_to_sq(&is, pos.sq, &capt_ri[k]);
    pos.sq[m] = pos.sq[n];
    int v = probe_wdl(&pos, -2, 2);
    kslice_bit_set(p[v + 2], idx);
  }
}

// Calculate aggregate bitmaps for subtables, one per loss/bloss/draw/cwin/win.
static void calc_sub_kslices(int stm)
{
  char done[16];
  sprintf(done, "sub/done_%c", "wb"[stm]);
  FILE *F = fopen(done, "rb");
  if (F) {
    for (int i = 0; i < 5; i++)
      file_read(&sub_cnt[stm][i], 8, F);
    fclose(F);
    return;
  }

  char phase[64];
  snprintf(phase, sizeof phase, "probing subtables for %s", side[stm]);

  char name[7][16];
  for (int i = 0; i < 7; i++) {
    strcat(strcpy(name[i], "sub/"), wdl_name[i]);
    create_dir(-1, stm, name[i]);
  }

  g_slice.stm = stm;

  for (int s = 0; s < 240; s++) {
    show_progress(phase, s, 240, false);

    uint64_t c[5];

    if (k16slice_sub_test_count(s, stm, name[6], -1, &c[0])) {
      for (int i = 0; i < 5; i++)
        k16slice_sub_test_count(s, stm, name[i], -1, &c[i]);
    } else {
      for (int i = 0; i < 5; i++)
        k16slice_sub_clear_addr(k16slice_sub_buf[i], stm);

      for (int r = 0; r < 16; r++) {
        g_slice.sq[0] = KK16Square[s][r][0];
        g_slice.sq[1] = KK16Square[s][r][1];
        work_r = r;

        if (is_broken(&g_slice))
          continue;

        for (int k = 0; k < ri.numsets; k++) {
          if ((g_set_pt[k] >> 3) != stm)
            continue;
          work_set = k;
          run_threaded(calc_sub_worker, &work_capt_dynamic[k]);
        }
      }

      for (int i = 0; i < 5; i++) {
        c[i] = k16slice_sub_count_addr(k16slice_sub_buf[i], stm);
        k16slice_sub_write_addr(k16slice_sub_buf[i], s, stm, name[i], c[i]);
      }

      k16slice_sub_or_addr(k16slice_sub_buf[0], k16slice_sub_buf[1], stm);
      k16slice_sub_or_addr(k16slice_sub_buf[0], k16slice_sub_buf[2], stm);
      k16slice_sub_write_addr(k16slice_sub_buf[0], s, stm, name[5],
          c[0] + c[1] + c[2]);

      // TODO: skip this is c[3] == 0.
      k16slice_sub_or_addr(k16slice_sub_buf[0], k16slice_sub_buf[3], stm);
      k16slice_sub_write_addr(k16slice_sub_buf[0], s, stm, name[6],
          c[0] + c[1] + c[2] + c[3]);
    }

    for (int i = 0; i < 5; i++)
      sub_cnt[stm][i] += c[i];
  }

  F = file_open_write(done);
  for (int i = 0; i < 5; i++)
    file_write(&sub_cnt[stm][i], 8, F);
  fclose(F);
  file_rename(done);

  show_progress(phase, 240, 240, true);
}

// Legality check for psub positions: the black pawn has been captured, and
// the old pawn square is occupied by the capturing white piece.
INLINE bool idx_state_legal_psub(const struct IdxState *is, Bitboard occ)
{
  int ksq = is->sq[WHITE];
  for (int i = 0; g_sets[BLACK][i] >= 0; i++) {
    int k = g_sets[BLACK][i];
    Bitboard b = non_king_piece_attacks(g_set_pt[k], ksq, occ);
    if (b & is->bb[k + 1])
      return false;
  }
  return true;
}

static void calc_psub_worker(struct ThreadData *thread)
{
  struct IdxState is;
  int k = work_set;
  int r = work_r;

  Position pos = g_pos;
  pos.sq[0] = g_slice.sq[0];
  pos.sq[1] = g_slice.sq[1];
  pos.sq[2] = g_slice.sq[2];
  pos.stm = BLACK;
  int n = --pos.num;
  int m = ri.last[k];
  pos.pt[2] = pos.pt[m];
  pos.pt[m] = pos.pt[n];

  uint8_t *restrict p[5];
  for (int i = 0; i < 5; i++)
    p[i] = k16slice_sub_buf[i] + psub_offset[k] + r * kslice_sub_alloc_size[k];

  Bitboard occ = idx_state_init(&is, thread->begin, g_slice.sq, &capt_ri[k],
      false);
  for (uint64_t idx = thread->begin, end = thread->end; idx < end;
      idx++, occ = idx_state_inc(&is, &capt_ri[k]))
  {
    if (!idx_state_legal_psub(&is, occ))
      continue;
    pos.occ = occ;
    idx_state_to_sq(&is, pos.sq, &capt_ri[k]);
    pos.sq[m] = pos.sq[n];
    int v = probe_wdl(&pos, -2, 2);
    kslice_bit_set(p[v + 2], idx);
  }
}

// Calculate aggregate bitmaps for subtables reached through a capture of
// the black pawn by a white piece (non-king).
static void calc_psub_kslices(void)
{
  char done[16];
  sprintf(done, "psub/done");
  FILE *F = fopen(done, "rb");
  if (F) {
    for (int i = 0; i < 5; i++)
      file_read(&psub_cnt[i], 8, F);
    fclose(F);
    return;
  }

  char name[7][16];
  for (int i = 0; i < 7; i++) {
    strcat(strcpy(name[i], "psub/"), wdl_name[i]);
    create_dir(-1, WHITE, name[i]);
  }

  g_slice.stm = BLACK;

  for (int s = 0; s < 240; s++) {
    show_progress("probing pawn subtables for black", s, 240, false);

    uint64_t c[5];

    if (k16slice_sub_test_count(s, WHITE, name[6], -1, &c[0])) {
      for (int i = 0; i < 5; i++)
        k16slice_sub_test_count(s, WHITE, name[i], -1, &c[i]);
    } else {
      for (int i = 0; i < 5; i++)
        k16slice_sub_clear_addr(k16slice_sub_buf[i], WHITE);

      for (int r = 0; r < 16; r++) {
        g_slice.sq[0] = KK16Square[s][r][0];
        g_slice.sq[1] = KK16Square[s][r][1];
        work_r = r;

        if (is_broken(&g_slice))
          continue;

        for (int i = 0; g_sets[WHITE][i] >= 0; i++) {
          int k = g_sets[WHITE][i];
          work_set = k;
          run_threaded(calc_psub_worker, &work_capt_dynamic[k]);
        }
      }

      for (int i = 0; i < 5; i++) {
        c[i] = k16slice_sub_count_addr(k16slice_sub_buf[i], WHITE);
        k16slice_sub_write_addr(k16slice_sub_buf[i], s, WHITE, name[i], c[i]);
      }

      k16slice_sub_or_addr(k16slice_sub_buf[0], k16slice_sub_buf[1], WHITE);
      k16slice_sub_or_addr(k16slice_sub_buf[0], k16slice_sub_buf[2], WHITE);
      k16slice_sub_write_addr(k16slice_sub_buf[0], s, WHITE, name[5],
          c[0] + c[1] + c[2]);

      k16slice_sub_or_addr(k16slice_sub_buf[0], k16slice_sub_buf[3], WHITE);
      k16slice_sub_write_addr(k16slice_sub_buf[0], s, WHITE, name[6],
          c[0] + c[1] + c[2] + c[3]);
    }

    for (int i = 0; i < 5; i++)
      psub_cnt[i] += c[i];
  }

  F = file_open_write(done);
  for (int i = 0; i < 5; i++)
    file_write(&psub_cnt[i], 8, F);
  fclose(F);
  file_rename(done);

  show_progress("probing pawn subtables for black", 240, 240, true);
//  for (int i = 0; i < 5; i++)
//    printf("psub_cnt[%d] = %lu\n", i, psub_cnt[i]);
}

static void predecessors_sub_worker(struct ThreadData *thread)
{
  struct IdxState is;
  int stm = g_slice.stm;
  int k = work_set;

  uint64_t *restrict p =
    (uint64_t *)kslice_sub_get_address(g_slice.sq, k) + (thread->begin >> 6);
  uint8_t *restrict const q = kslice_get_address(g_slice.sq);

  uint64_t last = thread->begin;
  idx_state_init(&is, last, g_slice.sq, &capt_ri[k], false);

  for (uint64_t idx = last, end = thread->end; idx < end; idx += 64) {
    uint64_t w = *p++;
    while (w) {
      uint64_t cur = idx + pop_lsb(&w);
      Bitboard occ = idx_state_add(&is, cur - last, &capt_ri[k]);
      last = cur;
      // Uncapture by king.
      mark_king_uncaptures(stm, k, occ, &is);
      // Uncapture by non-king pieces.
      mark_uncaptures(stm, KNIGHT, k, q, occ, &is);
      mark_uncaptures(stm, BISHOP, k, q, occ, &is);
      mark_uncaptures(stm, ROOK  , k, q, occ, &is);
      mark_uncaptures(stm, QUEEN , k, q, occ, &is);
    }
  }
}

// Uncapture pieces from the loaded subtable of stm^1 positions to stm
// positions.
static void predecessors_sub(int stm, int s)
{
  g_slice.stm = stm;
  for (int r = 0; r < 16; r++) {
    g_slice.sq[0] = KK16Square[s][r][0];
    g_slice.sq[1] = KK16Square[s][r][1];

    if (is_broken(&g_slice))
      continue;

    // Loop through the sets from which a piece is removed.
    for (int k = 0; k < ri.numsets; k++) {
      if ((g_set_pt[k] >> 3) == stm)
        continue;
      work_set = k;
      run_threaded(predecessors_sub_worker, &work_capt_dynamic[k]);
    }
  }
}

static void predecessors_psub_worker_tmpl(struct ThreadData *thread,
    const int pc)
{
  struct IdxState is;
  int k = g_piece_set[WHITE][pc];
  int psq = g_slice.sq[2];

  uint64_t *restrict p =
    (uint64_t *)kslice_psub_get_address(g_slice.sq, k) + (thread->begin >> 6);
  uint8_t *restrict const q = kslice_get_address(g_slice.sq);

  uint64_t last = thread->begin;
  idx_state_init(&is, last, g_slice.sq, &capt_ri[k], false);

  for (uint64_t idx = last, end = thread->end; idx < end; idx += 64) {
    uint64_t w = *p++;
    while (w) {
      uint64_t cur = idx + pop_lsb(&w);
      Bitboard occ = idx_state_add(&is, cur - last, &capt_ri[k]);
      last = cur;
      mark_uncapture_king_or_pawn(WHITE, pc, psq, q, occ, &is);
    }
  }
}

static void predecessors_psub_knight_worker(struct ThreadData *thread)
{
  predecessors_psub_worker_tmpl(thread, KNIGHT);
}

static void predecessors_psub_bishop_worker(struct ThreadData *thread)
{
  predecessors_psub_worker_tmpl(thread, BISHOP);
}

static void predecessors_psub_rook_worker(struct ThreadData *thread)
{
  predecessors_psub_worker_tmpl(thread, ROOK);
}

static void predecessors_psub_queen_worker(struct ThreadData *thread)
{
  predecessors_psub_worker_tmpl(thread, QUEEN);
}

// Uncapture the black pawn from the loaded subtable of btm positions
// to wtm positions.
static void predecessors_psub(int s)
{
  for (int r = 0; r < 16; r++) {
    g_slice.sq[0] = KK16Square[s][r][0];
    g_slice.sq[1] = KK16Square[s][r][1];

    if (is_broken(&g_slice))
      continue;

    // Loop through the sets from which a piece uncaptures the pawn.
    int k;
    if ((k = g_piece_set[WHITE][KNIGHT]) >= 0)
      run_threaded(predecessors_psub_knight_worker, &work_capt_dynamic[k]);
    if ((k = g_piece_set[WHITE][BISHOP]) >= 0)
      run_threaded(predecessors_psub_bishop_worker, &work_capt_dynamic[k]);
    if ((k = g_piece_set[WHITE][ROOK]) >= 0)
      run_threaded(predecessors_psub_rook_worker, &work_capt_dynamic[k]);
    if ((k = g_piece_set[WHITE][QUEEN]) >= 0)
      run_threaded(predecessors_psub_queen_worker, &work_capt_dynamic[k]);
  }
}

static void calc_king_captures_pawn_worker(struct ThreadData *thread)
{
  struct IdxState is;
  int lower = work_lower;
  int upper = work_upper;

  Position pos = g_pos;
  pos.sq[0] = g_slice.sq[0];
  pos.sq[1] = g_slice.sq[1];
  pos.sq[2] = g_slice.sq[2];
  pos.stm = g_slice.stm;
  int stm = pos.stm;

  uint8_t *restrict const p = kslice_get_address(g_slice.sq);

  Bitboard occ = idx_state_init(&is, thread->begin, g_slice.sq, &ri, false);
  for (uint64_t idx = thread->begin, end = thread->end; idx < end;
      idx++, occ = idx_state_inc(&is, &ri))
  {
    if (!idx_state_legal(&is, stm, occ))
      continue;
    idx_state_to_sq(&is, pos.sq, &ri);
    pos.occ = occ;
    if (do_capture(&pos, g_slice.sq[0], g_slice.sq[2], 0, 2)) {
      int v = -probe_wdl(&pos, -2, 2);
      if (v >= lower && v <= upper)
        kslice_bit_set(p, idx);
      undo_capture(&pos, g_slice.sq[0], g_slice.sq[2], 0, 2);
    }
  }
}

static void calc_king_captures_pawn(int s, int lower, int upper)
{
  work_lower = lower;
  work_upper = upper;
  g_slice.stm = WHITE;
  for (int r = 0; r < 16; r++) {
    g_slice.sq[0] = KK16Square[s][r][0];
    g_slice.sq[1] = KK16Square[s][r][1];

    if (is_broken(&g_slice))
      continue;

    if (!(king_attacks(g_slice.sq[0]) & bit(g_slice.sq[2])))
      continue;

    run_threaded(calc_king_captures_pawn_worker, &work_g_dynamic);
  }
}

static int wins_full[2][240];
static int wins_checked[2][240];
static uint64_t replay_cost[2][240];
static uint64_t write_cost[2][240];

static void add_wins(int s, int stm, int n, uint64_t cost)
{
  if ((replay_cost[stm][s] + cost) * 1.0 <= write_cost[stm][s])
    return;

  k16slice_or(-1, s);
  k16slice_write(-1, s, stm, "wins", n, nullptr);
  replay_cost[stm][s] = 0;

  k16slice_delete(s, stm, "wins", wins_full[stm][s]);
  for (int i = max(wins_full[stm][s] + 1, wins_checked[stm][s]); i < n; i++)
    k16slice_delete(s, stm, "wins", i);
  wins_full[stm][s] = wins_checked[stm][s] = n;
}

static void read_wins(int s, int slice, int stm, int n)
{
  if (n < 0)
    for (n = MAX_STATS / 2 - 4; n > 0; n--)
      if (g_stats[stm][stats_n(n)])
        break;

  if (wins_checked[stm][slice] < n) {
    int i;
    for (i = n; i > wins_checked[stm][slice]; i--)
      if (k16slice_test(slice, stm, "wins", i))
        break;
    if (i > wins_checked[stm][slice]) {
      for (int j = i - 1; j >= wins_full[stm][slice]; j--)
        k16slice_delete(slice, stm, "wins", j);
      wins_full[stm][slice] = i;
    }
    wins_checked[stm][slice] = n;
  }

  k16slice_read(s, slice, stm, "wins", wins_full[stm][slice]);
  write_cost[stm][slice] = k16slice_read_cost;
  k16slice_read_cost = 0;
  for (int i = wins_full[stm][slice] + 1; i <= n; i++)
    if (g_stats[stm][stats_n(i)])
      k16slice_read_or(s, slice, stm, "W", i);
  replay_cost[stm][slice] += k16slice_read_cost;
}

static void calc_capt(int stm, int wdl)
{
  char capt_name[32], pcapt_name[32], sub_name[32], psub_name[32];
  strcat(strcpy(capt_name , "capt/" ), wdl_name[2 + wdl]);
  strcat(strcpy(pcapt_name, "pcapt/"), wdl_name[2 + wdl]);
  strcat(strcpy(sub_name  , "sub/"  ), wdl_name[2 - wdl]);
  strcat(strcpy(psub_name , "psub/" ), wdl_name[2 - wdl]);

  char str[64];
  sprintf(str, "%s/%c/done", capt_name, "wb"[stm]);
  FILE *F = fopen(str, "rb");
  if (F) {
    file_read(&capt_cnt[stm][2 + wdl], 8, F);
    fclose(F);
    return;
  }

  char phase[64];
  snprintf(phase, sizeof phase, "calculating %s captures for %s",
      wdl_name[2 + wdl], side[stm]);

  bool partial = dir_exists(-1, stm, capt_name), done = partial;

  create_dir(-1, stm, capt_name);

  struct K16SliceIterator iter;
  uint64_t num[16], cnt = 0;
  int num_done = 0;

  k16slice_iter_init(&iter, stm);
  int s, s1;
  while (k16slice_iter_next(&iter, &s)) {
    show_progress(phase, num_done++, 240, false);

    while (k16slice_iter_in(&iter, &s1))
      if (partial && k16slice_test_count(s1, stm, capt_name, -1, num)) {
        cnt += sum16(num);
      } else {
        done = false;
        if (stm == WHITE)
          k16slice_clear(s1);
        else
          k16slice_read(s1, s1, BLACK, pcapt_name, -1);
      }

    if (!done) {
      k16slice_sub_read(s, s, stm ^ 1, sub_name);
      predecessors_sub(stm, s);

      if (stm == WHITE) {
        k16slice_sub_read(s, s, stm, psub_name);
        predecessors_psub(s);
        calc_king_captures_pawn(s, wdl, wdl);
      }
    }

    while (k16slice_iter_out(&iter, &s))
      if (!partial || !k16slice_test_count(s, stm, capt_name, -1, num)) {
        cnt += k16slice_count(s, num);
        k16slice_write(s, s, stm, capt_name, -1, num);
      }
  }

  capt_cnt[stm][2 + wdl] = cnt;

  F = file_open_write(str);
  file_write(&capt_cnt[stm][2 + wdl], 8, F);
  fclose(F);
  file_rename(str);

  snprintf(phase, sizeof phase, "%s %s captures: %lu", side[stm],
      wdl_name[2 + wdl], cnt);
  show_progress(phase, 240, 240, true);
//  printf("capt_cnt[%d][%d] = %lu\n", stm, 2 + wdl, cnt);
}

// Calculate positions with a capture or pawn move >= wdl (wdl = -1 or 0).
// Removing these positions from potential losses allows us to skip
// captures and pawn moves in check_successors().
static void calc_noloss(int stm, int wdl)
{
//  if (file_exists("capt/nobloss/done"))
//    return;

  const char *noloss = wdl == 0 ? "nobloss" : "noloss";

  char pawn_name[64], pcapt_name[64], sub_name[64], psub_name[64];
  strcat(strcpy(pawn_name , "pawn/" ), noloss);
  strcat(strcpy(pcapt_name, "pcapt/"), noloss);
  strcat(strcpy(sub_name  , "sub/"  ), noloss);
  strcat(strcpy(psub_name , "psub/" ), noloss);

  bool partial = dir_exists(-1, stm, noloss), done = partial;

  create_dir(-1, stm, noloss);

  uint64_t num[16];
  struct K16SliceIterator iter;
  int num_done = 0;

  char phase[64];
  snprintf(phase, sizeof phase, "calculating %s captures for %s",
      noloss, side[stm]);

  k16slice_iter_init(&iter, stm);
  int s, s1;
  while (k16slice_iter_next(&iter, &s)) {
    show_progress(phase, num_done++, 240, false);

    while (k16slice_iter_in(&iter, &s1))
      if (!partial || !k16slice_test_count(s1, stm, noloss, -1, num)) {
        done = false;
        if (stm == WHITE)
          k16slice_clear(s1);
        else {
          k16slice_read(s1, s1, BLACK, pcapt_name, -1);
          k16slice_read(-1, s1, BLACK, pawn_name, -1);
          k16slice_or(s1, -1);
        }
      }

    if (!done) {
      k16slice_sub_read(s, s, stm ^ 1, sub_name);
      predecessors_sub(stm, s);

      if (stm == WHITE) {
        k16slice_sub_read(s, s, stm, psub_name);
        predecessors_psub(s);
        calc_king_captures_pawn(s, wdl, 2);
      }
    }

    while (k16slice_iter_out(&iter, &s))
      if (!partial || !k16slice_test_count(s, stm, noloss, -1, num)) {
        // Add illegal positions to avoid having to check legality later. As
        // a bonus, add any faster wins already found to increase efficiency.
        read_wins(-1, s, stm, -1); // most recent "wins"
        k16slice_or(s, -1);
        k16slice_count(s, num);
        k16slice_write(s, s, stm, noloss, -1, num);
      }
  }

  snprintf(phase, sizeof phase, "%s %s captures done", side[stm], noloss);
  show_progress(phase, 240, 240, true);
}

static int probe_promotions(Position *pos, int v)
{
  pos->pt[2] = BQUEEN;
  for (int i = 0; i < 4; i++, pos->pt[2]--)
    if (v < 2)
      v = max(v, -probe_wdl(pos, -2, -v));
  return v;
}

static void calc_pawn_capts_worker(struct ThreadData *thread)
{
  struct IdxState is;

  Position pos = g_pos;
  pos.sq[0] = g_slice.sq[0];
  pos.sq[1] = g_slice.sq[1];
  pos.sq[2] = g_slice.sq[2];
  pos.stm = BLACK;

  size_t offset = k16offset(pos.sq);

  Bitboard occ = idx_state_init(&is, thread->begin, g_slice.sq, &ri, false);

  for (uint64_t idx = thread->begin, end = thread->end; idx < end;
      idx++, occ = idx_state_inc(&is, &ri))
  {
    Bitboard b = pawn_attacks(BLACK, pos.sq[2]) & occ;
    if (!b || !idx_state_legal(&is, BLACK, occ))
      continue;
    idx_state_to_sq(&is, pos.sq, &ri);
    pos.occ = occ;
    int v = -3;
    while (b) {
      int sq = pop_lsb(&b);
      int j = piece_idx(&pos, sq);
      if (pos.pt[j] & 0x08) continue;
      if (do_capture(&pos, g_slice.sq[2], sq, 2, j)) {
        v = sq < 8 ? probe_promotions(&pos, v)
                   : max(v, -probe_wdl(&pos, -2, -v));
        undo_capture(&pos, g_slice.sq[2], sq, 2, j);
      }
    }
    if (v > -3)
      kslice_bit_set(k16slice_buf[v + 2] + offset, idx);
  }
}

static void calc_pawn_capts(void)
{
  char done[16];
  sprintf(done, "pcapt/done");
  FILE *F = fopen(done, "rb");
  if (F) {
    for (int i = 0; i < 5; i++)
      file_read(&pcapt_cnt[i], 8, F);
    fclose(F);
    return;
  }

  char name[7][16];
  for (int i = 0; i < 7; i++) {
    strcat(strcpy(name[i], "pcapt/"), wdl_name[i]);
    create_dir(-1, BLACK, name[i]);
  }

  for (int s = 0; s < 240; s++) {
    show_progress("calculating pawn captures for white", s, 240, false);

    uint64_t c[5], num[16];

    if (k16slice_test_count(s, BLACK, name[6], -1, num)) {
      for (int i = 0; i < 5; i++) {
        k16slice_test_count(s, BLACK, name[i], -1, num);
        c[i] = sum16(num);
      }
    } else {
      for (int i = 0; i < 5; i++)
        k16slice_clear_addr(k16slice_buf[i]);

      for (int r = 0; r < 16; r++) {
        g_slice.sq[0] = KK16Square[s][r][0];
        g_slice.sq[1] = KK16Square[s][r][1];

        if (is_broken(&g_slice))
          continue;

        if (pawn_attacks(BLACK, g_slice.sq[2]) & bit(g_slice.sq[0]))
          continue;

        run_threaded(calc_pawn_capts_worker, &work_g_dynamic);
      }

      uint64_t num[16];

      for (int i = 0; i < 5; i++) {
        c[i] = k16slice_count_addr(k16slice_buf[i], num);
        k16slice_write_addr(k16slice_buf[i], s, BLACK, name[i], -1, num);
      }

      k16slice_or_addr(k16slice_buf[4], k16slice_buf[3]);
      k16slice_or_addr(k16slice_buf[4], k16slice_buf[2]);
      k16slice_write_addr(k16slice_buf[4], s, BLACK, name[5], -1, nullptr);

      k16slice_or_addr(k16slice_buf[4], k16slice_buf[1]);
      k16slice_write_addr(k16slice_buf[4], s, BLACK, name[6], -1, nullptr);
    }

    for (int i = 0; i < 5; i++)
      pcapt_cnt[i] += c[i];
  }

  show_progress("calculating pawn captures for white", 240, 240, true);
}

static void calc_pawn_prom_worker(struct ThreadData *thread)
{
  struct IdxState is;

  Position pos = g_pos;
  pos.sq[0] = g_slice.sq[0];
  pos.sq[1] = g_slice.sq[1];
  pos.sq[2] = g_slice.sq[2];
  pos.stm = WHITE;
  int s = pos.sq[2];

  uint8_t *restrict ilg = k16slice_buf[1] + k16offset(pos.sq);
  uint8_t *restrict p[5];
  for (int i = 0; i < 5; i++)
    p[i] = k16slice_buf[2 + i] + k16offset(pos.sq);

  uint64_t last = thread->begin;
  idx_state_init(&is, last, g_slice.sq, &ri, false);
  for (uint64_t idx = last, end = thread->end; idx < end; idx++) {
    if (kslice_bit_test(ilg, idx))
      continue;
    Bitboard occ = idx_state_add(&is, idx - last, &ri);
    last = idx;
    if (!(occ & bit(s - 8))) {
      pos.sq[2] = s - 8;
      occ ^= bit(s) ^ bit(s - 8);
      idx_state_to_sq(&is, pos.sq, &ri);
      pos.occ = occ;
      if (!opp_king_attacked(&pos)) {
        int v = probe_promotions(&pos, -2);
        kslice_bit_set(p[v + 2], idx);
      }
      pos.sq[2] = s;
    }
  }
}

static int merged_to_wdl[9] = {
  4, 3, 2, 1, 0, 3, 2, 1, 0
};

INLINE uint64_t pawn_push_idx(struct IdxState *is, int from, int to)
{
  Bitboard bb0 = is->bb[0];
  is->bb[0] = bb0 ^ bit(from) ^ bit(to);
  uint64_t idx = rank_bb(is->bb, &ri);
  is->bb[0] = bb0;
  return idx;
}

static void calc_pawn_push_worker(struct ThreadData *thread)
{
  struct IdxState is;
  int s = g_slice.sq[2];

  uint8_t *restrict ilg = k16slice_buf[1] + k16offset(g_slice.sq);
  uint8_t *restrict p[5];
  for (int i = 0; i < 5; i++)
    p[i] = k16slice_buf[2 + i] + k16offset(g_slice.sq);

  uint64_t last = thread->begin;
  idx_state_init(&is, last, g_slice.sq, &ri, false);
  for (uint64_t idx = last, end = thread->end; idx < end; idx++) {
    if (kslice_bit_test(ilg, idx))
      continue;
    Bitboard occ = idx_state_add(&is, idx - last, &ri);
    last = idx;
    if (!(occ & bit(s - 8))) {
      Bitboard occ2 = occ ^ bit(s) ^ bit(s - 8);
      if (idx_state_legal(&is, WHITE, occ2)) {
        uint64_t idx2 = pawn_push_idx(&is, s, s - 8);
        int v = merged_to_wdl[merged_table[idx2]];
        kslice_bit_set(p[v], idx);
      }
    }
  }
}

INLINE char fl(int sq)
{
  return 'a' + (sq & 7);
}

INLINE char rk(int sq)
{
  return '1' + (sq >> 3);
}

static void calc_pawn_push(void)
{
  char str[64];
  uint8_t dummy;

  if (g_slice.sq[2] - 8 == g_slice.sq[0] || g_slice.sq[2] - 8 == g_slice.sq[1])
    return;

  int q = PawnFlip[0][g_slice.sq[2] - 8];
  sprintf(str, "../%s/merged/wdl/w/%c%c%c%c", pawnstr[q], fl(g_slice.sq[0]),
      rk(g_slice.sq[0]), fl(g_slice.sq[1]), rk(g_slice.sq[1]));
  FILE *F = file_open_read(str);
  file_read(&dummy, 1, F);
  read_data(F, merged_table, kslice_size);
  fclose(F);

  run_threaded(calc_pawn_push_worker, &work_g_dynamic);
}

static void calc_pawn_double_push_worker(struct ThreadData *thread)
{
  struct IdxState is;
  int s = g_slice.sq[2];

  uint8_t *restrict ilg = k16slice_buf[1] + k16offset(g_slice.sq);
  uint8_t *restrict p[5];
  for (int i = 0; i < 5; i++)
    p[i] = k16slice_buf[2 + i] + k16offset(g_slice.sq);

  uint64_t last = thread->begin;
  idx_state_init(&is, last, g_slice.sq, &ri, false);
  for (uint64_t idx = last, end = thread->end; idx < end; idx++) {
    if (kslice_bit_test(ilg, idx))
      continue;
    Bitboard occ = idx_state_add(&is, idx - last, &ri);
    last = idx;
    if (!(occ & bit(s - 8))) {
      int v = -1;
      Bitboard occ8 = occ ^ bit(s) ^ bit(s - 8);
      if (idx_state_legal(&is, WHITE, occ8)) {
        uint64_t idx2 = pawn_push_idx(&is, s, s - 8);
        v = merged_to_wdl[merged_table[idx2]];
      }
      if (!(occ8 & bit(s - 16))) {
        Bitboard occ16 = occ ^ bit(s) ^ bit(s - 16);
        if (idx_state_legal(&is, WHITE, occ16)) {
          uint64_t idx2 = pawn_push_idx(&is, s, s - 16);
          v = max(v, merged_to_wdl[merged_table2[idx2]]);
        }
      }
      if (v >= 0)
        kslice_bit_set(p[v], idx);
    }
  }
}

static void calc_pawn_double_push(void)
{
  char str[64];
  uint8_t dummy;

  if (g_slice.sq[2] - 8 == g_slice.sq[0] || g_slice.sq[2] - 8 == g_slice.sq[1])
    return;

  int q = PawnFlip[0][g_slice.sq[2] - 8];
  sprintf(str, "../%s/merged/wdl/w/%c%c%c%c", pawnstr[q], fl(g_slice.sq[0]),
      rk(g_slice.sq[0]), fl(g_slice.sq[1]), rk(g_slice.sq[1]));
  FILE *F = file_open_read(str);
  file_read(&dummy, 1, F);
  read_data(F, merged_table, kslice_size);
  fclose(F);

  if (   g_slice.sq[2] - 16 == g_slice.sq[0]
      || g_slice.sq[2] - 16 == g_slice.sq[1])
  {
    run_threaded(calc_pawn_push_worker, &work_g_dynamic);
    return;
  }

  q = PawnFlip[0][g_slice.sq[2] - 16];
  sprintf(str, "../%s/merged/wdl/w/%c%c%c%c", pawnstr[q], fl(g_slice.sq[0]),
      rk(g_slice.sq[0]), fl(g_slice.sq[1]), rk(g_slice.sq[1]));
  F = file_open_read(str);
  file_read(&dummy, 1, F);
  read_data(F, merged_table2, kslice_size);
  fclose(F);

  run_threaded(calc_pawn_double_push_worker, &work_g_dynamic);
}

static void predecessors_worker(struct ThreadData *thread)
{
  struct IdxState is;
  int stm = g_slice.stm;

  uint64_t *restrict p = (uint64_t *)kslice_get_address_scratch(g_slice.sq);
  uint8_t *restrict const q = kslice_get_address(g_slice.sq);

  p += thread->begin >> 6;
  uint64_t last = thread->begin;
  idx_state_init(&is, last, g_slice.sq, &ri, false);
  for (uint64_t idx = last, end = thread->end; idx < end; idx += 64) {
    uint64_t w = *p++;
    while (w) {
      uint64_t cur = idx + pop_lsb(&w);
      Bitboard occ = idx_state_add(&is, cur - last, &ri);
      last = cur;
      mark_king_unmoves(stm, occ, &is);
      mark_unmoves(stm, KNIGHT, q, occ, &is);
      mark_unmoves(stm, BISHOP, q, occ, &is);
      mark_unmoves(stm, ROOK  , q, occ, &is);
      mark_unmoves(stm, QUEEN , q, occ, &is);
    }
  }
}

// Calculate stm predecessors of stm^1 positions in scratch.
// TODO: use k16slice_read_count[16] to skip empty kslices
static void predecessors(int stm, int s)
{
  g_slice.stm = stm;
  for (int r = 0; r < 16; r++) {
    g_slice.sq[0] = KK16Square[s][r][0];
    g_slice.sq[1] = KK16Square[s][r][1];

    if (is_broken(&g_slice))
      continue;

    run_threaded(predecessors_worker, &work_g_dynamic);
  }
}

static void check_successors_worker(struct ThreadData *thread)
{
  struct IdxState is;
  int stm = g_slice.stm;
  uint64_t cnt = 0;

  uint64_t *restrict p = (uint64_t *)kslice_get_address_scratch(g_slice.sq);
  uint8_t *restrict const q = kslice_get_address(g_slice.sq);

  p += thread->begin >> 6;
  uint64_t last = thread->begin;
  idx_state_init(&is, last, g_slice.sq, &ri, false);

  for (uint64_t idx = last, end = thread->end; idx < end; idx += 64, p++) {
    uint64_t todo = *p;
    if (!todo) continue;
    uint64_t kept = todo;
    while (todo) {
      unsigned bt = pop_lsb(&todo);
      uint64_t cur = idx + bt;
      Bitboard occ = idx_state_add(&is, cur - last, &ri);
      last = cur;
      if (   !check_moves(stm, QUEEN , q, occ, &is)
          || !check_moves(stm, ROOK  , q, occ, &is)
          || !check_moves(stm, BISHOP, q, occ, &is)
          || !check_moves(stm, KNIGHT, q, occ, &is)
          || !check_king_moves(stm, occ, &is))
        kept ^= bit(bt);
    }
    *p = kept;
    cnt += popcnt(kept);
  }

  thread->cnt += cnt;
}

// Verify stm positions as loss against stm^1 positions and return the
// number of verified losses.
static uint64_t check_successors(int stm, int s, uint64_t num[16])
{
  g_slice.stm = stm;
  for (int r = 0; r < 16; r++) {
    num[r] = 0;

    g_slice.sq[0] = KK16Square[s][r][0];
    g_slice.sq[1] = KK16Square[s][r][1];

    if (is_broken(&g_slice))
      continue;

    for (int t = 0; t < g_num_threads; t++)
      g_thread_data[t].cnt = 0;

    run_threaded(check_successors_worker, &work_g_dynamic);

    for (int t = 0; t < g_num_threads; t++)
      num[r] += g_thread_data[t].cnt;
  }

  uint64_t cnt = 0;
  for (int r = 0; r < 16; r++)
    cnt += num[r];

  return cnt;
}

static void calc_illegal_worker_tmpl(struct ThreadData *thread, const int pc)
{
  struct IdxState is;
  int stm = g_slice.stm;
  int k = g_piece_set[stm][pc];
  int ksq = g_slice.sq[stm ^ 1];

  uint8_t *restrict const p = k16slice_buf[stm] + k16offset(g_slice.sq);

  Bitboard occ = idx_state_init(&is, thread->begin, g_slice.sq, &capt_ri[k],
      false);

  for (uint64_t idx = thread->begin, end = thread->end; idx < end;
      idx++, occ = idx_state_inc(&is, &capt_ri[k]))
  {
    mark_uncapture_king_or_pawn(stm, pc, ksq, p, occ, &is);
  }
}

static void calc_illegal_knight_worker(struct ThreadData *thread)
{
  calc_illegal_worker_tmpl(thread, KNIGHT);
}

static void calc_illegal_bishop_worker(struct ThreadData *thread)
{
  calc_illegal_worker_tmpl(thread, BISHOP);
}

static void calc_illegal_rook_worker(struct ThreadData *thread)
{
  calc_illegal_worker_tmpl(thread, ROOK);
}

static void calc_illegal_queen_worker(struct ThreadData *thread)
{
  calc_illegal_worker_tmpl(thread, QUEEN);
}

static void calc_mate_worker(struct ThreadData *thread)
{
  struct IdxState is;

  uint64_t *restrict p0 = (uint64_t *)(k16slice_buf[0] + k16offset(g_slice.sq));
  uint64_t *restrict p1 = (uint64_t *)(k16slice_buf[1] + k16offset(g_slice.sq));
  uint64_t *restrict q0 = (uint64_t *)(k16slice_buf[2] + k16offset(g_slice.sq));
  uint64_t *restrict q1 = (uint64_t *)(k16slice_buf[3] + k16offset(g_slice.sq));

  uint64_t last = thread->begin;
  p0 += last >> 6;
  p1 += last >> 6;
  q0 += last >> 6;
  q1 += last >> 6;
  idx_state_init(&is, last, g_slice.sq, &ri, false);
  for (uint64_t idx = last, end = thread->end; idx < end;
      idx += 64, p0++, p1++, q0++, q1++)
  {
    uint64_t w = *p0 ^ *p1;
    if (!w) continue;
    uint64_t white = 0, black = 0;
    while (w) {
      unsigned bt = pop_lsb(&w);
      uint64_t cur = idx + bt;
      Bitboard occ = idx_state_add(&is, cur - last, &ri);
      last = cur;
      if (*p1 & bit(bt)) {
        if (idx_state_mate(&is, WHITE, occ))
          white |= bit(bt);
      } else {
        if (idx_state_mate(&is, BLACK, occ))
          black |= bit(bt);
      }
    }
    *q0 = white;
    *q1 = black;
  }
}

// Calc illegal, mate (L0) and pawn-push positions.
static void calc_illegal_and_mate_and_pawn_push(void)
{
  FILE *F = fopen("0/done", "rb");
  if (F) {
    for (int stm = 0; stm < 2; stm++) {
      file_read(&g_stats[stm][0], 8, F);
      file_read(&g_stats[stm][MAX_STATS - 1], 8, F);
    }
    file_read(&pawn_cnt, sizeof pawn_cnt, F);
    fclose(F);
    return;
  }

  bool partial = file_exists("0");

  for (int stm = 0; stm < 2; stm++) {
    create_dir(0, stm, "wins");
    create_dir(0, stm, "L");
  }

  char name[7][16];
  for (int i = 0; i < 7; i++) {
    if (i == 2) continue;
    strcat(strcpy(name[i], "pawn/"), wdl_name[i]);
    create_dir(-1, BLACK, name[i]);
  }

  uint64_t broken[2] = { 0 }, loss0[2] = { 0 }, num[16];

  for (int s = 0; s < 240; s++) {
    show_progress("calculating illegals, mates and pawn pushes", s, 240, false);
    if (partial) {
      char str[64];
      create_name(str, s, BLACK, name[6], -1);
      if (file_exists(str)) {
        for (int stm = 0; stm < 2; stm++) {
          if (k16slice_test_count(s, stm, "wins", 0, num))
            broken[stm] += sum16(num);
          if (k16slice_test_count(s, stm, "L", 0, num))
            loss0[stm] += sum16(num);
        }
        for (int i = 0; i < 5; i++)
          if (k16slice_test_count(s, BLACK, name[i], -1, num))
            pawn_cnt[i] += sum16(num);
        continue;
      }
    }

    k16slice_clear_addr(k16slice_buf[0]); // wtm illegal
    k16slice_clear_addr(k16slice_buf[1]); // btm illegal

    for (int r = 0; r < 16; r++) {
      g_slice.sq[0] = KK16Square[s][r][0];
      g_slice.sq[1] = KK16Square[s][r][1];

      if (is_broken(&g_slice))
        continue;

      g_slice.stm = WHITE;
      int k;
      if ((k = g_piece_set[WHITE][KNIGHT]) >= 0)
        run_threaded(calc_illegal_knight_worker, &work_capt_dynamic[k]);
      if ((k = g_piece_set[WHITE][BISHOP]) >= 0)
        run_threaded(calc_illegal_bishop_worker, &work_capt_dynamic[k]);
      if ((k = g_piece_set[WHITE][ROOK]) >= 0)
        run_threaded(calc_illegal_rook_worker, &work_capt_dynamic[k]);
      if ((k = g_piece_set[WHITE][QUEEN]) >= 0)
        run_threaded(calc_illegal_queen_worker, &work_capt_dynamic[k]);

      if (pawn_attacks(BLACK, g_slice.sq[2]) & bit(g_slice.sq[0])) {
        kslice_set_addr(k16slice_buf[1] + k16offset(g_slice.sq));
        continue;
      }

      g_slice.stm = BLACK;

      if ((k = g_piece_set[BLACK][KNIGHT]) >= 0)
        run_threaded(calc_illegal_knight_worker, &work_capt_dynamic[k]);
      if ((k = g_piece_set[BLACK][BISHOP]) >= 0)
        run_threaded(calc_illegal_bishop_worker, &work_capt_dynamic[k]);
      if ((k = g_piece_set[BLACK][ROOK]) >= 0)
        run_threaded(calc_illegal_rook_worker, &work_capt_dynamic[k]);
      if ((k = g_piece_set[BLACK][QUEEN]) >= 0)
        run_threaded(calc_illegal_queen_worker, &work_capt_dynamic[k]);
    }

    for (int stm = 0; stm < 2; stm++) {
      broken[stm] += k16slice_count_addr(k16slice_buf[stm], num);
      k16slice_write_addr(k16slice_buf[stm], s, stm, "wins", 0, num);
    }

    k16slice_clear_addr(k16slice_buf[2]); // wtm mate
    k16slice_clear_addr(k16slice_buf[3]); // btm mate

    for (int r = 0; r < 16; r++) {
      g_slice.sq[0] = KK16Square[s][r][0];
      g_slice.sq[1] = KK16Square[s][r][1];

      if (is_broken(&g_slice))
        continue;

      run_threaded(calc_mate_worker, &work_g_static);
    }

    for (int stm = 0; stm < 2; stm++) {
      loss0[stm] += k16slice_count_addr(k16slice_buf[2 + stm], num);
      k16slice_write_addr(k16slice_buf[2 + stm], s, stm, "L", 0, num);
    }

    for (int i = 0; i < 5; i++)
      k16slice_clear_addr(k16slice_buf[2 + i]);

    for (int r = 0; r < 16; r++) {
      g_slice.sq[0] = KK16Square[s][r][0];
      g_slice.sq[1] = KK16Square[s][r][1];

      if (is_broken(&g_slice))
        continue;

      if (g_slice.sq[2] < 16)
        run_threaded(calc_pawn_prom_worker, &work_g_dynamic);
      else if (g_slice.sq[2] < 48)
        calc_pawn_push();
      else
        calc_pawn_double_push();
    }

    for (int i = 0; i < 5; i++) {
      if (i == 2) continue;
      pawn_cnt[i] += k16slice_count_addr(k16slice_buf[2 + i], num);
      k16slice_write_addr(k16slice_buf[2 + i], s, BLACK, name[i], -1, num);
    }

    k16slice_or_addr(k16slice_buf[6], k16slice_buf[5]);
    k16slice_or_addr(k16slice_buf[6], k16slice_buf[4]);
    k16slice_write_addr(k16slice_buf[6], s, BLACK, name[5], -1, nullptr);
    k16slice_or_addr(k16slice_buf[6], k16slice_buf[3]);
    k16slice_write_addr(k16slice_buf[6], s, BLACK, name[6], -1, nullptr);
  }

  F = file_open_write("0/done");
  for (int stm = 0; stm < 2; stm++) {
    g_stats[stm][0] = broken[stm];
    g_stats[stm][MAX_STATS - 1] = loss0[stm];
    file_write(&g_stats[stm][0], 8, F);
    file_write(&g_stats[stm][MAX_STATS - 1], 8, F);
  }
  file_write(&pawn_cnt, sizeof pawn_cnt, F);
  fclose(F);
  file_rename("0/done");

  char phase[128];
  snprintf(phase, sizeof phase, "illegal: %lu, %lu; mate: %lu, %lu",
      broken[WHITE], broken[BLACK], loss0[WHITE], loss0[BLACK]);
  show_progress(phase, 240, 240, true);
//  for (int i = 0; i < 5; i++)
//    printf("pawn_cnt[%d] = %lu\n", i, pawn_cnt[i]);
}

// Calculate stm losses in n from stm^1 wins in n-1 (n > 1) or
// from stm^1 wins in sub tables reached through captures (n == 1).
static bool calc_L(int stm, int n, bool more_l)
{
  char str[64];
  sprintf(str, "%d/L/%c/done", n, "wb"[stm]);
  FILE *F = fopen(str, "rb");
  if (F) {
    file_read(&g_stats[stm][MAX_STATS - 1 - n], 8, F);
    fclose(F);
    return g_stats[stm][MAX_STATS - 1 - n] != 0;
  }

  char phase[64];

  struct K16SliceIterator iter;
  bool partial = true;
  int num_done = 0;

  if (dir_exists(n, stm, "L"))
    goto skip_X;

  partial = dir_exists(n, stm, "X");
  bool done = partial;
  uint64_t num[16];

  create_dir(n, stm, "X");

  snprintf(phase, sizeof phase, "%d/X/%c", n, "wb"[stm]);

  // Calculate potential losses in n = predecessors(W(n-1))
  k16slice_iter_init(&iter, stm);
  int s, s1;
  while (k16slice_iter_next(&iter, &s)) {
    show_progress(phase, num_done++, 240, false);

    bool pred = more_l && k16slice_test(s, stm ^ 1, "W", n - 1);

    if (pred || n == 1 || n == DRAW_RULE + 1) {
      while (k16slice_iter_in(&iter, &s1))
        if (!partial || !k16slice_test_count(s1, stm, "X", n, num)) {
          done = false;
          k16slice_clear(s1);
        }

      if (pred && !done) {
        k16slice_read(-1, s, stm ^ 1, "W", n - 1);
        predecessors(stm, s);
      }
    }

    while (k16slice_iter_out(&iter, &s))
      if (!partial || !k16slice_test_count(s, stm, "X", n, num)) {
        // Add any potential losses by capture or pawn push.
        if (n == 1) {
          k16slice_read(-1, s, stm, "capt/loss", -1);
          k16slice_or(s, -1);
          if (stm == BLACK) {
            k16slice_read(-1, s, stm, "pawn/loss", -1);
            k16slice_or(s, -1);
          }
        }
        else if (n == DRAW_RULE + 1) {
          k16slice_read(-1, s, stm , "capt/bloss", -1);
          k16slice_or(s, -1);
          if (stm == BLACK) {
            k16slice_read(-1, s, stm, "pawn/bloss", -1);
            k16slice_or(s, -1);
          }
        }

        // Remove positions with non-losing captures or pawn moves.
        k16slice_read(-1, s, stm, n <= DRAW_RULE ? "noloss" : "nobloss", -1);
        k16slice_andnot(s, -1);
        k16slice_clear_tail(s);
        k16slice_write(s, s, stm, "X", n, nullptr);

        if (n == 1) {
          k16slice_delete(s, stm, "capt/loss", -1);
          if (stm == BLACK)
            k16slice_delete(s, stm, "pawn/loss", -1);
        }
        else if (n == DRAW_RULE + 1 && stm == BLACK)
          k16slice_delete(s, stm, "pawn/bloss", -1);
      }
  }
  show_progress(phase, 240, 240, false);

  create_dir(n, stm, "L");
  partial = false;

skip_X:

  uint64_t cnt = 0;
  num_done = 0;
  snprintf(phase, sizeof phase, "\x1b[%sm%d/L/%c\x1b[0m",
      clr_L[(2 * n + stm) & 3], n, "wb"[stm]);

  // Verify potential losses.
  k16slice_iter_init(&iter, stm);
  while (k16slice_iter_next(&iter, &s)) {
    show_progress(phase, num_done++, 240, false);

    if (k16slice_test(s, stm , "X", n)) {
      while (k16slice_iter_in(&iter, &s1))
        read_wins(s1, s1, stm ^ 1, n - 1);

      k16slice_read(-1, s, stm, "X", n);
      cnt += check_successors(stm, s, num);
      k16slice_write(-1, s, stm, "L", n, num);
      if (0 && stm == BLACK && n == 4) {
//        list_positions(s, k16slice_buf[11]);
//        printf("s = %d, cnt = %lu\n", s, cnt);
      }
    } else if (partial && k16slice_test_count(s, stm, "L", n, num))
      cnt += sum16(num);

    while (k16slice_iter_out(&iter, &s))
      k16slice_delete(s, stm, "X", n);
  }

  g_stats[stm][MAX_STATS - 1 - n] = cnt;

  F = file_open_write(str);
  file_write(&g_stats[stm][MAX_STATS - 1 - n], 8, F);
  fclose(F);
  file_rename(str);

  snprintf(phase, sizeof phase, "\x1b[%sm%d/L/%c  %lu\x1b[0m",
      clr_L[(2 * n + stm) & 3], n, "wb"[stm], cnt);
  show_progress(phase, 240, 240, true);

  return cnt != 0;
}

static bool calc_W(int stm, int n, bool more_w)
{
  char str[64];
  sprintf(str, "%d/W/%c/done", n, "wb"[stm]);
  FILE *F = fopen(str, "rb");
  if (F) {
    if (n == 1) {
      file_read(&g_stats[stm][1], 8, F);
      if (stm == BLACK)
        file_read(&g_stats[stm][2], 8, F);
    }
    else if (n == DRAW_RULE + 1)
      file_read(&g_stats[stm][DRAW_RULE + 3], 8, F);
    file_read(&g_stats[stm][stats_n(n)], 8, F);
    fclose(F);
    return g_stats[stm][stats_n(n)] != 0;
  }

  struct K16SliceIterator iter;
  uint64_t cnt = 0, cnt_w = 0, cnt_pw = 0, num[16];
  int num_done = 0;

  bool partial = dir_exists(n, stm, "W"), done = partial;

  create_dir(n, stm, "W");
  create_dir(n, stm, "wins");

  char phase[64];
  snprintf(phase, sizeof phase, "\x1b[%sm%d/W/%c\x1b[0m",
      clr_W[(2 * n + stm) & 3], n, "wb"[stm]);

  // Calculate wins in n = predecessors(L(n-1))
  k16slice_iter_init(&iter, stm);
  int s, s1;
  while (k16slice_iter_next(&iter, &s)) {
    show_progress(phase, num_done++, 240, false);

    bool pred = more_w && k16slice_test(s, stm ^ 1, "L", n - 1);

    if (pred || n == 1 || n == DRAW_RULE + 1) {
      while (k16slice_iter_in(&iter, &s1))
        if (!partial || !k16slice_test_count(s1, stm, "W", n, num)) {
          done = false;
          if (n == 1) {
            // Add any wins by capture or pawn push.
            k16slice_read(s1, s1, stm, "capt/win", -1);
            // Remove illegal positions to count wins by capture.
            // FIXME: should probable do k16slice_read_andnot()
            k16slice_read(-1, s1, stm, "wins", 0);
            k16slice_andnot(s1, -1);
            cnt_w += k16slice_count(s1, num);
            if (stm == BLACK) {
              // Count wins by pawn push. These are all legal already.
              k16slice_read(-1, s1, stm, "pawn/win", -1);
              k16slice_or(s1, -1);
              cnt_pw += k16slice_count(s1, num);
            }
          }
          else if (n == DRAW_RULE + 1) {
            k16slice_read(s1, s1, stm, "capt/cwin", -1);
            if (stm == BLACK) {
              k16slice_read(-1, s1, stm, "pawn/cwin", -1);
              k16slice_or(s1, -1);
            }
            // FIXME: this might be wrong ("wins", 0)
            k16slice_read(-1, s1, stm, "wins", 0);
            k16slice_andnot(s1, -1);
            cnt_w += k16slice_count(s1, num);
          }
          else
            k16slice_clear(s1);
        }

      if (pred && !done) {
        k16slice_read(-1, s, stm ^ 1, "L", n - 1);
        predecessors(stm, s);
      }
    }

    while (k16slice_iter_out(&iter, &s))
      if (!partial || !k16slice_test_count(s, stm, "W", n, num)) {
        // Remove illegal positions and known faster wins.
        read_wins(-1, s, stm, n - 1);
        k16slice_andnot(s, -1);
        cnt += k16slice_count(s, num);
        uint64_t cost = k16slice_write(s, s, stm, "W", n, num);
        add_wins(s, stm, n, cost);

        if (0 && stm == WHITE && n == 1) {
//        list_positions(s, k16slice_get_address(s));
//        printf("s = %d, cnt = %lu\n", s, cnt);
        }
      }
  }

  F = file_open_write(str);
  if (n == 1) {
//    printf("capt_win_%c = %lu\n", "wb"[stm], cnt_w);
    g_stats[stm][1] = cnt_w;
    file_write(&g_stats[stm][1], 8, F);
    if (stm == BLACK) {
      g_stats[stm][2] = cnt_pw - cnt_w;
      file_write(&g_stats[stm][2], 8, F);
//      printf("pawn_win_%c = %lu\n", "wb"[stm], cnt_pw - cnt_w);
    }
  }
  else if (n == DRAW_RULE + 1) {
    g_stats[stm][DRAW_RULE + 3] = cnt_w;
    file_write(&g_stats[stm][DRAW_RULE + 3], 8, F);
//    printf("capt/pawn_cwin_%c = %lu\n", "wb"[stm], cnt_w);
  }
  g_stats[stm][stats_n(n)] = cnt;
  file_write(&g_stats[stm][stats_n(n)], 8, F);
  fclose(F);
  file_rename(str);

  snprintf(phase, sizeof phase, "\x1b[%sm%d/W/%c  %lu\x1b[0m",
      clr_W[(2 * n + stm) & 3], n, "wb"[stm], cnt);
  show_progress(phase, 240, 240, true);

  return cnt != 0;
}

// 0 = ILLEGAL
// 1 = CAPT_WIN
// 2 = PAWN_WIN
// 1+2 ... DRAW_RULE+2
// DRAW_RULE+3 = CAPT_CWIN
// DRAW_RULE+4 = PAWN_CWIN

void generate(void)
{
  FILE *F = fopen("generate_info", "rb");
  if (F) {
    read_data(F, &g_stats, sizeof g_stats);
    file_read(&capt_cnt, sizeof capt_cnt, F);
    file_read(&pawn_cnt, sizeof pawn_cnt, F);
    file_read(&max_iteration, sizeof max_iteration, F);
    fclose(F);
    printf("Table for pawn slice %s was already generated.\n",
        pawnstr[PawnFlip[0][g_slice.sq[2]]]);
    return;
  }

  printf("Generating table for pawn slice %s.\n",
      pawnstr[PawnFlip[0][g_slice.sq[2]]]);

  for (int i = 0; i < 7; i++)
    k16slice_buf[i] = alloc_k16slice();

  if (g_slice.sq[2] >= 16) {
    merged_table = alloc_huge(kslice_size);
    if (!merged_table)
      out_of_mem();
    if (g_slice.sq[2] >= 48) {
      merged_table2 = alloc_huge(kslice_size);
      if (!merged_table2)
        out_of_mem();
    }
  }

  memset(&pawn_cnt, 0, sizeof pawn_cnt);
  memset(&psub_cnt, 0, sizeof psub_cnt);
  memset(&sub_cnt, 0, sizeof sub_cnt);
  memset(&pcapt_cnt, 0, sizeof pcapt_cnt);
  memset(&capt_cnt, 0, sizeof capt_cnt);

  memset(&wins_full, 0, sizeof wins_full);
  memset(&wins_checked, 0, sizeof wins_checked);
  memset(&replay_cost, 0, sizeof replay_cost);

  calc_illegal_and_mate_and_pawn_push();

  if (merged_table) {
    free(merged_table);
    merged_table = nullptr;
    if (merged_table2) {
      free(merged_table2);
      merged_table2 = nullptr;
    }
  }

  kslice_alloc_buffers();

  // Separate psub files or part of sub??
  calc_sub_kslices(WHITE);
  calc_sub_kslices(BLACK);
  calc_psub_kslices();

  calc_pawn_capts();

  // Calculate positions with winning captures.
  calc_capt(WHITE, 2);
  calc_capt(BLACK, 2);

  // Calculate positions with losing captures.
  // (Only used to find loss-in-1 positions.)
  calc_capt(WHITE, -2);
  calc_capt(BLACK, -2);

  // Calculate noloss positions, i.e. positions with a capture or pawn move
  // that gives at least a blessed loss.
  calc_noloss(WHITE, -1);
  calc_noloss(BLACK, -1);

  bool more_ww = true, more_wb = true, more_lw = true, more_lb = true;
  bool more_wb_next, more_ww_next;

  int n;
  for (n = 1; n <= DRAW_RULE; n++) {
    more_wb_next = more_lw && calc_L(WHITE, n, more_lw);
    more_ww_next = more_lb && calc_L(BLACK, n, more_lb);

    more_lb = more_ww && calc_W(WHITE, n, more_ww);
    more_lw = more_wb && calc_W(BLACK, n, more_wb);

    more_wb = more_wb_next;
    more_ww = more_ww_next;
  }

  // Calculate positions with a capture that gives a cursed win.
  calc_capt(WHITE, 1);
  calc_capt(BLACK, 1);

  // Calculate positions with a capture that gives a blessed loss.
  calc_capt(WHITE, -1);
  calc_capt(BLACK, -1);

  // Calculate nobloss positions, i.e. positions with a capture or pawn move
  // that gives at least a draw.
  calc_noloss(WHITE, 0);
  calc_noloss(BLACK, 0);

  more_wb_next = calc_L(WHITE, n, more_lw);
  more_ww_next = calc_L(BLACK, n, more_lb);

  more_lb = calc_W(WHITE, n, more_ww);
  more_lw = calc_W(BLACK, n, more_wb);

  more_wb = more_wb_next;
  more_ww = more_ww_next;

  for (n++; more_ww || more_wb || more_lw || more_lb; n++) {
    more_wb_next = more_lw && calc_L(WHITE, n, more_lw);
    more_ww_next = more_lb && calc_L(BLACK, n, more_lb);

    more_lb = more_ww && calc_W(WHITE, n, more_ww);
    more_lw = more_wb && calc_W(BLACK, n, more_wb);

    more_wb = more_wb_next;
    more_ww = more_ww_next;
  }

  // Calculate positions with a drawing capture.
  // During WDL merging, these positions are set to partial don't care.
  calc_capt(WHITE, 0);
  calc_capt(BLACK, 0);

  max_iteration = n;

  uint64_t num_kslices = 0;
  for (int s = 0; s < 240; s++)
    for (int r = 0; r < 16; r++) {
      g_slice.sq[0] = KK16Square[s][r][0];
      g_slice.sq[1] = KK16Square[s][r][1];
      num_kslices += !is_broken(&g_slice);
    }

  for (int stm = 0; stm < 2; stm++) {
    // Remove some double counting.
    g_stats[stm][3] -= g_stats[stm][1] + g_stats[stm][2];
    g_stats[stm][DRAW_RULE + 5] -=
      g_stats[stm][DRAW_RULE + 3] - g_stats[stm][DRAW_RULE + 4];
    uint64_t total = 0;
    for (int i = 0; i < MAX_STATS; i++)
      total += g_stats[stm][i];
    g_stats[stm][MAX_STATS / 2 + 2] = num_kslices * kslice_size - total;
  }

  F = file_open_write("generate_info");
  write_data(F, &g_stats, sizeof g_stats);
  file_write(&capt_cnt, sizeof capt_cnt, F);
  file_write(&pawn_cnt, sizeof pawn_cnt, F);
  file_write(&max_iteration, sizeof max_iteration, F);
  fclose(F);
  file_rename("generate_info");

  kslice_free_buffers();
}

void delete_intermediate_slices(void)
{
  if (file_exists("0/wins")) {
    for (int stm = 0; stm < 2; stm++)
      for (int s = 0; s < 240; s++) {
        k16slice_delete(s, stm, "wins", wins_full[stm][s]);
        k16slice_delete(s, stm, "nobloss", -1);
        k16slice_delete(s, stm, "noloss", -1);
      }
    for (int n = max_iteration - 1; n >= 0; n--)
      delete_dir(n, "wins");
    delete_dir(-1, "nobloss");
    delete_dir(-1, "noloss");
  }

  if (file_exists("sub")) {
    for (int i = 0; i < 7; i++) {
      char name[64];
      sprintf(name, "sub/%s", wdl_name[i]);
      for (int stm = 0; stm < 2; stm++)
        for (int s = 0; s < 240; s++)
          k16slice_delete(s, stm, name, -1);
      delete_dir(-1, name);
      sprintf(name, "psub/%s", wdl_name[i]);
      for (int s = 0; s < 240; s++)
        k16slice_delete(s, WHITE, name, -1);
      delete_dir(-1, name);
      sprintf(name, "pcapt/%s", wdl_name[i]);
      for (int s = 0; s < 240; s++)
        k16slice_delete(s, BLACK, name, -1);
      delete_dir(-1, name);
    }
    unlink("sub/done_w");
    unlink("sub/done_b");
    rmdir("sub");
    unlink("psub/done");
    rmdir("psub");
    rmdir("pcapt");
  }

  if (file_exists("pawn/loss")) {
    for (int s = 0; s < 240; s++) {
      k16slice_delete(s, BLACK, "pawn/nobloss", -1);
      k16slice_delete(s, BLACK, "pawn/noloss", -1);
    }
    delete_dir(-1, "pawn/nobloss");
    delete_dir(-1, "pawn/noloss");
  }
}

// Cleanup the files that remain after merge()ing.
void cleanup_generation(void)
{
  char dir[64], file[64];
  for (int n = 0; n < max_iteration; n++) {
    sprintf(dir, "%d", n);
    if (!file_exists(dir)) continue;
    delete_dir(n, "L");
    delete_dir(n, "X");
    delete_dir(n, "W");
    delete_dir(n, "wins");
    sprintf(file, "%d/done", n);
    unlink(file);
    rmdir(dir);
  }
  for (int i = 0; i < 5; i++) {
    sprintf(dir, "capt/%s", wdl_name[i]);
    delete_dir(-1, dir);
  }
  delete_dir(-1, "capt");
  for (int i = 0; i < 5; i++) {
    sprintf(dir, "pawn/%s", wdl_name[i]);
    delete_dir(-1, dir);
  }
  delete_dir(-1, "pawn");
}
