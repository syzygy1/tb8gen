/*
  Copyright (c) 2026 Ronald de Man

  This file is distributed under the terms of the GNU GPL, version 2.
*/

#include <assert.h>
#include <errno.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#include "decompress.h"
#include "defs.h"
#include "index.h"
#include "kslicep.h"
#include "movegen.h"
#include "probe.h"
#include "stats.h"
#include "tb8gen.h"
#include "threads.h"
#include "types.h"
#include "util.h"
#include "hash/xxhash.h"

const char wdl_name[5][8] = { "loss", "bloss", "draw", "cwin", "win" };
const char side[2][8] = { "white", "black" };

static uint64_t capt_cnt[2][5], sub_cnt[2][5];
static uint64_t psub_cnt[5], pcapt_cnt[5], pawn_cnt[5];
static uint8_t *merged_table, *merged_table2;

static void k16slice_read_addr(void *p, int slice, int stm, const char *name,
    int n);
static void k16slice_read_path_addr(void *p, const char *str);

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

// Uncapture stm^1 king by an stm piece being added to set k.
INLINE void mark_uncapture_king_or_pawn(int stm, const int pc, int ksq,
    uint8_t *restrict const p, Bitboard occ, struct IdxState *is)
{
  int k = g_piece_set[stm][pc];
  if (k < 0) return;

  uint64_t idx0 = 0;
  for (int i = 0; i < k; i++)
    idx0 = idx0 * ri.factor[i] + is->sub[i];

  Bitboard b = piece_moves(pc, ksq, occ);
  while (b) {
    Bitboard from_bb = b & -b;
    is->bb[k + 1] ^= from_bb;
    uint64_t idx = rank_bb_from(is->bb, idx0, k, is->occ[k], &ri);
    kslice_bit_set_atomic(p, idx);
    is->bb[k + 1] ^= from_bb;
    b ^= from_bb;
  }
}

// Return true if one king move hits a set bit.
INLINE bool test_king_moves(int stm, Bitboard occ, struct IdxState *is)
{
  uint8_t ksq[2] = { is->sq[0], is->sq[1] };

  Bitboard b = king_attacks(ksq[stm]) & ~king_attacks(ksq[stm ^ 1]) & ~occ;
  while (b) {
    ksq[stm] = pop_lsb(&b);
    is->bb[0] = bit(ksq[0]) | bit(ksq[1]) | bit(is->sq[2]);
    uint64_t idx = rank_bb(is->bb, &ri);
    uint8_t *p = kslice_get_address(ksq);
    if (kslice_bit_test(p, idx)) {
      is->bb[0] = is->occ[0];
      return true;
    }
  }

  is->bb[0] = is->occ[0];
  return false;
}

INLINE bool test_moves(int stm, const int pc, uint8_t *restrict const p,
    Bitboard occ, struct IdxState *is)
{
  int k = g_piece_set[stm][pc];
  if (k < 0) return false;

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
      if (kslice_bit_test(p, idx)) {
        is->bb[k + 1] ^= from_bb;
        return true;
      }
      b ^= to_bb;
    }
    is->bb[k + 1] ^= from_bb;
  }

  return false;
}

static int work_slice, work_set, work_r;
static int work_lower, work_upper;

static struct Work work_g_dynamic, work_g_static;
static struct Work work_capt_dynamic[MAX_SETS];

static constexpr uint64_t VERIFY_MIN_CHUNK = 1ULL << 9;
static constexpr int VERIFY_DYNAMIC_FACTOR = 4;

void init_verification_work(void)
{
  work_init(&work_g_dynamic, kslice_size, 0x1ff, WORK_DYNAMIC,
      VERIFY_DYNAMIC_FACTOR, VERIFY_MIN_CHUNK);
  work_init(&work_g_static, kslice_size, 0x1ff, WORK_STATIC, 1,
      VERIFY_MIN_CHUNK);

  for (int k = 0; k < ri.numsets; k++)
    work_init(&work_capt_dynamic[k], capt_ri[k].sizes[0], 0x1ff, WORK_DYNAMIC,
        VERIFY_DYNAMIC_FACTOR, VERIFY_MIN_CHUNK);
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
    for (int i = 0; i < 5; i++)
      printf("sub_cnt[%d][%d] = %lu\n", stm, i, sub_cnt[stm][i]);
    return;
  }

  char phase[64];
  snprintf(phase, sizeof phase, "probing subtables for %s", side[stm]);

  char name[5][16];
  for (int i = 0; i < 5; i++) {
    strcat(strcpy(name[i], "sub/"), wdl_name[i]);
    create_dir(-1, stm, name[i]);
  }

  g_slice.stm = stm;

  for (int s = 0; s < 240; s++) {
    show_progress(phase, s, 240, false);

    uint64_t c[5];

    if (k16slice_sub_test_count(s, stm, name[4], -1, &c[0])) {
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
  for (int i = 0; i < 5; i++)
    printf("sub_cnt[%d][%d] = %lu\n", stm, i, sub_cnt[stm][i]);
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

  char name[5][16];
  for (int i = 0; i < 5; i++) {
    strcat(strcpy(name[i], "psub/"), wdl_name[i]);
    create_dir(-1, WHITE, name[i]);
  }

  g_slice.stm = BLACK;

  for (int s = 0; s < 240; s++) {
    show_progress("probing pawn subtables for black", s, 240, false);

    uint64_t c[5];

    if (k16slice_sub_test_count(s, WHITE, name[4], -1, &c[0])) {
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

// Uncapture from the loaded subtable of stm^1 positions to stm positions.
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

static void calc_capt(int stm, int wdl)
{
  if (!sub_cnt[stm ^ 1][2 - wdl])
    return;

  char capt_name[32], pcapt_name[32], sub_name[32], psub_name[32];
  strcat(strcpy(capt_name, "capt/" ), wdl_name[2 + wdl]);
  strcat(strcpy(pcapt_name, "pcapt/"), wdl_name[2 + wdl]);
  strcat(strcpy(sub_name , "sub/"  ), wdl_name[2 - wdl]);
  strcat(strcpy(psub_name , "psub/" ), wdl_name[2 - wdl]);

  char str[64];
  sprintf(str, "%s/%c/done", capt_name, "wb"[stm]);
  FILE *F = fopen(str, "rb");
  uint64_t cnt = 0;
  if (F) {
    file_read(&cnt, 8, F);
    capt_cnt[stm][2 + wdl] = cnt;
    fclose(F);
    goto finished;
  }

  char phase[64];
  snprintf(phase, sizeof phase, "calculating %s captures for %s",
      wdl_name[2 + wdl], side[stm]);

  bool partial = dir_exists(-1, stm, capt_name), done = partial;

  create_dir(-1, stm, capt_name);

  struct K16SliceIterator iter;
  uint64_t num[16];
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

finished:
  snprintf(phase, sizeof phase, "%s %s captures: %lu", side[stm],
      wdl_name[2 + wdl], cnt);
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

  char name[5][16];
  for (int i = 0; i < 5; i++) {
    strcat(strcpy(name[i], "pcapt/"), wdl_name[i]);
    create_dir(-1, BLACK, name[i]);
  }

  for (int s = 0; s < 240; s++) {
    show_progress("calculating pawn captures for white", s, 240, false);

    uint64_t c[5], num[16];

    if (k16slice_test_count(s, BLACK, name[4], -1, num)) {
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
    }

    for (int i = 0; i < 5; i++)
      pcapt_cnt[i] += c[i];
  }

  F = file_open_write(done);
  for (int i = 0; i < 5; i++)
    file_write(&pcapt_cnt[i], 8, F);
  fclose(F);
  file_rename(done);

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

#ifndef NDEBUG
INLINE uint64_t pawn_push_idx_check(struct IdxState *is, int from, int to)
{
  Bitboard bb0 = is->bb[0];
  is->bb[0] = bb0 ^ bit(from) ^ bit(to);
  uint64_t idx = rank_bb(is->bb, &ri);
  is->bb[0] = bb0;
  return idx;
}
#endif

static void calc_pawn_push_worker(struct ThreadData *thread)
{
  struct IdxState is;
  struct PawnIdxState pis;
  int s = g_slice.sq[2];

  uint8_t *restrict ilg = k16slice_buf[1] + k16offset(g_slice.sq);
  uint8_t *restrict p[5];
  for (int i = 0; i < 5; i++)
    p[i] = k16slice_buf[2 + i] + k16offset(g_slice.sq);

  uint64_t last = thread->begin;
  idx_state_init(&is, last, g_slice.sq, &ri, false);
  pawn_idx_init(&pis, bit(s) ^ bit(s - 8), &is, &ri);
  for (uint64_t idx = last, end = thread->end; idx < end; idx++) {
    if (kslice_bit_test(ilg, idx))
      continue;
    Bitboard occ = idx_state_pawn_add(&is, idx - last, &pis, nullptr, &ri);
    last = idx;
    if (occ & bit(s - 8)) continue;
    if (!idx_state_legal(&is, WHITE, occ ^ bit(s) ^ bit(s - 8))) continue;
    uint64_t idx2 = pawn_idx(&pis, &ri);
    assert(idx2 == pawn_push_idx_check(&is, s, s - 8));
    int v = merged_to_wdl[merged_table[idx2]];
    kslice_bit_set(p[v], idx);
  }
}

static void load_pawn_successors(uint8_t *table, int psq)
{
  memset(table, 8, kslice_size);
  if (psq == g_slice.sq[0] || psq == g_slice.sq[1]) {
    return;
  }

  uint8_t sq[3] = { g_slice.sq[0], g_slice.sq[1], psq };
  size_t offset = k16offset(sq);
  int s = kk_to_slice[_pext_u32(sq[0] + (sq[1] << 8),
      0b11011000110110)];
  int q = PawnFlip[0][psq];
  for (int i = 0; i < 5; i++) {
    char name[128], path[160];
    create_name(name, s, WHITE, wdl_name[i], -1);
    snprintf(path, sizeof path, "../%s/%s", pawnstr[q], name);
    k16slice_read_path_addr(k16slice_buf[7], path);
    for (uint64_t idx = 0; idx < kslice_size; idx++)
      if (kslice_bit_test(k16slice_buf[7] + offset, idx))
        table[idx] = i;
  }
}

static void calc_pawn_push(void)
{
  load_pawn_successors(merged_table, g_slice.sq[2] - 8);
  run_threaded(calc_pawn_push_worker, &work_g_dynamic);
}

static void calc_pawn_double_push_worker(struct ThreadData *thread)
{
  struct IdxState is;
  struct PawnIdxState pis8, pis16;
  int s = g_slice.sq[2];

  uint8_t *restrict ilg = k16slice_buf[1] + k16offset(g_slice.sq);
  uint8_t *restrict p[5];
  for (int i = 0; i < 5; i++)
    p[i] = k16slice_buf[2 + i] + k16offset(g_slice.sq);

  uint64_t last = thread->begin;
  idx_state_init(&is, last, g_slice.sq, &ri, false);
  pawn_idx_init(&pis8, bit(s) ^ bit(s - 8), &is, &ri);
  pawn_idx_init(&pis16, bit(s) ^ bit(s - 16), &is, &ri);
  for (uint64_t idx = last, end = thread->end; idx < end; idx++) {
    if (kslice_bit_test(ilg, idx))
      continue;
    Bitboard occ = idx_state_pawn_add(&is, idx - last, &pis8, &pis16, &ri);
    last = idx;
    if (occ & bit(s - 8)) continue;
    int v = -1;
    Bitboard occ8 = occ ^ bit(s) ^ bit(s - 8);
    if (idx_state_legal(&is, WHITE, occ8)) {
      uint64_t idx2 = pawn_idx(&pis8, &ri);
      assert(idx2 == pawn_push_idx_check(&is, s, s - 8));
      v = merged_to_wdl[merged_table[idx2]];
    }
    if (!(occ & bit(s - 16))) {
      Bitboard occ16 = occ ^ bit(s) ^ bit(s - 16);
      if (idx_state_legal(&is, WHITE, occ16)) {
        uint64_t idx2 = pawn_idx(&pis16, &ri);
        assert(idx2 == pawn_push_idx_check(&is, s, s - 16));
        v = max(v, merged_to_wdl[merged_table2[idx2]]);
      }
    }
    if (v >= 0)
      kslice_bit_set(p[v], idx);
  }
}

static void calc_pawn_double_push(void)
{
  load_pawn_successors(merged_table, g_slice.sq[2] - 8);
  load_pawn_successors(merged_table2, g_slice.sq[2] - 16);
  run_threaded(calc_pawn_double_push_worker, &work_g_dynamic);
}

static void calc_pawn_moves(void)
{
  char name[5][16];
  for (int i = 0; i < 5; i++) {
    strcat(strcpy(name[i], "pawn/"), wdl_name[i]);
    create_dir(-1, BLACK, name[i]);
  }

  uint64_t num[16];
  for (int s = 0; s < 240; s++) {
    show_progress("calculating pawn moves", s, 240, false);

    if (k16slice_test_count(s, BLACK, name[4], -1, num)) {
      for (int i = 0; i < 5; i++) {
        k16slice_test_count(s, BLACK, name[i], -1, num);
        pawn_cnt[i] += sum16(num);
      }
      continue;
    }

    k16slice_read_addr(k16slice_buf[1], s, BLACK, "illegal", -1);
    for (int i = 0; i < 5; i++)
      k16slice_clear_addr(k16slice_buf[2 + i]);

    for (int r = 0; r < 16; r++) {
      g_slice.sq[0] = KK16Square[s][r][0];
      g_slice.sq[1] = KK16Square[s][r][1];
      g_slice.stm = BLACK;

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
      pawn_cnt[i] += k16slice_count_addr(k16slice_buf[2 + i], num);
      k16slice_write_addr(k16slice_buf[2 + i], s, BLACK, name[i], -1, num);
    }
  }

  show_progress("calculating pawn moves", 240, 240, true);
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

void calc_illegal(int stm)
{
  bool partial = dir_exists(-1, stm, "illegal");

  create_dir(-1, stm, "illegal");

  g_slice.stm = stm;

  char phase[64];
  snprintf(phase, sizeof phase, "calculating %s illegal positions",
      side[stm]);

  uint64_t cnt = 0, num[16];

  for (int s = 0; s < 240; s++) {
    show_progress(phase, s, 240, false);

    if (partial && k16slice_test_count(s, stm, "illegal", -1, num)) {
      cnt += sum16(num);
      continue;
    }

    k16slice_clear_addr(k16slice_buf[0]);

    for (int r = 0; r < 16; r++) {
      g_slice.sq[0] = KK16Square[s][r][0];
      g_slice.sq[1] = KK16Square[s][r][1];

      if (is_broken(&g_slice))
        continue;

      if (    stm == BLACK
          && (pawn_attacks(BLACK, g_slice.sq[2]) & bit(g_slice.sq[0])))
      {
        k16slice_set_addr(k16slice_buf[0] + k16offset(g_slice.sq));
        continue;
      }

      int k;
      if ((k = g_piece_set[stm][KNIGHT]) >= 0)
        run_threaded(calc_illegal_knight_worker, &work_capt_dynamic[k]);
      if ((k = g_piece_set[stm][BISHOP]) >= 0)
        run_threaded(calc_illegal_bishop_worker, &work_capt_dynamic[k]);
      if ((k = g_piece_set[stm][ROOK]) >= 0)
        run_threaded(calc_illegal_rook_worker, &work_capt_dynamic[k]);
      if ((k = g_piece_set[stm][QUEEN]) >= 0)
        run_threaded(calc_illegal_queen_worker, &work_capt_dynamic[k]);
    }

    cnt += k16slice_count_addr(k16slice_buf[0], num);
    k16slice_write_addr(k16slice_buf[0], s, stm, "illegal", -1, num);
  }

  snprintf(phase, sizeof phase, "%lu %s illegal positions", cnt, side[stm]);
  show_progress(phase, 240, 240, true);
}

void calc_illegal_both(void)
{
  bool partial[2] = {
    dir_exists(-1, WHITE, "illegal"),
    dir_exists(-1, BLACK, "illegal")
  };

  create_dir(-1, WHITE, "illegal");
  create_dir(-1, BLACK, "illegal");

  uint64_t cnt[2] = { 0 }, num[16];

  for (int s = 0; s < 240; s++) {
    show_progress("calculating illegal positions", s, 240, false);

    if (   partial[WHITE] && partial[BLACK]
        && k16slice_test_count(s, WHITE, "illegal", -1, num))
    {
      cnt[WHITE] += sum16(num);
      if (k16slice_test_count(s, BLACK, "illegal", -1, num)) {
        cnt[BLACK] += sum16(num);
        continue;
      }
    }

    k16slice_clear_addr(k16slice_buf[WHITE]);
    k16slice_clear_addr(k16slice_buf[BLACK]);

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
        kslice_set_addr(k16slice_buf[BLACK] + k16offset(g_slice.sq));
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
      cnt[stm] += k16slice_count_addr(k16slice_buf[stm], num);
      k16slice_write_addr(k16slice_buf[stm], s, stm, "illegal", -1, num);
    }
  }

  char phase[128];
  snprintf(phase, sizeof phase, "%lu white, %lu black illegal positions",
      cnt[WHITE], cnt[BLACK]);
  show_progress(phase, 240, 240, true);
}

static mtx_t report_mutex;
static int num_fails;

static void report_fail(int s, uint64_t idx, Bitboard occ,
    const struct IdxState *is, int n)
{
  mtx_lock(&report_mutex);
  char fenstr[64];
  Position pos = g_pos;
  pos.sq[0] = g_slice.sq[0];
  pos.sq[1] = g_slice.sq[1];
  pos.sq[2] = g_slice.sq[2];
  pos.stm = g_slice.stm;
  pos.occ = occ;
  idx_state_to_sq(is, pos.sq, &ri);
  pos_to_fen(&pos, fenstr, false);
  int wdl = probe_wdl(&pos, -2, 2);
  printf("\nn = %d, slice = %d, idx = %lu, wdl = %d\n%s\n", n, s, idx, wdl,
      fenstr);
  num_fails++;
  if (num_fails == 10)
    exit(EXIT_FAILURE);
  mtx_unlock(&report_mutex);
}

static bool check_dtz_W101(struct IdxState *is, Bitboard occ)
{
  Position pos = g_pos;
  pos.sq[0] = g_slice.sq[0];
  pos.sq[1] = g_slice.sq[1];
  pos.sq[2] = g_slice.sq[2];
  pos.stm = g_slice.stm;
  pos.occ = occ;

  // First check that the position has dtz == 101.
  if (probe_dtz(&pos) != DRAW_RULE + 1)
    return false;

  // Now check that the best quiet move reaches dtz == -100.
  bool pos_ok = false;
  for (int i = 0; i < pos.num; i++) {
    if ((pos.pt[i] >> 3) != pos.stm) continue;
    int from = pos.sq[i];
    Bitboard b = piece_moves(pos.pt[i], from, pos.occ);
    while (b) {
      int to = pop_lsb(&b);
      if (do_move(&pos, from, to, i)) {
        bool capture_is_best;
        int dtz, wdl = probe_wdl_helper(&pos, &capture_is_best);
        if (wdl == -2 && !capture_is_best)
          dtz = probe_dtz_helper(&pos, wdl);
        undo_move(&pos, from, to, i);
        if (wdl == -2) {
          if (capture_is_best || dtz != -DRAW_RULE)
            return false;
          pos_ok = true;
        }
      }
    }
  }
  return pos_ok;
}

static bool check_dtz_L101(struct IdxState *is, Bitboard occ)
{
  Position pos = g_pos;
  pos.sq[0] = g_slice.sq[0];
  pos.sq[1] = g_slice.sq[1];
  pos.sq[2] = g_slice.sq[2];
  pos.stm = g_slice.stm;
  pos.occ = occ;

  // First check that the position has dtz == -101.
  if (probe_dtz(&pos) != -DRAW_RULE - 1)
    return false;

  // Now check that the best quiet move reaches dtz == 100.
  bool pos_ok = false;
  for (int i = 0; i < pos.num; i++) {
    if ((pos.pt[i] >> 3) != pos.stm) continue;
    int from = pos.sq[i];
    Bitboard b = piece_moves(pos.pt[i], from, pos.occ);
    while (b) {
      int to = pop_lsb(&b);
      if (do_move(&pos, from, to, i)) {
        bool capture_is_best;
        int dtz, wdl = probe_wdl_helper(&pos, &capture_is_best);
        if (wdl != 2)
          return false;
        if (!capture_is_best)
          dtz = probe_dtz_helper(&pos, wdl);
        undo_move(&pos, from, to, i);
        if (dtz == DRAW_RULE)
          pos_ok = true;
      }
    }
  }
  return pos_ok;
}

enum { CZ_REGULAR, CZ_CWIN };

static int work_slice;

INLINE void check_zero_worker_tmpl(struct ThreadData *thread, const int T)
{
  struct IdxState is;
  int stm = g_slice.stm;
  int s = work_slice;

  uint64_t *restrict p = (uint64_t *)kslice_get_address_scratch(g_slice.sq);
  uint8_t *restrict const q = kslice_get_address(g_slice.sq);

  p += thread->begin >> 6;
  uint64_t last = thread->begin;
  idx_state_init(&is, last, g_slice.sq, &ri, false);
  for (uint64_t idx = last, end = thread->end; idx < end; idx += 64, p++) {
    uint64_t w = *p;
    while (w) {
      uint64_t cur = idx + pop_lsb(&w);
      Bitboard occ = idx_state_add(&is, cur - last, &ri);
      last = cur;

      if (   !test_moves(stm, QUEEN , q, occ, &is)
          && !test_moves(stm, ROOK  , q, occ, &is)
          && !test_moves(stm, BISHOP, q, occ, &is)
          && !test_moves(stm, KNIGHT, q, occ, &is)
          && !test_king_moves(stm, occ, &is))
        continue;

      switch (T) {
      case CZ_REGULAR:
        report_fail(s, cur, occ, &is, 0);
        break;
      case CZ_CWIN:
        if (!check_dtz_W101(&is, occ))
          report_fail(s, cur, occ, &is, 0);
        break;
      }
    }
  }
}

static void check_zero_regular_worker(struct ThreadData *thread)
{
  check_zero_worker_tmpl(thread, CZ_REGULAR);
}

static void check_zero_cwin_worker(struct ThreadData *thread)
{
  check_zero_worker_tmpl(thread, CZ_CWIN);
}

static void check_zero(int stm, int s, const int T)
{
  work_slice = s;
  g_slice.stm = stm;
  for (int r = 0; r < 16; r++) {
    g_slice.sq[0] = KK16Square[s][r][0];
    g_slice.sq[1] = KK16Square[s][r][1];

    if (is_broken(&g_slice))
      continue;

    switch (T) {
      case CZ_REGULAR:
        run_threaded(check_zero_regular_worker, &work_g_dynamic);
        break;
      case CZ_CWIN:
        run_threaded(check_zero_cwin_worker, &work_g_dynamic);
        break;
    }
  }
}

enum { CO_REGULAR, CO_DRAW, CO_BLOSS, CO_LOSS };

INLINE void check_one_worker_tmpl(struct ThreadData *thread, const int T)
{
  struct IdxState is;
  int stm = g_slice.stm;
  int s = work_slice;

  uint64_t *restrict p = (uint64_t *)kslice_get_address_scratch(g_slice.sq);
  uint8_t *restrict const q = kslice_get_address(g_slice.sq);

  p += thread->begin >> 6;
  uint64_t last = thread->begin;
  idx_state_init(&is, last, g_slice.sq, &ri, false);
  for (uint64_t idx = last, end = thread->end; idx < end; idx += 64, p++) {
    uint64_t w = *p;
    while (w) {
      uint64_t cur = idx + pop_lsb(&w);
      Bitboard occ = idx_state_add(&is, cur - last, &ri);
      last = cur;

      if (   test_moves(stm, QUEEN , q, occ, &is)
          || test_moves(stm, ROOK  , q, occ, &is)
          || test_moves(stm, BISHOP, q, occ, &is)
          || test_moves(stm, KNIGHT, q, occ, &is)
          || test_king_moves(stm, occ, &is))
        continue;

      switch (T) {
      case CO_REGULAR:
        report_fail(s, cur, occ, &is, 1);
        break;
      case CO_DRAW:
        if (    idx_state_has_legal_moves(&is, stm, occ)
            || !idx_state_legal(&is, stm ^ 1, occ))
          report_fail(s, cur, occ, &is, 1);
        break;
      case CO_BLOSS:
        if (!check_dtz_L101(&is, occ))
          report_fail(s, cur, occ, &is, 1);
        break;
      case CO_LOSS:
        if (idx_state_legal(&is, stm ^ 1, occ))
          report_fail(s, cur, occ, &is, 1);
        break;
      default:
        abort();
      }
    }
  }
}

static void check_one_regular_worker(struct ThreadData *thread)
{
  check_one_worker_tmpl(thread, CO_REGULAR);
}

static void check_one_draw_worker(struct ThreadData *thread)
{
  check_one_worker_tmpl(thread, CO_DRAW);
}

static void check_one_bloss_worker(struct ThreadData *thread)
{
  check_one_worker_tmpl(thread, CO_BLOSS);
}

static void check_one_loss_worker(struct ThreadData *thread)
{
  check_one_worker_tmpl(thread, CO_LOSS);
}

static void check_one(int stm, int s, const int T)
{
  work_slice = s;
  g_slice.stm = stm;
  for (int r = 0; r < 16; r++) {
    g_slice.sq[0] = KK16Square[s][r][0];
    g_slice.sq[1] = KK16Square[s][r][1];

    if (is_broken(&g_slice))
      continue;

    switch (T) {
    case CO_REGULAR:
      run_threaded(check_one_regular_worker, &work_g_dynamic);
      break;
    case CO_DRAW:
      run_threaded(check_one_draw_worker, &work_g_dynamic);
      break;
    case CO_BLOSS:
      run_threaded(check_one_bloss_worker, &work_g_dynamic);
      break;
    case CO_LOSS:
      run_threaded(check_one_loss_worker, &work_g_dynamic);
      break;
    default:
      abort();
    }
  }
}

static void check_win_cwin(int stm)
{
  struct K16SliceIterator iter;

  char phase[64];
  snprintf(phase, sizeof phase, "check %s win, cwin-", side[stm]);
  int num_done = 0;

  k16slice_iter_init(&iter, stm);
  int s, s1;
  while (k16slice_iter_next(&iter, &s)) {
    show_progress(phase, num_done++, 240, false);

    while (k16slice_iter_in(&iter, &s1)) {
      k16slice_read(s1, s1, stm ^ 1, "loss", -1);
    }

    k16slice_read(-1, s, stm, "wins", -1);
    k16slice_read_andnot(-1, s, stm, "capt/win", -1);
    if (stm == BLACK)
      k16slice_read_andnot(-1, s, stm, "pawn/win", -1);
    check_one(stm, s, CO_REGULAR);

    k16slice_read(-1, s, stm, "cwin", -1);
    k16slice_read_andnot(-1, s, stm, "capt/cwin", -1);
    check_zero(stm, s, CZ_CWIN);

    while (k16slice_iter_out(&iter, &s));
  }
  show_progress(phase, 240, 240, true);
}

static void check_cwin_draw(int stm)
{
  struct K16SliceIterator iter;

  char phase[64];
  snprintf(phase, sizeof phase, "check %s cwin+, draw-", side[stm]);
  int num_done = 0;

  k16slice_iter_init(&iter, stm);
  int s, s1;
  while (k16slice_iter_next(&iter, &s)) {
    show_progress(phase, num_done++, 240, false);

    while (k16slice_iter_in(&iter, &s1)) {
      k16slice_read(s1, s1, stm ^ 1, "loss", -1);
      k16slice_read_or(s1, s1, stm ^ 1, "bloss", -1);
    }

    // We already know that no cwin position reaches "loss" (with the
    // exception of dtz==101 -> dtz==100). So it does not hurt to test against
    // "loss" and "bloss" here to ensure that one cwin position reaches
    // "bloss". As a bonus, no special check is needed for the 101->100 case.
    k16slice_read(-1, s, stm, "cwin", -1);
    k16slice_read_andnot(-1, s, stm, "capt/cwin", -1);
    if (stm == BLACK)
      k16slice_read_andnot(-1, s, stm, "pawn/cwin", -1);
    check_one(stm, s, CO_REGULAR);

    k16slice_read(-1, s, stm, "draw", -1);
    k16slice_read_andnot(-1, s, stm, "capt/draw", -1);
    check_zero(stm, s, CZ_REGULAR);

    while (k16slice_iter_out(&iter, &s));
  }
  show_progress(phase, 240, 240, true);
}

static void check_draw_bloss(int stm)
{
  struct K16SliceIterator iter;

  char phase[64];
  snprintf(phase, sizeof phase, "check %s draw+, bloss-", side[stm]);
  int num_done = 0;

  k16slice_iter_init(&iter, stm);
  int s, s1;
  while (k16slice_iter_next(&iter, &s)) {
    show_progress(phase, num_done++, 240, false);

    while (k16slice_iter_in(&iter, &s1)) {
      k16slice_read(s1, s1, stm ^ 1, "loss", -1);
      k16slice_read_or(s1, s1, stm ^ 1, "bloss", -1);
      k16slice_read_or(s1, s1, stm ^ 1, "draw", -1);
    }

    // Again, we already checked that drawn positions do not hit "loss" and
    // "bloss". So we can combine draw+ with bloss-.
    k16slice_read(-1, s, stm, "draw", -1);
    k16slice_read_andnot(-1, s, stm, "capt/draw", -1);
    if (stm == BLACK)
      k16slice_read_andnot(-1, s, stm, "pawn/draw", -1);
    check_one(stm, s, CO_DRAW);

    k16slice_read(-1, s, stm, "bloss", -1);
    k16slice_read_andnot(-1, s, stm, "capt/bloss", -1);
    check_zero(stm, s, CZ_REGULAR);

    while (k16slice_iter_out(&iter, &s));
  }
  show_progress(phase, 240, 240, true);
}

static void check_bloss_loss(int stm)
{
  struct K16SliceIterator iter;

  char phase[64];
  snprintf(phase, sizeof phase, "check %s bloss+, loss-", side[stm]);
  int num_done = 0;

  k16slice_iter_init(&iter, stm);
  int s, s1;
  while (k16slice_iter_next(&iter, &s)) {
    show_progress(phase, num_done++, 240, false);

    // Load "loss" | "bloss" | "draw" | "cwin".
    while (k16slice_iter_in(&iter, &s1)) {
      k16slice_read(s1, s1, stm ^ 1, "win", -1);
      k16slice_read_or(s1, s1, stm ^ 1, "illegal", -1);
      k16slice_not(s1);
    }

    k16slice_read(-1, s, stm, "bloss", -1);
    k16slice_read_andnot(-1, s, stm, "capt/bloss", -1);
    if (stm == BLACK)
      k16slice_read_andnot(-1, s, stm, "pawn/bloss", -1);
    check_one(stm, s, CO_BLOSS);

    k16slice_read(-1, s, stm, "loss", -1);
    check_zero(stm, s, CZ_REGULAR);

    while (k16slice_iter_out(&iter, &s));
  }
  show_progress(phase, 240, 240, true);
}

static void check_loss(int stm)
{
  struct K16SliceIterator iter;

  char phase[64];
  snprintf(phase, sizeof phase, "check %s loss+", side[stm]);
  int num_done = 0;

  k16slice_iter_init(&iter, stm);
  int s, s1;
  while (k16slice_iter_next(&iter, &s)) {
    show_progress(phase, num_done++, 240, false);

    while (k16slice_iter_in(&iter, &s1))
      k16slice_read(s1, s1, stm ^ 1, "win", -1);

    k16slice_read(-1, s, stm, "loss", -1);
    k16slice_read_andnot(-1, s, stm, "capt/loss", -1);
    if (stm == BLACK)
      k16slice_read_andnot(-1, s, stm, "pawn/loss", -1);
    check_one(stm, s, CO_LOSS);

    while (k16slice_iter_out(&iter, &s));
  }
  show_progress(phase, 240, 240, true);
}

static uint8_t *tb_table, *dtz_table;
static struct TbTable2 *table_info;
static struct RankInfo tb_ri[2];
static size_t work_k16_offset;
static bool work_pk_layout;
static bool flip[2], btm_side[2];

INLINE int current_k16slice(void)
{
  return kk_to_slice[_pext_u32(g_slice.sq[0] + (g_slice.sq[1] << 8),
      0b11011000110110)];
}

static void init_perm_rank_pawn(struct RankInfo *rank_ri,
    const struct TbTable2 *table, const struct RankInfo *dst_ri)
{
  *rank_ri = *table->ri;
  for (int k = 0; k < table->ri->numsets; k++) {
    if (table->first[k] == 0 || table->first[k] == 1) {
      int side = table->first[k];
      if (work_pk_layout)
        side = side == 0 ? g_slice.stm : g_slice.stm ^ 1;
      rank_ri->perm[k] = dst_ri->numsets + 1 + side;
      continue;
    }

    int j = 0;
    for (; j < dst_ri->numsets; j++)
      if (table->first[k] == dst_ri->first[j])
        break;
    assert(j < dst_ri->numsets);
    rank_ri->perm[k] = j + 1;
  }
}

static uint64_t depermute_pawn_idx(struct IdxState *is,
    const struct RankInfo *rank_ri)
{
  alignas(64) Bitboard bb[8];
  for (int k = 0; k <= ri.numsets; k++)
    bb[k] = is->bb[k];

  bb[0] = bit(g_slice.sq[2]);
  if (work_pk_layout)
    bb[0] |= bit(g_slice.sq[g_slice.stm]);
  bb[ri.numsets + 1] = bit(g_slice.sq[0]);
  bb[ri.numsets + 2] = bit(g_slice.sq[1]);

  return perm_rank_bb(bb, rank_ri);
}

static void table_squares_from_idx_state(uint8_t sq[MAX_PIECES],
    struct IdxState *is)
{
  int stm = g_slice.stm;
  uint8_t mirror = flip[stm] ? 0x38 : 0x00;

  sq[0] = g_slice.sq[flip[stm] ? BLACK : WHITE] ^ mirror;
  sq[1] = g_slice.sq[flip[stm] ? WHITE : BLACK] ^ mirror;
  sq[2] = g_slice.sq[2] ^ mirror;

  for (int k = 0; k < ri.numsets; k++) {
    int i = tb_ri[stm].first[k];
    Bitboard b = is->bb[k + 1];
    while (b)
      sq[i++] = pop_lsb(&b) ^ mirror;
  }

  if (sq[2] & 0x04)
    for (int i = 0; i < g_pos.num; i++)
      sq[i] ^= 0x07;
}

static uint64_t depermute_pawn_trivial_idx(struct IdxState *is,
    const struct TbTable2 *table)
{
  uint8_t sq[MAX_PIECES];
  int stm = g_slice.stm;
  table_squares_from_idx_state(sq, is);

  Bitboard occ = bit(sq[2]);
  if (work_pk_layout) {
    if (btm_side[stm])
      Swap(sq[0], sq[1]);
    occ |= bit(sq[0]);
  }

  return rank_trivial_from(sq, 0, occ, table->first, table->ri);
}

static int pawn_pk_tsq(int q)
{
  int stm = g_slice.stm;
  uint8_t mirror = flip[stm] ? 0x38 : 0x00;
  int p0 = g_slice.sq[flip[stm] ? BLACK : WHITE] ^ mirror;
  int p1 = g_slice.sq[flip[stm] ? WHITE : BLACK] ^ mirror;
  int psq = g_slice.sq[2] ^ mirror;

  if (psq & 0x04) {
    p0 ^= 0x07;
    p1 ^= 0x07;
    psq ^= 0x07;
  }
  if (btm_side[stm])
    Swap(p0, p1);

  assert(PawnFlip[0][psq] == q);
  return q * 63 + p0 - (p0 > psq);
}

INLINE void store_wdl(uint64_t idx, uint8_t *restrict tmp, int v)
{
  uint64_t *p[5];
  for (int i = 0; i < 5; i++)
    p[i] = (uint64_t *)(k16slice_buf[i] + work_k16_offset);

  tmp[idx & 63] = v;
  if ((idx & 63) == 63) {
#ifdef __AVX512BW__
    __m512i x = _mm512_load_si512((__m512i *)tmp);
    p[0][idx >> 6] = _mm512_cmpeq_epu8_mask(x, _mm512_set1_epi8(0));
    p[1][idx >> 6] = _mm512_cmpeq_epu8_mask(x, _mm512_set1_epi8(1));
    p[2][idx >> 6] = _mm512_cmpeq_epu8_mask(x, _mm512_set1_epi8(2));
    p[3][idx >> 6] = _mm512_cmpeq_epu8_mask(x, _mm512_set1_epi8(3));
    p[4][idx >> 6] = _mm512_cmpeq_epu8_mask(x, _mm512_set1_epi8(4));
#else
    uint64_t w[5] = { 0 };
    for (int k = 0; k < 64; k++)
      if (tmp[k] < 5)
        w[tmp[k]] |= 1ULL << k;
    p[0][idx >> 6] = w[0];
    p[1][idx >> 6] = w[1];
    p[2][idx >> 6] = w[2];
    p[3][idx >> 6] = w[3];
    p[4][idx >> 6] = w[4];
#endif
  }
}

#define NUM 16
static void depermute_wdl_pawn_worker(struct ThreadData *thread)
{
  alignas(64) uint8_t tmp[64];
  const uint8_t *restrict src = tb_table;
  const struct RankInfo *map_ri = &tb_ri[g_slice.stm];
  const struct RankInfo *dst_ri = &ri;
  struct RankInfo rank_ri;
  struct IdxState is;

  uint64_t idx_dec_buf[NUM];

  if (work_pk_layout)
    init_perm_rank_pawn(&rank_ri, table_info, map_ri);
  idx_state_init(&is, thread->begin, g_slice.sq, dst_ri, false);

  uint64_t idx = thread->begin, end = thread->end;
  int fill = 0;
  int head = 0;

  for (; fill < NUM && idx < end;
      fill++, idx++, idx_state_inc(&is, dst_ri)) {
    uint64_t idx_dec = depermute_pawn_trivial_idx(&is, table_info);
    __builtin_prefetch(src + idx_dec, 0, 3);
    idx_dec_buf[fill] = idx_dec;
  }

  for (; idx < end; idx++, idx_state_inc(&is, dst_ri)) {
    uint64_t idx_dec = depermute_pawn_trivial_idx(&is, table_info);
    __builtin_prefetch(src + idx_dec, 0, 3);
    store_wdl(idx - NUM, tmp, src[idx_dec_buf[head]]);
    idx_dec_buf[head] = idx_dec;
    head = (head + 1) & (NUM - 1);
  }

  for (uint64_t out = idx - fill; fill-- > 0; out++) {
    store_wdl(out, tmp, src[idx_dec_buf[head]]);
    head = (head + 1) & (NUM - 1);
  }
  for (; idx & 0x3f; idx++)
    store_wdl(idx, tmp, 0);
}

static uint64_t wdl_cnt[2][5];

static void k16slice_read_path_addr(void *p, const char *str)
{
  FILE *F = fopen(str, "rb");
  if (!F && errno != ENOENT) {
    fprintf(stderr, "Could not open %s for reading.\n", str);
    exit(EXIT_FAILURE);
  }
  if (!F) {
    k16slice_clear_addr(p);
    return;
  }

  struct stat st;
  fstat(fileno(F), &st);
  if (st.st_size > 0) {
    if (fseek(F, 16 * 8, SEEK_SET) != 0) {
      fprintf(stderr, "fseek() failed.\n");
      exit(EXIT_FAILURE);
    }
    read_data(F, p, k16slice_cache_lines << 6);
  } else
    k16slice_clear_addr(p);
  fclose(F);
}

static void k16slice_read_addr(void *p, int slice, int stm, const char *name,
    int n)
{
  char str[64];
  create_name(str, slice, stm, name, n);
  k16slice_read_path_addr(p, str);
}

static void k16slice_read_andnot_addr(void *p, int slice, int stm,
    const char *name, int n)
{
  char str[64];
  create_name(str, slice, stm, name, n);

  FILE *F = fopen(str, "rb");
  if (!F && errno != ENOENT) {
    fprintf(stderr, "Could not open %s for reading.\n", str);
    exit(EXIT_FAILURE);
  }
  if (!F)
    return;

  struct stat st;
  fstat(fileno(F), &st);
  if (st.st_size > 0) {
    uint8_t *tmp = alloc_k16slice();
    if (fseek(F, 16 * 8, SEEK_SET) != 0) {
      fprintf(stderr, "fseek() failed.\n");
      exit(EXIT_FAILURE);
    }
    read_data(F, tmp, k16slice_cache_lines << 6);
    uint64_t *restrict a = p;
    uint64_t *restrict b = (uint64_t *)tmp;
    for (size_t i = 0; i < k16slice_cache_lines << 3; i++)
      a[i] &= ~b[i];
    free(tmp);
  }
  fclose(F);
}

static void k16slice_abc_addr(void *p, void *q, void *r)
{
  uint64_t *restrict a0 = p;
  uint64_t *restrict b0 = q;
  uint64_t *restrict c0 = r;

  for (size_t i = 0; i < k16slice_cache_lines << 3; i++) {
    uint64_t a = a0[i];
    uint64_t b = b0[i];
    uint64_t c = c0[i];
    a0[i] = (a | c) & ~b;
    b0[i] = a | b | c;
  }
}

static void k16slice_andnot_addr(void *p, void *q)
{
  uint64_t *restrict a = p;
  uint64_t *restrict b = q;

  for (size_t i = 0; i < k16slice_cache_lines << 3; i++)
    a[i] &= ~b[i];
}

static void finish_wdl_kslice(int stm, int s)
{
  k16slice_read_addr(k16slice_buf[5], s, stm, "illegal", -1);
  k16slice_read_addr(k16slice_buf[6], s, stm, "capt/win", -1);
  k16slice_abc_addr(k16slice_buf[4], k16slice_buf[5], k16slice_buf[6]);
  k16slice_read_addr(k16slice_buf[6], s, stm, "capt/cwin", -1);
  k16slice_abc_addr(k16slice_buf[3], k16slice_buf[5], k16slice_buf[6]);
  k16slice_read_addr(k16slice_buf[6], s, stm, "capt/draw", -1);
  k16slice_abc_addr(k16slice_buf[2], k16slice_buf[5], k16slice_buf[6]);
  k16slice_read_addr(k16slice_buf[6], s, stm, "capt/bloss", -1);
  k16slice_abc_addr(k16slice_buf[1], k16slice_buf[5], k16slice_buf[6]);
  k16slice_andnot_addr(k16slice_buf[0], k16slice_buf[5]);

  uint64_t num[16];
  for (int i = 0; i < 5; i++) {
    wdl_cnt[stm][i] += k16slice_count_addr(k16slice_buf[i], num);
    k16slice_write_addr(k16slice_buf[i], s, stm, wdl_name[i], -1, num);
  }
}

static void decompress_kslice_wdl_p(struct Tbase *tb, int stm, int q)
{
  int psq = InvPawnFlip[0][q];
  int t = !symmetric ? 2 * q + btm_side[stm] : q;
  if (!tb->table[t])
    tb->table[t] = init_new_table(tb, g_pos.num, WDL, t, q);
  struct TbTable2 *table = tb->table[t];

  int const_value = -1;
  if (!table->precomp) {
    struct TbTableConst *tbl = (struct TbTableConst *)table;
    const_value = tbl->const_value;
  } else {
    decompress_table(table, tb_table, 63 * 62 * kslice_size);
    table_info = table;
  }

  g_slice.sq[2] = psq;
  g_slice.stm = stm;
  work_pk_layout = false;
  for (int s = 0; s < 240; s++) {
    for (int i = 0; i < 5; i++)
      k16slice_clear_addr(k16slice_buf[i]);

    for (int r = 0; r < 16; r++) {
      g_slice.sq[0] = KK16Square[s][r][0];
      g_slice.sq[1] = KK16Square[s][r][1];

      if (is_broken(&g_slice))
        continue;
      work_k16_offset = k16offset(g_slice.sq);
      if (const_value >= 0) {
        kslice_set_addr(k16slice_buf[const_value] + work_k16_offset);
        for (int i = 0; i < 5; i++)
          if (i != const_value)
            memset(k16slice_buf[i] + work_k16_offset, 0, kslice_alloc_size);
      } else
        run_threaded(depermute_wdl_pawn_worker, &work_g_static);
    }
    finish_wdl_kslice(stm, s);
  }
}

static void fill_wdl_pk_kslice(struct Tbase *tb, int stm, int q,
    int ksq)
{
  int psq = InvPawnFlip[0][q];
  int tsq = q * 63 + ksq - (ksq > psq);
  int t = !symmetric ? 2 * tsq + btm_side[stm] : tsq;
  if (!tb->table[t])
    tb->table[t] = init_new_table(tb, g_pos.num, WDL, t, tsq);
  struct TbTable2 *table = tb->table[t];

  int const_value = -1;
  if (!table->precomp) {
    struct TbTableConst *tbl = (struct TbTableConst *)table;
    const_value = tbl->const_value;
  } else {
    decompress_table(table, tb_table, 62 * kslice_size);
    table_info = table;
  }

  g_slice.sq[2] = psq;
  g_slice.sq[stm] = ksq;
  g_slice.stm = stm;
  work_pk_layout = true;
  work_k16_offset = k16offset(g_slice.sq);
  if (const_value >= 0) {
    kslice_set_addr(k16slice_buf[const_value] + work_k16_offset);
    for (int i = 0; i < 5; i++)
      if (i != const_value)
        memset(k16slice_buf[i] + work_k16_offset, 0, kslice_alloc_size);
  } else
    run_threaded(depermute_wdl_pawn_worker, &work_g_static);
}

static void depermute_dtz_pawn_worker(struct ThreadData *thread)
{
  const uint8_t *restrict src = tb_table;
  uint8_t *restrict dst = dtz_table;
  const struct RankInfo *map_ri = &tb_ri[g_slice.stm];
  const struct RankInfo *dst_ri = &ri;
  struct RankInfo rank_ri;
  struct IdxState is;

  uint64_t idx_dec_buf[NUM];

  if (work_pk_layout)
    init_perm_rank_pawn(&rank_ri, table_info, map_ri);
  idx_state_init(&is, thread->begin, g_slice.sq, dst_ri, false);

  uint64_t idx = thread->begin, end = thread->end;
  int fill = 0;
  int head = 0;

  for (; fill < NUM && idx < end;
      fill++, idx++, idx_state_inc(&is, dst_ri)) {
    uint64_t idx_dec = depermute_pawn_trivial_idx(&is, table_info);
    __builtin_prefetch(src + idx_dec, 0, 3);
    idx_dec_buf[fill] = idx_dec;
  }

  for (; idx < end; idx++, idx_state_inc(&is, dst_ri)) {
    uint64_t idx_dec = depermute_pawn_trivial_idx(&is, table_info);
    __builtin_prefetch(src + idx_dec, 0, 3);
    dst[idx - NUM] = src[idx_dec_buf[head]];
    idx_dec_buf[head] = idx_dec;
    head = (head + 1) & (NUM - 1);
  }

  for (uint64_t out = idx - fill; fill-- > 0; out++) {
    dst[out] = src[idx_dec_buf[head]];
    head = (head + 1) & (NUM - 1);
  }
}

static bool dtz_stored[2][5];
static uint16_t freq_map[MAX_VAL];
static uint64_t (*dtz_freq)[2][2][MAX_VAL + DRAW_RULE];
static int win_loss;
static struct DtzTable2 unmapped_dtz_table;

static void update_stats_worker(struct ThreadData *thread)
{
  const uint16_t *restrict map = freq_map;
  uint64_t *restrict freq = dtz_freq[thread->thread_id][g_slice.stm][win_loss];
  const uint64_t *restrict p =
      (uint64_t *)(k16slice_buf[0] + work_k16_offset);
  uint8_t *restrict q = dtz_table;
  p += thread->begin >> 6;
  for (uint64_t idx = thread->begin, end = thread->end; idx < end; idx += 64) {
    uint64_t w = *p++;
    while (w) {
      uint64_t cur = idx + pop_lsb(&w);
      freq[map[q[cur]]]++;
    }
  }
}

static void update_dtz_stats(struct DtzTable2 *table, int stm, int s)
{
  static constexpr int WdlToMap[5] = { 1, 3, 0, 2, 0 };
  static constexpr int offset[4] = { 0, 0, 100, 100 };

  int bound[4] = { 0 };
  switch (table->mapped) {
  case 1:
    const uint8_t *data = table->map_dtz;
    for (int i = 0; i < 4; i++) {
      bound[i] = data[0];
      data += 1 + data[0];
    }
    break;
  case 2:
    const uint16_t *data16 = table->map_dtz16;
    for (int i = 0; i < 4; i++) {
      bound[i] = data16[0];
      data16 += 1 + data16[0];
    }
    break;
  }

  for (int i = 0; i < 5; i++) {
    if (!dtz_stored[stm][i])
      continue;
    k16slice_read_addr(k16slice_buf[0], s, stm, wdl_name[i], -1);
    if (i == 3 || i == 4) {
      char capt_name[32];
      strcat(strcpy(capt_name, "capt/"), wdl_name[i]);
      k16slice_read_andnot_addr(k16slice_buf[0], s, stm, capt_name, -1);
      if (stm == BLACK) {
        char pawn_name[32];
        strcat(strcpy(pawn_name, "pawn/"), wdl_name[i]);
        k16slice_read_andnot_addr(k16slice_buf[0], s, stm, pawn_name, -1);
      }
    }
    int m = WdlToMap[i];
    win_loss = i > 2 ? 0 : 1;
    memset(freq_map, 0, sizeof freq_map);
    switch (table->mapped) {
      case 0:
        for (int j = 0; j < MAX_VAL; j++)
          freq_map[j] = offset[m] + j;
        break;
      case 1:
        for (int j = 0; j < bound[m]; j++)
          freq_map[j] = offset[m] + table->map_dtz[table->map_dtz_idx[m] + j];
        break;
      case 2:
        for (int j = 0; j < bound[m]; j++)
          freq_map[j] =
              offset[m] + table->map_dtz16[table->map_dtz_idx[m] + j];
        break;
    }
    run_threaded(update_stats_worker, &work_g_static);
  }
}

static void decompress_kslice_dtz_p(struct Tbase *tb, int stm, int q)
{
  int psq = InvPawnFlip[0][q];
  int t = (tb->dist_format & TWO_SIDED) ? 2 * q + btm_side[stm] : q;
  if (!tb->table[t])
    tb->table[t] = init_new_table(tb, g_pos.num, DTZ, t, q);
  struct DtzTable2 *table = tb->table[t];
  if (table == (struct DtzTable2 *)1)
    return;

  int const_value = -1;
  if (!table->precomp) {
    struct TbTableConst *tbl = (struct TbTableConst *)table;
    const_value = tbl->const_value;
    table = &unmapped_dtz_table;
  } else {
    decompress_table((struct TbTable2 *)table, tb_table,
        63 * 62 * kslice_size);
    table_info = (struct TbTable2 *)table;
  }

  g_slice.sq[2] = psq;
  g_slice.stm = stm;
  work_pk_layout = false;
  for (int k1 = 0; k1 < 64; k1++) {
    if (k1 == psq) continue;
    g_slice.sq[0] = k1;
    for (int k2 = 0; k2 < 64; k2++) {
      if (k2 == k1 || k2 == psq) continue;
      g_slice.sq[1] = k2;

      if (is_broken(&g_slice))
        continue;
      int s = current_k16slice();
      work_k16_offset = k16offset(g_slice.sq);
      if (const_value >= 0)
        memset(dtz_table, const_value, kslice_size);
      else
        run_threaded(depermute_dtz_pawn_worker, &work_g_static);

      update_dtz_stats(table, stm, s);
    }
  }
}

static void decompress_kslice_dtz_pk(struct Tbase *tb, int stm, int q,
    int ksq)
{
  int psq = InvPawnFlip[0][q];
  int tsq = q * 63 + ksq - (ksq > psq);
  int t = (tb->dist_format & TWO_SIDED) ? 2 * tsq + btm_side[stm] : tsq;
  if (!tb->table[t])
    tb->table[t] = init_new_table(tb, g_pos.num, DTZ, t, tsq);
  struct DtzTable2 *table = tb->table[t];
  if (table == (struct DtzTable2 *)1)
    return;

  int const_value = -1;
  if (!table->precomp) {
    struct TbTableConst *tbl = (struct TbTableConst *)table;
    const_value = tbl->const_value;
    table = &unmapped_dtz_table;
  } else {
    decompress_table((struct TbTable2 *)table, tb_table, 62 * kslice_size);
    table_info = (struct TbTable2 *)table;
  }

  g_slice.sq[2] = psq;
  g_slice.sq[stm] = ksq;
  g_slice.stm = stm;
  work_pk_layout = true;
  for (int k2 = 0; k2 < 64; k2++) {
    if (k2 == ksq || k2 == psq) continue;
    g_slice.sq[stm ^ 1] = k2;

    if (is_broken(&g_slice))
      continue;
    int s = current_k16slice();
    work_k16_offset = k16offset(g_slice.sq);
    if (const_value >= 0)
      memset(dtz_table, const_value, kslice_size);
    else
      run_threaded(depermute_dtz_pawn_worker, &work_g_static);

    update_dtz_stats(table, stm, s);
  }
}

void verify(void)
{
  // Initialize and mmap() tablebase file.
  struct TbEntry entry = {
    .has_pawns = true,
    .symmetric = symmetric,
    .num = g_pos.num
  };

  char name[64];
  if (g_tablename[0] != '/')
    snprintf(name, sizeof name, "../%s", g_tablename);
  else
    strcpy(name, g_tablename);

  struct Tbase *wdl_tb = init_tbase(&entry, name, WDL, false);
  if (!wdl_tb) {
    fprintf(stderr, "Could not open %s.rtbw.\n", g_tablename);
    exit(EXIT_FAILURE);
  }

  if (wdl_tb->layout != LT_PAWN_P && wdl_tb->layout != LT_PAWN_PK) {
    fprintf(stderr, "Layout type %d is currently not supported.\n",
        wdl_tb->layout);
    exit(EXIT_FAILURE);
  }

  struct Tbase *dtz_tb = init_tbase(&entry, name, DTZ, false);
  if (!dtz_tb) {
    fprintf(stderr, "Could not open %s.rtbz.\n", g_tablename);
    exit(EXIT_FAILURE);
  }

  if (dtz_tb->layout != LT_PAWN_P && dtz_tb->layout != LT_PAWN_PK) {
    fprintf(stderr, "Layout type %d is currently not supported.\n",
        dtz_tb->layout);
    exit(EXIT_FAILURE);
  }

  uint64_t *p;
  bool table_is_ok;

  mtx_init(&report_mutex, mtx_plain);

  size_t wdl_table_size = wdl_tb->layout == LT_PAWN_P
                        ? 63 * 62 * kslice_size
                        : 62 * kslice_size;
  uint8_t *wdl_tb_table = alloc_huge((wdl_table_size + 63) & ~0x3f);
  if (!wdl_tb_table)
    out_of_mem();

  size_t dtz_tb_table_size = dtz_tb->layout == LT_PAWN_P
                           ? 63 * 62 * kslice_size
                           : 62 * kslice_size;
  uint8_t *dtz_tb_table = alloc_huge((dtz_tb_table_size + 63) & ~0x3f);
  dtz_table = alloc_huge((kslice_size + 63) & ~0x3f);
  if (!dtz_tb_table || !dtz_table)
    out_of_mem();

  bool one_sided = dtz_tb->dist_format & WTM_OR_BTM;
  int one_sided_stm = (dtz_tb->dist_format & WTM_ONLY) ? WHITE : BLACK;
  bool wins_only = dtz_tb->dist_format & WIN_ONLY;
  bool has_stm[2] = {
    !one_sided || one_sided_stm == WHITE,
    !symmetric && (!one_sided || one_sided_stm == BLACK)
  };

  memset(dtz_stored, 0, sizeof dtz_stored);
  printf("DTZ format: ");
  if (one_sided) {
    printf("%s to move only\n", side[one_sided_stm]);
    dtz_stored[one_sided_stm][0] = dtz_stored[one_sided_stm][1] = true;
    dtz_stored[one_sided_stm][3] = dtz_stored[one_sided_stm][4] = true;
  }
  else if (wins_only) {
    printf("wins only\n");
    dtz_stored[WHITE][3] = dtz_stored[WHITE][4] = true;
    dtz_stored[BLACK][3] = dtz_stored[BLACK][4] = true;
  }
  else {
    printf("losses only\n");
    dtz_stored[WHITE][0] = dtz_stored[WHITE][1] = true;
    dtz_stored[BLACK][0] = dtz_stored[BLACK][1] = true;
  }

  dtz_freq = alloc_aligned(g_num_threads * sizeof *dtz_freq, 64);
  memset(dtz_freq, 0, g_num_threads * sizeof *dtz_freq);

  for (int q = 0; q < 24; q++) {
    g_slice.sq[2] = InvPawnFlip[0][q];
    make_dir(pawnstr[q]);
    change_dir(pawnstr[q]);

    for (int i = 1; i < 12; i++)
      if (!k16slice_buf[i])
        k16slice_buf[i] = alloc_k16slice();

    calc_sub_kslices(WHITE);
    calc_sub_kslices(BLACK);
    calc_psub_kslices();

    calc_pawn_capts();

    calc_illegal_both();

    calc_capt(WHITE, 2);
    calc_capt(BLACK, 2);
    calc_capt(WHITE, 1);
    calc_capt(BLACK, 1);
    calc_capt(WHITE, 0);
    calc_capt(BLACK, 0);
    calc_capt(WHITE, -1);
    calc_capt(BLACK, -1);
    calc_capt(WHITE, -2);
    calc_capt(BLACK, -2);

    // Release memory for K-slice bitmaps we don't need for now.
    for (int i = 7; i < 12; i++)
      if (k16slice_buf[i]) {
        free(k16slice_buf[i]);
        k16slice_buf[i] = nullptr;
      }

    tb_table = wdl_tb_table;

    for (int stm = 0; stm < 2; stm++) {
      char phase[64];
      snprintf(phase, sizeof phase, "loading %ctm WDL slices", "wb"[stm]);

      for (int i = 0; i < 5; i++)
        create_dir(-1, stm, wdl_name[i]);

      tb_ri[stm] = ri;
      g_slice.stm = stm;
      flip[stm] = !symmetric ? !wdl_tb->flipped : stm != WHITE;
      btm_side[stm] = (stm == WHITE) == flip[stm];

      for (int k = 0; k < ri.numsets; k++) {
        int pt = g_set_pt[k] ^ (flip[stm] << 3);
        for (int j = 0; j < g_pos.num; j++)
          if (wdl_tb->pt[j] == pt) {
            tb_ri[stm].first[k] = j;
            break;
          }
      }

      if (wdl_tb->layout == LT_PAWN_P) {
        show_progress(phase, 0, 1, false);
        decompress_kslice_wdl_p(wdl_tb, stm, q);
        show_progress(phase, 1, 1, true);
      } else {
        int num_done = 0;
        int psq = InvPawnFlip[0][q];
        for (int s = 0; s < 240; s++) {
          show_progress(phase, num_done++, 240, false);
          for (int i = 0; i < 5; i++)
            k16slice_clear_addr(k16slice_buf[i]);

          for (int r = 0; r < 16; r++) {
            g_slice.sq[2] = psq;
            g_slice.sq[0] = KK16Square[s][r][0];
            g_slice.sq[1] = KK16Square[s][r][1];
            g_slice.stm = stm;
            if (is_broken(&g_slice))
              continue;
            int ksq = g_slice.sq[stm];
            fill_wdl_pk_kslice(wdl_tb, stm, q, ksq);
          }
          finish_wdl_kslice(stm, s);
        }
        show_progress(phase, 240, 240, true);
      }

      for (int i = 4; i >= 0; i--)
        printf("wdl_cnt[%d][%d] = %lu\n", stm, i, wdl_cnt[stm][i]);
    }

    // Allocate memory for the K-slice bitmaps we released earlier.
    for (int i = 7; i < 12; i++)
      k16slice_buf[i] = alloc_k16slice();

    g_slice.sq[2] = InvPawnFlip[0][q];
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
    calc_pawn_moves();

    if (merged_table) {
      free(merged_table);
      merged_table = nullptr;
      if (merged_table2) {
        free(merged_table2);
        merged_table2 = nullptr;
      }
    }

    check_win_cwin(WHITE);
    check_cwin_draw(WHITE);
    check_draw_bloss(WHITE);
    check_bloss_loss(WHITE);
    check_loss(WHITE);

    check_win_cwin(BLACK);
    check_cwin_draw(BLACK);
    check_draw_bloss(BLACK);
    check_bloss_loss(BLACK);
    check_loss(BLACK);

    for (int i = 1; i < 12; i++)
      if (k16slice_buf[i]) {
        free(k16slice_buf[i]);
        k16slice_buf[i] = nullptr;
      }

    tb_table = dtz_tb_table;
    for (int stm = 0; stm < 2; stm++) {
      tb_ri[stm] = ri;
      g_slice.stm = stm;
      flip[stm] = !symmetric ? !dtz_tb->flipped : stm != WHITE;
      btm_side[stm] = (stm == WHITE) == flip[stm];

      for (int k = 0; k < ri.numsets; k++) {
        int pt = g_set_pt[k] ^ (flip[stm] << 3);
        for (int j = 0; j < g_pos.num; j++)
          if (dtz_tb->pt[j] == pt) {
            tb_ri[stm].first[k] = j;
            break;
          }
      }
    }

    int num_total = (dtz_tb->layout == LT_PAWN_P ? 1 : 63)
      * (has_stm[WHITE] + has_stm[BLACK]);
    int num_done = 0;
    if (dtz_tb->layout == LT_PAWN_P) {
      for (int stm = 0; stm < 2; stm++) {
        if (!has_stm[stm]) continue;
        show_progress("loading DTZ slices", num_done++, num_total, false);
        decompress_kslice_dtz_p(dtz_tb, stm, q);
      }
    } else {
      int psq = InvPawnFlip[0][q];
      for (int ksq = 0; ksq < 64; ksq++) {
        if (ksq == psq) continue;
        for (int stm = 0; stm < 2; stm++) {
          if (!has_stm[stm]) continue;
          show_progress("loading DTZ slices", num_done++, num_total, false);
          decompress_kslice_dtz_pk(dtz_tb, stm, q, ksq);
        }
      }
    }
    show_progress("loading DTZ slices", num_total, num_total, true);

    tb_table = wdl_tb_table;
    change_dir("..");
  }

  XXH128_hash_t wdl_checksum = XXH3_128bits(wdl_cnt, sizeof wdl_cnt);
  p = (uint64_t *)((uint8_t *)wdl_tb->data + 4);
  table_is_ok = memcmp(&wdl_checksum, p, 16) == 0;
  if (table_is_ok)
    printf("\x1b[92mWDL counts checksum is OK.\x1b[0m\n");
  else
    printf("\x1b[91mWDL counts checksum does not match.\x1b[0m\n");

  if (table_is_ok && num_fails == 0)
    printf("\x1b[92mWDL table is OK.\x1b[0m\n");
  else
    printf("\x1b[91mWDL table is not OK.\x1b[0m\n");

  free(wdl_tb_table);
  free(dtz_tb_table);
  free(dtz_table);
  kslice_free_buffers();
  mtx_destroy(&report_mutex);

  uint64_t freq[2][2][MAX_VAL + DRAW_RULE] = { 0 };
  for (int t = 0; t < g_num_threads; t++)
    for (int stm = 0; stm < 2; stm++)
      for (int wl = 0; wl < 2; wl++)
        for (int i = 0; i < MAX_VAL + DRAW_RULE; i++)
          freq[stm][wl][i] += dtz_freq[t][stm][wl][i];

  if (symmetric)
    memcpy(freq[1], freq[0], sizeof freq[1]);

  XXH128_hash_t dtz_checksum = freq_to_hash(MAX_VAL + DRAW_RULE,
      freq[0][0], freq[0][1], freq[1][0], freq[1][1]);

  p = (uint64_t *)((uint8_t *)dtz_tb->data + 20);
  table_is_ok = memcmp(&dtz_checksum, p, 16) == 0;
  if (table_is_ok)
    printf("\x1b[92mDTZ stored values checksum is OK.\x1b[0m\n");
  else
    printf("\x1b[91mDTZ stored values checksum does not match.\x1b[0m\n");

#if 0
  for (int stm = 0; stm < 2; stm++)
    for (int wl = 0; wl < 2; wl++)
      for (int i = 0; i < MAX_VAL + DRAW_RULE; i++)
        if (freq[stm][wl][i])
          printf("freq[%d][%d][%d] = %lu\n", stm, wl, i, freq[stm][wl][i]);
#endif

  strcat(name, ".txt");
  FILE *F = fopen(name, "rb");
  if (!F) {
    printf("%s not found.\n", name);
    exit(EXIT_SUCCESS);
  }

  printf("reading %s.txt\n", g_tablename);

  uint64_t stats[2][2][MAX_STATS / 2] = { 0 };
  char line[128];
  int stm = -1;

  while (fgets(line, sizeof line, F)) {
    if (strcmp(line, "White to move:\n") == 0) {
      stm = WHITE;
      continue;
    }

    if (strcmp(line, "Black to move:\n") == 0) {
      stm = BLACK;
      continue;
    }

    if (stm < 0)
      continue;

    uint64_t n, dummy;
    unsigned ply;
    char what[8];

    if (sscanf(line, "%lu positions %7s in %u ply.", &n, what, &ply) != 3) {
      if (sscanf(line, "%lu (%lu) positions %7s in %u ply.", &n, &dummy, what,
            &ply) != 4)
        continue;
    }

    int wl;
    if (strcmp(what, "win") == 0)
      wl = 0;
    else if (strcmp(what, "lose") == 0)
      wl = 1;
    else
      continue;

    if (ply >= MAX_STATS / 2)
      break;

    stats[stm][wl][ply] = n;
  }
  fclose(F);

  dtz_checksum = freq_to_hash(MAX_STATS / 2, stats[0][0], stats[0][1],
      stats[1][0], stats[1][1]);
  table_is_ok = memcmp(&dtz_checksum, p - 2, 16) == 0;
  if (table_is_ok)
    printf("\x1b[92mDTZ full stats checksum is OK.\x1b[0m\n");
  else
    printf("\x1b[91mDTZ full stats checksum does not match.\x1b[0m\n");
}
