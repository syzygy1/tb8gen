/*
  Copyright (c) 2026 Ronald de Man

  This file is distributed under the terms of the GNU GPL, version 2.
*/

#include <inttypes.h>
#include <stdint.h>
#include <stdlib.h>

#include "decompress.h"
#include "defs.h"
#include "index.h"
#include "kslice.h"
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

static uint64_t sub_cnt[2][5];

// Dummy definitions.
bool one_sided, wins_only;
int one_sided_stm;

INLINE void mark_king_uncaptures(int stm, int k, Bitboard occ,
    struct IdxState2 *is)
{
  alignas(64) Bitboard bb[8];
  uint8_t ksq[2] = { is->sq[0], is->sq[1] };

  // The stm king uncaptures, so we need to place a piece where the king was.
  is->bb[k + 1] |= bit(ksq[stm]);

  Bitboard b = king_attacks(ksq[stm]) & ~king_attacks(ksq[stm ^ 1]) & ~occ;
  while (b) {
    ksq[stm] = pop_lsb(&b);
    is->bb[0] = bit(ksq[0]) | bit(ksq[1]);
    int s = KKMap[ksq[0]][ksq[1]];
    int t = KK_transform[ksq[0]][ksq[1]];
    uint64_t idx;
    if (!t) {
      idx =  s < 441 ? rank_bb(is->bb, &ri) : rank_bb_ref(is->bb, &ri);
    } else {
      transform_set_bb(t, is->bb, bb);
      idx =  s < 441 ? rank_bb(bb, &ri) : rank_bb_ref(bb, &ri);
    }
    uint8_t *p = kslice_get_address(s);
    kslice_bit_set_atomic(p, idx);
  }

  // Restore the set occupancy bitboards.
  is->bb[0] = is->occ[0];
  is->bb[k + 1] ^= bit(is->sq[stm]);
}

// Uncapture a piece in set k by a piece in set j.
INLINE void mark_uncaptures(int stm, const int pc, int k,
    uint8_t *restrict const p, Bitboard occ, struct IdxState2 *is,
    const bool ref)
{
  int j = g_piece_set[stm][pc];
  if (j < 0) return;

  uint64_t idx0 = 0;
  int l;
  if (!ref) {
    l = min(j, k);
    for (int i = 0; i < l; i++)
      idx0 = idx0 * ri.factor[i] + is->sub[i];
  }

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
      uint64_t idx = ref ? rank_bb_ref(is->bb, &ri)
        : rank_bb_from(is->bb, idx0, l, is->occ[l], &ri);
      kslice_bit_set_atomic(p, idx);
      is->bb[j + 1] ^= from_bb;
      b ^= from_bb;
    }
    is->bb[j + 1] ^= to_bb;
    is->bb[k + 1] ^= to_bb;
  }
}

// Uncapture stm^1 king by an stm piece being added to set k.
INLINE void mark_uncapture_king(int stm, const int pc, int ksq,
    uint8_t *restrict const p, Bitboard occ, struct IdxState2 *is,
    const bool ref)
{
  int k = g_piece_set[stm][pc];
  if (k < 0) return;

  uint64_t idx0 = 0;
  if (!ref) {
    for (int i = 0; i < k; i++)
      idx0 = idx0 * ri.factor[i] + is->sub[i];
  }

  Bitboard b = piece_moves(pc, ksq, occ);
  while (b) {
    Bitboard from_bb = b & -b;
    is->bb[k + 1] ^= from_bb;
    uint64_t idx = ref ? rank_bb_ref(is->bb, &ri)
      : rank_bb_from(is->bb, idx0, k, is->occ[k], &ri);
    kslice_bit_set_atomic(p, idx);
    is->bb[k + 1] ^= from_bb;
    b ^= from_bb;
  }
}

// Return true if one move hits a set bit.
INLINE bool test_king_moves(int stm, Bitboard occ, struct IdxState2 *is)
{
  alignas(64) Bitboard bb[8];
  uint8_t ksq[2] = { is->sq[0], is->sq[1] };

  Bitboard b = king_attacks(ksq[stm]) & ~king_attacks(ksq[stm ^ 1]) & ~occ;
  while (b) {
    ksq[stm] = pop_lsb(&b);
    is->bb[0] = bit(ksq[0]) | bit(ksq[1]);
    int s = KKMap[ksq[0]][ksq[1]];
    int t = KK_transform[ksq[0]][ksq[1]];
    uint64_t idx;
    if (!t) {
      idx =  s < 441 ? rank_bb(is->bb, &ri) : rank_bb_ref(is->bb, &ri);
    } else {
      transform_set_bb(t, is->bb, bb);
      idx =  s < 441 ? rank_bb(bb, &ri) : rank_bb_ref(bb, &ri);
    }
    uint8_t *p = kslice_get_address(s);
    if (kslice_bit_test(p, idx)) {
      is->bb[0] = is->occ[0];
      return true;
    }
  }

  is->bb[0] = is->occ[0];
  return false;
}

INLINE bool test_moves(int stm, const int pc, uint8_t *restrict const p,
    Bitboard occ, struct IdxState2 *is, const bool ref)
{
  int k = g_piece_set[stm][pc];
  if (k < 0) return false;

  uint64_t idx0 = 0;
  if (!ref) {
    for (int i = 0; i < k; i++)
      idx0 = idx0 * ri.factor[i] + is->sub[i];
  }

  Bitboard bb = is->bb[k + 1];
  while (bb) {
    Bitboard from_bb = bb & -bb;
    is->bb[k + 1] ^= from_bb;
    bb ^= from_bb;
    Bitboard b = piece_moves(pc, lsb(from_bb), occ);
    while (b) {
      Bitboard to_bb = b & -b;
      is->bb[k + 1] ^= to_bb;
      uint64_t idx = ref ? rank_bb_ref(is->bb, &ri)
        : rank_bb_from(is->bb, idx0, k, is->occ[k], &ri);
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

static int work_slice, work_set;

static struct Work work_g_dynamic[2], work_g_static[2];
static struct Work work_capt_dynamic[MAX_SETS];

static constexpr uint64_t VERIFY_MIN_CHUNK = 1ULL << 9;
static constexpr int VERIFY_DYNAMIC_FACTOR = 4;

void init_verification_work(void)
{
  work_init(&work_g_dynamic[0], kslice_sizes[0], 0x1ff, WORK_DYNAMIC,
      VERIFY_DYNAMIC_FACTOR, VERIFY_MIN_CHUNK);
  work_init(&work_g_dynamic[1], kslice_sizes[1], 0x1ff, WORK_DYNAMIC,
      VERIFY_DYNAMIC_FACTOR, VERIFY_MIN_CHUNK);
  work_init(&work_g_static[0], kslice_sizes[0], 0x1ff, WORK_STATIC, 1,
      VERIFY_MIN_CHUNK);
  work_init(&work_g_static[1], kslice_sizes[1], 0x1ff, WORK_STATIC, 1,
      VERIFY_MIN_CHUNK);

  for (int k = 0; k < ri.numsets; k++)
    work_init(&work_capt_dynamic[k], capt_ri[k].sizes[0], 0x1ff, WORK_DYNAMIC,
        VERIFY_DYNAMIC_FACTOR, VERIFY_MIN_CHUNK);
}

static void calc_sub_worker(struct ThreadData *thread)
{
  struct IdxState2 is;
  int k = work_set;
  int stm = g_slice.stm;

  Position pos = g_pos;
  pos.sq[0] = g_slice.sq[0];
  pos.sq[1] = g_slice.sq[1];
  pos.stm = stm;
  int m = ri.last[k];
  int n = --pos.num;
  pos.pt[m] = pos.pt[n];

  uint8_t *restrict p[5];
  for (int i = 0; i < 5; i++)
    p[i] = kslice_sub_buf[i] + sub_offset[k];

  Bitboard occ = idx_state2_init(&is, thread->begin, g_slice.sq, &capt_ri[k],
      false);
  for (uint64_t idx = thread->begin, end = thread->end; idx < end;
      idx++, occ = idx_state2_inc(&is, &capt_ri[k]))
  {
    if (!idx_state2_legal(&is, stm, occ))
      continue;
    pos.occ = occ;
    idx_state2_to_sq(&is, pos.sq, &capt_ri[k]);
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

  for (int s = 0; s < 462; s++) {
    show_progress(phase, s, 462, false);

    uint64_t c[5];

    if (kslice_test_count(s, stm, name[4], -1, &c[4])) {
      for (int i = 0; i < 4; i++)
        kslice_test_count(s, stm, name[i], -1, &c[i]);
    } else {
      for (int i = 0; i < 5; i++)
        kslice_sub_clear_addr(kslice_sub_buf[i], stm);

      for (int t = 0; t < g_num_threads; t++)
        g_thread_data[t].cnt = 0;

      g_slice.sq[0] = KKSquare[s][0];
      g_slice.sq[1] = KKSquare[s][1];

      for (int k = 0; k < ri.numsets; k++) {
        if ((g_set_pt[k] >> 3) != stm)
          continue;
        work_set = k;
        run_threaded(calc_sub_worker, &work_capt_dynamic[k]);
      }

      for (int i = 0; i < 5; i++) {
        c[i] = kslice_sub_count_addr(kslice_sub_buf[i], stm);
        kslice_sub_write_addr(kslice_sub_buf[i], s, stm, name[i], c[i]);
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

  show_progress(phase, 462, 462, true);
  for (int i = 0; i < 5; i++)
    printf("sub_cnt[%d][%d] = %lu\n", stm, i, sub_cnt[stm][i]);
}

static void predecessors_sub_worker(struct ThreadData *thread)
{
  struct IdxState2 is;
  int stm = g_slice.stm;
  int k = work_set;
  int s = work_slice;

  uint64_t *restrict p =
    (uint64_t *)kslice_sub_get_address(s, k) + (thread->begin >> 6);
  uint8_t *restrict const q = kslice_get_address(s);

  uint64_t last = thread->begin;
  idx_state2_init(&is, last, g_slice.sq, &capt_ri[k], false);

  if (s < 441) {
    for (uint64_t idx = last, end = thread->end; idx < end; idx += 64) {
      uint64_t w = *p++;
      while (w) {
        uint64_t cur = idx + pop_lsb(&w);
        Bitboard occ = idx_state2_add(&is, cur - last, &capt_ri[k]);
        last = cur;
        // Uncapture by king.
        mark_king_uncaptures(stm, k, occ, &is);
        // Uncapture by non-king pieces.
        mark_uncaptures(stm, KNIGHT, k, q, occ, &is, false);
        mark_uncaptures(stm, BISHOP, k, q, occ, &is, false);
        mark_uncaptures(stm, ROOK  , k, q, occ, &is, false);
        mark_uncaptures(stm, QUEEN , k, q, occ, &is, false);
      }
    }
  } else {
    for (uint64_t idx = last, end = thread->end; idx < end; idx += 64) {
      uint64_t w = *p++;
      while (w) {
        uint64_t cur = idx + pop_lsb(&w);
        Bitboard occ = idx_state2_add(&is, cur - last, &capt_ri[k]);
        last = cur;
        // Uncapture by king.
        mark_king_uncaptures(stm, k, occ, &is);
        // Uncapture by non-king pieces.
        mark_uncaptures(stm, KNIGHT, k, q, occ, &is, true);
        mark_uncaptures(stm, BISHOP, k, q, occ, &is, true);
        mark_uncaptures(stm, ROOK  , k, q, occ, &is, true);
        mark_uncaptures(stm, QUEEN , k, q, occ, &is, true);
      }
    }
  }
}

// Uncapture from the loaded subtable of stm^1 positions to stm positions.
static void predecessors_sub(int stm, int s)
{
  work_slice = s;

  g_slice.stm = stm;
  g_slice.sq[0] = KKSquare[s][0];
  g_slice.sq[1] = KKSquare[s][1];

  // Loop through the sets from which a piece is removed.
  for (int k = 0; k < ri.numsets; k++) {
    if ((g_set_pt[k] >> 3) == stm)
      continue;
    work_set = k;
    run_threaded(predecessors_sub_worker, &work_capt_dynamic[k]);
  }
}

static void calc_capt(int stm, int wdl)
{
  if (!sub_cnt[stm ^ 1][2 - wdl])
    return;

  uint64_t cnt = 0;

  char capt_name[32], sub_name[32];
  strcat(strcpy(capt_name, "capt/"), wdl_name[2 + wdl]);
  strcat(strcpy(sub_name , "sub/" ), wdl_name[2 - wdl]);

  char str[64];
  sprintf(str, "%s/%c/done", capt_name, "wb"[stm]);
  FILE *F = fopen(str, "rb");
  if (F) {
    file_read(&cnt, 8, F);
    fclose(F);
    goto finished;
  }

  char phase[64];
  snprintf(phase, sizeof phase, "calculating %s captures for %s",
      wdl_name[2 + wdl], side[stm]);

  bool partial = dir_exists(-1, stm, capt_name), done = partial;

  create_dir(-1, stm, capt_name);

  struct KSliceIterator iter;
  uint64_t num;
  int num_done = 0;

  kslice_iter_init(&iter, stm);
  int s, s1;
  while (kslice_iter_next(&iter, &s)) {
    show_progress(phase, num_done++, 462, false);

    if (kslice_test(s, stm ^ 1, sub_name, -1)) {
      while (kslice_iter_in(&iter, &s1)) {
        if (partial && kslice_test_count(s1, stm, capt_name, -1, &num)) {
          cnt += num;
        } else {
          done = false;
          kslice_clear(s1);
        }
      }

      if (!done) {
        kslice_sub_read(s, s, stm ^ 1, sub_name);
        predecessors_sub(stm, s);
      }
    }

    while (kslice_iter_out(&iter, &s)) {
      if (!partial || !kslice_test_count(s, stm, capt_name, -1, &num)) {
        cnt += num = kslice_count(s);
        kslice_write(s, s, stm, capt_name, -1, num);
      }
    }
  }

  F = file_open_write(str);
  file_write(&cnt, 8, F);
  fclose(F);
  file_rename(str);

finished:
  snprintf(phase, sizeof phase, "%s %s captures: %lu", side[stm],
      wdl_name[2 + wdl], cnt);
  show_progress(phase, 462, 462, true);
}

static void calc_illegal_worker_tmpl(struct ThreadData *thread, const int pc,
    const bool ref)
{
  struct IdxState2 is;
  int stm = g_slice.stm;
  int k = g_piece_set[stm][pc];
  int ksq = g_slice.sq[stm ^ 1];

  uint8_t *restrict const p = kslice_buf[0];

  Bitboard occ = idx_state2_init(&is, thread->begin, g_slice.sq, &capt_ri[k],
      false);

  for (uint64_t idx = thread->begin, end = thread->end; idx < end;
      idx++, occ = idx_state2_inc(&is, &capt_ri[k]))
  {
    mark_uncapture_king(stm, pc, ksq, p, occ, &is, ref);
  }
}

static void calc_illegal_knight_worker(struct ThreadData *thread)
{
  calc_illegal_worker_tmpl(thread, KNIGHT, false);
}

static void calc_illegal_bishop_worker(struct ThreadData *thread)
{
  calc_illegal_worker_tmpl(thread, BISHOP, false);
}

static void calc_illegal_rook_worker(struct ThreadData *thread)
{
  calc_illegal_worker_tmpl(thread, ROOK, false);
}

static void calc_illegal_queen_worker(struct ThreadData *thread)
{
  calc_illegal_worker_tmpl(thread, QUEEN, false);
}

static void calc_illegal_knight_ref_worker(struct ThreadData *thread)
{
  calc_illegal_worker_tmpl(thread, KNIGHT, true);
}

static void calc_illegal_bishop_ref_worker(struct ThreadData *thread)
{
  calc_illegal_worker_tmpl(thread, BISHOP, true);
}

static void calc_illegal_rook_ref_worker(struct ThreadData *thread)
{
  calc_illegal_worker_tmpl(thread, ROOK, true);
}

static void calc_illegal_queen_ref_worker(struct ThreadData *thread)
{
  calc_illegal_worker_tmpl(thread, QUEEN, true);
}

void calc_illegal(int stm)
{
  bool partial = dir_exists(-1, stm, "illegal");

  create_dir(-1, stm, "illegal");

  g_slice.stm = stm;

  char phase[64];
  snprintf(phase, sizeof phase, "calculating %s illegal positions",
      side[stm]);

  uint64_t cnt = 0, num;

  for (int s = 0; s < 462; s++) {
    show_progress(phase, s, 462, false);

    if (partial && kslice_test_count(s, stm, "illegal", -1, &num)) {
      cnt += num;
      continue;
    }

    g_slice.sq[0] = KKSquare[s][0];
    g_slice.sq[1] = KKSquare[s][1];

    kslice_clear_addr(kslice_buf[0], s);

    int k;
    if (s < 441) {
      if ((k = g_piece_set[stm][KNIGHT]) >= 0)
        run_threaded(calc_illegal_knight_worker, &work_capt_dynamic[k]);
      if ((k = g_piece_set[stm][BISHOP]) >= 0)
        run_threaded(calc_illegal_bishop_worker, &work_capt_dynamic[k]);
      if ((k = g_piece_set[stm][ROOK]) >= 0)
        run_threaded(calc_illegal_rook_worker, &work_capt_dynamic[k]);
      if ((k = g_piece_set[stm][QUEEN]) >= 0)
        run_threaded(calc_illegal_queen_worker, &work_capt_dynamic[k]);
    } else {
      if ((k = g_piece_set[stm][KNIGHT]) >= 0)
        run_threaded(calc_illegal_knight_ref_worker, &work_capt_dynamic[k]);
      if ((k = g_piece_set[stm][BISHOP]) >= 0)
        run_threaded(calc_illegal_bishop_ref_worker, &work_capt_dynamic[k]);
      if ((k = g_piece_set[stm][ROOK]) >= 0)
        run_threaded(calc_illegal_rook_ref_worker, &work_capt_dynamic[k]);
      if ((k = g_piece_set[stm][QUEEN]) >= 0)
        run_threaded(calc_illegal_queen_ref_worker, &work_capt_dynamic[k]);
    }

    cnt += num = kslice_count_addr(kslice_buf[0], s);
    kslice_write_addr(kslice_buf[0], s, stm, "illegal", -1, num);
  }

  snprintf(phase, sizeof phase, "%lu %s illegal positions", cnt, side[stm]);
  show_progress(phase, 462, 462, true);
}

static mtx_t report_mutex;
static int num_fails;

static void report_fail(int s, uint64_t idx, Bitboard occ,
    const struct IdxState2 *is)
{
  mtx_lock(&report_mutex);
  char fenstr[64];
  Position pos = g_pos;
  pos.sq[0] = g_slice.sq[0];
  pos.sq[1] = g_slice.sq[1];
  pos.stm = g_slice.stm;
  pos.occ = occ;
  idx_state2_to_sq(is, pos.sq, &ri);
  pos_to_fen(&pos, fenstr, false);
  int wdl = probe_wdl(&pos, -2, 2);
  printf("\nslice = %d, idx = %lu, wdl = %d\n%s\n", s, idx, wdl, fenstr);
  num_fails++;
  if (num_fails == 10)
    exit(EXIT_FAILURE);
  mtx_unlock(&report_mutex);
}

static bool check_dtz_W101(struct IdxState2 *is, Bitboard occ)
{
  Position pos = g_pos;
  pos.sq[0] = g_slice.sq[0];
  pos.sq[1] = g_slice.sq[1];
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

static bool check_dtz_L101(struct IdxState2 *is, Bitboard occ)
{
  Position pos = g_pos;
  pos.sq[0] = g_slice.sq[0];
  pos.sq[1] = g_slice.sq[1];
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

INLINE void check_zero_worker_tmpl(struct ThreadData *thread, const int T,
    const bool ref)
{
  struct IdxState2 is;
  int stm = g_slice.stm;
  int s = work_slice;

  uint64_t *restrict p = (uint64_t *)kslice_get_address(-1);
  uint8_t *restrict const q = kslice_get_address(s);

  p += thread->begin >> 6;
  uint64_t last = thread->begin;
  idx_state2_init(&is, last, g_slice.sq, &ri, ref);
  for (uint64_t idx = last, end = thread->end; idx < end; idx += 64, p++) {
    uint64_t w = *p;
    while (w) {
      uint64_t cur = idx + pop_lsb(&w);
      Bitboard occ = ref ? unrank_bb_ref(cur, is.bb, &ri)
        : idx_state2_add(&is, cur - last, &ri);
      last = cur;

      if (   !test_moves(stm, QUEEN , q, occ, &is, ref)
          && !test_moves(stm, ROOK  , q, occ, &is, ref)
          && !test_moves(stm, BISHOP, q, occ, &is, ref)
          && !test_moves(stm, KNIGHT, q, occ, &is, ref)
          && !test_king_moves(stm, occ, &is))
        continue;

      switch (T) {
      case CZ_REGULAR:
        report_fail(s, cur, occ, &is);
        break;
      case CZ_CWIN:
        if (!check_dtz_W101(&is, occ))
          report_fail(s, cur, occ, &is);
        break;
      }
    }
  }
}

static void check_zero_regular_worker(struct ThreadData *thread)
{
  check_zero_worker_tmpl(thread, CZ_REGULAR, false);
}

static void check_zero_cwin_worker(struct ThreadData *thread)
{
  check_zero_worker_tmpl(thread, CZ_CWIN, false);
}

static void check_zero_ref_regular_worker(struct ThreadData *thread)
{
  check_zero_worker_tmpl(thread, CZ_REGULAR, true);
}

static void check_zero_ref_cwin_worker(struct ThreadData *thread)
{
  check_zero_worker_tmpl(thread, CZ_CWIN, true);
}

static void check_zero(int stm, int s, const int T)
{
  work_slice = s;
  g_slice.stm = stm;
  g_slice.sq[0] = KKSquare[s][0];
  g_slice.sq[1] = KKSquare[s][1];

  switch (T) {
  case CZ_REGULAR:
    if (s < 441)
      run_threaded(check_zero_regular_worker, &work_g_dynamic[0]);
    else
      run_threaded(check_zero_ref_regular_worker, &work_g_dynamic[1]);
    break;
  case CZ_CWIN:
    if (s < 441)
      run_threaded(check_zero_cwin_worker, &work_g_dynamic[0]);
    else
      run_threaded(check_zero_ref_cwin_worker, &work_g_dynamic[1]);
    break;
  }
}

enum { CO_REGULAR, CO_DRAW, CO_BLOSS, CO_LOSS };

INLINE void check_one_worker_tmpl(struct ThreadData *thread, const int T,
    const bool ref)
{
  struct IdxState2 is;
  int stm = g_slice.stm;
  int s = work_slice;

  uint64_t *restrict p = (uint64_t *)kslice_get_address(-1);
  uint8_t *restrict const q = kslice_get_address(s);

  p += thread->begin >> 6;
  uint64_t last = thread->begin;
  idx_state2_init(&is, last, g_slice.sq, &ri, ref);
  for (uint64_t idx = last, end = thread->end; idx < end; idx += 64, p++) {
    uint64_t w = *p;
    while (w) {
      uint64_t cur = idx + pop_lsb(&w);
      Bitboard occ = ref ? unrank_bb_ref(cur, is.bb, &ri)
          : idx_state2_add(&is, cur - last, &ri);
      last = cur;

      if (   test_moves(stm, QUEEN , q, occ, &is, ref)
          || test_moves(stm, ROOK  , q, occ, &is, ref)
          || test_moves(stm, BISHOP, q, occ, &is, ref)
          || test_moves(stm, KNIGHT, q, occ, &is, ref)
          || test_king_moves(stm, occ, &is))
        continue;

      switch (T) {
      case CO_REGULAR:
        report_fail(s, cur, occ, &is);
        break;
      case CO_DRAW:
        if (    idx_state2_has_legal_moves(&is, stm, occ)
            || !idx_state2_legal(&is, stm ^ 1, occ))
          report_fail(s, cur, occ, &is);
        break;
      case CO_BLOSS:
        if (!check_dtz_L101(&is, occ))
          report_fail(s, cur, occ, &is);
        break;
      case CO_LOSS:
        if (idx_state2_legal(&is, stm ^ 1, occ))
          report_fail(s, cur, occ, &is);
        break;
      default:
        abort();
      }
    }
  }
}

static void check_one_regular_worker(struct ThreadData *thread)
{
  check_one_worker_tmpl(thread, CO_REGULAR, false);
}

static void check_one_draw_worker(struct ThreadData *thread)
{
  check_one_worker_tmpl(thread, CO_DRAW, false);
}

static void check_one_bloss_worker(struct ThreadData *thread)
{
  check_one_worker_tmpl(thread, CO_BLOSS, false);
}

static void check_one_loss_worker(struct ThreadData *thread)
{
  check_one_worker_tmpl(thread, CO_LOSS, false);
}

static void check_one_regular_ref_worker(struct ThreadData *thread)
{
  check_one_worker_tmpl(thread, CO_REGULAR, true);
}

static void check_one_draw_ref_worker(struct ThreadData *thread)
{
  check_one_worker_tmpl(thread, CO_DRAW, true);
}

static void check_one_bloss_ref_worker(struct ThreadData *thread)
{
  check_one_worker_tmpl(thread, CO_BLOSS, true);
}

static void check_one_loss_ref_worker(struct ThreadData *thread)
{
  check_one_worker_tmpl(thread, CO_LOSS, true);
}

static void check_one(int stm, int s, const int T)
{
  work_slice = s;
  g_slice.stm = stm;
  g_slice.sq[0] = KKSquare[s][0];
  g_slice.sq[1] = KKSquare[s][1];

  switch (T) {
  case CO_REGULAR:
    if (s < 441)
      run_threaded(check_one_regular_worker, &work_g_dynamic[0]);
    else
      run_threaded(check_one_regular_ref_worker, &work_g_dynamic[1]);
    break;
  case CO_DRAW:
    if (s < 441)
      run_threaded(check_one_draw_worker, &work_g_dynamic[0]);
    else
      run_threaded(check_one_draw_ref_worker, &work_g_dynamic[1]);
    break;
  case CO_BLOSS:
    if (s < 441)
      run_threaded(check_one_bloss_worker, &work_g_dynamic[0]);
    else
      run_threaded(check_one_bloss_ref_worker, &work_g_dynamic[1]);
    break;
  case CO_LOSS:
    if (s < 441)
      run_threaded(check_one_loss_worker, &work_g_dynamic[0]);
    else
      run_threaded(check_one_loss_ref_worker, &work_g_dynamic[1]);
    break;
  default:
    abort();
  }
}

static void check_win_cwin(int stm)
{
  struct KSliceIterator iter;

  char phase[64];
  snprintf(phase, sizeof phase, "check %s win, cwin-", side[stm]);
  int num_done = 0;

  kslice_iter_init(&iter, stm);
  int s, s1;
  while (kslice_iter_next(&iter, &s)) {
    show_progress(phase, num_done++, 462, false);

    while (kslice_iter_in(&iter, &s1))
      kslice_read(s1, s1, stm ^ 1, "loss", -1);

    kslice_read(-1, s, stm, "wins", -1);
    kslice_read_andnot(-1, s, stm, "capt/win", -1);
    check_one(stm, s, CO_REGULAR);

    kslice_read(-1, s, stm, "cwin", -1);
    kslice_read_andnot(-1, s, stm, "capt/cwin", -1);
    check_zero(stm, s, CZ_CWIN);

    while (kslice_iter_out(&iter, &s));
  }
  show_progress(phase, 462, 462, true);
}

static void check_cwin_draw(int stm)
{
  struct KSliceIterator iter;

  char phase[64];
  snprintf(phase, sizeof phase, "check %s cwin+, draw-", side[stm]);
  int num_done = 0;

  kslice_iter_init(&iter, stm);
  int s, s1;
  while (kslice_iter_next(&iter, &s)) {
    show_progress(phase, num_done++, 462, false);

    while (kslice_iter_in(&iter, &s1)) {
      kslice_read(s1, s1, stm ^ 1, "loss", -1);
      kslice_read_or(s1, s1, stm ^ 1, "bloss", -1);
    }

    // We already know that no cwin position reaches "loss" (with the
    // exception of dtz==101 -> dtz==100). So it does not hurt to test against
    // "loss" and "bloss" here to ensure that one cwin position reaches
    // "bloss". As a bonus, no special check is needed for the 101->100 case.
    kslice_read(-1, s, stm, "cwin", -1);
    kslice_read_andnot(-1, s, stm, "capt/cwin", -1);
    check_one(stm, s, CO_REGULAR);

    kslice_read(-1, s, stm, "draw", -1);
    kslice_read_andnot(-1, s, stm, "capt/draw", -1);
    check_zero(stm, s, CZ_REGULAR);

    while (kslice_iter_out(&iter, &s));
  }
  show_progress(phase, 462, 462, true);
}

static void check_draw_bloss(int stm)
{
  struct KSliceIterator iter;

  char phase[64];
  snprintf(phase, sizeof phase, "check %s draw+, bloss-", side[stm]);
  int num_done = 0;

  kslice_iter_init(&iter, stm);
  int s, s1;
  while (kslice_iter_next(&iter, &s)) {
    show_progress(phase, num_done++, 462, false);

    while (kslice_iter_in(&iter, &s1)) {
      kslice_read(s1, s1, stm ^ 1, "loss", -1);
      kslice_read_or(s1, s1, stm ^ 1, "bloss", -1);
      kslice_read_or(s1, s1, stm ^ 1, "draw", -1);
    }

    // Again, we already checked that drawn positions do not hit "loss" and
    // "bloss". So we can combine draw+ with bloss-.
    kslice_read(-1, s, stm, "draw", -1);
    kslice_read_andnot(-1, s, stm, "capt/draw", -1);
    check_one(stm, s, CO_DRAW);

    kslice_read(-1, s, stm, "bloss", -1);
    kslice_read_andnot(-1, s, stm, "capt/bloss", -1);
    check_zero(stm, s, CZ_REGULAR);

    while (kslice_iter_out(&iter, &s));
  }
  show_progress(phase, 462, 462, true);
}

static void check_bloss_loss(int stm)
{
  struct KSliceIterator iter;

  char phase[64];
  snprintf(phase, sizeof phase, "check %s bloss+, loss-", side[stm]);
  int num_done = 0;

  kslice_iter_init(&iter, stm);
  int s, s1;
  while (kslice_iter_next(&iter, &s)) {
    show_progress(phase, num_done++, 462, false);

    // Load "loss" | "bloss" | "draw" | "cwin".
    while (kslice_iter_in(&iter, &s1)) {
      kslice_read(s1, s1, stm ^ 1, "win", -1);
      kslice_read_or(s1, s1, stm ^ 1, "illegal", -1);
      kslice_not(s1);
    }

    kslice_read(-1, s, stm, "bloss", -1);
    kslice_read_andnot(-1, s, stm, "capt/bloss", -1);
    check_one(stm, s, CO_BLOSS);

    kslice_read(-1, s, stm, "loss", -1);
    check_zero(stm, s, CZ_REGULAR);

    while (kslice_iter_out(&iter, &s));
  }
  show_progress(phase, 462, 462, true);
}

static void check_loss(int stm)
{
  struct KSliceIterator iter;

  char phase[64];
  snprintf(phase, sizeof phase, "check %s loss+", side[stm]);
  int num_done = 0;

  kslice_iter_init(&iter, stm);
  int s, s1;
  while (kslice_iter_next(&iter, &s)) {
    show_progress(phase, num_done++, 462, false);

    while (kslice_iter_in(&iter, &s1))
      kslice_read(s1, s1, stm ^ 1, "win", -1);

    kslice_read(-1, s, stm, "loss", -1);
    kslice_read_andnot(-1, s, stm, "capt/loss", -1);
    check_one(stm, s, CO_LOSS);

    while (kslice_iter_out(&iter, &s));
  }
  show_progress(phase, 462, 462, true);
}

static uint8_t *tb_table, *dtz_table;
static struct TbTable2 *table_info;
static struct RankInfo tb_ri[2];

static void init_perm_rank_ri(struct RankInfo *rank_ri,
    const struct TbTable2 *table, const struct RankInfo *dst_ri)
{
  *rank_ri = *table->ri;
  for (int k = 0; k < table->ri->numsets; k++) {
    int j = 0;
    for (; j < dst_ri->numsets; j++)
      if (table->first[k] == dst_ri->first[j])
        break;
    assert(j < dst_ri->numsets);
    rank_ri->perm[k] = j + 1;
  }
}

static uint64_t depermute_462_idx(struct IdxState2 *is,
    const struct RankInfo *rank_ri, int t)
{
  alignas(64) Bitboard bb[8];
  const Bitboard *set_bb = is->bb;
  if (t) {
    transform_set_bb(t, is->bb, bb);
    set_bb = bb;
  }
  return perm_rank_bb(set_bb, rank_ri);
}

static uint64_t depermute_462_ref_idx(struct IdxState2 *is,
    const struct RankInfo *rank_ri, int t)
{
  alignas(64) Bitboard bb[8];
  const Bitboard *set_bb = is->bb;
  if (t) {
    transform_set_bb(t, is->bb, bb);
    set_bb = bb;
  }
  return perm_rank_bb_ref(set_bb, rank_ri);
}

INLINE void store_wdl(uint64_t idx, uint8_t *restrict tmp, int v)
{
  uint64_t **const p = (uint64_t **)kslice_buf;
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
static void depermute_wdl_462_worker(struct ThreadData *thread)
{
  alignas(64) uint8_t tmp[64];
  uint8_t *restrict src = tb_table;
  const struct RankInfo *dst_ri = &tb_ri[g_slice.stm];
  struct RankInfo rank_ri;
  struct IdxState2 is;

  uint64_t idx_dec_buf[NUM];

  init_perm_rank_ri(&rank_ri, table_info, dst_ri);
  int t = KK_transform[g_slice.sq[0]][g_slice.sq[1]];
  idx_state2_init(&is, thread->begin, g_slice.sq, dst_ri, false);

  uint64_t idx = thread->begin, end = thread->end;
  int fill = 0;
  int head = 0; // Next buffered element to consume.

  // Fill pipeline.
  for (; fill < NUM && idx < end;
      fill++, idx++, idx_state2_inc(&is, dst_ri))
  {
    uint64_t idx_dec = depermute_462_idx(&is, &rank_ri, t);
    __builtin_prefetch(src + idx_dec, 0, 3);
    idx_dec_buf[fill] = idx_dec;
  }

  // Steady-state pipeline.
  for (; idx < end; idx++, idx_state2_inc(&is, dst_ri)) {
    uint64_t idx_dec = depermute_462_idx(&is, &rank_ri, t);
    __builtin_prefetch(src + idx_dec, 0, 3);
    store_wdl(idx - NUM, tmp, src[idx_dec_buf[head]]);
    idx_dec_buf[head] = idx_dec;
    head = (head + 1) & (NUM - 1);
  }

  // Drain pipeline.
  for (uint64_t out = idx - fill; fill-- > 0; out++) {
    store_wdl(out, tmp, src[idx_dec_buf[head]]);
    head = (head + 1) & (NUM - 1);
  }
  for (; idx & 0x3f; idx++)
    store_wdl(idx, tmp, 0);
}

static void depermute_wdl_462_ref_worker(struct ThreadData *thread)
{
  alignas(64) uint8_t tmp[64];
  uint8_t *restrict src = tb_table;
  const struct RankInfo *dst_ri = &tb_ri[g_slice.stm];
  struct RankInfo rank_ri;
  struct IdxState2 is;

  uint64_t idx_dec_buf[NUM];

  init_perm_rank_ri(&rank_ri, table_info, dst_ri);
  int t = KK_transform[g_slice.sq[0]][g_slice.sq[1]];
  is.bb[0] = bit(g_slice.sq[0]) | bit(g_slice.sq[1]);

  uint64_t idx = thread->begin, end = thread->end;
  int fill = 0;
  int head = 0; // Next buffered element to consume.

  // Fill pipeline.
  for (; fill < NUM && idx < end; fill++, idx++) {
    unrank_bb_ref(idx, is.bb, dst_ri);
    uint64_t idx_dec = depermute_462_ref_idx(&is, &rank_ri, t);
    __builtin_prefetch(src + idx_dec, 0, 3);
    idx_dec_buf[fill] = idx_dec;
  }

  // Steady-state pipeline.
  for (; idx < end; idx++) {
    unrank_bb_ref(idx, is.bb, dst_ri);
    uint64_t idx_dec = depermute_462_ref_idx(&is, &rank_ri, t);
    __builtin_prefetch(src + idx_dec, 0, 3);
    store_wdl(idx - NUM, tmp, src[idx_dec_buf[head]]);
    idx_dec_buf[head] = idx_dec;
    head = (head + 1) & (NUM - 1);
  }

  // Drain pipeline.
  for (uint64_t out = idx - fill; fill-- > 0; out++) {
    store_wdl(out, tmp, src[idx_dec_buf[head]]);
    head = (head + 1) & (NUM - 1);
  }
  for (; idx & 0x3f; idx++)
    store_wdl(idx, tmp, 0);
}

static bool flip[2], btm_side[2];

static uint64_t wdl_cnt[2][5];

static void decompress_kslice_wdl_462(struct Tbase *tb, int stm, int s)
{
  char str[64];
  create_name(str, s, stm, wdl_name[4], -1);
  if (file_exists(str)) {
    uint64_t num;
    for (int i = 0; i < 5; i++)
      if (kslice_test_count(s, stm, wdl_name[i], -1, &num))
        wdl_cnt[stm][i] += num;
    return;
  }

  int tsq = KKMap[g_slice.sq[0]][g_slice.sq[1]];
  int t = !symmetric ? 2 * tsq + btm_side[stm] : tsq;
  if (!tb->table[t])
    tb->table[t] = init_new_table(tb, g_pos.num, WDL, t, tsq);
  struct TbTable2 *table = tb->table[t];

  if (!table->precomp) {
    struct TbTableConst *tbl = (struct TbTableConst *)table;
    int v = tbl->const_value;
    kslice_set_addr(kslice_buf[v], s);
    for (int i = 0; i < 5; i++) {
      if (i == v) continue;
      kslice_clear_addr(kslice_buf[i], s);
    }
  }
  else {
    decompress_table(table, tb_table, kslice_sizes[s >= 441]);

    table_info = table;

    if (s < 441)
      run_threaded(depermute_wdl_462_worker, &work_g_static[0]);
    else
      run_threaded(depermute_wdl_462_ref_worker, &work_g_static[1]);
  }

  // Reconstruct proper win/cwin/draw/bloss/loss slices of legal positions.
  kslice_read_addr(kslice_buf[5], s, stm, "illegal", -1);
  kslice_read_addr(kslice_buf[6], s, stm, "capt/win", -1);
  kslice_abc_addr(kslice_buf[4], kslice_buf[5], kslice_buf[6], s);
  kslice_read_addr(kslice_buf[6], s, stm, "capt/cwin", -1);
  kslice_abc_addr(kslice_buf[3], kslice_buf[5], kslice_buf[6], s);
  kslice_read_addr(kslice_buf[6], s, stm, "capt/draw", -1);
  kslice_abc_addr(kslice_buf[2], kslice_buf[5], kslice_buf[6], s);
  kslice_read_addr(kslice_buf[6], s, stm, "capt/bloss", -1);
  kslice_abc_addr(kslice_buf[1], kslice_buf[5], kslice_buf[6], s);
  kslice_andnot_addr(kslice_buf[0], kslice_buf[5], s);

  // Count and save win/cwin/draw/bloss/loss slices.
  uint64_t num;
  for (int i = 0; i < 5; i++) {
    wdl_cnt[stm][i] += num = kslice_count_addr(kslice_buf[i], s);
    kslice_write_addr(kslice_buf[i], s, stm, wdl_name[i], -1, num);
  }
}

static uint8_t Offset10[64];

static void init_perm_rank_10(int8_t *restrict perm,
    const struct TbTable2 *table, const struct RankInfo *dst_ri)
{
  const struct RankInfo10 *ri = table->ri_10;
  for (int k = 0; k < ri->numsets; k++) {
    if (ri->mult[k] == 0) {
      perm[k] = -1;
      continue;
    }
    int j = 0;
    for (; j < dst_ri->numsets; j++)
      if (table->first[k] == dst_ri->first[j])
        break;
    assert(j < dst_ri->numsets);
    perm[k] = j + 1;
  }
}

static uint64_t rank_10_bb(const Bitboard *restrict set_bb,
    const struct TbTable2 *table, const int8_t *restrict perm, int k2sq)
{
  const struct RankInfo10 *ri = table->ri_10;
  Bitboard occ = set_bb[0];
  uint64_t idx = 0;
  for (int k = 0; k < ri->numsets; k++) {
    size_t s = 0;
    int m = ri->mult[k];
    if (m == 0)
      s = Offset10[k2sq];
    else {
      Bitboard bb = set_bb[perm[k]];
      Bitboard b = _pext_u64(bb, ~occ);
      occ |= bb;
      for (int j = 1; b; j++)
        s += Binomial[j][pop_lsb(&b)];
    }
    idx = idx * ri->factor[k] + s;
  }
  return idx;
}

static void init_10_ksq(uint8_t *restrict sq)
{
  int s = KKMap[g_slice.sq[0]][g_slice.sq[1]];
  sq[0] = KKSquare[s][0];
  sq[1] = KKSquare[s][1];
  if (flip[g_slice.stm] != btm_side[g_slice.stm])
    Swap(sq[0], sq[1]);
}

static uint64_t depermute_10_idx(struct IdxState2 *is,
    const struct TbTable2 *table, const int8_t *restrict perm, int t, int k2sq)
{
  alignas(64) Bitboard bb[8];
  const Bitboard *set_bb = is->bb;
  if (t) {
    transform_set_bb(t, is->bb, bb);
    set_bb = bb;
  }
  return rank_10_bb(set_bb, table, perm, k2sq);
}

INLINE int k2sq_10(void)
{
  return g_slice.sq[g_slice.stm ^ 1];
}

// Do the inverse of permute10.c/_tmpl.c
static void depermute_wdl_10_worker(struct ThreadData *thread)
{
  alignas(64) uint8_t tmp[64];
  uint8_t sq[MAX_PIECES];
  const uint8_t *restrict src = tb_table;
  const struct TbTable2 *table = table_info;
  const struct RankInfo *dst_ri = &tb_ri[g_slice.stm];
  struct IdxState2 is;
  int8_t perm[MAX_SETS + 1];

  uint64_t idx_dec_buf[NUM];

  init_10_ksq(sq);
  init_perm_rank_10(perm, table, dst_ri);
  int t = KK_transform[sq[0]][sq[1]];
  int k2sq = k2sq_10();

  idx_state2_init(&is, thread->begin, sq, dst_ri, false);

  uint64_t idx = thread->begin, end = thread->end;
  int fill = 0;
  int head = 0; // Next buffered element to consume.

  // We loop through the K-slice bitmaps.
  // So we'll need a separate function for depermute_wdl_10_ref_worker().

  // Fill pipeline.
  for (; fill < NUM && idx < end;
      fill++, idx++, idx_state2_inc(&is, dst_ri))
  {
    // Now basically just probe the wdl_10 table slice.
    uint64_t idx_dec = depermute_10_idx(&is, table, perm, t, k2sq);
    __builtin_prefetch(src + idx_dec, 0, 3);
    idx_dec_buf[fill] = idx_dec;
  }

  // Steady-state pipeline.
  for (; idx < end; idx++, idx_state2_inc(&is, dst_ri)) {
    uint64_t idx_dec = depermute_10_idx(&is, table, perm, t, k2sq);
    __builtin_prefetch(src + idx_dec, 0, 3);
    store_wdl(idx - NUM, tmp, src[idx_dec_buf[head]]);
    idx_dec_buf[head] = idx_dec;
    head = (head + 1) & (NUM - 1);
  }

  // Drain pipeline.
  for (uint64_t out = idx - fill; fill-- > 0; out++) {
    store_wdl(out, tmp, src[idx_dec_buf[head]]);
    head = (head + 1) & (NUM - 1);
  }
  for (; idx & 0x3f; idx++)
    store_wdl(idx, tmp, 0);
}

static void depermute_wdl_10_ref_worker(struct ThreadData *thread)
{
  alignas(64) uint8_t tmp[64];
  uint8_t sq[MAX_PIECES];
  uint8_t *restrict src = tb_table;
  const struct TbTable2 *table = table_info;
  const struct RankInfo *dst_ri = &tb_ri[g_slice.stm];
  struct IdxState2 is;
  int8_t perm[MAX_SETS + 1];

  uint64_t idx_dec_buf[NUM];

  init_10_ksq(sq);
  init_perm_rank_10(perm, table, dst_ri);
  int t = KK_transform[sq[0]][sq[1]];
  is.bb[0] = bit(sq[0]) | bit(sq[1]);
  int k2sq = k2sq_10();

  uint64_t idx = thread->begin, end = thread->end;
  int fill = 0;
  int head = 0; // Next buffered element to consume.

  // Fill pipeline.
  for (; fill < NUM && idx < end; fill++, idx++) {
    unrank_bb_ref(idx, is.bb, dst_ri);
    uint64_t idx_dec = depermute_10_idx(&is, table, perm, t, k2sq);
    __builtin_prefetch(src + idx_dec, 0, 3);
    idx_dec_buf[fill] = idx_dec;
  }

  // Steady-state pipeline.
  for (; idx < end; idx++) {
    unrank_bb_ref(idx, is.bb, dst_ri);
    uint64_t idx_dec = depermute_10_idx(&is, table, perm, t, k2sq);
    __builtin_prefetch(src + idx_dec, 0, 3);
    store_wdl(idx - NUM, tmp, src[idx_dec_buf[head]]);
    idx_dec_buf[head] = idx_dec;
    head = (head + 1) & (NUM - 1);
  }

  // Drain pipeline.
  for (uint64_t out = idx - fill; fill-- > 0; out++) {
    store_wdl(out, tmp, src[idx_dec_buf[head]]);
    head = (head + 1) & (NUM - 1);
  }
  for (; idx & 0x3f; idx++)
    store_wdl(idx, tmp, 0);
}

static constexpr uint8_t knum[] = { 58, 58, 58, 55, 55, 55, 33, 30, 30, 30 };
static constexpr uint8_t MAX_KNUM = 58;

static void decompress_kslice_wdl_10(struct Tbase *tb, int stm, int k)
{
  int t = !symmetric ? 2 * k + btm_side[stm] : k;
  if (!tb->table[t])
    tb->table[t] = init_new_table(tb, g_pos.num, WDL, t, k);
  struct TbTable2 *table = tb->table[t];

  for (int l = 0, n = 0; l < 64; l++) {
    if (KKIdx[k][l] < 0) continue;
    Offset10[l] = n++;
  }

  int const_value = -1;
  if (!table->precomp) {
    struct TbTableConst *tbl = (struct TbTableConst *)table;
    const_value = tbl->const_value;
  }
  else {
    decompress_table(table, tb_table, knum[k] * kslice_sizes[0]);
    table_info = table;
  }

  // Now decompress K-slice by K-slice inside the big slice.
  // Note that everything is mixed up in the big slice.
  for (int l = 0; l < 64; l++) {
    if (KKIdx[k][l] < 0)
      continue;

    g_slice.sq[stm] = InvTriangle[k];
    g_slice.sq[stm ^ 1] = l;
    g_slice.stm = stm;

    int s = KKMap[g_slice.sq[0]][g_slice.sq[1]];

    if (const_value >= 0) {
      kslice_set_addr(kslice_buf[const_value], s);
      for (int i = 0; i < 5; i++) {
        if (i == const_value) continue;
        kslice_clear_addr(kslice_buf[i], s);
      }
    } else {
      if (s < 441)
        run_threaded(depermute_wdl_10_worker, &work_g_static[0]);
      else
        run_threaded(depermute_wdl_10_ref_worker, &work_g_static[1]);
    }

    // Reconstruct proper win/cwin/draw/bloss/loss slices of legal positions.
    kslice_read_addr(kslice_buf[5], s, stm, "illegal", -1);
    kslice_read_addr(kslice_buf[6], s, stm, "capt/win", -1);
    kslice_abc_addr(kslice_buf[4], kslice_buf[5], kslice_buf[6], s);
    kslice_read_addr(kslice_buf[6], s, stm, "capt/cwin", -1);
    kslice_abc_addr(kslice_buf[3], kslice_buf[5], kslice_buf[6], s);
    kslice_read_addr(kslice_buf[6], s, stm, "capt/draw", -1);
    kslice_abc_addr(kslice_buf[2], kslice_buf[5], kslice_buf[6], s);
    kslice_read_addr(kslice_buf[6], s, stm, "capt/bloss", -1);
    kslice_abc_addr(kslice_buf[1], kslice_buf[5], kslice_buf[6], s);
    kslice_andnot_addr(kslice_buf[0], kslice_buf[5], s);

    // Count and save win/cwin/draw/bloss/loss slices.
    uint64_t num;
    for (int i = 0; i < 5; i++) {
      wdl_cnt[stm][i] += num = kslice_count_addr(kslice_buf[i], s);
      kslice_write_addr(kslice_buf[i], s, stm, wdl_name[i], -1, num);
    }
  }
}

void depermute_dtz_462_worker(struct ThreadData *thread)
{
  uint8_t *restrict src = tb_table;
  uint8_t *restrict dst = dtz_table;
  const struct RankInfo *dst_ri = &tb_ri[g_slice.stm];
  struct RankInfo rank_ri;
  struct IdxState2 is;

  uint64_t idx_dec_buf[NUM];

  init_perm_rank_ri(&rank_ri, table_info, dst_ri);
  int t = KK_transform[g_slice.sq[0]][g_slice.sq[1]];
  idx_state2_init(&is, thread->begin, g_slice.sq, dst_ri, false);

  uint64_t idx = thread->begin, end = thread->end;
  int fill = 0;
  int head = 0; // Next buffered element to consume.

  // Fill pipeline.
  for (; fill < NUM && idx < end;
      fill++, idx++, idx_state2_inc(&is, dst_ri))
  {
    uint64_t idx_dec = depermute_462_idx(&is, &rank_ri, t);
    __builtin_prefetch(src + idx_dec, 0, 3);
    idx_dec_buf[fill] = idx_dec;
  }

  // Steady-state pipeline.
  for (; idx < end; idx++, idx_state2_inc(&is, dst_ri)) {
    uint64_t idx_dec = depermute_462_idx(&is, &rank_ri, t);
    __builtin_prefetch(src + idx_dec, 0, 3);
    dst[idx - NUM] = src[idx_dec_buf[head]];
    idx_dec_buf[head] = idx_dec;
    head = (head + 1) & (NUM - 1);
  }

  // Drain pipeline.
  for (uint64_t out = idx - fill; fill-- > 0; out++) {
    dst[out] = src[idx_dec_buf[head]];
    head = (head + 1) & (NUM - 1);
  }
}

void depermute_dtz_462_ref_worker(struct ThreadData *thread)
{
  uint8_t *restrict src = tb_table;
  uint8_t *restrict dst = dtz_table;
  const struct RankInfo *dst_ri = &tb_ri[g_slice.stm];
  struct RankInfo rank_ri;
  struct IdxState2 is;

  uint64_t idx_dec_buf[NUM];

  init_perm_rank_ri(&rank_ri, table_info, dst_ri);
  int t = KK_transform[g_slice.sq[0]][g_slice.sq[1]];
  is.bb[0] = bit(g_slice.sq[0]) | bit(g_slice.sq[1]);

  uint64_t idx = thread->begin, end = thread->end;
  int fill = 0;
  int head = 0; // Next buffered element to consume.

  // Fill pipeline.
  for (; fill < NUM && idx < end; fill++, idx++) {
    unrank_bb_ref(idx, is.bb, dst_ri);
    uint64_t idx_dec = depermute_462_ref_idx(&is, &rank_ri, t);
    __builtin_prefetch(src + idx_dec, 0, 3);
    idx_dec_buf[fill] = idx_dec;
  }

  // Steady-state pipeline.
  for (; idx < end; idx++) {
    unrank_bb_ref(idx, is.bb, dst_ri);
    uint64_t idx_dec = depermute_462_ref_idx(&is, &rank_ri, t);
    __builtin_prefetch(src + idx_dec, 0, 3);
    dst[idx - NUM] = src[idx_dec_buf[head]];
    idx_dec_buf[head] = idx_dec;
    head = (head + 1) & (NUM - 1);
  }

  // Drain pipeline.
  for (uint64_t out = idx - fill; fill-- > 0; out++) {
    dst[out] = src[idx_dec_buf[head]];
    head = (head + 1) & (NUM - 1);
  }
}

void depermute_dtz_10_worker(struct ThreadData *thread)
{
  uint8_t sq[MAX_PIECES];
  const uint8_t *restrict src = tb_table;
  uint8_t *restrict dst = dtz_table;
  const struct TbTable2 *table = table_info;
  const struct RankInfo *dst_ri = &tb_ri[g_slice.stm];
  struct IdxState2 is;
  int8_t perm[MAX_SETS + 1];

  uint64_t idx_dec_buf[NUM];

  init_10_ksq(sq);
  init_perm_rank_10(perm, table, dst_ri);
  int t = KK_transform[sq[0]][sq[1]];
  int k2sq = k2sq_10();

  idx_state2_init(&is, thread->begin, sq, dst_ri, false);

  uint64_t idx = thread->begin, end = thread->end;
  int fill = 0;
  int head = 0; // Next buffered element to consume.

  // Fill pipeline.
  for (; fill < NUM && idx < end;
      fill++, idx++, idx_state2_inc(&is, dst_ri))
  {
    uint64_t idx_dec = depermute_10_idx(&is, table, perm, t, k2sq);
    __builtin_prefetch(src + idx_dec, 0, 3);
    idx_dec_buf[fill] = idx_dec;
  }

  // Steady-state pipeline.
  for (; idx < end; idx++, idx_state2_inc(&is, dst_ri)) {
    uint64_t idx_dec = depermute_10_idx(&is, table, perm, t, k2sq);
    __builtin_prefetch(src + idx_dec, 0, 3);
    dst[idx - NUM] = src[idx_dec_buf[head]];
    idx_dec_buf[head] = idx_dec;
    head = (head + 1) & (NUM - 1);
  }

  // Drain pipeline.
  for (uint64_t out = idx - fill; fill-- > 0; out++) {
    dst[out] = src[idx_dec_buf[head]];
    head = (head + 1) & (NUM - 1);
  }
}

void depermute_dtz_10_ref_worker(struct ThreadData *thread)
{
  uint8_t sq[MAX_PIECES];
  const uint8_t *restrict src = tb_table;
  uint8_t *restrict dst = dtz_table;
  const struct TbTable2 *table = table_info;
  const struct RankInfo *dst_ri = &tb_ri[g_slice.stm];
  struct IdxState2 is;
  int8_t perm[MAX_SETS + 1];

  uint64_t idx_dec_buf[NUM];

  init_10_ksq(sq);
  init_perm_rank_10(perm, table, dst_ri);
  int t = KK_transform[sq[0]][sq[1]];
  is.bb[0] = bit(sq[0]) | bit(sq[1]);
  int k2sq = k2sq_10();

  uint64_t idx = thread->begin, end = thread->end;
  int fill = 0;
  int head = 0; // Next buffered element to consume.

  // Fill pipeline.
  for (; fill < NUM && idx < end; fill++, idx++) {
    unrank_bb_ref(idx, is.bb, dst_ri);
    uint64_t idx_dec = depermute_10_idx(&is, table, perm, t, k2sq);
    __builtin_prefetch(src + idx_dec, 0, 3);
    idx_dec_buf[fill] = idx_dec;
  }

  // Steady-state pipeline.
  for (; idx < end; idx++) {
    unrank_bb_ref(idx, is.bb, dst_ri);
    uint64_t idx_dec = depermute_10_idx(&is, table, perm, t, k2sq);
    __builtin_prefetch(src + idx_dec, 0, 3);
    dst[idx - NUM] = src[idx_dec_buf[head]];
    idx_dec_buf[head] = idx_dec;
    head = (head + 1) & (NUM - 1);
  }

  // Drain pipeline.
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

void update_stats_worker(struct ThreadData *thread)
{
  const uint16_t *restrict map = freq_map;
  uint64_t *restrict freq = dtz_freq[thread->thread_id][g_slice.stm][win_loss];
  const uint64_t *restrict p = (uint64_t *)kslice_buf[0];
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

  int bound[4];
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
    kslice_read_addr(kslice_buf[0], s, stm, wdl_name[i], -1);
    if (i == 3 || i == 4) {
      // Remove CAPT_WIN and CAPT_CWIN
      char capt_name[32];
      strcat(strcpy(capt_name, "capt/"), wdl_name[i]);
      kslice_read_andnot_addr(kslice_buf[0], s, stm, capt_name, -1);
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
          freq_map[j] = offset[m] + table->map_dtz16[table->map_dtz_idx[m] + j];
        break;
    }
    run_threaded(update_stats_worker, &work_g_static[s >= 441]);
  }
}

void decompress_kslice_dtz_462(struct Tbase *tb, int stm, int s)
{
  int tsq = KKMap[g_slice.sq[0]][g_slice.sq[1]];
  int t = (tb->dist_format & TWO_SIDED) ? 2 * tsq + btm_side[stm] : tsq;
  if (!tb->table[t])
    tb->table[t] = init_new_table(tb, g_pos.num, DTZ, t, tsq);
  struct DtzTable2 *table = tb->table[t];

  if (!table->precomp) {
    struct TbTableConst *tbl = (struct TbTableConst *)table;
    int v = tbl->const_value;
    memset(dtz_table, v, kslice_sizes[s >= 441]);
    table = &unmapped_dtz_table;
  }
  else {
    decompress_table((struct TbTable2 *)table, tb_table, kslice_sizes[s >= 441]);
    table_info = (struct TbTable2 *)table;
    if (s < 441)
      run_threaded(depermute_dtz_462_worker, &work_g_static[0]);
    else
      run_threaded(depermute_dtz_462_ref_worker, &work_g_static[1]);
  }

  update_dtz_stats(table, stm, s);
}

void decompress_kslice_dtz_10(struct Tbase *tb, int stm, int k)
{
  int t = (tb->dist_format & TWO_SIDED) ? 2 * k + btm_side[stm] : k;
  if (!tb->table[t])
    tb->table[t] = init_new_table(tb, g_pos.num, DTZ, t, k);
  struct DtzTable2 *table = tb->table[t];

  for (int l = 0, n = 0; l < 64; l++) {
    if (KKIdx[k][l] < 0) continue;
    Offset10[l] = n++;
  }

  int const_value = -1;
  if (!table->precomp) {
    struct TbTableConst *tbl = (struct TbTableConst *)table;
    const_value = tbl->const_value;
    table = &unmapped_dtz_table;
  }
  else {
    decompress_table((struct TbTable2 *)table, tb_table,
        knum[k] * kslice_sizes[0]);
    table_info = (struct TbTable2 *)table;
  }

  for (int l = 0; l < 64; l++) {
    if (KKIdx[k][l] < 0)
      continue;

    g_slice.sq[stm] = InvTriangle[k];
    g_slice.sq[stm ^ 1] = l;
    g_slice.stm = stm;

    int s = KKMap[g_slice.sq[0]][g_slice.sq[1]];
    if (const_value >= 0)
      memset(dtz_table, const_value, kslice_sizes[s >= 441]);
    else if (s < 441)
      run_threaded(depermute_dtz_10_worker, &work_g_static[0]);
    else
      run_threaded(depermute_dtz_10_ref_worker, &work_g_static[1]);

    update_dtz_stats(table, stm, s);
  }
}

void verify(void)
{
  // Initialize and mmap() tablebase file.
  struct TbEntry entry = {
    .has_pawns = false,
    .symmetric = symmetric,
    .num = g_pos.num
  };

  char name[64];
  if (g_tablename[0] != '/')
    snprintf(name, sizeof name, "../%s", g_tablename);
  else
    strcpy(name, g_tablename);

  struct Tbase *tb = init_tbase(&entry, name, WDL, false);
  if (!tb) {
    fprintf(stderr, "Could not open %s.rtbw.\n", g_tablename);
    exit(EXIT_FAILURE);
  }

  if (tb->layout != LT_PIECE_K && tb->layout != LT_PIECE_KK) {
    fprintf(stderr, "Layout type %d is currently not supported.\n", tb->layout);
    exit(EXIT_FAILURE);
  }

  calc_sub_kslices(WHITE);
  calc_sub_kslices(BLACK);

  calc_illegal(WHITE);
  calc_illegal(BLACK);

  calc_capt(WHITE, 2);
  calc_capt(BLACK, 2);
  calc_capt(WHITE, 1);
  calc_capt(BLACK, 1);
  calc_capt(WHITE, 0);
  calc_capt(BLACK, 0);
  calc_capt(WHITE, -1);
  calc_capt(BLACK, -1);
  calc_capt(WHITE, -2);
  if (!symmetric)
    calc_capt(BLACK, -2);

  // Release memory for K-slice bitmaps we don't need for now.
  for (int i = 7; i < 20; i++)
    if (kslice_buf[i]) {
      free(kslice_buf[i]);
      kslice_buf[i] = nullptr;
    }

  // Allocate memory for decompressed WDL data.
  size_t wdl_table_size = tb->layout == LT_PIECE_K
                        ? MAX_KNUM * kslice_sizes[0]
                        : kslice_sizes[0];
  tb_table = alloc_huge((wdl_table_size + 63) & ~0x3f);
  if (!tb_table)
    out_of_mem();

  for (int stm = 0; stm < 2; stm++) {
    char phase[64];
    snprintf(phase, sizeof phase, "loading %ctm WDL slices", "wb"[stm]);

    for (int i = 0; i < 5; i++)
      create_dir(-1, stm, wdl_name[i]);

    if (tb->layout == LT_PIECE_KK) {

      tb_ri[stm] = ri;
      g_slice.stm = stm;
      flip[stm] = !symmetric ? !tb->flipped : stm != WHITE;
      btm_side[stm] = (stm == WHITE) == flip[stm];

      for (int k = 0; k < ri.numsets; k++) {
        int pt = g_set_pt[k] ^ (flip[stm] << 3);
        for (int j = 0; j < g_pos.num; j++)
          if (tb->pt[j] == pt) {
            tb_ri[stm].first[k] = j;
            break;
          }
      }

      for (int s = 0; s < 462; s++) {
        show_progress(phase, s, 462, false);

        g_slice.sq[0] = KKSquare[s][0];
        g_slice.sq[1] = KKSquare[s][1];

        if (flip[stm])
          Swap(g_slice.sq[0], g_slice.sq[1]);

        decompress_kslice_wdl_462(tb, stm, s);
      }
      show_progress(phase, 462, 462, true);

    } else { /* tb_layout == LT_PIECE_K */

      tb_ri[stm] = ri;
      g_slice.stm = stm;
      flip[stm] = !symmetric ? !tb->flipped : stm != WHITE;
      btm_side[stm] = (stm == WHITE) == flip[stm];

      for (int k = 0; k < ri.numsets; k++) {
        int pt = g_set_pt[k] ^ (flip[stm] << 3);
        for (int j = 0; j < g_pos.num; j++)
          if (tb->pt[j] == pt) {
            tb_ri[stm].first[k] = j;
            break;
          }
      }

      for (int k = 0; k < 10; k++) {
        show_progress(phase, k, 10, false);

        g_slice.sq[stm] = InvTriangle[k];
        decompress_kslice_wdl_10(tb, stm, k);
      }
      show_progress(phase, 10, 10, true);

    }

    for (int i = 4; i >= 0; i--)
      printf("wdl_cnt[%d][%d] = %lu\n", stm, i, wdl_cnt[stm][i]);
  }

  free(tb_table);

  // Allocate memory for the K-slice bitmaps we released earlier.
  for (int i = 7; i < 20; i++)
    kslice_buf[i] = kslice_alloc();

  XXH128_hash_t wdl_checksum = XXH3_128bits(wdl_cnt, sizeof wdl_cnt);
  uint64_t *p = (uint64_t *)((uint8_t *)tb->data + 4);
  bool table_is_ok = memcmp(&wdl_checksum, p, 16) == 0;
  if (table_is_ok)
    printf("\x1b[92mWDL counts checksum is OK.\x1b[0m\n");
  else
    printf("\x1b[91mWDL counts checksum does not match.\x1b[0m\n");

  mtx_init(&report_mutex, mtx_plain);

  check_win_cwin(WHITE);
  check_cwin_draw(WHITE);
  check_draw_bloss(WHITE);
  check_bloss_loss(WHITE);
  check_loss(WHITE);

  if (!symmetric) {
    check_win_cwin(BLACK);
    check_cwin_draw(BLACK);
    check_draw_bloss(BLACK);
    check_bloss_loss(BLACK);
    check_loss(BLACK);
  }

  if (table_is_ok && num_fails == 0)
    printf("\x1b[92mWDL table is OK.\x1b[0m\n");
  else
    printf("\x1b[91mWDL table is not OK.\x1b[0m\n");

  // Release memory for K-slice bitmaps we don't need for now.
  for (int i = 1; i < 20; i++)
    if (kslice_buf[i]) {
      free(kslice_buf[i]);
      kslice_buf[i] = nullptr;
    }

  tb = init_tbase(&entry, name, DTZ, false);
  if (!tb) {
    fprintf(stderr, "Could not open %s.rtbz.\n", g_tablename);
    exit(EXIT_FAILURE);
  }

  if (tb->layout != LT_PIECE_K && tb->layout != LT_PIECE_KK) {
    fprintf(stderr, "Layout type %d is currently not supported.\n", tb->layout);
    exit(EXIT_FAILURE);
  }

  // Allocate memory for decompressed and depermuted DTZ data
  // (1 byte per position).
  size_t dtz_tb_table_size = tb->layout == LT_PIECE_K
                           ? MAX_KNUM * kslice_sizes[0]
                           : kslice_sizes[0];
  tb_table = alloc_huge((dtz_tb_table_size + 63) & ~0x3f);
  dtz_table = alloc_huge((kslice_size + 63) & ~0x3f);
  if (!tb_table || !dtz_table)
    out_of_mem();

  bool one_sided = tb->dist_format & WTM_OR_BTM;
  int one_sided_stm = (tb->dist_format & WTM_ONLY) ? WHITE : BLACK;
  bool wins_only = tb->dist_format & WIN_ONLY;
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

  if (tb->layout == LT_PIECE_KK) {

    for (int stm = 0; stm < 2; stm++) {
      tb_ri[stm] = ri;
      g_slice.stm = stm;
      flip[stm] = !symmetric ? !tb->flipped : stm != WHITE;
      btm_side[stm] = (stm == WHITE) == flip[stm];

      for (int k = 0; k < ri.numsets; k++) {
        int pt = g_set_pt[k] ^ (flip[stm] << 3);
        for (int j = 0; j < g_pos.num; j++)
          if (tb->pt[j] == pt) {
            tb_ri[stm].first[k] = j;
            break;
          }
      }
    }

    for (int s = 0; s < 462; s++) {
      show_progress("loading DTZ slices", s, 462, false);

      for (int stm = 0; stm < 2; stm++) {
        if (!has_stm[stm]) continue;

        g_slice.sq[0] = KKSquare[s][0];
        g_slice.sq[1] = KKSquare[s][1];
        g_slice.stm = stm;

        if (flip[stm])
          Swap(g_slice.sq[0], g_slice.sq[1]);

        decompress_kslice_dtz_462(tb, stm, s);
      }
    }
    show_progress("loading DTZ slices", 462, 462, true);
  }
  else if (tb->layout == LT_PIECE_K) {
    for (int stm = 0; stm < 2; stm++) {
      tb_ri[stm] = ri;
      g_slice.stm = stm;
      flip[stm] = !symmetric ? !tb->flipped : stm != WHITE;
      btm_side[stm] = (stm == WHITE) == flip[stm];

      for (int k = 0; k < ri.numsets; k++) {
        int pt = g_set_pt[k] ^ (flip[stm] << 3);
        for (int j = 0; j < g_pos.num; j++)
          if (tb->pt[j] == pt) {
            tb_ri[stm].first[k] = j;
            break;
          }
      }
    }

    int num_done = 0;
    int num_total = 10 * (has_stm[WHITE] + has_stm[BLACK]);
    for (int k = 0; k < 10; k++) {
      for (int stm = 0; stm < 2; stm++) {
        if (!has_stm[stm]) continue;

        show_progress("loading DTZ slices", num_done++, num_total, false);

        g_slice.sq[stm] = InvTriangle[k];
        g_slice.stm = stm;

        decompress_kslice_dtz_10(tb, stm, k);
      }
    }
    show_progress("loading DTZ slices", num_total, num_total, true);
  }

  free(tb_table);
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

  p = (uint64_t *)((uint8_t *)tb->data + 20);
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
