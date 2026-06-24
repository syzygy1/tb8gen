/*
  Copyright (c) 2026 Ronald de Man

  This file is distributed under the terms of the GNU GPL, version 2.
*/

#include <inttypes.h>
#include <stdint.h>
#include <stdlib.h>
#include <immintrin.h>

#include "defs.h"
#include "generate.h"
#include "index.h"
#include "kslice.h"
#include "movegen.h"
#include "probe.h"
#include "tb8gen.h"
#include "threads.h"
#include "types.h"
#include "util.h"

const char wdl_name[5][8] = { "loss", "bloss", "draw", "cwin", "win" };
const char side[2][8] = { "white", "black" };
const char clr_L[4][4] = { "31", "32", "33", "34" };
const char clr_W[4][4] = { "94", "93", "92", "91" };

uint64_t sub_cnt[2][5];
int max_iteration;

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

INLINE void mark_king_unmoves(int stm, Bitboard occ, struct IdxState2 *is)
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
    kslice_bit_set_atomic(p, idx);
  }

  is->bb[0] = is->occ[0];
}

INLINE void mark_unmoves(int stm, const int pc, uint8_t *restrict const p,
    Bitboard occ, struct IdxState2 *is, const bool ref)
{
  int k = g_piece_set[stm][pc];
  if (k < 0) return;

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
      kslice_bit_set_atomic(p, idx);
      is->bb[k + 1] ^= to_bb;
      b ^= to_bb;
    }
    is->bb[k + 1] ^= from_bb;
  }
}

// Uncapture stm king by a piece being added to set k.
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

INLINE bool check_king_moves(int stm, Bitboard occ, struct IdxState2 *is)
{
  alignas(64) Bitboard bb[8];
  uint8_t ksq[2] = { is->sq[0], is->sq[1] };

  Bitboard b = king_attacks(ksq[stm]) & ~king_attacks(ksq[stm ^ 1]);

#if 1
  Bitboard attacks = b & occ & ~is->bb[0];
  while (attacks) {
    Bitboard to_bb = attacks & -attacks;
    attacks ^= to_bb;
    int j = 0;
    while (!(to_bb & is->bb[j + 1]))
      j++;
    if ((g_set_pt[j] >> 3) == stm) continue;
    ksq[stm] = lsb(to_bb);
    is->bb[0] = bit(ksq[0]) | bit(ksq[1]);
    is->bb[j + 1] ^= to_bb;
    int s = KKMap[ksq[0]][ksq[1]];
    int t = KK_transform[ksq[0]][ksq[1]];
    uint64_t idx;
    if (!t) {
      idx = rank_bb(is->bb, &capt_ri[j]);
    } else {
      transform_set_bb(t, is->bb, bb);
      idx = rank_bb(bb, &capt_ri[j]);
    }
    uint8_t *p = kslice_sub_get_address(s, j);
    is->bb[j + 1] ^= to_bb;
    if (!kslice_bit_test(p, idx)) {
      is->bb[0] = is->occ[0];
      return false;
    }
  }
#endif

  b &= ~occ;
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
    if (!kslice_bit_test(p, idx)) {
      is->bb[0] = is->occ[0];
      return false;
    }
  }

  is->bb[0] = is->occ[0];
  return true;
}

INLINE bool check_moves(int stm, const int pc, int s,
    uint8_t *restrict const p, Bitboard occ, struct IdxState2 *is,
    const bool ref)
{
  int k = g_piece_set[stm][pc];
  if (k < 0) return true;

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

    Bitboard b = piece_attacks(pc, lsb(from_bb), occ);

#if 1
    Bitboard attacks = b & occ & ~is->bb[0];
    while (attacks) {
      Bitboard to_bb = attacks & -attacks;
      attacks ^= to_bb;
      int j = 0;
      while (!(to_bb & is->bb[j + 1]))
        j++;
      if (!((g_set_pt[j] ^ g_set_pt[k]) & 8)) continue;
      is->bb[k + 1] ^= to_bb;
      is->bb[j + 1] ^= to_bb;
      uint64_t idx = rank_bb(is->bb, &capt_ri[j]);
      is->bb[k + 1] ^= to_bb;
      is->bb[j + 1] ^= to_bb;
      uint8_t *restrict q = kslice_sub_get_address(s, j);
      if (!kslice_bit_test(q, idx)) {
        is->bb[k + 1] ^= from_bb;
        return false;
      }
    }
#endif

    b &= ~occ;
    while (b) {
      Bitboard to_bb = b & -b;
      is->bb[k + 1] ^= to_bb;
      uint64_t idx = ref ? rank_bb_ref(is->bb, &ri)
        : rank_bb_from(is->bb, idx0, k, is->occ[k], &ri);
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

static int work_slice, work_set;

static struct Work work_g_dynamic[2], work_g_static[2];
static struct Work work_capt_dynamic[MAX_SETS];

static constexpr uint64_t GENERATE_MIN_CHUNK = 1ULL << 9;
static constexpr int GENERATE_DYNAMIC_FACTOR = 4;

void init_generation_work(void)
{
  work_init(&work_g_dynamic[0], kslice_sizes[0], 0x1ff, WORK_DYNAMIC,
      GENERATE_DYNAMIC_FACTOR, GENERATE_MIN_CHUNK);
  work_init(&work_g_dynamic[1], kslice_sizes[1], 0x1ff, WORK_DYNAMIC,
      GENERATE_DYNAMIC_FACTOR, GENERATE_MIN_CHUNK);
  work_init(&work_g_static[0], kslice_sizes[0], 0x1ff, WORK_STATIC, 1,
      GENERATE_MIN_CHUNK);
  work_init(&work_g_static[1], kslice_sizes[1], 0x1ff, WORK_STATIC, 1,
      GENERATE_MIN_CHUNK);

  for (int k = 0; k < ri.numsets; k++)
    work_init(&work_capt_dynamic[k], capt_ri[k].sizes[0], 0x1ff, WORK_DYNAMIC,
        GENERATE_DYNAMIC_FACTOR, GENERATE_MIN_CHUNK);
}

static void calc_sub_worker(struct ThreadData *thread)
{
  struct IdxState2 is;
  Position pos = g_pos;
  int k = work_set;
  int m = ri.last[k];
  int n = --pos.num;
  uint64_t cnt = 0;

  pos.pt[m] = pos.pt[n];
  uint8_t *restrict p[5];
  for (int i = 0; i < 5; i++)
    p[i] = kslice_sub_buf[i] + sub_offset[k];

  Bitboard occ = idx_state2_init(&is, thread->begin, pos.sq, &capt_ri[k]);
  for (uint64_t idx = thread->begin, end = thread->end; idx < end;
      idx++, occ = idx_state2_inc(&is, &capt_ri[k]))
  {
    if (!idx_state2_legal(&is, pos.stm, occ)) {
      // Include illegal positions in sub_cwin and sub_win.
      kslice_bit_set(p[3], idx);
      kslice_bit_set(p[4], idx);
      cnt++;
    } else {
      pos.occ = occ;
      idx_state2_to_sq(&is, pos.sq, &capt_ri[k]);
      pos.sq[m] = pos.sq[n];
      int v = probe_wdl(&pos, -2, 2);
      kslice_bit_set(p[v + 2], idx);
      // Add sub_win to sub_cwin.
      if (v == 2)
        kslice_bit_set(p[3], idx);
    }
  }

  thread->cnt += cnt;
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

  char name[5][16];
  for (int i = 0; i < 5; i++) {
    strcat(strcpy(name[i], "sub/"), wdl_name[i]);
    create_dir(-1, stm, name[i]);
  }

  g_pos.stm = stm;

  for (int s = 0; s < 462; s++) {
    show_progress(phase, s, 462, false);

    uint64_t c[5], cnt_ilgl = 0;

    if (kslice_test_count(s, stm, name[4], -1, &c[4])) {
      for (int i = 0; i < 4; i++)
        kslice_test_count(s, stm, name[i], -1, &c[i]);
    } else {
      for (int i = 0; i < 5; i++)
        kslice_sub_clear_addr(kslice_sub_buf[i], stm);

      for (int t = 0; t < g_num_threads; t++)
        g_thread_data[t].cnt = 0;

      g_pos.sq[0] = KKSquare[s][0];
      g_pos.sq[1] = KKSquare[s][1];

      for (int k = 0; k < ri.numsets; k++) {
        if ((g_pos.pt[ri.first[k]] >> 3) != stm)
          continue;
        work_set = k;
        run_threaded(calc_sub_worker, &work_capt_dynamic[k]);
      }

      for (int t = 0; t < g_num_threads; t++)
        cnt_ilgl += g_thread_data[t].cnt;

      for (int i = 0; i < 5; i++) {
        c[i] = kslice_sub_count_addr(kslice_sub_buf[i], stm);
        kslice_sub_write_addr(kslice_sub_buf[i], s, stm, name[i], c[i]);
      }
    }

    sub_cnt[stm][0] += c[0];
    sub_cnt[stm][1] += c[1];
    sub_cnt[stm][2] += c[2];
    sub_cnt[stm][3] += c[3] - c[4];
    sub_cnt[stm][4] += c[4] - cnt_ilgl;
  }

  F = file_open_write(done);
  for (int i = 0; i < 5; i++)
    file_write(&sub_cnt[stm][i], 8, F);
  fclose(F);
  file_rename(done);

  show_progress(phase, 462, 462, true);
//  for (int i = 0; i < 5; i++)
//    printf("sub_cnt[%d][%d] = %lu\n", stm, i, sub_cnt[stm][i]);
}

static bool work_legality;

static void predecessors_sub_worker(struct ThreadData *thread)
{
  struct IdxState2 is;
  Position pos = g_pos;
  int stm = pos.stm;
  pos.stm ^= 1;
  int k = work_set;
  int s = work_slice;
  bool legality = work_legality;

  uint64_t *restrict p =
    (uint64_t *)kslice_sub_get_address(s, k) + (thread->begin >> 6);
  uint8_t *restrict const q = kslice_get_address(s);

  uint64_t last = thread->begin;
  idx_state2_init(&is, last, pos.sq, &capt_ri[k]);

  if (s < 441) {
    for (uint64_t idx = last, end = thread->end; idx < end; idx += 64) {
      uint64_t w = *p++;
      while (w) {
        uint64_t cur = idx + pop_lsb(&w);
        Bitboard occ = idx_state2_add(&is, cur - last, &capt_ri[k]);
        last = cur;
        if (legality && !idx_state2_legal(&is, pos.stm, occ))
          continue;
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
        if (legality && !idx_state2_legal(&is, pos.stm, occ))
          continue;
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
static void predecessors_sub(int stm, int s, bool legality)
{
  work_slice = s;
  work_legality = legality;

  g_pos.stm = stm;
  g_pos.sq[0] = KKSquare[s][0];
  g_pos.sq[1] = KKSquare[s][1];

  // Loop through the sets from which a piece is removed.
  for (int k = 0; k < ri.numsets; k++) {
    if ((g_set_pt[k] >> 3) == stm)
      continue;
    work_set = k;
    run_threaded(predecessors_sub_worker, &work_capt_dynamic[k]);
  }
}

static int wins_full[2][462];
static int wins_checked[2][462];
static uint64_t replay_cost[2][462];
static uint64_t write_cost[2][462];

static void add_wins(int s, int stm, int n, uint64_t cost)
{
  if ((replay_cost[stm][s] + cost) * 1.0 <= write_cost[stm][s])
    return;

  kslice_or(-1, s);
  kslice_write(-1, s, stm, "wins", n, UINT64_MAX);
  replay_cost[stm][s] = 0;

  kslice_delete(s, stm, "wins", wins_full[stm][s]);
  for (int i = max(wins_full[stm][s] + 1, wins_checked[stm][s]); i < n; i++)
    kslice_delete(s, stm, "wins", i);
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
      if (kslice_test(slice, stm, "wins", i))
        break;
    if (i > wins_checked[stm][slice]) {
      for (int j = i - 1; j >= wins_full[stm][slice]; j--)
        kslice_delete(slice, stm, "wins", j);
      wins_full[stm][slice] = i;
    }
    wins_checked[stm][slice] = n;
  }

  kslice_read(s, slice, stm, "wins", wins_full[stm][slice]);
  write_cost[stm][slice] = kslice_read_cost;
  kslice_read_cost = 0;
  for (int i = wins_full[stm][slice] + 1; i <= n; i++)
    if (g_stats[stm][stats_n(i)])
      kslice_read_or(s, slice, stm, "W", i);
  replay_cost[stm][slice] += kslice_read_cost;
}

static void calc_capt(int stm, int wdl, int n)
{
  if (!sub_cnt[stm ^ 1][2 - wdl])
    return;

  char capt_name[32], sub_name[32];
  strcat(strcpy(capt_name, "capt/"), wdl_name[2 + wdl]);
  strcat(strcpy(sub_name , "sub/" ), wdl_name[2 - wdl]);

  char str[64];
  sprintf(str, "%s/%c/done", capt_name, "wb"[stm]);
  FILE *F = fopen(str, "rb");
  if (F) {
    file_read(&g_stats[stm][n], 8, F);
    fclose(F);
    return;
  }

  char phase[64];
  snprintf(phase, sizeof phase, "calculating %s captures for %s",
      wdl_name[2 + wdl], side[stm]);

  bool partial = dir_exists(-1, stm, capt_name), done = partial;

  create_dir(-1, stm, capt_name);

  struct KSliceIterator iter;
  uint64_t num, cnt = 0;
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
        predecessors_sub(stm, s, false);
      }
    }

    while (kslice_iter_out(&iter, &s)) {
      if (!partial || !kslice_test_count(s, stm, capt_name, -1, &num)) {
        read_wins(-1, s, stm, -1); // most recent "wins".
        kslice_andnot(s, -1); // Remove illegal positions from capt/win.
        cnt += num = kslice_count(s);
        kslice_write(s, s, stm, capt_name, -1, num);
      }
    }
  }

  g_stats[stm][n] = cnt;

  F = file_open_write(str);
  file_write(&g_stats[stm][n], 8, F);
  fclose(F);
  file_rename(str);

  snprintf(phase, sizeof phase, "%s %s captures: %lu", side[stm],
      wdl_name[2 + wdl], cnt);
  show_progress(phase, 462, 462, true);
}

static void calc_capt_bloss(int stm)
{
  if (!sub_cnt[stm ^ 1][3])
    return;

  char str[64];
  sprintf(str, "capt/bloss/%c/done", "wb"[stm]);
  if (file_exists(str))
    return;

  bool partial = dir_exists(-1, stm, "capt/bloss"), done = partial;

  create_dir(-1, stm, "capt/bloss");

  struct KSliceIterator iter;
  uint64_t dummy;
  int num_done = 0;

  char phase[64];
  snprintf(phase, sizeof phase, "calculating bloss captures for %s", side[stm]);

  kslice_iter_init(&iter, stm);
  int s, s1;
  while (kslice_iter_next(&iter, &s)) {
    show_progress(phase, num_done++, 462, false);

    while (kslice_iter_in(&iter, &s1))
      if (!partial || !kslice_test_count(s1, stm, "capt/bloss", -1, &dummy)) {
        done = false;
        kslice_clear(s1);
      }

    if (!done) {
      // Subtract sub_win from sub_cwin.
      kslice_sub_read(-1, s, stm ^ 1, "sub/win");
      kslice_sub_read(s, s, stm ^ 1, "sub/cwin");
      kslice_sub_andnot(s, -1, stm ^ 1);
      predecessors_sub(stm, s, false);
    }

    while (kslice_iter_out(&iter, &s))
      if (!partial || !kslice_test_count(s, stm, "capt/bloss", -1, &dummy)) {
        uint64_t num = kslice_count(s);
        if (num)
          kslice_write(s, s, stm, "capt/bloss", -1, num);
      }
  }

  create_empty(str);

  snprintf(phase, sizeof phase, "%s bloss captures done", side[stm]);
  show_progress(phase, 462, 462, true);
}

static void predecessors_worker(struct ThreadData *thread)
{
  struct IdxState2 is;
  Position pos = g_pos;
  int stm = pos.stm;
  int s = work_slice;

  uint64_t *restrict p = (uint64_t *)kslice_get_address(-1);
  uint8_t *restrict const q = kslice_get_address(s);

  p += thread->begin >> 6;
  uint64_t last = thread->begin;
  idx_state2_init(&is, last, pos.sq, &ri);
  for (uint64_t idx = last, end = thread->end; idx < end; idx += 64) {
    uint64_t w = *p++;
    while (w) {
      uint64_t cur = idx + pop_lsb(&w);
      Bitboard occ = idx_state2_add(&is, cur - last, &ri);
      last = cur;
      mark_king_unmoves(stm, occ, &is);
      mark_unmoves(stm, KNIGHT, q, occ, &is, false);
      mark_unmoves(stm, BISHOP, q, occ, &is, false);
      mark_unmoves(stm, ROOK  , q, occ, &is, false);
      mark_unmoves(stm, QUEEN , q, occ, &is, false);
    }
  }
}

static void predecessors_ref_worker(struct ThreadData *thread)
{
  struct IdxState2 is;
  Position pos = g_pos;
  int stm = pos.stm;
  int s = work_slice;

  is.sq[0] = pos.sq[0];
  is.sq[1] = pos.sq[1];
  is.occ[0] = is.bb[0] = bit(pos.sq[0]) | bit(pos.sq[1]);

  uint64_t *restrict p = (uint64_t *)kslice_get_address(-1);
  uint8_t *restrict const q = kslice_get_address(s);

  p += thread->begin >> 6;
  uint64_t last = thread->begin;
  for (uint64_t idx = last, end = thread->end; idx < end; idx += 64) {
    uint64_t w = *p++;
    while (w) {
      uint64_t cur = idx + pop_lsb(&w);
      Bitboard occ = unrank_bb_ref(cur, is.bb, &ri);
      mark_king_unmoves(stm, occ, &is);
      mark_unmoves(stm, KNIGHT, q, occ, &is, true);
      mark_unmoves(stm, BISHOP, q, occ, &is, true);
      mark_unmoves(stm, ROOK  , q, occ, &is, true);
      mark_unmoves(stm, QUEEN , q, occ, &is, true);
    }
  }
}

// Calculate stm predecessors of stm^1 positions in scratch.
static void predecessors(int stm, int s)
{
  work_slice = s;
  g_pos.stm = stm;
  g_pos.sq[0] = KKSquare[s][0];
  g_pos.sq[1] = KKSquare[s][1];

  if (s < 441)
    run_threaded(predecessors_worker, &work_g_dynamic[0]);
  else
    run_threaded(predecessors_ref_worker, &work_g_dynamic[1]);
}

static void check_successors_worker(struct ThreadData *thread)
{
  struct IdxState2 is;
  Position pos = g_pos;
  int stm = pos.stm;
  int s = work_slice;
  uint64_t cnt = 0;

  uint64_t *restrict p = (uint64_t *)kslice_get_address(-1);
  uint8_t *restrict const q = kslice_get_address(s);

  p += thread->begin >> 6;
  uint64_t last = thread->begin;
  idx_state2_init(&is, last, pos.sq, &ri);
  for (uint64_t idx = last, end = thread->end; idx < end; idx += 64, p++) {
    uint64_t w = *p;
    if (!w) continue;
    uint64_t w2 = w;
    while (w) {
      unsigned bt = pop_lsb(&w);
      uint64_t cur = idx + bt;
      Bitboard occ = idx_state2_add(&is, cur - last, &ri);
      last = cur;
      // Legality check not necessary if we already removed illegal positions.
      // Currently, we need to test.
      if (   !idx_state2_legal(&is, stm, occ)
          || !check_moves(stm, QUEEN , s, q, occ, &is, false)
          || !check_moves(stm, ROOK  , s, q, occ, &is, false)
          || !check_moves(stm, BISHOP, s, q, occ, &is, false)
          || !check_moves(stm, KNIGHT, s, q, occ, &is, false)
          || !check_king_moves(stm, occ, &is))
        w2 ^= bit(bt);
    }
    *p = w2;
    cnt += popcnt(w2);
  }

  thread->cnt += cnt;
}

static void check_successors_ref_worker(struct ThreadData *thread)
{
  struct IdxState2 is;
  Position pos = g_pos;
  int stm = pos.stm;
  int s = work_slice;
  uint64_t cnt = 0;

  is.sq[0] = pos.sq[0];
  is.sq[1] = pos.sq[1];
  is.occ[0] = is.bb[0] = bit(pos.sq[0]) | bit(pos.sq[1]);

  uint64_t *restrict p = (uint64_t *)kslice_get_address(-1);
  uint8_t *restrict const q = kslice_get_address(s);

  p += thread->begin >> 6;
  uint64_t last = thread->begin;
  for (uint64_t idx = last, end = thread->end; idx < end; idx += 64, p++) {
    uint64_t w = *p;
    if (!w) continue;
    uint64_t w2 = w;
    while (w) {
      unsigned bt = pop_lsb(&w);
      uint64_t cur = idx + bt;
      Bitboard occ = unrank_bb_ref(cur, is.bb, &ri);
      // Legality check not necessary if we already removed illegal positions.
      // Currently, we need to test.
      if (   !idx_state2_legal(&is, stm, occ)
          || !check_moves(stm, QUEEN , s, q, occ, &is, true)
          || !check_moves(stm, ROOK  , s, q, occ, &is, true)
          || !check_moves(stm, BISHOP, s, q, occ, &is, true)
          || !check_moves(stm, KNIGHT, s, q, occ, &is, true)
          || !check_king_moves(stm, occ, &is))
        w2 ^= bit(bt);
    }
    *p = w2;
    cnt += popcnt(w2);
  }

  thread->cnt += cnt;
}

// Verify stm positions as loss against stm^1 positions and return the
// number of verified losses.
static uint64_t check_successors(int stm, int s)
{
  work_slice = s;
  g_pos.stm = stm;
  g_pos.sq[0] = KKSquare[s][0];
  g_pos.sq[1] = KKSquare[s][1];

  for (int t = 0; t < g_num_threads; t++)
    g_thread_data[t].cnt = 0;

  if (s < 441)
    run_threaded(check_successors_worker, &work_g_dynamic[0]);
  else
    run_threaded(check_successors_ref_worker, &work_g_dynamic[1]);

  uint64_t cnt = 0;
  for (int t = 0; t < g_num_threads; t++)
    cnt += g_thread_data[t].cnt;

  return cnt;
}

static void calc_illegal_worker_tmpl(struct ThreadData *thread, const int pc,
    const bool ref)
{
  struct IdxState2 is;
  int stm = g_pos.stm;
  int k = g_piece_set[stm][pc];
  int ksq = g_pos.sq[stm ^ 1];

  uint8_t *restrict const p = kslice_buf[stm];

  Bitboard occ = idx_state2_init(&is, thread->begin, g_pos.sq, &capt_ri[k]);

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

static void calc_mate_worker(struct ThreadData *thread)
{
  struct IdxState2 is;
  Position pos = g_pos;

  uint64_t *restrict p0 = (uint64_t *)kslice_buf[0];
  uint64_t *restrict p1 = (uint64_t *)kslice_buf[1];
  uint64_t *restrict q0 = (uint64_t *)kslice_buf[2];
  uint64_t *restrict q1 = (uint64_t *)kslice_buf[3];

  uint64_t last = thread->begin;
  p0 += last >> 6;
  p1 += last >> 6;
  q0 += last >> 6;
  q1 += last >> 6;
  idx_state2_init(&is, last, pos.sq, &ri);
  for (uint64_t idx = last, end = thread->end; idx < end;
      idx += 64, p0++, p1++, q0++, q1++)
  {
    uint64_t w = *p0 ^ *p1;
    if (!w) continue;
    uint64_t white = 0, black = 0;
    while (w) {
      unsigned bt = pop_lsb(&w);
      uint64_t cur = idx + bt;
      Bitboard occ = idx_state2_add(&is, cur - last, &ri);
      last = cur;
      if (*p1 & bit(bt)) {
        if (idx_state2_mate(&is, WHITE, occ))
          white |= bit(bt);
      } else {
        if (idx_state2_mate(&is, BLACK, occ))
          black |= bit(bt);
      }
    }
    *q0 = white;
    *q1 = black;
  }
}

static void calc_mate_ref_worker(struct ThreadData *thread)
{
  struct IdxState2 is;
  Position pos = g_pos;

  is.sq[0] = pos.sq[0];
  is.sq[1] = pos.sq[1];
  is.occ[0] = is.bb[0] = bit(pos.sq[0]) | bit(pos.sq[1]);

  uint64_t *restrict p0 = (uint64_t *)kslice_buf[0];
  uint64_t *restrict p1 = (uint64_t *)kslice_buf[1];
  uint64_t *restrict q0 = (uint64_t *)kslice_buf[2];
  uint64_t *restrict q1 = (uint64_t *)kslice_buf[3];

  uint64_t last = thread->begin;
  p0 += last >> 6;
  p1 += last >> 6;
  q0 += last >> 6;
  q1 += last >> 6;
  for (uint64_t idx = last, end = thread->end; idx < end;
      idx += 64, p0++, p1++, q0++, q1++)
  {
    uint64_t w = *p0 ^ *p1;
    if (!w) continue;
    uint64_t white = 0, black = 0;
    while (w) {
      unsigned bt = pop_lsb(&w);
      uint64_t cur = idx + bt;
      Bitboard occ = unrank_bb_ref(cur, is.bb, &ri);
      if (*p1 & bit(bt)) {
        if (idx_state2_mate(&is, WHITE, occ))
          white |= bit(bt);
      } else {
        if (idx_state2_mate(&is, BLACK, occ))
          black |= bit(bt);
      }
    }
    *q0 = white;
    *q1 = black;
  }
}

// Calc illegal and mate (L0) positions.
static void calc_illegal_and_mate(void)
{
  FILE *F = fopen("0/done", "rb");
  if (F) {
    for (int stm = 0; stm < 2; stm++) {
      file_read(&g_stats[stm][0], 8, F);
      file_read(&g_stats[stm][MAX_STATS - 1], 8, F);
    }
    fclose(F);
    return;
  }

  bool partial = file_exists("0");

  for (int stm = 0; stm < 2; stm++) {
    create_dir(0, stm, "wins");
    create_dir(0, stm, "L");
  }

  uint64_t broken[2] = { 0 }, loss0[2] = { 0 }, num;

  for (int s = 0; s < 462; s++) {
    show_progress("calculating illegal and mate positions", s, 462, false);

    if (partial) {
      char str[64];
      create_name(str, s, BLACK, "L", 0);
      if (file_exists(str)) {
        for (int stm = 0; stm < 2; stm++) {
          if (kslice_test_count(s, stm, "wins", 0, &num))
            broken[stm] += num;
          if (kslice_test_count(s, stm, "L", 0, &num))
            loss0[stm] += num;
        }
        continue;
      }
    }

    g_pos.sq[0] = KKSquare[s][0];
    g_pos.sq[1] = KKSquare[s][1];

    kslice_clear_addr(kslice_buf[0], s); // wtm illegal
    kslice_clear_addr(kslice_buf[1], s); // btm illegal

    for (int stm = 0; stm < 2; stm++) {
      int k;
      g_pos.stm = stm;
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
    }

    for (int stm = 0; stm < 2; stm++) {
      broken[stm] += num = kslice_count_addr(kslice_buf[stm], s);
      kslice_write_addr(kslice_buf[stm], s, stm, "wins", 0, num);
    }

    kslice_clear_addr(kslice_buf[2], s); // wtm mate
    kslice_clear_addr(kslice_buf[3], s); // btm mate

    if (s < 441)
      run_threaded(calc_mate_worker, &work_g_static[0]);
    else
      run_threaded(calc_mate_ref_worker, &work_g_static[1]);

    for (int stm = 0; stm < 2; stm++) {
      loss0[stm] += num = kslice_count_addr(kslice_buf[2 + stm], s);
      kslice_write_addr(kslice_buf[2 + stm], s, stm, "L", 0, num);
    }
  }

  F = file_open_write("0/done");
  for (int stm = 0; stm < 2; stm++) {
    g_stats[stm][0] = broken[stm];
    g_stats[stm][MAX_STATS - 1] = loss0[stm];
    file_write(&g_stats[stm][0], 8, F);
    file_write(&g_stats[stm][MAX_STATS - 1], 8, F);
  }
  fclose(F);
  file_rename("0/done");

  char phase[128];
  snprintf(phase, sizeof phase, "illegal: %lu, %lu; mate: %lu, %lu",
      broken[WHITE], broken[BLACK], loss0[WHITE], loss0[BLACK]);
  show_progress(phase, 462, 462, true);
}

// Calculate stm losses in n from stm^1 wins in n-1 (n > 1) or
// from stm^1 wins in subtables reached through captures (n == 1).
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

  struct KSliceIterator iter;
  bool partial = true;
  int num_done = 0;

  if (dir_exists(n, stm, "L"))
    goto skip_X;

  partial = dir_exists(n, stm, "X");
  bool done = partial;
  uint64_t num;

  create_dir(n, stm, "X");

  snprintf(phase, sizeof phase, "%d/X/%c", n, "wb"[stm]);

  // Calculate potential losses in n = predecessors(W(n-1))
  kslice_iter_init(&iter, stm);
  int s, s1;
  while (kslice_iter_next(&iter, &s)) {
    show_progress(phase, num_done++, 462, false);

    bool pred_sub =   (n == 1 && sub_cnt[stm ^ 1][4])
                   || (n == DRAW_RULE + 1 && sub_cnt[stm ^ 1][3]);
    bool pred = more_l && kslice_test(s, stm ^ 1, "W", n - 1);

    if (pred_sub || pred) {
      while (kslice_iter_in(&iter, &s1))
        if (!partial || !kslice_test_count(s1, stm, "X", n, &num)) {
          done = false;
          kslice_clear(s1);
        }

      if (pred_sub && !done) {
        if (n == 1) {
          kslice_sub_read(s, s, stm ^ 1, "sub/win");
          predecessors_sub(stm, s, true);
        } else {
          // We must subtract sub/win from sub/cwin here.
          kslice_sub_read(-1, s, stm ^ 1, "sub/win");
          kslice_sub_read(s, s, stm ^ 1, "sub/cwin");
          kslice_sub_andnot(s, -1, stm ^ 1);
          predecessors_sub(stm, s, false);
        }
      }

      if (pred && !done) {
        kslice_read(-1, s, stm ^ 1, "W", n - 1);
        predecessors(stm, s);
      }
    }

#if 0
      // If there are many predecessors, it might be more efficient to
      // remove illegal positions and positions with non-losing captures
      // here.
      kslice_read(-1, s, stm, n <= DRAW_RULE ? "noloss" : "nobloss", -1);
#endif
    while (kslice_iter_out(&iter, &s))
      if (!partial || !kslice_test_count(s, stm, "X", n, &num)) {
        kslice_clear_tail(s, s);
        kslice_write(s, s, stm, "X", n, UINT64_MAX);
      }
  }
  show_progress(phase, 462, 462, true);

  create_dir(n, stm, "L");
  partial = false;

skip_X:

  uint64_t cnt = 0;
  num_done = 0;
  snprintf(phase, sizeof phase, "\x1b[%sm%d/L/%c\x1b[0m",
      clr_L[(2 * n + stm) & 3], n, "wb"[stm]);

  // Verify potential losses.
  kslice_iter_init(&iter, stm);
  while (kslice_iter_next(&iter, &s)) {
    show_progress(phase, num_done++, 462, false);

    if (kslice_test(s, stm, "X", n)) {
      while (kslice_iter_in(&iter, &s1)) {
        read_wins(s1, s1, stm ^ 1, n - 1);
        // If there are very few predecessors, it might be more efficient to
        // directly probe_wdl() their captures.
        kslice_sub_read(s1, s1, stm ^ 1,
            n <= DRAW_RULE || !sub_cnt[stm ^ 1][3] ? "sub/win" : "sub/cwin");
      }

      kslice_read(-1, s, stm, "X", n);
      cnt += num = check_successors(stm, s);
      kslice_write(-1, s, stm, "L", n, num);
    } else if (partial && kslice_test_count(s, stm, "L", n, &num))
      cnt += num;

    while (kslice_iter_out(&iter, &s))
      kslice_delete(s, stm, "X", n);
  }

  g_stats[stm][MAX_STATS - 1 - n] = cnt;

  F = file_open_write(str);
  file_write(&g_stats[stm][MAX_STATS - 1 - n], 8, F);
  fclose(F);
  file_rename(str);

  snprintf(phase, sizeof phase, "\x1b[%sm%d/L/%c  %lu\x1b[0m",
      clr_L[(2 * n + stm) & 3], n, "wb"[stm], cnt);
  show_progress(phase, 462, 462, true);

  return cnt != 0;
}

static bool calc_W(int stm, int n, bool more_w)
{
  char str[64];
  sprintf(str, "%d/W/%c/done", n, "wb"[stm]);
  FILE *F = fopen(str, "rb");
  if (F) {
    file_read(&g_stats[stm][stats_n(n)], 8, F);
    fclose(F);
    return g_stats[stm][stats_n(n)] != 0;
  }

  char phase[64];

  struct KSliceIterator iter;
  uint64_t cnt = 0, num;
  int num_done = 0;

  bool partial = dir_exists(n, stm, "W"), done = partial;

  create_dir(n, stm, "W");
  create_dir(n, stm, "wins");

  snprintf(phase, sizeof phase, "\x1b[%sm%d/W/%c\x1b[0m",
      clr_W[(2 * n + stm) & 3], n, "wb"[stm]);

  // Calculate wins in n = predecessors(L(n-1))
  kslice_iter_init(&iter, stm);
  int s, s1;
  while (kslice_iter_next(&iter, &s)) {
    show_progress(phase, num_done++, 462, false);

    bool pred_sub =   (n == 1             && sub_cnt[stm ^ 1][0])
                   || (n == DRAW_RULE + 1 && sub_cnt[stm ^ 1][1]);
    bool pred = more_w && kslice_test(s, stm ^ 1, "L", n - 1);

    if (pred_sub || pred) {
      while (kslice_iter_in(&iter, &s1))
        if (partial && kslice_test_count(s1, stm, "W", n, &num)) {
          cnt += num;
        } else {
          done = false;
          kslice_clear(s1);
        }

      if (n == DRAW_RULE + 1 && sub_cnt[stm ^ 1][1] && !done) {
        kslice_sub_read(s, s, stm ^ 1, "sub/bloss");
        predecessors_sub(stm, s, false);
      }

      if (pred && !done) {
        kslice_read(-1, s, stm ^ 1, "L", n - 1);
        predecessors(stm, s);
      }
    }

    while (kslice_iter_out(&iter, &s)) {
      if (!partial || !kslice_test_count(s, stm, "W", n, &num)) {
        // We should probably do the same for n == DRAW_RULE + 1.
        if (n == 1 && sub_cnt[stm ^ 1][0]) {
          kslice_read(-1, s, stm, "capt/win", -1);
          kslice_or(s, -1);
        }
        // Remove illegal positions and known faster wins.
        read_wins(-1, s, stm, n - 1);
        kslice_andnot(s, -1);
        num = kslice_count(s);
        uint64_t cost = kslice_write(s, s, stm, "W", n, num);
        add_wins(s, stm, n, cost);
        cnt += num;
      }
    }
  }

  g_stats[stm][stats_n(n)] = cnt;

  F = file_open_write(str);
  file_write(&g_stats[stm][stats_n(n)], 8, F);
  fclose(F);
  file_rename(str);

  snprintf(phase, sizeof phase, "\x1b[%sm%d/W/%c  %lu\x1b[0m",
      clr_W[(2 * n + stm) & 3], n, "wb"[stm], cnt);
  show_progress(phase, 462, 462, true);

  return cnt != 0;
}

void generate(void)
{
  FILE *F = fopen("generate_info", "rb");
  if (F) {
    read_data(F, &g_stats, sizeof g_stats);
    file_read(&sub_cnt, sizeof sub_cnt, F);
    file_read(&max_iteration, sizeof max_iteration, F);
    fclose(F);
    printf("Skipped generation phase.\n");
    return;
  }

  // Calculate kslices for positions reached through a capture.
  calc_sub_kslices(WHITE);
  calc_sub_kslices(BLACK);

  calc_capt_bloss(WHITE);
  calc_capt_bloss(BLACK);

  memset(&wins_full, 0, sizeof wins_full);
  memset(&wins_checked, 0, sizeof wins_checked);
  memset(&replay_cost, 0, sizeof replay_cost);

  calc_illegal_and_mate();

  // CAPT_WIN
  calc_capt(WHITE, 2, 1);
  calc_capt(BLACK, 2, 1);

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

  // CAPT_CWIN
  calc_capt(WHITE, 1, 3 + DRAW_RULE);
  calc_capt(BLACK, 1, 3 + DRAW_RULE);

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

  // CAPT_DRAW
  calc_capt(WHITE, 0, MAX_STATS / 2 + 1);
  calc_capt(BLACK, 0, MAX_STATS / 2 + 1);

  max_iteration = n;

  // Remove some double counting.
  for (int stm = 0; stm < 2; stm++) {
    g_stats[stm][3] -= g_stats[stm][1];
    g_stats[stm][4 + DRAW_RULE] -= g_stats[stm][3 + DRAW_RULE];
    uint64_t tot = 0;
    for (int i = 0; i < MAX_STATS; i++)
      tot += g_stats[stm][i];
    g_stats[stm][MAX_STATS / 2 + 2] = 462 * kslice_size - tot;
  }

  F = file_open_write("generate_info");
  write_data(F, &g_stats, sizeof g_stats);
  file_write(&sub_cnt, sizeof sub_cnt, F);
  file_write(&max_iteration, sizeof max_iteration, F);
  fclose(F);
  file_rename("generate_info");
}

void delete_intermediate_slices(void)
{
  if (file_exists("0/wins")) {
    for (int stm = 0; stm < 2; stm++)
      for (int s = 0; s < 462; s++)
        kslice_delete(s, stm, "wins", wins_full[stm][s]);
    for (int n = max_iteration - 1; n >= 0; n--)
      delete_dir(n, "wins");
  }

  if (file_exists("sub")) {
    for (int i = 0; i < 5; i++) {
      char name[64];
      sprintf(name, "sub/%s", wdl_name[i]);
      for (int stm = 0; stm < 2; stm++)
        for (int s = 0; s < 462; s++)
          kslice_delete(s, stm, name, -1);
      delete_dir(-1, name);
    }
    unlink("sub/done_w");
    unlink("sub/done_b");
    rmdir("sub");
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
  rmdir("capt/bloss/done");
  for (int i = 0; i < 5; i++) {
    sprintf(dir, "capt/%s", wdl_name[i]);
    delete_dir(-1, dir);
  }
  delete_dir(-1, "capt");
}
