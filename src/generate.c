/*
  Copyright (c) 2026 Ronald de Man

  This file is distributed under the terms of the GNU GPL, version 2.
*/

#include <inttypes.h>
#include <stdint.h>

#include "defs.h"
#include "index.h"
#include "kslice.h"
#include "movegen.h"
#include "probe.h"
#include "tb8gen.h"
#include "threads.h"
#include "types.h"

uint64_t sub_cnt[2][5];
int max_iteration;

INLINE void mark_king_unmoves(int stm, Bitboard occ, uint8_t *restrict sq,
    int s)
{
  uint8_t sq2[MAX_PIECES];
  Bitboard b = king_attacks(sq[stm]) & ~(occ | king_attacks(sq[stm ^ 1]));
  while (b) {
    sq[stm] = pop_lsb(&b);
    normalize(sq, sq2);
    int s2 = KKMap[sq2[0]][sq2[1]];
    uint8_t *p = kslice_get_address(s2);
    kslice_bit_set_atomic(p, sq_to_idx(sq2));
    if (s < 441 && s2 >= 441) {
      mirror_diagonal(sq2);
      kslice_bit_set_atomic(p, sq_to_idx(sq2));
    }
  }
}

INLINE void mark_unmoves(int k, uint8_t *restrict const p, Bitboard occ,
    uint8_t *restrict sq)
{
  uint8_t sq2[MAX_PIECES];
  Bitboard b = non_king_piece_moves(g_pos.pt[k], sq[k], occ);
  while (b) {
    sq[k] = pop_lsb(&b);
    for (int i = 0; i < MAX_PIECES; i++)
      sq2[i] = sq[i];
    uint64_t idx = sq_to_idx(sq2);
    kslice_bit_set_atomic(p, idx);
  }
}

static int work_slice, work_set;

static void calc_sub_worker(struct ThreadData *thread)
{
  Position pos = g_pos;
  uint32_t sub[MAX_SETS];
  int k = work_set;
  int m = ii.last[k];
  int n = --pos.num;
  uint64_t cnt = 0;

  pos.pt[m] = pos.pt[n];
  uint8_t *restrict p[5];
  for (int i = 0; i < 5; i++) {
    p[i] = kslice_sub_buf[i] + sub_offset[k];
    memset(p[i] + ((thread->begin + 7) >> 3), 0x00,
        (thread->end - thread->begin + 7) >> 3);
  }

  idx_to_sq_init(thread->begin, sub, &capt_ii[k]);
  for (uint64_t idx = thread->begin, end = thread->end; idx < end;
      idx++, idx_to_sq_inc(sub, &capt_ii[k]))
  {
    pos.occ = capt_idx_to_sq(sub, pos.sq, k);
    pos.sq[m] = pos.sq[n];
    if (opp_king_attacked(&pos)) {
      // Include illegal positions in sub_cwin and sub_win.
      kslice_bit_set(p[3], idx);
      kslice_bit_set(p[4], idx);
      cnt++;
    } else {
      int v = probe_wdl(&pos, -2, 2);
      kslice_bit_set(p[v + 2], idx);
      // add sub_win to sub_cwin
      if (v == 2)
        kslice_bit_set(p[3], idx);
    }
  }

  thread->cnt += cnt;
}

// Calculate aggregate bitmaps for subtables, one per loss/bloss/draw/cwin/win.
static void calc_sub_kslices(int stm)
{
  g_pos.stm = stm;


  for (int s = 0; s < 462; s++) {
    for (int t = 0; t < g_num_threads; t++)
      g_thread_data[t].cnt = 0;

    g_pos.sq[0] = KKSquare[s][0];
    g_pos.sq[1] = KKSquare[s][1];

    for (int k = 0; k < ii.numsets; k++) {
      if ((g_pos.pt[ii.first[k]] >> 3) != stm)
        continue;
      work_set = k;
      run_threaded(calc_sub_worker, work_capt[k], 0);
    }

    uint64_t cnt_ilgl = 0, cnt_l, cnt_bl, cnt_d, cnt_cw, cnt_w;

    for (int t = 0; t < g_num_threads; t++)
      cnt_ilgl += g_thread_data[t].cnt;

    cnt_l = kslice_sub_count_addr(kslice_sub_buf[0], stm);
    if (cnt_l) {
      create_dir(-1, stm, "sub/loss");
      kslice_sub_write_addr(kslice_sub_buf[0], s, stm, "sub/loss", cnt_l);
    }

    cnt_bl = kslice_sub_count_addr(kslice_sub_buf[1], stm);
    if (cnt_bl) {
      create_dir(-1, stm, "sub/bloss");
      kslice_sub_write_addr(kslice_sub_buf[1], s, stm, "sub/bloss", cnt_bl);
    }

    cnt_d = kslice_sub_count_addr(kslice_sub_buf[2], stm);
    if (cnt_d) {
      create_dir(-1, stm, "sub/draw");
      kslice_sub_write_addr(kslice_sub_buf[2], s, stm, "sub/draw", cnt_d);
    }

    cnt_cw = kslice_sub_count_addr(kslice_sub_buf[3], stm);
    if (cnt_cw) {
      create_dir(-1, stm, "sub/cwin");
      kslice_sub_write_addr(kslice_sub_buf[3], s, stm, "sub/cwin", cnt_cw);
    }

    cnt_w = kslice_sub_count_addr(kslice_sub_buf[4], stm);
    create_dir(-1, stm, "sub/win");
    kslice_sub_write_addr(kslice_sub_buf[4], s, stm, "sub/win", cnt_w);

    sub_cnt[stm][0] += cnt_l;
    sub_cnt[stm][1] += cnt_bl;
    sub_cnt[stm][2] += cnt_d;
    sub_cnt[stm][3] += cnt_cw - cnt_w;
    sub_cnt[stm][4] += cnt_w - cnt_ilgl;
  }
}

static bool work_legality;

static void predecessors_sub_worker(struct ThreadData *thread)
{
  Position pos = g_pos;
  uint32_t sub[MAX_SETS];
  int stm = pos.stm;
  pos.stm ^= 1;
  int n = --pos.num;
  int k = work_set;
  int s = work_slice;
  bool legality = work_legality;

  int m = ii.last[k];
  pos.pt[m] = pos.pt[n];

  uint64_t *restrict p = (uint64_t *)kslice_sub_get_address(s, k);
  uint8_t *restrict const q = kslice_get_address(s);

  uint64_t last = thread->begin;
  p += thread->begin >> 6;
  idx_to_sq_init(last, sub, &capt_ii[k]);
  for (uint64_t idx = last, end = thread->end; idx < end; idx += 64) {
    uint64_t w = *p++;
    while (w) {
      unsigned bt = pop_lsb(&w);
      if (idx + bt >= end) break;
      idx_to_sq_add(idx + bt - last, sub, &capt_ii[k]);
      last = idx + bt;
      Bitboard occ = capt_idx_to_sq(sub, pos.sq, k);
      if (legality) {
        pos.occ = occ;
        pos.sq[m] = pos.sq[n];
        if (opp_king_attacked(&pos))
          continue;
      }
      // Uncapture by king.
      pos.sq[m] = pos.sq[stm];
      mark_king_unmoves(stm, occ, pos.sq, s);
      pos.sq[stm] = pos.sq[m];
      // Uncapture by non-king pieces.
      for (int i = 1; pos.pcs[stm][i] >= 0; i++) {
        int j = pos.pcs[stm][i];
        pos.sq[m] = pos.sq[j];
        mark_unmoves(j, q, occ, pos.sq);
        pos.sq[j] = pos.sq[m];
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
  for (int k = 0; k < ii.numsets; k++) {
    int m = ii.last[k];
    if ((g_pos.pt[m] >> 3) == stm)
      continue;
    work_set = k;
    run_threaded(predecessors_sub_worker, work_capt[k], 0);
  }
}

static void calc_capt(int stm)
{
  struct KSliceIterator iter;

  if (sub_cnt[stm ^ 1][0]) {
    create_dir(-1, stm, "capt/win");

    kslice_iter_init(&iter, stm);
    int s, s1;
    while (kslice_iter_next(&iter, &s)) {

      while (kslice_iter_in(&iter, &s1))
        kslice_clear(s1);

      kslice_sub_read(s, s, stm ^ 1, "sub/loss");
      predecessors_sub(stm, s, false);

      while (kslice_iter_out(&iter, &s))
        kslice_write(s, s, stm, "capt/win", -1, UINT64_MAX);
    }
  }

  if (sub_cnt[stm ^ 1][1]) {
    create_dir(-1, stm, "capt/cwin");

    kslice_iter_init(&iter, stm);
    int s, s1;
    while (kslice_iter_next(&iter, &s)) {

      while (kslice_iter_in(&iter, &s1))
        kslice_clear(s1);

      kslice_sub_read(s, s, stm ^ 1, "sub/bloss");
      predecessors_sub(stm, s, false);

      while (kslice_iter_out(&iter, &s))
        kslice_write(s, s, stm, "capt/cwin", -1, UINT64_MAX);
    }
  }

  if (sub_cnt[stm ^ 1][2]) {
    create_dir(-1, stm, "capt/draw");

    kslice_iter_init(&iter, stm);
    int s, s1;
    while (kslice_iter_next(&iter, &s)) {

      while (kslice_iter_in(&iter, &s1))
        kslice_clear(s1);

      kslice_sub_read(s, s, stm ^ 1, "sub/draw");
      predecessors_sub(stm, s, false);

      while (kslice_iter_out(&iter, &s))
        kslice_write(s, s, stm, "capt/draw", -1, UINT64_MAX);
    }
  }

  if (sub_cnt[stm ^ 1][3]) {
    create_dir(-1, stm, "capt/bloss");

    kslice_iter_init(&iter, stm);
    int s, s1;
    while (kslice_iter_next(&iter, &s)) {

      while (kslice_iter_in(&iter, &s1))
        kslice_clear(s1);

      kslice_sub_read(s, s, stm ^ 1, "sub/cwin");
      predecessors_sub(stm, s, true);

      while (kslice_iter_out(&iter, &s))
        kslice_write(s, s, stm, "capt/bloss", -1, UINT64_MAX);
    }
  }
}

static void predecessors_worker(struct ThreadData *thread)
{
  uint32_t sub[MAX_SETS];
  Position pos = g_pos;
  int stm = pos.stm;
  int s = work_slice;

  uint64_t *restrict p = (uint64_t *)kslice_get_address(-1);
  uint8_t *restrict const q = kslice_get_address(s);

  p += thread->begin >> 6;
  uint64_t last = thread->begin;
  idx_to_sq_init(last, sub, &ii);
  for (uint64_t idx = last, end = thread->end; idx < end; idx += 64) {
    uint64_t w = *p++;
    while (w) {
      unsigned bt = pop_lsb(&w);
      idx_to_sq_add(idx + bt - last, sub, &ii);
      last = idx + bt;
      if (last >= end) break;  // we can remove this check later if safe
      Bitboard occ = idx_to_sq(sub, pos.sq);
      uint8_t tmp = pos.sq[stm];
      mark_king_unmoves(stm, occ, pos.sq, s);
      pos.sq[stm] = tmp;
      for (int i = 1; pos.pcs[stm][i] >= 0; i++) {
        int j = pos.pcs[stm][i];
        uint8_t tmp = pos.sq[j];
        mark_unmoves(j, q, occ, pos.sq);
        pos.sq[j] = tmp;
      }
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

  run_threaded(predecessors_worker, work_g, 0);
}

INLINE int get_idx(uint8_t *restrict sq, int s)
{
  for (int i = 0; ; i++)
    if (sq[i] == s)
      return i;
  unreachable();
}

INLINE bool check_king_moves(int stm, Bitboard occ, uint8_t *restrict sq)
{
  uint8_t sq2[MAX_PIECES];

  Bitboard b = king_attacks(sq[stm]) & ~king_attacks(sq[stm ^ 1]);
#if 1
  Bitboard attacks = b & occ;
  while (attacks) {
    int to = pop_lsb(&attacks);
    int j = get_idx(sq, to);
    if ((g_pos.pt[j] >> 3) == stm) continue;
    sq[stm] = to;
    normalize(sq, sq2);
    int l = pc_to_set[j];
    sq2[j] = sq2[ii.last[l]];
    int s2 = KKMap[sq2[0]][sq2[1]];
    uint8_t *p = kslice_sub_get_address(s2, l);
    if (!kslice_bit_test(p, capt_sq_to_idx(sq2, l)))
      return false;
  }
#endif

  b &= ~occ;
  while (b) {
    sq[stm] = pop_lsb(&b);
    normalize(sq, sq2);
    int s2 = KKMap[sq2[0]][sq2[1]];
    uint8_t *p = kslice_get_address(s2);
    if (!kslice_bit_test(p, sq_to_idx(sq2)))
      return false;
  } 

  return true;
}

INLINE bool check_moves(int k, int s, uint8_t *restrict const p, Bitboard occ,
    uint8_t *restrict sq)
{
  uint8_t sq2[MAX_PIECES];

  Bitboard b = non_king_piece_attacks(g_pos.pt[k], sq[k], occ);

#if 1
  Bitboard attacks = b & occ;
  while (attacks) {
    int to = pop_lsb(&attacks);
    int j = get_idx(sq, to);
    if (!((g_pos.pt[k] ^ g_pos.pt[j]) & 8)) continue;
    for (int i = 0; i < MAX_PIECES; i++)
      sq2[i] = sq[i];
    int l = pc_to_set[j];
    sq2[k] = to;
    sq2[j] = sq2[ii.last[l]];
    uint8_t *restrict q = kslice_sub_get_address(s, l);
    if (!kslice_bit_test(q, capt_sq_to_idx(sq2, l)))
      return false;
  }
#endif

  b &= ~occ;
  while (b) {
    sq[k] = pop_lsb(&b);
    for (int i = 0; i < MAX_PIECES; i++)
      sq2[i] = sq[i];
    if (!kslice_bit_test(p, sq_to_idx(sq2)))
      return false;
  }

  return true;
}

static void check_successors_worker(struct ThreadData *thread)
{
  uint32_t sub[MAX_SETS];
  Position pos = g_pos;
  int stm = pos.stm;
  int s = work_slice;
  uint64_t cnt = 0;

  uint64_t *restrict p = (uint64_t *)kslice_get_address(-1);
  uint8_t *restrict const q = kslice_get_address(s);

  p += thread->begin >> 6;
  uint64_t last = thread->begin;
  idx_to_sq_init(last, sub, &ii);
  for (uint64_t idx = last, end = thread->end; idx < end; idx += 64, p++) {
    uint64_t w = *p;
    if (!w) continue;
    uint64_t w2 = w;
    while (w) {
      unsigned bt = pop_lsb(&w);
      idx_to_sq_add(idx + bt - last, sub, &ii);
      last = idx + bt;
      if (last >= end) break; // we can remove this check later if safe
      Bitboard occ = pos.occ = idx_to_sq(sub, pos.sq);
      // Legality check not necessary if we already remove illegal positions.
      // Currently, we need to test.
      if (opp_king_attacked(&pos))
        goto clear_bit;
      for (int i = 1; pos.pcs[stm][i] >= 0; i++) {
        int j = pos.pcs[stm][i];
        uint8_t tmp = pos.sq[j];
        bool v = check_moves(j, s, q, occ, pos.sq);
        pos.sq[j] = tmp;
        if (!v) goto clear_bit;
      }
      uint8_t tmp = pos.sq[stm];
      bool v = check_king_moves(stm, occ, pos.sq);
      pos.sq[stm] = tmp;
      if (v) continue;
clear_bit:
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

  run_threaded(check_successors_worker, work_g, 0);

  uint64_t cnt = 0;
  for (int t = 0; t < g_num_threads; t++)
    cnt += g_thread_data[t].cnt;

  return cnt;
}

static void calc_illegal_worker(struct ThreadData *thread)
{
  uint32_t sub[MAX_SETS];
  Position pos = g_pos;
  int k = work_set;
  int m = ii.last[k];
  int stm = g_pos.pt[m] >> 3;
  int king_sq = pos.sq[stm ^ 1];

  uint8_t *restrict const p = kslice_buf[stm];

  idx_to_sq_init(thread->begin, sub, &capt_ii[k]);

  for (uint64_t idx = thread->begin, end = thread->end; idx < end;
      idx++, idx_to_sq_inc(sub, &capt_ii[k]))
  {
    Bitboard occ = capt_idx_to_sq(sub, pos.sq, k);
    pos.sq[m] = king_sq;
    mark_unmoves(m, p, occ, pos.sq);
  }
}

static void calc_mate_worker(struct ThreadData *thread)
{
  uint32_t sub[MAX_SETS];
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
  idx_to_sq_init(last, sub, &ii);
  for (uint64_t idx = last, end = thread->end; idx < end;
      idx += 64, p0++, p1++, q0++, q1++)
  {
    uint64_t w = *p0 ^ *p1;
    if (!w) continue;
    uint64_t white = 0, black = 0;
    while (w) {
      unsigned bt = pop_lsb(&w);
      idx_to_sq_add(idx + bt - last, sub, &ii);
      last = idx + bt;
      if (last >= end) break;
      pos.occ = idx_to_sq(sub, pos.sq);
      if (*p1 & bit(bt)) {
        pos.stm = WHITE;
        if (!has_legal_moves(&pos) && !has_legal_caps(&pos))
          white |= bit(bt);
      } else {
        pos.stm = BLACK;
        if (!has_legal_moves(&pos) && !has_legal_caps(&pos))
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
  uint64_t broken_w = 0, broken_b = 0, loss0_w = 0, loss0_b = 0, num;

  for (int s = 0; s < 462; s++) {
    g_pos.sq[0] = KKSquare[s][0];
    g_pos.sq[1] = KKSquare[s][1];

    kslice_clear_addr(kslice_buf[0]); // wtm illegal
    kslice_clear_addr(kslice_buf[1]); // btm illegal

    for (int k = 0; k < ii.numsets; k++) {
      work_set = k;
      run_threaded(calc_illegal_worker, work_capt[k], 0);
    }

    broken_w += num = kslice_count_addr(kslice_buf[0]);
    kslice_write_addr(kslice_buf[0], s, WHITE, "wins", 0, num);

    broken_b += num = kslice_count_addr(kslice_buf[1]);
    kslice_write_addr(kslice_buf[1], s, BLACK, "wins", 0, num);

    kslice_clear_addr(kslice_buf[2]); // wtm mate
    kslice_clear_addr(kslice_buf[3]); // btm mate

    run_threaded(calc_mate_worker, work_g, 0);

    loss0_w += num = kslice_count_addr(kslice_buf[2]);
    kslice_write_addr(kslice_buf[2], s, WHITE, "L", 0, num);

    loss0_b += num = kslice_count_addr(kslice_buf[3]);
    kslice_write_addr(kslice_buf[3], s, BLACK, "L", 0, num);
  }

  printf("broken_w = %lu\n", broken_w);
  printf("broken_b = %lu\n", broken_b);
  printf("l0_w = %lu\n", loss0_w);
  printf("l0_b = %lu\n", loss0_b);
}

// Calculate stm losses in n from stm^1 wins in n-1 (n > 1) or
// from stm^1 wins in sub tables reached through captures (n == 1).
static bool calc_L(int stm, int n, bool more_l)
{
  struct KSliceIterator iter;
  uint64_t cnt = 0, num;

  create_dir(n, stm, "X");

  // Calculate potential losses in n = predecessors(W(n-1))
  kslice_iter_init(&iter, stm);
  int s, s1;
  while (kslice_iter_next(&iter, &s)) {
    bool pred_sub =   (n == 1 && sub_cnt[stm ^ 1][4])
                   || (n == DRAW_RULE + 1 && sub_cnt[stm ^ 1][3]);
    bool pred = more_l && kslice_test(s, stm ^ 1, "W", n - 1);

    if (pred_sub || pred) {
      while (kslice_iter_in(&iter, &s1))
        kslice_clear(s1);

      if (pred_sub) {
        if (n == 1) {
          kslice_sub_read(s, s, stm ^ 1, "sub/win");
          predecessors_sub(stm, s, true);
        } else {
          // We must subtract sub/win from sub/cwin here.
          kslice_sub_read(-1, s, stm ^ 1, "sub/win");
          kslice_sub_read(s, s, stm ^ 1, "sub/cwin");
          kslice_sub_and_not(s, -1, stm ^ 1);
          predecessors_sub(stm, s, false);
        }
      }

      if (pred) {
        kslice_read(-1, s, stm ^ 1, "W", n - 1);
        predecessors(stm, s);
      }
    }

    while (kslice_iter_out(&iter, &s)) {
#if 0
        // If there are many predecessors, it might be more efficient to
        // filter them with this method.
        //
        // "wins" include illegal positions, so this also removes illegal
        // positions from X, which means we can remove the legality check
        // from check_successors().
        kslice_read(-1, s, stm, "wins", 0);
        kslice_and_not(s, -1);
        //
        // Removing positions with a drawing capture means we don't need to
        // test captures in check_successors().
        kslice_read(-1, s, stm, n <= DRAW_RULE ? "capt_bloss" : "capt_draw", 0);
        kslice_and_not(s, -1);
#endif
      kslice_write(s, s, stm, "X", n, UINT64_MAX); // FIXME
    }
  }

  create_dir(n, stm, "L");

  // Verify potential losses.
  kslice_iter_init(&iter, stm);
  while (kslice_iter_next(&iter, &s)) {

    if (kslice_test(s, stm, "X", n)) {
      while (kslice_iter_in(&iter, &s1)) {
        kslice_read(s1, s1, stm ^ 1, "wins", 0);
        // If there are very few predecessors, it might be more efficient to
        // directly probe_wdl() their captures.
        kslice_sub_read(s1, s1, stm ^ 1,
            n <= DRAW_RULE || !sub_cnt[stm ^ 1][cnt] ? "sub/win" : "sub/cwin");
      }

      kslice_read(-1, s, stm, "X", n);
      cnt += num = check_successors(stm, s);
      kslice_write(-1, s, stm, "L", n, num);
    }
    kslice_delete(s, stm, "X", n);

    while (kslice_iter_out(&iter, &s));
  }

  printf("l%d_%c = %lu\n", n, "wb"[stm], cnt);
  return cnt != 0;
}

static bool calc_W(int stm, int n, bool more_w)
{
  struct KSliceIterator iter;
  uint64_t cnt = 0, num;

  create_dir(n, stm, "W");

  // Calculate wins in n = predecessors(L(n-1))
  kslice_iter_init(&iter, stm);
  int s, s1;
  while (kslice_iter_next(&iter, &s)) {
    bool pred_sub =   (n == 1 && sub_cnt[stm ^ 1][0])
                   || (n == DRAW_RULE + 1 && sub_cnt[stm ^ 1][1]);
    bool pred = more_w && kslice_test(s, stm ^ 1, "L", n - 1);

    if (pred_sub || pred) {
      while (kslice_iter_in(&iter, &s1))
        kslice_clear(s1);

      if (pred_sub) {
        if (n == 1) {
          kslice_sub_read(s, s, stm ^ 1, "sub/loss");
          predecessors_sub(stm, s, false);
        } else {
          kslice_sub_read(s, s, stm ^ 1, "sub/bloss");
          predecessors_sub(stm, s, false);
        }
      }

      if (pred) {
        kslice_read(-1, s, stm ^ 1, "L", n - 1);
        predecessors(stm, s);
      }
    }

    while (kslice_iter_out(&iter, &s)) {
#if 0
      if (n == 1) {
        kslice_read(-1, s, stm, "capt_win", 0);
        kslice_or(s, -1);
      } else if (n == DRAW_RULE + 1) {
        kslice_read(-1, s, stm, "capt_cwin", 0);
        kslice_or(s, -1);
      }
#endif
      // We are currently removing illegal positions and known faster wins.
      // We could check easily for illegal positions.
      // But we MUST remove W_in_<=(N-1) and thus load "wins".
      // Again, we could work with deltas to avoid having to write "wins"
      // each time.
      kslice_read(-1, s, stm, "wins", 0);
      kslice_and_not(s, -1);
      cnt += num = kslice_count(s);
      kslice_write(s, s, stm, "W", n, num);
      kslice_or(-1, s);
      kslice_write(-1, s, stm, "wins", 0, UINT64_MAX);
    }
  }

  printf("w%d_%c = %lu\n", n, "wb"[stm], cnt);
  return cnt != 0;
}

void generate(void)
{
  // Calculate kslices for positions reached through a capture.
  make_dir("sub");
  calc_sub_kslices(WHITE);
  calc_sub_kslices(BLACK);

  calc_capt(WHITE);
  calc_capt(BLACK);

  create_dir(0, WHITE, "wins");
  create_dir(0, BLACK, "wins");
  create_dir(0, WHITE, "L");
  create_dir(0, BLACK, "L");
  calc_illegal_and_mate();

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

  max_iteration = n;
}
