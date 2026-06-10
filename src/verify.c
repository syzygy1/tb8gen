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
#include "tb8gen.h"
#include "threads.h"
#include "types.h"
#include "util.h"

const char wdl_name[5][8] = { "loss", "bloss", "draw", "cwin", "win" };
const char side[2][8] = { "white", "black" };
const char clr_L[4][4] = { "31", "32", "33", "34" };
const char clr_W[4][4] = { "94", "93", "92", "91" };

static uint64_t sub_cnt[2][5];

INLINE void mark_king_unmoves(int stm, Bitboard occ, uint8_t *restrict sq)
{
  uint8_t sq2[MAX_PIECES];
  uint8_t tmp = sq[stm];
  Bitboard b = king_attacks(sq[stm]) & ~(occ | king_attacks(sq[stm ^ 1]));
  while (b) {
    sq[stm] = pop_lsb(&b);
    normalize2(sq, sq2);
    int s2 = KKMap[sq2[0]][sq2[1]];
    uint8_t *p = kslice_get_address(s2);
    uint64_t idx = s2 < 441 ? sq_to_idx(sq2) : sq_to_idx_ref(sq2);
    kslice_bit_set_atomic(p, idx);
  }
  sq[stm] = tmp;
}

INLINE void mark_unmoves(int k, uint8_t *restrict const p, Bitboard occ,
    const uint8_t *restrict sq)
{
  uint8_t sq2[MAX_PIECES];
  Bitboard b = non_king_piece_moves(g_pos.pt[k], sq[k], occ);
  while (b) {
    memcpy(sq2, sq, sizeof sq2);
    sq2[k] = pop_lsb(&b);
    uint64_t idx = sq_to_idx(sq2);
    kslice_bit_set_atomic(p, idx);
  }
}

INLINE void mark_unmoves_ref(int k, uint8_t *restrict const p, Bitboard occ,
    const uint8_t *restrict sq)
{
  uint8_t sq2[MAX_PIECES];
  Bitboard b = non_king_piece_moves(g_pos.pt[k], sq[k], occ);
  while (b) {
    memcpy(sq2, sq, sizeof sq2);
    sq2[k] = pop_lsb(&b);
    uint64_t idx = sq_to_idx_ref(sq2);
    kslice_bit_set_atomic(p, idx);
  }
}

// Return true if one move hits a set bit.
INLINE bool test_king_moves(int stm, Bitboard occ, uint8_t *restrict sq)
{
  uint8_t sq2[MAX_PIECES];

  Bitboard b = king_attacks(sq[stm]) & ~king_attacks(sq[stm ^ 1]) & ~occ;

  while (b) {
    sq[stm] = pop_lsb(&b);
    normalize2(sq, sq2);
    int s = KKMap[sq2[0]][sq2[1]];
    uint64_t idx = s < 441 ? sq_to_idx(sq2) : sq_to_idx_ref(sq2);
    uint8_t *p = kslice_get_address(s);
    if (kslice_bit_test(p, idx))
      return true;
  } 

  return false;
}

INLINE bool test_moves(int k, int s, uint8_t *restrict const p, Bitboard occ,
    const uint8_t *restrict sq)
{
  uint8_t sq2[MAX_PIECES];

  Bitboard b = non_king_piece_attacks(g_pos.pt[k], sq[k], occ) & ~occ;

  while (b) {
    memcpy(sq2, sq, sizeof sq2);
    sq2[k] = pop_lsb(&b);
    uint64_t idx = sq_to_idx(sq2);
    if (kslice_bit_test(p, idx))
      return true;
  }

  return false;
}

INLINE bool test_moves_ref(int k, int s, uint8_t *restrict const p,
    Bitboard occ, const uint8_t *restrict sq)
{
  uint8_t sq2[MAX_PIECES];

  Bitboard b = non_king_piece_attacks(g_pos.pt[k], sq[k], occ) & ~occ;

  while (b) {
    memcpy(sq2, sq, sizeof sq2);
    sq2[k] = pop_lsb(&b);
    uint64_t idx = sq_to_idx_ref(sq2);
    if (kslice_bit_test(p, idx))
      return true;
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
  struct IdxState is;
  Position pos = g_pos;
  int k = work_set;
  int m = ri.last[k];
  int n = --pos.num;

  pos.pt[m] = pos.pt[n];
  uint8_t *restrict p[5];
  for (int i = 0; i < 5; i++)
    p[i] = kslice_sub_buf[i] + sub_offset[k];

  idx_state_init(&is, thread->begin, pos.sq, &capt_ri[k]);
  for (uint64_t idx = thread->begin, end = thread->end; idx < end;
      idx++, idx_state_inc(&is, &capt_ri[k]))
  {
    pos.occ = idx_state_to_sq(&is, pos.sq, &capt_ri[k]);
    pos.sq[m] = pos.sq[n];
    if (!opp_king_attacked(&pos)) {
      int v = probe_wdl(&pos, -2, 2);
      kslice_bit_set(p[v + 2], idx);
    }
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
    sub_cnt[stm][3] += c[3];
    sub_cnt[stm][4] += c[4];
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
  struct IdxState is;
  Position pos = g_pos;
  int stm = pos.stm;
  pos.stm ^= 1;
  int n = --pos.num;
  int k = work_set;
  int s = work_slice;

  int m = ri.last[k];
  pos.pt[m] = pos.pt[n];

  uint64_t *restrict p =
    (uint64_t *)kslice_sub_get_address(s, k) + (thread->begin >> 6);
  uint8_t *restrict const q = kslice_get_address(s);

  uint64_t last = thread->begin;
  idx_state_init(&is, last, pos.sq, &capt_ri[k]);
  idx_state_to_sq(&is, pos.sq, &capt_ri[k]);

  for (uint64_t idx = last, end = thread->end; idx < end; idx += 64) {
    uint64_t w = *p++;
    while (w) {
      uint64_t cur = idx + pop_lsb(&w);
      idx_state_add(&is, cur - last, &capt_ri[k]);
      last = cur;
      Bitboard occ = idx_state_to_sq(&is, pos.sq, &capt_ri[k]);
      // Uncapture by king.
      pos.sq[m] = pos.sq[stm];
      mark_king_unmoves(stm, occ, pos.sq);
      // Uncapture by non-king pieces.
      for (int i = 0; pos.pcs[stm][i] >= 0; i++) {
        int j = pos.pcs[stm][i];
        pos.sq[m] = pos.sq[j];
        if (s < 441)
          mark_unmoves(j, q, occ, pos.sq);
        else
          mark_unmoves_ref(j, q, occ, pos.sq);
      }
    }
  }
}

// Uncapture from the loaded subtable of stm^1 positions to stm positions.
static void predecessors_sub(int stm, int s)
{
  work_slice = s;

  g_pos.stm = stm;
  g_pos.sq[0] = KKSquare[s][0];
  g_pos.sq[1] = KKSquare[s][1];

  // Loop through the sets from which a piece is removed.
  for (int k = 0; k < ri.numsets; k++) {
    int m = ri.last[k];
    if ((g_pos.pt[m] >> 3) == stm)
      continue;
    work_set = k;
    run_threaded(predecessors_sub_worker, &work_capt_dynamic[k]);
  }
}

static void calc_capt(int stm, int wdl)
{
  if (!sub_cnt[stm ^ 1][2 - wdl])
    return;

  char capt_name[32], sub_name[32];
  strcat(strcpy(capt_name, "capt/"), wdl_name[2 + wdl]);
  strcat(strcpy(sub_name , "sub/" ), wdl_name[2 - wdl]);

  char str[64];
  sprintf(str, "%s/%c/done", capt_name, "wb"[stm]);
  if (file_exists(str))
    return;
#if 0
  FILE *F = fopen(str, "rb");
  if (F) {
    file_read(&g_stats[stm][n], 8, F);
    fclose(F);
    return;
  }
#endif

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

  create_empty(str);

  snprintf(phase, sizeof phase, "%s %s captures: %lu", side[stm],
      wdl_name[2 + wdl], cnt);
  show_progress(phase, 462, 462, true);
}

static void calc_illegal_worker(struct ThreadData *thread)
{
  struct IdxState is;
  Position pos = g_pos;
  int k = work_set;
  int m = ri.last[k];
  int stm = pos.stm;
  int king_sq = pos.sq[stm ^ 1];

  uint8_t *restrict const p = kslice_buf[0];

  idx_state_init(&is, thread->begin, pos.sq, &capt_ri[k]);

  for (uint64_t idx = thread->begin, end = thread->end; idx < end;
      idx++, idx_state_inc(&is, &capt_ri[k]))
  {
    Bitboard occ = idx_state_to_sq(&is, pos.sq, &capt_ri[k]);
    pos.sq[m] = king_sq;
    mark_unmoves(m, p, occ, pos.sq);
  }
}

static void calc_illegal_ref_worker(struct ThreadData *thread)
{
  struct IdxState is;
  Position pos = g_pos;
  int k = work_set;
  int m = ri.last[k];
  int stm = pos.stm;
  int king_sq = pos.sq[stm ^ 1];

  uint8_t *restrict const p = kslice_buf[0];

  idx_state_init(&is, thread->begin, pos.sq, &capt_ri[k]);

  for (uint64_t idx = thread->begin, end = thread->end; idx < end;
      idx++, idx_state_inc(&is, &capt_ri[k]))
  {
    Bitboard occ = idx_state_to_sq(&is, pos.sq, &capt_ri[k]);
    pos.sq[m] = king_sq;
    mark_unmoves_ref(m, p, occ, pos.sq);
  }
}

void calc_illegal(int stm)
{
  bool partial = dir_exists(-1, stm, "illegal");

  create_dir(-1, stm, "illegal");

  g_pos.stm = stm;

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

    g_pos.sq[0] = KKSquare[s][0];
    g_pos.sq[1] = KKSquare[s][1];

    kslice_clear_addr(kslice_buf[0], s);

    for (int k = 0; k < ri.numsets; k++) {
      if ((g_pos.pt[ri.first[k]] >> 3) != stm)
        continue;
      work_set = k;
      if (s < 441)
        run_threaded(calc_illegal_worker, &work_capt_dynamic[k]);
      else
        run_threaded(calc_illegal_ref_worker, &work_capt_dynamic[k]);
    }

    cnt += num = kslice_count_addr(kslice_buf[0], s);
    kslice_write_addr(kslice_buf[0], s, stm, "illegal", -1, num);
  }

  snprintf(phase, sizeof phase, "%lu %s illegal positions", cnt, side[stm]);
  show_progress(phase, 462, 462, true);
}

static mtx_t report_mutex;
static int num_fails;

static void report_fail(int s, uint64_t idx, Position *pos)
{
  mtx_lock(&report_mutex);
  char fenstr[64];
  pos_to_fen(pos, fenstr, false);
  int wdl = probe_wdl(pos, -2, 2);
  printf("\nslice = %d, idx = %lu, wdl = %d\n%s\n", s, idx, wdl, fenstr);
  num_fails++;
  if (num_fails == 10)
    exit(EXIT_FAILURE);
  mtx_unlock(&report_mutex);
}

static int work_slice;
static bool check_stalemate;

static void check_zero_worker(struct ThreadData *thread)
{
  struct IdxState is;
  Position pos = g_pos;
  int stm = pos.stm;
  int s = work_slice;

  uint64_t *restrict p = (uint64_t *)kslice_get_address(-1);
  uint8_t *restrict const q = kslice_get_address(s);

  p += thread->begin >> 6;
  uint64_t last = thread->begin;
  idx_state_init(&is, last, pos.sq, &ri);
  idx_state_to_sq(&is, pos.sq, &ri);
  for (uint64_t idx = last, end = thread->end; idx < end; idx += 64, p++) {
    uint64_t w = *p;
    while (w) {
      uint64_t cur = idx + pop_lsb(&w);
      idx_state_add(&is, cur - last, &ri);
      last = cur;
      Bitboard occ = pos.occ = idx_state_to_sq(&is, pos.sq, &ri);
      for (int i = 0; pos.pcs[stm][i] >= 0; i++) {
        int j = pos.pcs[stm][i];
        if (test_moves(j, s, q, occ, pos.sq))
          goto fail;
      }
      uint8_t tmp = pos.sq[stm];
      bool v = test_king_moves(stm, occ, pos.sq);
      pos.sq[stm] = tmp;
      if (v || (check_stalemate && !has_legal_moves(&pos))) {
fail:
        report_fail(s, cur, &pos);
      }
    }
  }
}

static void check_zero_ref_worker(struct ThreadData *thread)
{
  Position pos = g_pos;
  int stm = pos.stm;
  int s = work_slice;
  Bitboard kings = bit(pos.sq[0]) | bit(pos.sq[1]);

  uint64_t *restrict p = (uint64_t *)kslice_get_address(-1);
  uint8_t *restrict const q = kslice_get_address(s);

  p += thread->begin >> 6;
  uint64_t last = thread->begin;
  for (uint64_t idx = last, end = thread->end; idx < end; idx += 64, p++) {
    uint64_t w = *p;
    while (w) {
      unsigned bt = pop_lsb(&w);
      uint64_t cur = idx + bt;
      Bitboard occ = pos.occ = unrank_reflection(cur, pos.sq, kings, &ri);
      for (int i = 0; pos.pcs[stm][i] >= 0; i++) {
        int j = pos.pcs[stm][i];
        if (test_moves_ref(j, s, q, occ, pos.sq))
          goto fail;
      }
      uint8_t tmp = pos.sq[stm];
      bool v = test_king_moves(stm, occ, pos.sq);
      pos.sq[stm] = tmp;
      if (v || (check_stalemate && !has_legal_moves(&pos))) {
fail:
        report_fail(s, cur, &pos);
      }
    }
  }
}

static void check_zero(int stm, int s)
{
  work_slice = s;
  g_pos.stm = stm;
  g_pos.sq[0] = KKSquare[s][0];
  g_pos.sq[1] = KKSquare[s][1];

  if (s < 441)
    run_threaded(check_zero_worker, &work_g_dynamic[0]);
  else
    run_threaded(check_zero_ref_worker, &work_g_dynamic[1]);
}

static void check_one_worker(struct ThreadData *thread)
{
  struct IdxState is;
  Position pos = g_pos;
  int stm = pos.stm;
  int s = work_slice;

  uint64_t *restrict p = (uint64_t *)kslice_get_address(-1);
  uint8_t *restrict const q = kslice_get_address(s);

  p += thread->begin >> 6;
  uint64_t last = thread->begin;
  idx_state_init(&is, last, pos.sq, &ri);
  idx_state_to_sq(&is, pos.sq, &ri);
  for (uint64_t idx = last, end = thread->end; idx < end; idx += 64, p++) {
    uint64_t w = *p;
    while (w) {
      uint64_t cur = idx + pop_lsb(&w);
      idx_state_add(&is, cur - last, &ri);
      last = cur;
      Bitboard occ = idx_state_to_sq(&is, pos.sq, &ri);
      for (int i = 0; pos.pcs[stm][i] >= 0; i++) {
        int j = pos.pcs[stm][i];
        if (test_moves(j, s, q, occ, pos.sq))
          goto next;
      }
      uint8_t tmp = pos.sq[stm];
      bool v = test_king_moves(stm, occ, pos.sq);
      pos.sq[stm] = tmp;
      if (!v) {
        pos.occ = occ;
        if (!check_stalemate || has_legal_moves(&pos))
          report_fail(s, cur, &pos);
      }
    }
next:
  }
}

static void check_one_ref_worker(struct ThreadData *thread)
{
  Position pos = g_pos;
  int stm = pos.stm;
  int s = work_slice;
  Bitboard kings = bit(pos.sq[0]) | bit(pos.sq[1]);

  uint64_t *restrict p = (uint64_t *)kslice_get_address(-1);
  uint8_t *restrict const q = kslice_get_address(s);

  p += thread->begin >> 6;
  uint64_t last = thread->begin;
  for (uint64_t idx = last, end = thread->end; idx < end; idx += 64, p++) {
    uint64_t w = *p;
    while (w) {
      unsigned bt = pop_lsb(&w);
      uint64_t cur = idx + bt;
      Bitboard occ = unrank_reflection(cur, pos.sq, kings, &ri);
      for (int i = 0; pos.pcs[stm][i] >= 0; i++) {
        int j = pos.pcs[stm][i];
        if (test_moves_ref(j, s, q, occ, pos.sq))
          goto next;
      }
      uint8_t tmp = pos.sq[stm];
      bool v = test_king_moves(stm, occ, pos.sq);
      pos.sq[stm] = tmp;
      if (!v) {
        pos.occ = occ;
        if (!check_stalemate || has_legal_moves(&pos))
          report_fail(s, cur, &pos);
      }
    }
next:
  }
}

static void check_one(int stm, int s)
{
  work_slice = s;
  g_pos.stm = stm;
  g_pos.sq[0] = KKSquare[s][0];
  g_pos.sq[1] = KKSquare[s][1];

  if (s < 441)
    run_threaded(check_one_worker, &work_g_dynamic[0]);
  else
    run_threaded(check_one_ref_worker, &work_g_dynamic[1]);
}

static void check_win(int stm)
{
  struct KSliceIterator iter;

  char phase[64];
  snprintf(phase, sizeof phase, "check %s win, cwin-", side[stm]);
  int num_done = 0;
  check_stalemate = false;

  kslice_iter_init(&iter, stm);
  int s, s1;
  while (kslice_iter_next(&iter, &s)) {
    show_progress(phase, num_done++, 462, false);

    while (kslice_iter_in(&iter, &s1))
      kslice_read(s1, s1, stm ^ 1, "loss", -1);

    kslice_read(-1, s, stm, "wins", -1);
    kslice_read_andnot(-1, s, stm, "capt/win", -1);
    check_one(stm, s);

    kslice_read(-1, s, stm, "cwin", -1);
    kslice_read_andnot(-1, s, stm, "capt/cwin", -1);
    check_zero(stm, s);

    while (kslice_iter_out(&iter, &s));
  }
  show_progress(phase, 462, 462, true);
}

static void check_cwin(int stm)
{
  struct KSliceIterator iter;

  char phase[64];
  snprintf(phase, sizeof phase, "check %s cwin+", side[stm]);
  int num_done = 0;
  check_stalemate = false;

  kslice_iter_init(&iter, stm);
  int s, s1;
  while (kslice_iter_next(&iter, &s)) {
    show_progress(phase, num_done++, 462, false);

    while (kslice_iter_in(&iter, &s1))
      kslice_read(s1, s1, stm ^ 1, "bloss", -1);

    kslice_read(-1, s, stm, "cwin", -1);
    kslice_read_andnot(-1, s, stm, "capt/cwin", -1);
    check_one(stm, s);

    while (kslice_iter_out(&iter, &s));
  }
  show_progress(phase, 462, 462, true);
}

static void check_draw(int stm)
{
  struct KSliceIterator iter;

  char phase[64];
  snprintf(phase, sizeof phase, "check %s draw-", side[stm]);
  int num_done = 0;
  check_stalemate = false;

  kslice_iter_init(&iter, stm);
  int s, s1;
  while (kslice_iter_next(&iter, &s)) {
    show_progress(phase, num_done++, 462, false);

    while (kslice_iter_in(&iter, &s1)) {
      kslice_read(s1, s1, stm ^ 1, "loss", -1);
      kslice_read_or(s1, s1, stm ^ 1, "bloss", -1);
    }

    kslice_read(-1, s, stm, "draw", -1);
    kslice_read_andnot(-1, s, stm, "capt/draw", -1);
    check_zero(stm, s);

    while (kslice_iter_out(&iter, &s));
  }
  show_progress(phase, 462, 462, true);

  snprintf(phase, sizeof phase, "check %s draw+", side[stm]);
  num_done = 0;
  check_stalemate = true;
  kslice_iter_init(&iter, stm);
  while (kslice_iter_next(&iter, &s)) {
    show_progress(phase, num_done++, 462, false);

    while (kslice_iter_in(&iter, &s1))
      kslice_read(s1, s1, stm ^ 1, "draw", -1);

    kslice_read(-1, s, stm, "draw", -1);
    kslice_read_andnot(-1, s, stm, "capt/draw", -1);
    check_one(stm, s);

    while (kslice_iter_out(&iter, &s));
  }
  show_progress(phase, 462, 462, true);
}

static void check_bloss(int stm)
{
  struct KSliceIterator iter;

  char phase[64];
  snprintf(phase, sizeof phase, "check %s bloss-", side[stm]);
  int num_done = 0;
  check_stalemate = false;

  kslice_iter_init(&iter, stm);
  int s, s1;
  while (kslice_iter_next(&iter, &s)) {
    show_progress(phase, num_done++, 462, false);

    while (kslice_iter_in(&iter, &s1)) {
      kslice_read(s1, s1, stm ^ 1, "loss", -1);
      kslice_read_or(s1, s1, stm ^ 1, "bloss", -1);
      kslice_read_or(s1, s1, stm ^ 1, "draw", -1);
    }

    kslice_read(-1, s, stm, "bloss", -1);
    kslice_read_andnot(-1, s, stm, "capt/bloss", -1);
    check_zero(stm, s);

    while (kslice_iter_out(&iter, &s));
  }
  show_progress(phase, 462, 462, true);

  snprintf(phase, sizeof phase, "check %s bloss+", side[stm]);
  num_done = 0;
  kslice_iter_init(&iter, stm);
  while (kslice_iter_next(&iter, &s)) {
    show_progress(phase, num_done++, 462, false);

    while (kslice_iter_in(&iter, &s1))
      kslice_read(s1, s1, stm ^ 1, "cwin", -1);

    kslice_read(-1, s, stm, "bloss", -1);
    kslice_read_andnot(-1, s, stm, "capt/bloss", -1);
    check_one(stm, s);

    while (kslice_iter_out(&iter, &s));
  }
  show_progress(phase, 462, 462, true);
}

static void check_loss(int stm)
{
  struct KSliceIterator iter;

  char phase[64];
  snprintf(phase, sizeof phase, "check %s loss", side[stm]);
  int num_done = 0;
  check_stalemate = true;

  kslice_iter_init(&iter, stm);
  int s, s1;
  while (kslice_iter_next(&iter, &s)) {
    show_progress(phase, num_done++, 462, false);

    while (kslice_iter_in(&iter, &s1)) {
      kslice_read(s1, s1, stm ^ 1, "win", -1);
      kslice_read_or(s1, s1, stm ^ 1, "illegal", -1);
      kslice_not(s1);
    }

    kslice_read(-1, s, stm, "bloss", -1);
    kslice_read_andnot(-1, s, stm, "capt/bloss", -1);
    check_zero(stm, s);

    while (kslice_iter_out(&iter, &s));
  }
  show_progress(phase, 462, 462, true);
}

uint8_t *tb_table;
struct TbTable2 *table_info;
struct RankInfo tb_ri[2];

INLINE void store_wdl(uint64_t idx, uint8_t *restrict tmp, int v)
{
  uint64_t **const p = (uint64_t **)kslice_buf;
  tmp[idx & 63] = v;
  if ((idx & 63) == 63) {
#ifdef __AVX512BW__
    __m512i x = _mm512_load_si512(tmp);
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
void depermute_worker(struct ThreadData *thread)
{
  alignas(64) uint8_t tmp[64];
  uint8_t sq[MAX_PIECES];
  uint8_t *restrict src = tb_table;
  struct TbTable2 *table = table_info;
  struct RankInfo *dst_ri = &tb_ri[g_pos.stm];
  struct IdxState is;

  uint64_t idx_dec_buf[NUM];

  sq[0] = g_pos.sq[0];
  sq[1] = g_pos.sq[1];
  Bitboard occ = bit(sq[0]) | bit(sq[1]);

  idx_state_init(&is, thread->begin, sq, dst_ri);

  uint64_t idx = thread->begin, end = thread->end;
  int fill = 0;
  int head = 0; // Next buffered element to consume.

  // Fill pipeline.
  for (; fill < NUM && idx < end; fill++, idx++, idx_state_inc(&is, dst_ri)) {
    idx_state_to_sq(&is, sq, dst_ri);
    uint64_t idx_dec = rank_trivial_from(sq, 0, occ, table->first, table->ri);
    __builtin_prefetch(src + idx_dec, 0, 3);
    idx_dec_buf[fill] = idx_dec;
  }

  // Steady-state pipeline.
  for (; idx < end; idx++, idx_state_inc(&is, dst_ri)) {
    idx_state_to_sq(&is, sq, dst_ri);
    uint64_t idx_dec = rank_trivial_from(sq, 0, occ, table->first, table->ri);
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

void depermute_ref_worker(struct ThreadData *thread)
{
  alignas(64) uint8_t tmp[64];
  uint8_t sq[MAX_PIECES];
  uint8_t *restrict src = tb_table;
  struct TbTable2 *table = table_info;
  struct RankInfo *dst_ri = &tb_ri[g_pos.stm];

  uint64_t idx_dec_buf[NUM];

  sq[0] = g_pos.sq[0];
  sq[1] = g_pos.sq[1];
  Bitboard occ = bit(sq[0]) | bit(sq[1]);

  uint64_t idx = thread->begin, end = thread->end;
  int fill = 0;
  int head = 0; // Next buffered element to consume.

  // Fill pipeline.
  for (; fill < NUM && idx < end; fill++, idx++) {
    unrank_reflection(idx, sq, occ, dst_ri);
    uint64_t idx_dec = rank_reflection(sq, occ, table->first, table->ri);
    __builtin_prefetch(src + idx_dec, 0, 3);
    idx_dec_buf[fill] = idx_dec;
  }

  // Steady-state pipeline.
  for (; idx < end; idx++) {
    unrank_reflection(idx, sq, occ, dst_ri);
    uint64_t idx_dec = rank_reflection(sq, occ, table->first, table->ri);
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

bool flip[2], btm_side[2];

uint64_t wdl_cnt[2][5];

void decompress_kslice(struct Tbase *tb, int stm, int s)
{
  Position pos = g_pos;

  int tsq = KKMap[pos.sq[0]][pos.sq[1]];
  int t = !symmetric ? 2 * tsq + btm_side[stm] : tsq;
  if (!tb->table[t])
    tb->table[t] = init_new_table(tb, pos.num, WDL, t, tsq);
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
      run_threaded(depermute_worker, &work_g_static[0]);
    else
      run_threaded(depermute_ref_worker, &work_g_static[1]);
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

void verify(void)
{
  // Allocate memory for decompressed WDL K-slice.
  tb_table = alloc_huge((kslice_size + 63) & ~0x3f);
  if (!tb_table)
    out_of_mem();

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
    fprintf(stderr, "Could not open %s.\n", g_tablename);
    exit(EXIT_FAILURE);
  }

  if (tb->layout != LT_PIECE_KK) {
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

  for (int stm = 0; stm < 2; stm++) {
    char phase[64];
    snprintf(phase, sizeof phase, "loading %ctm slices", "wb"[stm]);

    for (int i = 0; i < 5; i++)
      create_dir(-1, stm, wdl_name[i]);

    tb_ri[stm] = ri;
    g_pos.stm = stm;
    flip[stm] = !symmetric ? !tb->flipped : stm != WHITE;
    btm_side[stm] = (stm == WHITE) == flip[stm];

    for (int k = 0; k < ri.numsets; k++) {
      int pt = g_pos.pt[ri.first[k]] ^ (flip[stm] << 3);
      for (int j = 0; j < g_pos.num; j++)
        if (tb->pt[j] == pt) {
          tb_ri[stm].first[k] = j;
          break;
        }
    }

    for (int s = 0; s < 462; s++) {
      show_progress(phase, s, 462, false);

      g_pos.sq[0] = KKSquare[s][0];
      g_pos.sq[1] = KKSquare[s][1];
      g_pos.stm = stm;

      if (flip[stm])
        Swap(g_pos.sq[0], g_pos.sq[1]);

      decompress_kslice(tb, stm, s);
    }
    show_progress(phase, 462, 462, true);

    for (int i = 0; i < 5; i++)
      printf("wdl_cnt[%d][%d] = %lu\n", stm, i, wdl_cnt[stm][i]);
  }

  free(tb_table);

  mtx_init(&report_mutex, mtx_plain);

  check_win(WHITE);
  check_cwin(WHITE);
  check_draw(WHITE);
  check_bloss(WHITE);
  check_loss(WHITE);

  check_win(BLACK);
  check_cwin(BLACK);
  check_draw(BLACK);
  check_bloss(BLACK);
  check_loss(BLACK);

  if (num_fails == 0)
    printf("WDL table is OK.\n");

  mtx_destroy(&report_mutex);
}
