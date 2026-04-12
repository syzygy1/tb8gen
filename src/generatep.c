/*
  Copyright (c) 2026 Ronald de Man

  This file is distributed under the terms of the GNU GPL, version 2.
*/

#include <inttypes.h>
#include <stdint.h>

#include "defs.h"
#include "generatep.h"
#include "indexp.h"
#include "kslicep.h"
#include "movegen.h"
#include "probe.h"
#include "tb8genp.h"
#include "threads.h"
#include "types.h"

uint64_t sub_cnt[2][5];
uint64_t psub_cnt[5];
uint64_t pawn_cnt[5];
int max_iteration;
static uint8_t *merged_table, *merged_table2;

void list_positions(int s, uint8_t *restrict q)
{
  uint32_t sub[MAX_SETS];
  Position pos = g_pos;

  for (int r = 0; r < 16; r++) {
    pos.sq[0] = KK16Square[s][r][0];
    pos.sq[1] = KK16Square[s][r][1];

    uint8_t *restrict p = q + k16offset(pos.sq);

    idx_to_sq_init(0, sub, &ii);
    for (uint64_t idx = 0; idx < kslice_size; idx++, idx_to_sq_inc(sub, &ii)) {
      idx_to_sq(sub, pos.sq);
      if (kslice_bit_test(p, idx)) {
        printf("%d %d %d %d %d\n", pos.sq[0], pos.sq[1], pos.sq[3],
            pos.sq[4], pos.sq[2]);
      }
    }
  }
}

INLINE void mark_king_unmoves(int stm, Bitboard occ, const uint8_t *restrict sq)
{
  uint8_t sq2[MAX_PIECES];
  memcpy(sq2, sq, sizeof sq2);
  Bitboard b = king_attacks(sq[stm]) & ~(occ | king_attacks(sq[stm ^ 1]));
  while (b) {
    sq2[stm] = pop_lsb(&b);
    uint8_t *p = kslice_get_address(sq2);
    kslice_bit_set_atomic(p, sq_to_idx(sq2));
  }
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

static int work_set, work_r, work_lower, work_upper;

static void calc_sub_worker(struct ThreadData *thread)
{
  Position pos = g_pos;
  uint32_t sub[MAX_SETS];
  int k = work_set;
  int r = work_r;
  int m = ii.last[k];
  int n = --pos.num;

  pos.pt[m] = pos.pt[n];
  uint8_t *restrict p[5];
  for (int i = 0; i < 5; i++)
    p[i] = k16slice_sub_buf[i] + sub_offset[k] + r * kslice_sub_alloc_size[k];

  idx_to_sq_init(thread->begin, sub, &capt_ii[k]);
  for (uint64_t idx = thread->begin, end = thread->end; idx < end;
      idx++, idx_to_sq_inc(sub, &capt_ii[k]))
  {
    pos.occ = capt_idx_to_sq(sub, pos.sq, k);
    pos.sq[m] = pos.sq[n];
    if (opp_king_attacked(&pos))
      continue;
    int v = probe_wdl(&pos, -2, 2);
    kslice_bit_set(p[v + 2], idx);
  }
}

// Calculate aggregate bitmaps for subtables, one per loss/bloss/draw/cwin/win.
static void calc_sub_kslices(int stm)
{
  create_dir(-1, stm, "sub/loss");
  create_dir(-1, stm, "sub/bloss");
  create_dir(-1, stm, "sub/draw");
  create_dir(-1, stm, "sub/cwin");
  create_dir(-1, stm, "sub/win");
  create_dir(-1, stm, "sub/noloss");
  create_dir(-1, stm, "sub/nobloss");

  g_pos.stm = stm;

  for (int s = 0; s < 240; s++) {
    for (int i = 0; i < 5; i++)
      k16slice_sub_clear_addr(k16slice_sub_buf[i], stm);

    for (int t = 0; t < g_num_threads; t++)
      g_thread_data[t].cnt = 0;

    for (int r = 0; r < 16; r++) {
      g_pos.sq[0] = KK16Square[s][r][0];
      g_pos.sq[1] = KK16Square[s][r][1];
      work_r = r;

      if (is_broken(&g_pos))
        continue;

      for (int k = 0; k < ii.numsets; k++) {
        if ((g_pos.pt[ii.first[k]] >> 3) != stm)
          continue;
        work_set = k;
        run_threaded(calc_sub_worker, work_capt[k], 0);
      }
    }

    uint64_t cnt_l, cnt_bl, cnt_d, cnt_cw, cnt_w;

    cnt_l = k16slice_sub_count_addr(k16slice_sub_buf[0], stm);
    k16slice_sub_write_addr(k16slice_sub_buf[0], s, stm, "sub/loss", cnt_l);

    cnt_bl = k16slice_sub_count_addr(k16slice_sub_buf[1], stm);
    k16slice_sub_write_addr(k16slice_sub_buf[1], s, stm, "sub/bloss", cnt_bl);

    cnt_d = k16slice_sub_count_addr(k16slice_sub_buf[2], stm);
    k16slice_sub_write_addr(k16slice_sub_buf[2], s, stm, "sub/draw", cnt_d);

    cnt_cw = k16slice_sub_count_addr(k16slice_sub_buf[3], stm);
    k16slice_sub_write_addr(k16slice_sub_buf[3], s, stm, "sub/cwin", cnt_cw);

    cnt_w = k16slice_sub_count_addr(k16slice_sub_buf[4], stm);
    k16slice_sub_write_addr(k16slice_sub_buf[4], s, stm, "sub/win", cnt_w);

    sub_cnt[stm][0] += cnt_l;
    sub_cnt[stm][1] += cnt_bl;
    sub_cnt[stm][2] += cnt_d;
    sub_cnt[stm][3] += cnt_cw;
    sub_cnt[stm][4] += cnt_w;

    k16slice_sub_or_addr(k16slice_sub_buf[0], k16slice_sub_buf[1], stm);
    k16slice_sub_or_addr(k16slice_sub_buf[0], k16slice_sub_buf[2], stm);
    k16slice_sub_write_addr(k16slice_sub_buf[0], s, stm, "sub/nobloss",
        UINT64_MAX);

    k16slice_sub_or_addr(k16slice_sub_buf[0], k16slice_sub_buf[3], stm);
    k16slice_sub_write_addr(k16slice_sub_buf[0], s, stm, "sub/noloss",
        UINT64_MAX);
  }
}

static void calc_psub_worker(struct ThreadData *thread)
{
  Position pos = g_pos;
  uint32_t sub[MAX_SETS];
  int k = work_set;
  int r = work_r;
  int m = ii.last[k];
  int n = --pos.num;

  pos.pt[2] = pos.pt[m];
  pos.pt[m] = pos.pt[n];
  uint8_t *restrict p[5];
  for (int i = 0; i < 5; i++)
    p[i] = k16slice_sub_buf[i] + psub_offset[k] + r * kslice_sub_alloc_size[k];

  idx_to_sq_init(thread->begin, sub, &capt_ii[k]);
  for (uint64_t idx = thread->begin, end = thread->end; idx < end;
      idx++, idx_to_sq_inc(sub, &capt_ii[k]))
  {
    pos.occ = capt_idx_to_sq(sub, pos.sq, k);
    pos.sq[m] = pos.sq[n];
    if (opp_king_attacked(&pos))
      continue;
    int v = probe_wdl(&pos, -2, 2);
    kslice_bit_set(p[v + 2], idx);
  }
}

// Calculate aggregate bitmaps for subtables reached through a capture of
// the black pawn by a white piece (non-king).
static void calc_psub_kslices(void)
{
  create_dir(-1, WHITE, "psub/loss");
  create_dir(-1, WHITE, "psub/bloss");
  create_dir(-1, WHITE, "psub/draw");
  create_dir(-1, WHITE, "psub/cwin");
  create_dir(-1, WHITE, "psub/win");
  create_dir(-1, WHITE, "psub/noloss");
  create_dir(-1, WHITE, "psub/nobloss");

  g_pos.stm = BLACK;

  for (int s = 0; s < 240; s++) {
    for (int i = 0; i < 5; i++)
      k16slice_sub_clear_addr(k16slice_sub_buf[i], WHITE);

    for (int t = 0; t < g_num_threads; t++)
      g_thread_data[t].cnt = 0;

    for (int r = 0; r < 16; r++) {
      g_pos.sq[0] = KK16Square[s][r][0];
      g_pos.sq[1] = KK16Square[s][r][1];
      work_r = r;

      if (is_broken(&g_pos))
        continue;

      for (int k = 0; k < ii.numsets; k++) {
        if ((g_pos.pt[ii.first[k]] >> 3) != WHITE)
          continue;
        work_set = k;
        run_threaded(calc_psub_worker, work_capt[k], 0);
      }
    }

    uint64_t cnt_l, cnt_bl, cnt_d, cnt_cw, cnt_w;

    cnt_l = k16slice_sub_count_addr(k16slice_sub_buf[0], WHITE);
    k16slice_sub_write_addr(k16slice_sub_buf[0], s, WHITE, "psub/loss", cnt_l);

    cnt_bl = k16slice_sub_count_addr(k16slice_sub_buf[1], WHITE);
    k16slice_sub_write_addr(k16slice_sub_buf[1], s, WHITE, "psub/bloss",
        cnt_bl);

    cnt_d = k16slice_sub_count_addr(k16slice_sub_buf[2], WHITE);
    k16slice_sub_write_addr(k16slice_sub_buf[2], s, WHITE, "psub/draw", cnt_d);

    cnt_cw = k16slice_sub_count_addr(k16slice_sub_buf[3], WHITE);
    k16slice_sub_write_addr(k16slice_sub_buf[3], s, WHITE, "psub/cwin", cnt_cw);

    cnt_w = k16slice_sub_count_addr(k16slice_sub_buf[4], WHITE);
    k16slice_sub_write_addr(k16slice_sub_buf[4], s, WHITE, "psub/win", cnt_w);

    psub_cnt[0] += cnt_l;
    psub_cnt[1] += cnt_bl;
    psub_cnt[2] += cnt_d;
    psub_cnt[3] += cnt_cw;
    psub_cnt[4] += cnt_w;
  
    k16slice_sub_or_addr(k16slice_sub_buf[0], k16slice_sub_buf[1], WHITE);
    k16slice_sub_or_addr(k16slice_sub_buf[0], k16slice_sub_buf[2], WHITE);
    k16slice_sub_write_addr(k16slice_sub_buf[0], s, WHITE, "psub/nobloss",
        UINT64_MAX);

    k16slice_sub_or_addr(k16slice_sub_buf[0], k16slice_sub_buf[3], WHITE);
    k16slice_sub_write_addr(k16slice_sub_buf[0], s, WHITE, "psub/noloss",
        UINT64_MAX);
  }

  for (int i = 0; i < 5; i++)
    printf("psub_cnt[%d] = %lu\n", i, psub_cnt[i]);
}

static void predecessors_sub_worker(struct ThreadData *thread)
{
  Position pos = g_pos;
  uint32_t sub[MAX_SETS];
  int stm = pos.stm;
  int k = work_set;
  int m = ii.last[k];

  uint64_t *restrict p =
    (uint64_t *)kslice_sub_get_address(pos.sq, k) + (thread->begin >> 6);
  uint8_t *restrict const q = kslice_get_address(pos.sq);

  uint64_t last = thread->begin;
  idx_to_sq_init(last, sub, &capt_ii[k]);

  for (uint64_t idx = last, end = thread->end; idx < end; idx += 64) {
    uint64_t w = *p++;
    while (w) {
      uint64_t cur = idx + pop_lsb(&w);
      idx_to_sq_add(cur - last, sub, &capt_ii[k]);
      last = cur;
      Bitboard occ = capt_idx_to_sq(sub, pos.sq, k);
      // Uncapture by king.
      pos.sq[m] = pos.sq[stm];
      mark_king_unmoves(stm, occ, pos.sq);
      // Uncapture by non-king pieces.
      for (int i = 0; pos.pcs[stm][i] >= 0; i++) {
        int j = pos.pcs[stm][i];
        pos.sq[m] = pos.sq[j];
        mark_unmoves(j, q, occ, pos.sq);
      }
    }
  }
}

// Uncapture pieces from the loaded subtable of stm^1 positions to stm
// positions.
static void predecessors_sub(int stm, int s)
{
  g_pos.stm = stm;
  for (int r = 0; r < 16; r++) {
    g_pos.sq[0] = KK16Square[s][r][0];
    g_pos.sq[1] = KK16Square[s][r][1];

    if (is_broken(&g_pos))
      continue;

    // Loop through the sets from which a piece is removed.
    for (int k = 0; k < ii.numsets; k++) {
      int m = ii.last[k];
      if ((g_pos.pt[m] >> 3) == stm)
        continue;
      work_set = k;
      run_threaded(predecessors_sub_worker, work_capt[k], 0);
    }
  }
}

static void predecessors_psub_worker(struct ThreadData *thread)
{
  Position pos = g_pos;
  uint32_t sub[MAX_SETS];
  int k = work_set;
  int m = ii.last[k];

  uint64_t *restrict p =
    (uint64_t *)kslice_psub_get_address(pos.sq, k) + (thread->begin >> 6);
  uint8_t *restrict const q = kslice_get_address(pos.sq);

  uint64_t last = thread->begin;
  idx_to_sq_init(last, sub, &capt_ii[k]);

  for (uint64_t idx = last, end = thread->end; idx < end; idx += 64) {
    uint64_t w = *p++;
    while (w) {
      uint64_t cur = idx + pop_lsb(&w);
      idx_to_sq_add(cur - last, sub, &capt_ii[k]);
      last = cur;
      Bitboard occ = capt_idx_to_sq(sub, pos.sq, k);
      pos.sq[m] = pos.sq[2];
      mark_unmoves(m, q, occ, pos.sq);
    }
  }
}

// Uncapture the black pawn from the loaded subtable of btm positions
// to wtm positions.
static void predecessors_psub(int s)
{
  for (int r = 0; r < 16; r++) {
    g_pos.sq[0] = KK16Square[s][r][0];
    g_pos.sq[1] = KK16Square[s][r][1];

    if (is_broken(&g_pos))
      continue;

    // Loop through the sets from which a piece uncaptures the pawn.
    for (int k = 0; k < ii.numsets; k++) {
      int m = ii.last[k];
      if ((g_pos.pt[m] >> 3) != WHITE)
        continue;
      work_set = k;
      run_threaded(predecessors_psub_worker, work_capt[k], 0);
    }
  }
}

static void calc_king_captures_pawn_worker(struct ThreadData *thread)
{
  uint32_t sub[MAX_SETS];
  Position pos = g_pos;
  int lower = work_lower;
  int upper = work_upper;

  uint8_t *restrict const p = kslice_get_address(pos.sq);

  idx_to_sq_init(thread->begin, sub, &ii);
  for (uint64_t idx = thread->begin, end = thread->end; idx < end;
      idx++, idx_to_sq_inc(sub, &ii))
  {
    pos.occ = idx_to_sq(sub, pos.sq);
    if (opp_king_attacked(&pos))
      continue;
    if (do_capture(&pos, g_pos.sq[0], g_pos.sq[2], 0, 2)) {
      int v = -probe_wdl(&pos, -2, 2);
      if (v >= lower && v <= upper)
        kslice_bit_set(p, idx);
      undo_capture(&pos, g_pos.sq[0], g_pos.sq[2], 0, 2);
    }
  }
}

static void calc_king_captures_pawn(int s, int lower, int upper)
{
  work_lower = lower;
  work_upper = upper;
  g_pos.stm = WHITE;
  for (int r = 0; r < 16; r++) {
    g_pos.sq[0] = KK16Square[s][r][0];
    g_pos.sq[1] = KK16Square[s][r][1];

    if (is_broken(&g_pos))
      continue;

    if (!(king_attacks(g_pos.sq[0]) & bit(g_pos.sq[2])))
      continue;

    run_threaded(calc_king_captures_pawn_worker, work_g, 0);
  }
}

const char *wdl_name[5] = { "loss", "bloss", "draw", "cwin", "win" };

static void calc_capt(int stm, int wdl)
{
  struct K16SliceIterator iter;
  uint64_t num[16], cnt = 0;

  char capt_name[128], pcapt_name[128], sub_name[128], psub_name[128];
  strcat(strcpy(capt_name , "capt/" ), wdl_name[2 + wdl]);
  strcat(strcpy(pcapt_name, "pcapt/"), wdl_name[2 + wdl]);
  strcat(strcpy(sub_name  , "sub/"  ), wdl_name[2 - wdl]);
  strcat(strcpy(psub_name , "psub/" ), wdl_name[2 - wdl]);

  if (true || sub_cnt[stm ^ 1][2 - wdl]) {
    create_dir(-1, stm, capt_name);

    k16slice_iter_init(&iter, stm);
    int s, s1;
    while (k16slice_iter_next(&iter, &s)) {

      if (true || k16slice_test(s, stm ^ 1, sub_name, -1)) {
        while (k16slice_iter_in(&iter, &s1)) {
          if (stm == WHITE)
            k16slice_clear(s1);
          else
            k16slice_read(s1, s1, BLACK, pcapt_name, -1);
        }

        k16slice_sub_read(s, s, stm ^ 1, sub_name);
        predecessors_sub(stm, s);

        if (stm == WHITE) {
          k16slice_sub_read(s, s, stm, psub_name);
          predecessors_psub(s);
          calc_king_captures_pawn(s, wdl, wdl);
        }
      }

      while (k16slice_iter_out(&iter, &s)) {
        cnt += k16slice_count(s, num);
        k16slice_write(s, s, stm, capt_name, -1, num);
      }
    }
  }

//  g_stats[stm][n] = cnt;
  printf("capt_%s_%c = %lu\n", wdl_name[2 + wdl], "wb"[stm], cnt);
}

// Calculate positions with a capture or pawn move >= wdl (wdl = -1 or 0).
// Removing these positions from potential losses allows us to skip
// captures and pawn moves in check_successors().
static void calc_noloss(int stm, int wdl)
{
  uint64_t num[16];
  struct K16SliceIterator iter;

  const char *noloss = wdl == 0 ? "nobloss" : "noloss";

  create_dir(-1, stm, noloss);

  char pawn_name[128], pcapt_name[128], sub_name[128], psub_name[128];
  strcat(strcpy(pawn_name , "pawn/" ), noloss);
  strcat(strcpy(pcapt_name, "pcapt/"), noloss);
  strcat(strcpy(sub_name  , "sub/"  ), noloss);
  strcat(strcpy(psub_name , "psub/" ), noloss);

  k16slice_iter_init(&iter, stm);
  int s, s1;
  while (k16slice_iter_next(&iter, &s)) {

    while (k16slice_iter_in(&iter, &s1)) {
      k16slice_clear(s1);
      if (stm == WHITE)
        k16slice_clear(s1);
      else {
        k16slice_read(s1, s1, BLACK, pcapt_name, -1);
        k16slice_read(-1, s1, BLACK, pawn_name, -1);
        k16slice_or(s1, -1);
      }
    }

    k16slice_sub_read(s, s, stm ^ 1, sub_name);
    predecessors_sub(stm, s);

    if (stm == WHITE) {
      k16slice_sub_read(s, s, stm, psub_name);
      predecessors_psub(s);
      calc_king_captures_pawn(s, wdl, 2);
    }

    while (k16slice_iter_out(&iter, &s)) {
      // Add illegal positions to avoid having to check legality later.
      // As a bonus, add any faster wins already found to increase efficiency.
      k16slice_read(-1, s, stm, "wins", 0);
      k16slice_or(s, -1);
      k16slice_count(s, num);
      k16slice_write(s, s, stm, noloss, -1, num);
    }
  }
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
  uint32_t sub[MAX_SETS];
  Position pos = g_pos;
  pos.stm = BLACK;

  size_t offset = k16offset(pos.sq);

  idx_to_sq_init(thread->begin, sub, &ii);

  for (uint64_t idx = thread->begin, end = thread->end; idx < end;
      idx++, idx_to_sq_inc(sub, &ii))
  {
    pos.occ = idx_to_sq(sub, pos.sq);
    Bitboard b = pawn_attacks(BLACK, pos.sq[2]) & pos.occ;
    if (!b || opp_king_attacked(&pos))
      continue;
    int v = -3;
    while (b) {
      int sq = pop_lsb(&b);
      int j = piece_idx(&pos, sq);
      if (pos.pt[j] & 0x08) continue;
      if (do_capture(&pos, g_pos.sq[2], sq, 2, j)) {
        v = sq < 8 ? probe_promotions(&pos, v)
                   : max(v, -probe_wdl(&pos, -2, -v));
        undo_capture(&pos, g_pos.sq[2], sq, 2, j);
      }
    }
    if (v > -3)
      kslice_bit_set(k16slice_buf[v + 2] + offset, idx);
  }
}

static void calc_pawn_capts(void)
{
  create_dir(-1, BLACK, "pcapt/loss");
  create_dir(-1, BLACK, "pcapt/bloss");
  create_dir(-1, BLACK, "pcapt/draw");
  create_dir(-1, BLACK, "pcapt/cwin");
  create_dir(-1, BLACK, "pcapt/win");
  create_dir(-1, BLACK, "pcapt/noloss");
  create_dir(-1, BLACK, "pcapt/nobloss");

  for (int s = 0; s < 240; s++) {
    for (int i = 0; i < 5; i++)
      k16slice_clear_addr(k16slice_buf[i]);

    for (int r = 0; r < 16; r++) {
      g_pos.sq[0] = KK16Square[s][r][0];
      g_pos.sq[1] = KK16Square[s][r][1];

      if (is_broken(&g_pos))
        continue;

      if (pawn_attacks(BLACK, g_pos.sq[2]) & bit(g_pos.sq[0]))
        continue;

      run_threaded(calc_pawn_capts_worker, work_g, 0);
    }

    uint64_t cnt_l, cnt_bl, cnt_d, cnt_cw, cnt_w;
    uint64_t num[16];

    cnt_l = k16slice_count_addr(k16slice_buf[0], num);
    k16slice_write_addr(k16slice_buf[0], s, BLACK, "pcapt/loss", -1, num);

    cnt_bl = k16slice_count_addr(k16slice_buf[1], num);
    k16slice_write_addr(k16slice_buf[1], s, BLACK, "pcapt/bloss", -1, num);

    cnt_d = k16slice_count_addr(k16slice_buf[2], num);
    k16slice_write_addr(k16slice_buf[2], s, BLACK, "pcapt/draw", -1, num);

    cnt_cw = k16slice_count_addr(k16slice_buf[3], num);
    k16slice_write_addr(k16slice_buf[3], s, BLACK, "pcapt/cwin", -1, num);

    cnt_w = k16slice_count_addr(k16slice_buf[4], num);
    k16slice_write_addr(k16slice_buf[4], s, BLACK, "pcapt/win", -1, num);

    psub_cnt[0] += cnt_l;
    psub_cnt[1] += cnt_bl;
    psub_cnt[2] += cnt_d;
    psub_cnt[3] += cnt_cw;
    psub_cnt[4] += cnt_w;

    k16slice_or_addr(k16slice_buf[4], k16slice_buf[3]);
    k16slice_or_addr(k16slice_buf[4], k16slice_buf[2]);
    k16slice_write_addr(k16slice_buf[4], s, BLACK, "pcapt/nobloss", -1,
        nullptr);

    k16slice_or_addr(k16slice_buf[4], k16slice_buf[1]);
    k16slice_write_addr(k16slice_buf[4], s, BLACK, "pcapt/noloss", -1, nullptr);
  }
}

static void calc_pawn_prom_worker(struct ThreadData *thread)
{
  uint32_t sub[MAX_SETS];
  Position pos = g_pos;
  pos.stm = WHITE;
  int s = g_pos.sq[2];

  uint8_t *restrict ilg = k16slice_buf[1] + k16offset(pos.sq);
  uint8_t *restrict p[5];
  for (int i = 0; i < 5; i++)
    p[i] = k16slice_buf[2 + i] + k16offset(pos.sq);

  idx_to_sq_init(thread->begin, sub, &ii);
  for (uint64_t idx = thread->begin, end = thread->end; idx < end;
      idx++, idx_to_sq_inc(sub, &ii))
  {
    if (kslice_bit_test(ilg, idx))
      continue;
    pos.occ = idx_to_sq(sub, pos.sq);
    if (!(pos.occ & bit(s - 8))) {
      pos.sq[2] = s - 8;
      pos.occ ^= bit(s) ^ bit(s - 8);
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

static void calc_pawn_push_worker(struct ThreadData *thread)
{
  uint32_t sub[MAX_SETS];
  Position pos = g_pos;
  pos.stm = WHITE;
  int s = g_pos.sq[2];

  uint8_t *restrict ilg = k16slice_buf[1] + k16offset(pos.sq);
  uint8_t *restrict p[5];
  for (int i = 0; i < 5; i++)
    p[i] = k16slice_buf[2 + i] + k16offset(pos.sq);

  idx_to_sq_init(thread->begin, sub, &ii);
  for (uint64_t idx = thread->begin, end = thread->end; idx < end;
      idx++, idx_to_sq_inc(sub, &ii))
  {
    if (kslice_bit_test(ilg, idx))
      continue;
    pos.occ = idx_to_sq(sub, pos.sq);
    if (!(pos.occ & bit(s - 8))) {
      pos.sq[2] = s - 8;
      pos.occ ^= bit(s) ^ bit(s - 8);
      if (!opp_king_attacked(&pos)) {
        uint64_t idx2 = sq_to_idx(pos.sq);
        int v = merged_to_wdl[merged_table[idx2]];
        kslice_bit_set(p[v], idx);
      }
      pos.sq[2] = s;
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
  char str[128];

  if (g_pos.sq[2] - 8 == g_pos.sq[0] || g_pos.sq[2] - 8 == g_pos.sq[1])
    return;

  sprintf(str, "../%c%c/merged/wdl/w/%c%c%c%c", fl(g_pos.sq[2]),
      rk(g_pos.sq[2] - 8), fl(g_pos.sq[0]), rk(g_pos.sq[0]),
      fl(g_pos.sq[1]), rk(g_pos.sq[1]));
  FILE *F = fopen(str, "rb");
  if (!F) {
    fprintf(stderr, "Could not open %s.\n", str);
    exit(EXIT_FAILURE);
  }
  read_data(F, merged_table, kslice_size);
  fclose(F);

  run_threaded(calc_pawn_push_worker, work_g, 0);
}

static void calc_pawn_double_push_worker(struct ThreadData *thread)
{
  uint32_t sub[MAX_SETS];
  Position pos = g_pos;
  pos.stm = WHITE;
  int s = g_pos.sq[2];

  uint8_t *restrict ilg = k16slice_buf[1] + k16offset(pos.sq);
  uint8_t *restrict p[5];
  for (int i = 0; i < 5; i++)
    p[i] = k16slice_buf[2 + i] + k16offset(pos.sq);

  idx_to_sq_init(thread->begin, sub, &ii);
  for (uint64_t idx = thread->begin, end = thread->end; idx < end;
      idx++, idx_to_sq_inc(sub, &ii))
  {
    if (kslice_bit_test(ilg, idx))
      continue;
    Bitboard occ = idx_to_sq(sub, pos.sq);
    if (!(occ & bit(s - 8))) {
      pos.sq[2] = s - 8;
      pos.occ ^= bit(s) ^ bit(s - 8);
      if (!opp_king_attacked(&pos)) {
        uint64_t idx2 = sq_to_idx(pos.sq);
        int v = merged_to_wdl[merged_table[idx2]];
        if (!(occ & bit(s - 16))) {
          pos.sq[2] = s - 16;
          pos.occ ^= bit(s - 8) ^ bit(s - 16);
          if (!opp_king_attacked(&pos)) {
            idx2 = sq_to_idx(pos.sq);
            v = max(v, merged_to_wdl[merged_table2[idx2]]);
          }
        }
        kslice_bit_set(p[v], idx);
      }
      pos.sq[2] = s;
    }
  }
}

static void calc_pawn_double_push(void)
{
  char str[128];

  if (g_pos.sq[2] - 8 == g_pos.sq[0] || g_pos.sq[2] - 8 == g_pos.sq[1])
    return;

  sprintf(str, "../%c%c/merged/wdl/w/%c%c%c%c", fl(g_pos.sq[2]),
      rk(g_pos.sq[2] - 8), fl(g_pos.sq[0]), rk(g_pos.sq[0]),
      fl(g_pos.sq[1]), rk(g_pos.sq[1]));
  FILE *F = fopen(str, "rb");
  if (!F) {
    fprintf(stderr, "Could not open %s.\n", str);
    exit(EXIT_FAILURE);
  }
  read_data(F, merged_table, kslice_size);
  fclose(F);

  if (g_pos.sq[2] - 16 == g_pos.sq[0] || g_pos.sq[2] - 16 == g_pos.sq[1])
    return;

  sprintf(str, "../%c%c/merged/wdl/w/%c%c%c%c", fl(g_pos.sq[2]),
      rk(g_pos.sq[2] - 16), fl(g_pos.sq[0]), rk(g_pos.sq[0]),
      fl(g_pos.sq[1]), rk(g_pos.sq[1]));
  F = fopen(str, "rb");
  if (!F) {
    fprintf(stderr, "Could not open %s.\n", str);
    exit(EXIT_FAILURE);
  }
  read_data(F, merged_table2, kslice_size);
  fclose(F);

  run_threaded(calc_pawn_double_push_worker, work_g, 0);
}

static void predecessors_worker(struct ThreadData *thread)
{
  uint32_t sub[MAX_SETS];
  Position pos = g_pos;
  int stm = pos.stm;

  uint64_t *restrict p = (uint64_t *)kslice_get_address_scratch(pos.sq);
  uint8_t *restrict const q = kslice_get_address(pos.sq);

  p += thread->begin >> 6;
  uint64_t last = thread->begin;
  idx_to_sq_init(last, sub, &ii);
  for (uint64_t idx = last, end = thread->end; idx < end; idx += 64) {
    uint64_t w = *p++;
    while (w) {
      uint64_t cur = idx + pop_lsb(&w);
      idx_to_sq_add(cur - last, sub, &ii);
      last = cur;
      Bitboard occ = idx_to_sq(sub, pos.sq);
      mark_king_unmoves(stm, occ, pos.sq);
      for (int i = 0; pos.pcs[stm][i] >= 0; i++) {
        int j = pos.pcs[stm][i];
        mark_unmoves(j, q, occ, pos.sq);
      }
    }
  }
}

// Calculate stm predecessors of stm^1 positions in scratch.
// TODO: use k16slice_read_count[16] to skip empty kslices
static void predecessors(int stm, int s)
{
  g_pos.stm = stm;
  for (int r = 0; r < 16; r++) {
    g_pos.sq[0] = KK16Square[s][r][0];
    g_pos.sq[1] = KK16Square[s][r][1];

    if (is_broken(&g_pos))
      continue;

    run_threaded(predecessors_worker, work_g, 0);
  }
}

INLINE int get_idx(const uint8_t *restrict sq, int s)
{
  for (int i = 0; ; i++)
    if (sq[i] == s)
      return i;
  unreachable();
}

INLINE bool check_king_moves(int stm, Bitboard occ, const uint8_t *restrict sq)
{
  uint8_t sq2[MAX_PIECES];
  memcpy(sq2, sq, sizeof sq2);

  Bitboard b = king_attacks(sq[stm]) & ~(occ | king_attacks(sq[stm ^ 1]));
  while (b) {
    sq2[stm] = pop_lsb(&b);
    uint8_t *p = kslice_get_address(sq2);
    if (!kslice_bit_test(p, sq_to_idx(sq2)))
      return false;
  } 

  return true;
}

INLINE bool check_moves(int k, const uint8_t *restrict const p, Bitboard occ,
    const uint8_t *restrict sq)
{
  uint8_t sq2[MAX_PIECES];

  Bitboard b = non_king_piece_attacks(g_pos.pt[k], sq[k], occ) & ~occ;
  while (b) {
    memcpy(sq2, sq, sizeof sq2);
    sq2[k] = pop_lsb(&b);
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
  uint64_t cnt = 0;

  uint64_t *restrict p = (uint64_t *)kslice_get_address_scratch(pos.sq);
  uint8_t *restrict const q = kslice_get_address(pos.sq);

  p += thread->begin >> 6;
  uint64_t last = thread->begin;
  idx_to_sq_init(last, sub, &ii);

  for (uint64_t idx = last, end = thread->end; idx < end; idx += 64, p++) {
    uint64_t todo = *p;
    if (!todo) continue;
    uint64_t kept = todo;
    while (todo) {
      unsigned bt = pop_lsb(&todo);
      uint64_t cur = idx + bt;
      idx_to_sq_add(cur - last, sub, &ii);
      last = cur;
      Bitboard occ = pos.occ = idx_to_sq(sub, pos.sq);
      for (int i = 0; pos.pcs[stm][i] >= 0; i++) {
        int j = pos.pcs[stm][i];
        if (!check_moves(j, q, occ, pos.sq))
          goto clear_bit;
      }
      if (check_king_moves(stm, occ, pos.sq))
        continue;
clear_bit:
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
  g_pos.stm = stm;
  for (int r = 0; r < 16; r++) {
    num[r] = 0;

    g_pos.sq[0] = KK16Square[s][r][0];
    g_pos.sq[1] = KK16Square[s][r][1];

    if (is_broken(&g_pos))
      continue;

    for (int t = 0; t < g_num_threads; t++)
      g_thread_data[t].cnt = 0;

    run_threaded(check_successors_worker, work_g, 0);

    for (int t = 0; t < g_num_threads; t++)
      num[r] += g_thread_data[t].cnt;
  }

  uint64_t cnt = 0;
  for (int r = 0; r < 16; r++)
    cnt += num[r];

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

  uint8_t *restrict const p = k16slice_buf[stm] + k16offset(pos.sq);

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

  uint64_t *restrict p0 = (uint64_t *)(k16slice_buf[0] + k16offset(pos.sq));
  uint64_t *restrict p1 = (uint64_t *)(k16slice_buf[1] + k16offset(pos.sq));
  uint64_t *restrict q0 = (uint64_t *)(k16slice_buf[2] + k16offset(pos.sq));
  uint64_t *restrict q1 = (uint64_t *)(k16slice_buf[3] + k16offset(pos.sq));

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
      uint64_t cur = idx + bt;
      idx_to_sq_add(cur - last, sub, &ii);
      last = cur;
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

// Calc illegal, mate (L0) and pawn-push positions.
static void calc_illegal_and_mate_and_pawn_push(void)
{
  uint64_t broken_w = 0, broken_b = 0, loss0_w = 0, loss0_b = 0, num[16];

  create_dir(0, WHITE, "wins");
  create_dir(0, BLACK, "wins");
  create_dir(0, WHITE, "L");
  create_dir(0, BLACK, "L");
  create_dir(-1, BLACK, "pawn/loss");
  create_dir(-1, BLACK, "pawn/bloss");
  create_dir(-1, BLACK, "pawn/draw");
  create_dir(-1, BLACK, "pawn/cwin");
  create_dir(-1, BLACK, "pawn/win");
  create_dir(-1, BLACK, "pawn/noloss");
  create_dir(-1, BLACK, "pawn/nobloss");

  for (int s = 0; s < 240; s++) {
    k16slice_clear_addr(k16slice_buf[0]); // wtm illegal
    k16slice_clear_addr(k16slice_buf[1]); // btm illegal

    for (int r = 0; r < 16; r++) {
      g_pos.sq[0] = KK16Square[s][r][0];
      g_pos.sq[1] = KK16Square[s][r][1];

      if (is_broken(&g_pos))
        continue;

      bool btm_illegal = pawn_attacks(BLACK, g_pos.sq[2]) & bit(g_pos.sq[0]);
      if (btm_illegal)
        kslice_set_addr(k16slice_buf[1] + k16offset(g_pos.sq));

      for (int k = 0; k < ii.numsets; k++) {
        if (btm_illegal && (g_pos.pt[ii.first[k]] & 0x08))
          continue;
        work_set = k;
        run_threaded(calc_illegal_worker, work_capt[k], 0);
      }
    }

    broken_w += k16slice_count_addr(k16slice_buf[0], num);
    k16slice_write_addr(k16slice_buf[0], s, WHITE, "wins", 0, num);

    broken_b += k16slice_count_addr(k16slice_buf[1], num);
    k16slice_write_addr(k16slice_buf[1], s, BLACK, "wins", 0, num);

    k16slice_clear_addr(k16slice_buf[2]); // wtm mate
    k16slice_clear_addr(k16slice_buf[3]); // btm mate

    for (int r = 0; r < 16; r++) {
      g_pos.sq[0] = KK16Square[s][r][0];
      g_pos.sq[1] = KK16Square[s][r][1];

      if (is_broken(&g_pos))
        continue;

      run_threaded(calc_mate_worker, work_g, 0);
    }

    loss0_w += k16slice_count_addr(k16slice_buf[2], num);
    k16slice_write_addr(k16slice_buf[2], s, WHITE, "L", 0, num);

    loss0_b += k16slice_count_addr(k16slice_buf[3], num);
    k16slice_write_addr(k16slice_buf[3], s, BLACK, "L", 0, num);

    for (int i = 0; i < 5; i++)
      k16slice_clear_addr(k16slice_buf[2 + i]);

    for (int r = 0; r < 16; r++) {
      g_pos.sq[0] = KK16Square[s][r][0];
      g_pos.sq[1] = KK16Square[s][r][1];

      if (is_broken(&g_pos))
        continue;

      if (g_pos.sq[2] < 16)
        run_threaded(calc_pawn_prom_worker, work_g, 0);
      else if (g_pos.sq[2] < 48)
        calc_pawn_push();
      else
        calc_pawn_double_push();
    }

    pawn_cnt[0] += k16slice_count_addr(k16slice_buf[2], num);
    k16slice_write_addr(k16slice_buf[2], s, BLACK, "pawn/loss", -1, num);
    pawn_cnt[1] += k16slice_count_addr(k16slice_buf[3], num);
    k16slice_write_addr(k16slice_buf[3], s, BLACK, "pawn/bloss", -1, num);
    pawn_cnt[2] += k16slice_count_addr(k16slice_buf[4], num);
    k16slice_write_addr(k16slice_buf[4], s, BLACK, "pawn/draw", -1, num);
    pawn_cnt[3] += k16slice_count_addr(k16slice_buf[5], num);
    k16slice_write_addr(k16slice_buf[5], s, BLACK, "pawn/cwin", -1, num);
    pawn_cnt[4] += k16slice_count_addr(k16slice_buf[6], num);
    k16slice_write_addr(k16slice_buf[6], s, BLACK, "pawn/win", -1, num);

    k16slice_or_addr(k16slice_buf[6], k16slice_buf[5]);
    k16slice_or_addr(k16slice_buf[6], k16slice_buf[4]);
    k16slice_write_addr(k16slice_buf[6], s, BLACK, "pawn/nobloss", -1, nullptr);
    k16slice_or_addr(k16slice_buf[6], k16slice_buf[3]);
    k16slice_write_addr(k16slice_buf[6], s, BLACK, "pawn/noloss", -1, nullptr);
  }

  g_stats[WHITE][0] = broken_w;
  g_stats[BLACK][0] = broken_b;
  g_stats[WHITE][MAX_STATS - 1] = loss0_w;
  g_stats[BLACK][MAX_STATS - 1] = loss0_b;

  printf("broken_w = %lu\n", broken_w);
  printf("broken_b = %lu\n", broken_b);
  printf("l0_w = %lu\n", loss0_w);
  printf("l0_b = %lu\n", loss0_b);
  for (int i = 0; i < 5; i++)
    printf("pawn_cnt[%d] = %lu\n", i, pawn_cnt[i]);
}

// Calculate stm losses in n from stm^1 wins in n-1 (n > 1) or
// from stm^1 wins in sub tables reached through captures (n == 1).
static bool calc_L(int stm, int n, bool more_l)
{
  struct K16SliceIterator iter;
  uint64_t cnt = 0, num[16];

  create_dir(n, stm, "X");

  // Calculate potential losses in n = predecessors(W(n-1))
  k16slice_iter_init(&iter, stm);
  int s, s1;
  while (k16slice_iter_next(&iter, &s)) {
    bool pred = more_l && k16slice_test(s, stm ^ 1, "W", n - 1);

    if (pred || n == 1 || n == DRAW_RULE + 1) {
      while (k16slice_iter_in(&iter, &s1))
        k16slice_clear(s1);

      if (pred) {
        k16slice_read(-1, s, stm ^ 1, "W", n - 1);
        predecessors(stm, s);
      }
    }

    while (k16slice_iter_out(&iter, &s)) {
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
      k16slice_and_not(s, -1);
      k16slice_clear_tail(s);
      k16slice_write(s, s, stm, "X", n, nullptr); // FIXME
    }
  }

  create_dir(n, stm, "L");

  // Verify potential losses.
  k16slice_iter_init(&iter, stm);
  while (k16slice_iter_next(&iter, &s)) {

    if (k16slice_test(s, stm , "X", n)) {
      while (k16slice_iter_in(&iter, &s1))
        k16slice_read(s1, s1, stm ^ 1, "wins", 0);

      k16slice_read(-1, s, stm, "X", n);
      cnt += check_successors(stm, s, num);
      k16slice_write(-1, s, stm, "L", n, num);
      if (0 && stm == BLACK && n == 4) {
        list_positions(s, k16slice_buf[11]);
//        printf("s = %d, cnt = %lu\n", s, cnt);
      }
    }
    k16slice_delete(s, stm, "X", n);

    while (k16slice_iter_out(&iter, &s));
  }

  g_stats[stm][MAX_STATS - 1 - n] = cnt;
  printf("l%d_%c = %lu\n", n, "wb"[stm], cnt);
  return cnt != 0;
}

static bool calc_W(int stm, int n, bool more_w)
{
  struct K16SliceIterator iter;
  uint64_t cnt = 0, cnt_w = 0, cnt_pw = 0, num[16];

  create_dir(n, stm, "W");

  // Calculate wins in n = predecessors(L(n-1))
  k16slice_iter_init(&iter, stm);
  int s, s1;
  while (k16slice_iter_next(&iter, &s)) {
    bool pred = more_w && k16slice_test(s, stm ^ 1, "L", n - 1);

    if (pred || n == 1 || n == DRAW_RULE + 1) {
      while (k16slice_iter_in(&iter, &s1)) {
        if (n == 1) {
          // Add any wins by capture or pawn push.
          k16slice_read(s1, s1, stm, "capt/win", -1);
          // Remove illegal positions to count wins by capture.
          k16slice_read(-1, s1, stm, "wins", 0);
          k16slice_and_not(s1, -1);
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
          k16slice_read(-1, s1, stm, "wins", 0);
          k16slice_and_not(s1, -1);
          cnt_w += k16slice_count(s1, num);
        }
        else
          k16slice_clear(s1);
      }

      if (pred) {
        k16slice_read(-1, s, stm ^ 1, "L", n - 1);
        predecessors(stm, s);
      }
    }

    while (k16slice_iter_out(&iter, &s)) {
      // Remove illegal positions and known faster wins.
      k16slice_read(-1, s, stm, "wins", 0);
      k16slice_and_not(s, -1);
      cnt += k16slice_count(s, num);
      k16slice_write(s, s, stm, "W", n, num);

      if (0 && stm == WHITE && n == 1) {
        list_positions(s, k16slice_get_address(s));
//        printf("s = %d, cnt = %lu\n", s, cnt);
      }

      // TODO: do not update "wins" each iteration but work with small deltas.
      k16slice_or(-1, s);
      k16slice_write(-1, s, stm, "wins", 0, nullptr);
    }
  }

  if (n == 1) {
    printf("capt_win_%c = %lu\n", "wb"[stm], cnt_w);
    g_stats[stm][1] = cnt_w;
    if (stm == BLACK) {
      g_stats[stm][2] = cnt_pw - cnt_w;
      printf("pawn_win_%c = %lu\n", "wb"[stm], cnt_pw - cnt_w);
    }
  }
  else if (n == DRAW_RULE + 1) {
    printf("capt_cwin_%c = %lu\n", "wb"[stm], cnt_w);
    g_stats[stm][DRAW_RULE + 3] = cnt_w;
    if (stm == BLACK) {
      printf("pawn_cwin_%c = %lu\n", "wb"[stm], cnt_pw - cnt_w);
      g_stats[stm][DRAW_RULE + 4] = cnt_pw - cnt_w;
    }
  }

  g_stats[stm][2 + n + 2 * (n > DRAW_RULE)] = cnt;
  printf("w%d_%c = %lu\n", n, "wb"[stm], cnt);
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
  printf("Generating table for pawn slice %c%c.\n", 'a' + (g_pos.sq[2] & 7),
      '1' + (g_pos.sq[2] >> 3));

  for (int i = 0; i < 7; i++)
    k16slice_buf[i] = alloc_k16slice();

  if (g_pos.sq[2] >= 16) {
    merged_table = alloc_huge(kslice_size);
    if (!merged_table)
      out_of_mem();
    if (g_pos.sq[2] >= 48) {
      merged_table2 = alloc_huge(kslice_size);
      if (!merged_table2)
        out_of_mem();
    }
  }

  for (int i = 0; i < 5; i++)
    pawn_cnt[i] = psub_cnt[i] = sub_cnt[WHITE][i] = sub_cnt[BLACK][i] = 0;

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

  kslice_free_buffers();
}
