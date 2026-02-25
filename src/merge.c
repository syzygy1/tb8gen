/*
  Copyright (c) 2026 Ronald de Man

  This file is distributed under the terms of the GNU GPL, version 2.
*/

#include <stdint.h>
#include <stdio.h>

#include "defs.h"
#include "index.h"
#include "generate.h"
#include "kslice.h"
#include "movegen.h"
#include "tb8gen.h"
#include "threads.h"
#include "types.h"

static uint8_t *restrict merge_table;
static int merge_n;
static int work_set, work_slice;

static void merge_draw_worker(struct ThreadData *thread)
{
  uint8_t *restrict q = merge_table;
  memset(q + thread->begin, 129, thread->end - thread->begin);
}

static void merge_worker(struct ThreadData *thread)
{
  int n = merge_n;
  uint64_t *restrict p = (uint64_t *)kslice_get_address(-1);
  uint8_t *restrict q = merge_table;
  p += thread->begin >> 6;
  for (uint64_t idx = thread->begin, end = thread->end; idx < end; idx += 64) {
    uint64_t w = *p++;
    if (!w) continue;
    while (w)
      q[idx + pop_lsb(&w)] = n;
  }
}

static void merge_capt_worker(struct ThreadData *thread)
{
  int n = merge_n;
  uint64_t *restrict p = (uint64_t *)kslice_get_address(-1);
  uint8_t *restrict q = merge_table;
  p += thread->begin >> 6;
  for (uint64_t idx = thread->begin, end = thread->end; idx < end; idx += 64) {
    uint64_t w = *p++;
    if (!w) continue;
    while (w) {
      unsigned bt = pop_lsb(&w);
      q[idx + bt] = min(q[idx + bt], n);
    }
  }
}

INLINE void merge_mark_unmoves(int k, uint8_t *restrict const p, Bitboard occ,
    uint8_t *restrict sq)
{
  uint8_t sq2[MAX_PIECES];
  Bitboard b = non_king_piece_moves(g_pos.pt[k], sq[k], occ);
  while (b) {
    sq[k] = pop_lsb(&b);
    for (int i = 0; i < MAX_PIECES; i++)
      sq2[i] = sq[i];
    uint64_t idx = sq_to_idx(sq2);
    p[idx] = 0;
  }
}

static void merge_illegal_worker(struct ThreadData *thread)
{
  uint32_t sub[MAX_SETS];
  Position pos = g_pos;
  int k = work_set;
  int m = ii.last[k];
  int stm = pos.stm;
  int king_sq = pos.sq[stm ^ 1];

  uint8_t *restrict const p = merge_table;

  idx_to_sq_init(thread->begin, sub, &capt_ii[k]);

  for (uint64_t idx = thread->begin, end = thread->end; idx < end;
      idx++, idx_to_sq_inc(sub, &capt_ii[k]))
  {
    Bitboard occ = capt_idx_to_sq(sub, pos.sq, k);
    pos.sq[m] = king_sq;
    merge_mark_unmoves(m, p, occ, pos.sq);
  }
}

static uint64_t stats[256];
alignas(64) static uint64_t thread_stats[32][256];

static void merge_statistics_worker(struct ThreadData *thread)
{
  int s = work_slice;
  int t = thread->thread_id;

  uint8_t *restrict const p = merge_table;

  if (s < 441) {
    for (uint64_t idx = thread->begin, end = thread->end; idx < end; idx++)
      thread_stats[t][p[idx]] += 2;
  } else {
    uint32_t sub[MAX_SETS];
    Position pos = g_pos;
    idx_to_sq_init(thread->begin, sub, &ii);

    for (uint64_t idx = thread->begin, end = thread->end; idx < end;
        idx++, idx_to_sq_inc(sub, &ii))
    {
      idx_to_sq(sub, pos.sq);
      mirror_diagonal(pos.sq);
      uint64_t idx2 = sq_to_idx(pos.sq);
      thread_stats[t][p[idx]] += idx == idx2 ? 2 : 1;
    }
  }
}

// ILLEGAL     = 0
// CAPT_WIN    = 1
// wins          2..101 -> win-in-1 to win-in-100
// CAPT_CWIN   = 102
// wins          103-127  -> win-in-101 to-win-in-125
// CAPT_DRAW   = 128
// unknown       129 -> draw
// CAPT_CLOSS  = x
// losses        130..154 -> loss-in-125 to loss-in-101
// losses        155..254 -> loss-in-100 to loss-in-1
// losses        255 -> mate

// FIXME: check that max_iteration <= 126 (-> max dtz = 125)
static void merge_slices(int stm, int s)
{
  // DRAW
  run_threaded(merge_draw_worker, work_g, 0);

  g_pos.stm = stm;
  g_pos.sq[0] = KKSquare[s][0];
  g_pos.sq[1] = KKSquare[s][1];

  // ILLEGAL
  for (int k = 0; k < ii.numsets; k++) {
    if ((g_pos.pt[ii.first[k]] & 8) != stm)
      continue;
    run_threaded(merge_illegal_worker, work_capt[k], 0);
  }

  // losses
  for (int n = 0; n < max_iteration; n++)
    if (kslice_test(s, stm, "L", n)) {
      kslice_read(-1, s, stm, "L", n);
      merge_n = 255 - n;
      run_threaded(merge_worker, work_g, 0);
    }

  // wins
  for (int n = 1; n < max_iteration; n++)
    if (kslice_test(s, stm, "W", n)) {
      kslice_read(-1, s, stm, "W", n);
      merge_n = n <= DRAW_RULE ? 1 + n : 2 + n;
      run_threaded(merge_worker, work_g, 0);
    }

  // CAPT_WIN
  if (sub_cnt[stm ^ 1][0]) {
    kslice_read(-1, s, stm, "capt/win", -1);
    merge_n = 1;
    run_threaded(merge_capt_worker, work_g, 0);
  }

  // CAPT_CWIN
  if (sub_cnt[stm ^ 1][1]) {
    kslice_read(-1, s, stm, "capt/cwin", -1);
    merge_n = 102;
    run_threaded(merge_capt_worker, work_g, 0);
  }

  // CAPT_DRAW
  if (sub_cnt[stm ^ 1][2]) {
    kslice_read(-1, s, stm, "capt/draw", -1);
    merge_n = 128;
    run_threaded(merge_capt_worker, work_g, 0);
  }

  memset(stats, 0, 256 * sizeof(uint64_t));
  for (int t = 0; t < g_num_threads; t++)
    memset(thread_stats[t], 0, 256 * sizeof(uint64_t));

  work_slice = s;
  run_threaded(merge_statistics_worker, work_g, 0);

  for (int t = 0; t < g_num_threads; t++)
    for (int i = 0; i < 256; i++)
      stats[i] += thread_stats[t][i];

  char str[128];
  create_name(str, s, stm, "stats", -1);
  FILE *F = fopen(str, "wb");
  file_write(stats, sizeof(stats), F);
  fclose(F);

#if 0
  // CAPT_BLOSS
  if (sub_cnt[stm ^ 1][3]) {
    kslice_read(-1, s, stm, "capt/bloss", -1);
  }
#endif

  // merge in the following order:
  // - CAPT_WIN -> do not overrule ILLEGAL
  // - CAPT_CWIN -> do not overrule wins
  // - CAPT_DRAW -> do not overrule cursed wins
  // - CAPT_BLOSS -> difficult!
  //   - for wdl: priority over all other bloss
  //   - for dtz: ignore -> keep dtz blessed loss value
  //   -> 2 versions of each bloss value?
  //      or merge later? i.e. first output dtz, then merge and output wdl.
}

void collect_stats(int stm)
{
  char str[128];
  uint64_t tmp[256];
  int i;

  memset(stats, 0, sizeof(stats));

  for (int s = 0; s < 462; s++) {
    create_name(str, s, stm, "stats", -1);
    FILE *F = fopen(str, "rb");
    file_read(tmp, sizeof(tmp), F);
    fclose(F);
    for (int i = 0; i < 256; i++)
      stats[i] += tmp[i];
  }
  for (int i = 0; i < 256; i++)
    stats[i] >>= 1;

  printf("\n%s to move:\n\n", stm == WHITE ? "White" : "Black");

  if (stats[1] + stats[2])
    printf("%lu (%lu) positions win in %d ply.\n", stats[1] + stats[2],
        stats[1], 1);
  for (i = 2; i <= DRAW_RULE; i++)
    if (stats[1 + i])
      printf("%lu positions win in %d ply.\n", stats[1 + i], i);
  if (stats[1 + i] + stats[2 + i])
    printf("%lu (%lu) positions win in %d ply.\n",
        (stats[1 + i] + stats[2 + i]) >> 1, stats[1 + i], i);
  for (i++; i < 126; i++)
    if (stats[2 + i])
      printf("%lu positions win in %d ply.\n", stats[1 + i], i);
  printf("\n");

  uint64_t tot = 0;
  for (i = 1; i <= 1 + DRAW_RULE; i++)
    tot += stats[i];
  printf("%lu positions are wins.\n", tot);
  tot = 0;
  for (; i < 128; i++)
    tot += stats[i];
  if (tot)
    printf("%lu positions are cursed wins.\n", tot);
  if (stats[128] + stats[129])
    printf("%lu (%lu) positions are draws.\n", stats[128] + stats[129],
        stats[128]);
  tot = 0;
  for (i = DRAW_RULE + 1; i < 126; i++)
    tot += stats[255 - i];
  if (tot)
    printf("%lu positions are blessed losses.\n", tot);
  tot = 0;
  for (i = 0; i <= DRAW_RULE; i++)
    tot += stats[255 - i];
  printf("%lu positions are losses.\n\n", tot);

  for (i = 0; i < 126; i++)
    if (stats[255 - i])
      printf("%lu positions lose in %d ply.\n", stats[255 - i], i);
}

void merge(void)
{
  merge_table = alloc_huge(kslice_size);
  if (!merge_table)
    out_of_mem();

  create_dir(-1, WHITE, "stats");
  create_dir(-1, BLACK, "stats");
  for (int s = 0; s < 462; s++) {
    merge_slices(WHITE, s);
    merge_slices(BLACK, s);
  }
}
