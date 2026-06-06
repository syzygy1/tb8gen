/*
  Copyright (c) 2026 Ronald de Man

  This file is distributed under the terms of the GNU GPL, version 2.
*/

#include <fcntl.h>
#include <limits.h>
#include <stdatomic.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>

#include "defs.h"
#include "index.h"
#include "generate.h"
#include "kslice.h"
#include "merge.h"
#include "movegen.h"
#include "stats.h"
#include "tb8gen.h"
#include "threads.h"
#include "types.h"
#include "util.h"

struct MergeInfo mi;

static void *merge_table, *merge_w;
static int merge_n;
static int work_set;
static bool include_wins, include_losses;
static struct Work work_g_merge_dynamic[2], work_g_merge_static[2];
static struct Work work_capt_merge_dynamic[MAX_SETS];
static int merge_num_active_threads;
static int merge_active_threads[MAX_THREADS];

static constexpr uint64_t MERGE_MIN_DYNAMIC_CHUNK = 1ULL << 18;
static constexpr int MERGE_DYNAMIC_FACTOR = 4;

alignas(64) static uint64_t thread_stats[MAX_THREADS][MAX_STATS];

static void set_merge_active_threads(struct Work *work)
{
  if (work->units <= 1) {
    merge_num_active_threads = 1;
    merge_active_threads[0] = g_num_threads - 1;
  } else {
    merge_num_active_threads = g_num_threads;
    for (int t = 0; t < g_num_threads; t++)
      merge_active_threads[t] = t;
  }
}

void init_merge_work(void)
{
  work_init(&work_g_merge_dynamic[0], kslice_sizes[0], 0x1ff, WORK_DYNAMIC,
      MERGE_DYNAMIC_FACTOR, MERGE_MIN_DYNAMIC_CHUNK);
  work_init(&work_g_merge_dynamic[1], kslice_sizes[1], 0x1ff, WORK_DYNAMIC,
      MERGE_DYNAMIC_FACTOR, MERGE_MIN_DYNAMIC_CHUNK);
  work_init(&work_g_merge_static[0], kslice_sizes[0], 0x1ff, WORK_STATIC, 1,
      MERGE_MIN_DYNAMIC_CHUNK);
  work_init(&work_g_merge_static[1], kslice_sizes[1], 0x1ff, WORK_STATIC, 1,
      MERGE_MIN_DYNAMIC_CHUNK);
  for (int k = 0; k < ri.numsets; k++)
    work_init(&work_capt_merge_dynamic[k], capt_ri[k].sizes[0], 0x1ff,
        WORK_DYNAMIC, MERGE_DYNAMIC_FACTOR, MERGE_MIN_DYNAMIC_CHUNK);
}

static _Atomic uint64_t found_idx;
static _Atomic bool found_capt_bloss;

static void find_position_worker(struct ThreadData *thread)
{
  uint64_t cur = atomic_load_explicit(&found_idx, memory_order_relaxed);
  if (thread->begin >= cur)
    return;

  uint64_t *restrict p =
    (uint64_t *)kslice_get_address(-1) + (thread->begin >> 6);

  uint64_t idx, end;
  for (idx = thread->begin, end = thread->end; idx < end; idx += 64) {
    uint64_t w = *p++;
    if (w) {
      idx += pop_lsb(&w);
      break;
    }
  }
  if (idx >= end)
    return;

  uint64_t old = atomic_load_explicit(&found_idx, memory_order_relaxed);
  while (idx < old) {
    if (atomic_compare_exchange_weak_explicit(&found_idx, &old, idx,
          memory_order_relaxed, memory_order_relaxed))
      break;
  }
}

static void find_position(int s, int stm, bool loss, bool cursed)
{
  atomic_store_explicit(&found_idx, UINT64_MAX, memory_order_relaxed);

  run_threaded(find_position_worker, &work_g_merge_dynamic[s >= 441]);

  uint64_t idx = atomic_load_explicit(&found_idx, memory_order_relaxed);
  if (idx >= kslice_sizes[s >= 441])
    return;

  Position pos = g_pos;
  pos.stm = stm ^ loss;
  if (s < 441) {
    struct IdxState is;
    idx_state_init(&is, idx, pos.sq, &ri);
    idx_state_to_sq(&is, pos.sq, &ri);
  } else {
    Bitboard occ = bit(pos.sq[0]) | bit(pos.sq[1]);
    unrank_reflection(idx, pos.sq, occ, &ri);
  }
  pos_to_fen(&pos, mf.fen[stm][cursed], false);
  mf.found[stm][cursed] = true;

  FILE *F = file_open_write("maxfens");
  file_write(&mf, sizeof mf, F);
  fclose(F);
  file_rename("maxfens");
}

#define T u8
#define MAX 256
#include "merge_tmpl.c"
#undef MAX
#undef T

#define T u16
#define MAX MAX_STATS
#include "merge_tmpl.c"
#undef MAX
#undef T

void delete_bitmaps(int stm, int s)
{
  for (int n = 0; n < max_iteration; n++) {
    kslice_delete(s, stm, "L", n);
    kslice_delete(s, stm, "W", n);
    kslice_delete(s, stm, "wins", n);
  }
  kslice_delete(s, stm, "capt/win", -1);
  kslice_delete(s, stm, "capt/cwin", -1);
  kslice_delete(s, stm, "capt/draw", -1);
  kslice_delete(s, stm, "capt/bloss", -1);
}

void merge(int stm)
{
  char str[64];
  sprintf(str, "merge_info.%c", "wb"[stm]);
  if (file_exists(str))
    return;

  char phase[32];
  snprintf(phase, sizeof phase, "merging %s slices",
      stm == WHITE ? "white" : "black");

  uint64_t *stats = g_stats[stm];

  create_dir(-1, stm, "stats");
  if (    (one_sided && one_sided_stm == stm)
      || (!one_sided && (!symmetric || stm == WHITE)))
    create_dir(-1, stm, "merged/dtz");
  if (!symmetric || stm == WHITE)
    create_dir(-1, stm, "merged/wdl");

  // Determine whether we need to store wins, losses or both in the merged
  // in-RAM table.
  bool wins, losses;
  if (one_sided)
    wins = losses = stm == one_sided_stm;
  else {
    wins = wins_only;
    losses = !wins_only;
  }

  // Count the number of distinct values to see if we can fit them all
  // in one byte.
  int win_vals = 0, cwin_vals = 0, bloss_vals = 0, loss_vals = 0;
  for (int i = 2; i <= DRAW_RULE + 1; i++)
    win_vals += (stats[i] != 0);
  for (int i = DRAW_RULE + 3; i < MAX_STATS / 2; i++)
    cwin_vals += (stats[i] != 0);
  for (int i = 0; i <= DRAW_RULE; i++)
    loss_vals += (stats[MAX_STATS - 1 - i] != 0);
  for (int i = DRAW_RULE + 1; i < MAX_STATS / 2 - 2; i++)
    bloss_vals += (stats[MAX_STATS - 1 - i] != 0);

  mi.v_wdl[0] = loss_vals;
  mi.v_wdl[1] = bloss_vals;
  mi.v_wdl[2] = stats[MAX_STATS / 2 + 1];
  mi.v_wdl[3] = cwin_vals;
  mi.v_wdl[4] = win_vals;

  bool dc[4] = {
    sub_cnt[stm ^ 1][3], stats[MAX_STATS / 2], stats[DRAW_RULE + 2], true
  };

  int i, j;
  for (i = 0; i < 4; i++)
    if (dc[i]) break;
  for (j = 0; j < 5; j++)
    if (mi.v_wdl[j]) break;
  if (j > i + 1)
    mi.v_wdl[0] = true;

  int special = 1 + (stats[1] != 0)
                  + (stats[DRAW_RULE + 2] != 0)
                  + (stats[MAX_STATS / 2] != 0)
                  + (stats[MAX_STATS / 2 + 1] != 0);

  int wins_red = (win_vals != 0) + (cwin_vals != 0);
  int losses_red = (loss_vals != 0) + (bloss_vals != 0);

  int tot_vals =  special + (wins ? win_vals + cwin_vals : wins_red)
                + (losses ? loss_vals + bloss_vals : losses_red);

  printf("tot_vals = %d (%d)\n", tot_vals,
      special + win_vals + cwin_vals + loss_vals + bloss_vals);

  include_wins = wins;
  include_losses = losses;

  mi.wide = tot_vals > 256;

  if (!mi.wide) {
    // One byte suffices.

    int n = init_merge_value_map_u8(stats);
    if (n > 255) {
      fprintf(stderr, "Internal error.\n");
      exit(EXIT_FAILURE);
    }

    merge_table = alloc_huge(sizeof(u8) * kslice_size);
    if (!merge_table)
      out_of_mem();

    for (int s = 0; s < 462; s++) {
      show_progress(phase, s, 462, false);
      merge_bitmaps_u8(stm, s);
      delete_bitmaps(stm, s);
    }
    show_progress(phase, 462, 462, true);

    free(merge_table);

  } else {

    init_merge_value_map_u16(stats);

    merge_table = alloc_huge(sizeof(u16) * kslice_size);
    if (!merge_table)
      out_of_mem();

    for (int s = 0; s < 462; s++) {
      show_progress(phase, s, 462, false);
      merge_bitmaps_u16(stm, s);
      delete_bitmaps(stm, s);
    }
    show_progress(phase, 462, 462, true);

    free(merge_table);

  }

  FILE *F = file_open_write(str);
  write_data(F, &mi, sizeof mi);
  fclose(F);
  file_rename(str);
}
