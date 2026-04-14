/*
  Copyright (c) 2026 Ronald de Man

  This file is distributed under the terms of the GNU GPL, version 2.
*/

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "defs.h"
#include "indexp.h"
#include "generatep.h"
#include "kslicep.h"
#include "merge.h"
#include "movegen.h"
#include "tb8genp.h"
#include "threads.h"
#include "types.h"

struct MergeInfo mi;

static void *merge_table;
static int merge_n, merge_r;
static void *merge_w;
static int work_set;
static bool include_wins, include_losses;

alignas(64) static uint64_t thread_stats[MAX_THREADS][MAX_STATS];

static void stat_count(uint64_t stats[16][MAX_STATS], int n)
{
  for (int r = 0; r < 16; r++)
    stats[r][n] = k16slice_read_count[r];
}

#define T u8
#define MAX 256
#include "mergep_tmpl.c"
#undef MAX
#undef T

#define T u16
#define MAX MAX_STATS
#include "mergep_tmpl.c"
#undef MAX
#undef T

void merge(int stm)
{
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
  for (int i = 3; i <= DRAW_RULE + 2; i++)
    win_vals += (stats[i] != 0);
  for (int i = DRAW_RULE + 5; i < MAX_STATS / 2 + 1; i++)
    cwin_vals += (stats[i] != 0);
  for (int i = 0; i <= DRAW_RULE; i++)
    loss_vals += (stats[MAX_STATS - 1 - i] != 0);
  for (int i = DRAW_RULE + 1; i < MAX_STATS / 2 - 3; i++)
    bloss_vals += (stats[MAX_STATS - 1 - i] != 0);

  mi.v_wdl[0] = loss_vals;
  mi.v_wdl[1] = bloss_vals;
  mi.v_wdl[2] = stats[MAX_STATS / 2 + 2];
  mi.v_wdl[3] = stats[DRAW_RULE + 4] || cwin_vals;
  mi.v_wdl[4] = stats[2] || win_vals;

  bool dc[4] = {
    capt_cnt[stm][1], capt_cnt[stm][2], stats[DRAW_RULE + 3], true
  };

  int i, j;
  for (i = 0; i < 4; i++)
    if (dc[i]) break;
  for (j = 0; j < 5; j++)
    if (mi.v_wdl[j]) break;
  if (j > i + 1)
    mi.v_wdl[0] = true;

  int special = 1 + (stats[1] != 0) + (stats[2] != 0)
                  + (stats[DRAW_RULE + 3] != 0) + (stats[DRAW_RULE + 4] != 0)
                  + (capt_cnt[stm][2] != 0)
                  + (stats[MAX_STATS / 2 + 2] != 0);

  int wins_red = (win_vals != 0) + (cwin_vals != 0);
  int losses_red = (loss_vals != 0) + (bloss_vals != 0);

  int tot_vals =  special + (wins ? win_vals + cwin_vals : wins_red)
                + (losses ? loss_vals + bloss_vals : losses_red);

  printf("tot_vals = %d (%d)\n", tot_vals,
      special + win_vals + cwin_vals + loss_vals + bloss_vals);

  mi.wide = tot_vals > 256;

  if (!mi.wide) {
    // One byte suffices.

    // Include more if it fits. This slightly speeds up counting statistics.
    include_wins = wins;
    include_losses = losses;
#if 0
    if (special + win_vals + cwin_vals + bloss_vals + loss_vals <= 256)
      include_wins = include_losses = true;
    else if (!wins && !losses && special + win_vals + cwin_vals <= 256)
      include_wins = true;
    else if (!wins && !losses && special + loss_vals + bloss_vals <= 256)
      include_losses = true;
#endif

    // Create the corresponding mapping from u16 to u8.
    int n = 0;
    mi.v_u8[0] = n;
    n += (stats[1] != 0);
    mi.v_u8[1] = n;
    n += (stats[2] != 0);
    mi.v_u8[2] = n;
    if (include_wins) {
      for (int i = 3; i <= DRAW_RULE + 2; i++) {
        n += (stats[i] != 0);
        mi.v_u8[i] = n;
      }
      n += (capt_cnt[stm][3] != 0);
      mi.v_u8[DRAW_RULE + 3] = n;
      n += (stm == BLACK && pawn_cnt[3] != 0);
      mi.v_u8[DRAW_RULE + 4] = n;
      for (int i = DRAW_RULE + 5; i < MAX_STATS / 2 + 1; i++) {
        n += (stats[i] != 0);
        mi.v_u8[i] = n;
      }
    } else {
      n += (win_vals != 0);
      for (int i = 3; i <= DRAW_RULE + 2; i++)
        mi.v_u8[i] = n;
      n += (capt_cnt[stm][3] != 0);
      mi.v_u8[DRAW_RULE + 3] = n;
      n += (stm == BLACK && pawn_cnt[3] != 0);
      mi.v_u8[DRAW_RULE + 4] = n;
      n += (cwin_vals != 0);
      for (int i = DRAW_RULE + 5; i < MAX_STATS / 2 + 1; i++)
        mi.v_u8[i] = n;
    }
    n += (capt_cnt[stm][2] != 0);
    mi.v_u8[MAX_STATS / 2 + 1] = n;
    n += (stats[MAX_STATS / 2 + 2] != 0);
    mi.v_u8[MAX_STATS / 2 + 2] = n;
    if (include_losses) {
      for (int i = MAX_STATS / 2 - 4; i >= 0; i--) {
        n += (stats[MAX_STATS - 1 - i] != 0);
        mi.v_u8[MAX_STATS - 1 - i] = n;
      }
    } else {
      n += (bloss_vals != 0);
      for (int i = MAX_STATS / 2 - 4; i >= DRAW_RULE + 1; i--)
        mi.v_u8[MAX_STATS - 1 - i] = n;
      n += (loss_vals != 0);
      for (int i = DRAW_RULE; i >= 0; i--)
        mi.v_u8[MAX_STATS - 1 - i] = n;
    }
    assert(n <= 255);

    for (int i = 0, j = -1; i < MAX_STATS; i++)
      if (mi.v_u8[i] != j)
        mi.v_inv_u8[j = mi.v_u8[i]] = i;

    merge_table = alloc_huge(sizeof(u8) * 8 * k16slice_alloc_size);
    if (!merge_table)
      out_of_mem();

    for (int s = 0; s < 240; s++)
      merge_bitmaps_u8(stm, s);

    free(merge_table);

  } else {

    // We need to use u16. This makes the mapping part straightfoward.
    include_wins = include_losses = true;
    for (int i = 0; i < MAX_STATS; i++)
      mi.v_u16[i] = mi.v_inv_u16[i] = i;

    merge_table = alloc_huge(sizeof(u16) * 8 * k16slice_alloc_size);
    if (!merge_table)
      out_of_mem();

    for (int s = 0; s < 240; s++)
      merge_bitmaps_u16(stm, s);

    free(merge_table);

  }

  char str[128], tmp[128];
  sprintf(str, "merge_info.%c", "wb"[stm]);
  strcat(strcpy(tmp, str), ".tmp");
  FILE *F = fopen(tmp, "wb");
  if (!F) {
    fprintf(stderr, "Could not open %s.\n", tmp);
    exit(EXIT_FAILURE);
  }
  file_write(&mi, sizeof mi, F);
  fclose(F);
  rename(tmp, str);
}
