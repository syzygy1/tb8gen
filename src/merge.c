/*
  Copyright (c) 2026 Ronald de Man

  This file is distributed under the terms of the GNU GPL, version 2.
*/

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "defs.h"
#include "index.h"
#include "generate.h"
#include "kslice.h"
#include "merge.h"
#include "movegen.h"
#include "tb8gen.h"
#include "threads.h"
#include "types.h"

static void *merge_table;
static int merge_n;
static int work_set, work_slice;

static uint64_t g_stats[2][MAX_STATS];
alignas(64) static uint64_t thread_stats[32][MAX_STATS];

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

void collect_stats(int stm)
{
  char str[128];
  uint64_t tmp[MAX_STATS];
  int i;

  int mx = max_iteration > 126 ? MAX_STATS : 256;
  uint64_t *stats = g_stats[stm];

  memset(stats, 0, sizeof(g_stats[stm]));

  for (int s = 0; s < 462; s++) {
    create_name(str, s, stm, "stats", -1);
    FILE *F = fopen(str, "rb");
    read_data(F, (void *)tmp, mx * sizeof(uint64_t));
    fclose(F);
    for (int i = 0; i < mx / 2; i++)
      stats[i] += tmp[i];
    stats[MAX_STATS / 2] += tmp[mx / 2];
    stats[MAX_STATS / 2 + 1] += tmp[mx / 2 + 1];
    for (int i = 0; i < mx / 2 - 2; i++)
      stats[MAX_STATS - 1 - i] += tmp[mx - 1 - i];
  }
  for (int i = 0; i < MAX_STATS; i++)
    stats[i] >>= 1;

  printf("\n%s to move:\n\n", stm == WHITE ? "White" : "Black");

  if (stats[1] + stats[2])
    printf("%lu (%lu) positions win in %d ply.\n", stats[1] + stats[2],
        stats[1], 1);
  for (i = 2; i <= DRAW_RULE; i++)
    if (stats[1 + i])
      printf("%lu positions win in %d ply.\n", stats[1 + i], i);
  if (stats[1 + i] + stats[2 + i])
    printf("%lu (%lu) positions win in %d ply.\n", stats[1 + i] + stats[2 + i],
        stats[1 + i], i);
  for (i += 2; i < MAX_STATS / 2 - 2; i++)
    if (stats[2 + i])
      printf("%lu positions win in %d ply.\n", stats[2 + i], i);
  printf("\n");

  uint64_t tot = 0;
  for (i = 1; i <= 1 + DRAW_RULE; i++)
    tot += stats[i];
  printf("%lu positions are wins.\n", tot);
  tot = 0;
  for (; i < MAX_STATS / 2; i++)
    tot += stats[i];
  if (tot)
    printf("%lu positions are cursed wins.\n", tot);
  if (stats[MAX_STATS / 2] + stats[MAX_STATS / 2 + 1])
    printf("%lu (%lu) positions are draws.\n",
        stats[MAX_STATS / 2] + stats[MAX_STATS / 2 + 1], stats[MAX_STATS / 2]);
  tot = 0;
  for (i = DRAW_RULE + 1; i < MAX_STATS / 2 - 2; i++)
    tot += stats[MAX_STATS - 1 - i];
  if (tot)
    printf("%lu positions are blessed losses.\n", tot);
  tot = 0;
  for (i = 0; i <= DRAW_RULE; i++)
    tot += stats[MAX_STATS - 1 - i];
  printf("%lu positions are losses.\n\n", tot);

  for (i = 0; i < MAX_STATS / 2 - 2; i++)
    if (stats[MAX_STATS - 1 - i])
      printf("%lu positions lose in %d ply.\n", stats[MAX_STATS - 1 - i], i);

  tot = 0;
  for (int i = 0; i < MAX_STATS; i++)
    tot += stats[i];
  printf("\n%lu positions out of %lu are illegal.\n", stats[0], tot);
}

// calculate DTZ entropy
static double entropy_helper(uint64_t *stats, uint64_t removed)
{
  for (int i = 0; i < MAX_STATS / 2; i++)
    for (int j = i + 1; j < MAX_STATS / 2; j++)
      if (stats[i] < stats[j])
        Swap(stats[i], stats[j]);
  for (int i = MAX_STATS / 2; i < MAX_STATS; i++)
    for (int j = i + 1; j < MAX_STATS; j++)
      if (stats[i] < stats[j])
        Swap(stats[i], stats[j]);
  for (int i = 0; i < MAX_STATS / 2; i++)
    stats[i] += stats[MAX_STATS / 2 + i];
  stats[0] += removed;

  uint64_t tot = 0;
  for (int i = 0; i < MAX_STATS / 2; i++)
    tot += stats[i];
  double entropy = 0;
  for (int i = 0; i < MAX_STATS / 2 && stats[i]; i++) {
    double p = (double)stats[i] / tot;
    entropy += -p * log(p);
  }
  entropy /= log(2.0);
  entropy = entropy * (double)tot / 8.0;

  return entropy;
}

double entropy_one_sided(int stm)
{
  uint64_t stats[MAX_STATS];
  memcpy(stats, g_stats[stm], sizeof(stats));

  uint64_t tot = stats[0] + stats[1] + stats[DRAW_RULE + 2]
    + stats[MAX_STATS / 2] + stats[MAX_STATS / 2 + 1];
  stats[0] = stats[1] = stats[DRAW_RULE + 2] = stats[MAX_STATS / 2]
    = stats[MAX_STATS / 2 + 1] = 0;

  return entropy_helper(stats, tot);
}

double entropy_loss_only(int stm)
{
  uint64_t stats[MAX_STATS];
  memcpy(stats, g_stats[stm], sizeof(stats));

  uint64_t tot = 0;
  for (int i = 0; i < MAX_STATS / 2 + 2; i++) {
    tot += stats[i];
    stats[i] = 0;
  }

  return entropy_helper(stats, tot);
}

double entropy_win_only(int stm)
{
  uint64_t stats[MAX_STATS];
  memcpy(stats, g_stats[stm], sizeof(stats));

  uint64_t tot = stats[0] + stats[1] + stats[DRAW_RULE + 2]
    + stats[MAX_STATS / 2] + stats[MAX_STATS / 2 + 1];
  stats[0] = stats[1] = stats[DRAW_RULE + 2] = stats[MAX_STATS / 2]
    = stats[MAX_STATS / 2 + 1] = 0;
  for (int i = MAX_STATS / 2 + 2; i < MAX_STATS; i++) {
    tot += stats[i];
    stats[i] = 0;
  }

  return entropy_helper(stats, tot);
}

void merge(int merge_type)
{
  bool wide = max_iteration > 126;

  merge_table = alloc_huge((1 + wide) * kslice_size);
  if (!merge_table)
    out_of_mem();

  if (merge_type == MERGE_SAVE) {
    create_dir(-1, WHITE, "merged");
    create_dir(-1, BLACK, "merged");
  }

  create_dir(-1, WHITE, "stats");
  create_dir(-1, BLACK, "stats");

  if (!wide)
    for (int s = 0; s < 462; s++) {
      merge_bitmaps_u8(WHITE, s);
      merge_bitmaps_u8(BLACK, s);
    }
  else
    for (int s = 0; s < 462; s++) {
      merge_bitmaps_u16(WHITE, s);
      merge_bitmaps_u16(BLACK, s);
    }
}
