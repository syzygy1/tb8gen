/*
  Copyright (c) 2026 Ronald de Man

  This file is distributed under the terms of the GNU GPL, version 2.
*/

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "compress.h"
#include "defs.h"
#include "generatep.h"
#include "kslicep.h"
#include "movegen.h"
#include "stats.h"
#include "tb8genp.h"
#include "types.h"
#include "util.h"

uint64_t g_stats[2][MAX_STATS];

void collect_stats(int stm)
{
  char str[64];
  uint64_t tmp[MAX_STATS];

  uint64_t *stats = g_stats[stm];

  memset(stats, 0, sizeof g_stats[stm]);

  for (int s = 0; s < 240; s++)
    for (int r = 0; r < 16; r++) {
      g_pos.sq[0] = KK16Square[s][r][0];
      g_pos.sq[1] = KK16Square[s][r][1];

      if (is_broken(&g_pos))
        continue;

      create_name_sq(str, g_pos.sq[0], g_pos.sq[1], stm, "stats", -1);
      FILE *F = file_open_read(str);
      read_data(F, tmp, sizeof tmp);
      fclose(F);
      for (int i = 0; i < MAX_STATS; i++)
        stats[i] += tmp[i];
    }
}

void print_stats(FILE *F, int stm)
{
  uint64_t *stats = g_stats[stm ^ flipped];

  fprintf(F, "\n%s to move:\n\n", stm == WHITE ? "White" : "Black");

  if (stats[1] + stats[2] + stats[3])
    fprintf(F, "%lu (%lu,%lu) positions win in %d ply.\n",
        stats[1] + stats[2] + stats[3], stats[1], stats[2], 1);
  for (int i = 4; i <= DRAW_RULE + 2; i++)
    if (stats[i])
      fprintf(F, "%lu positions win in %d ply.\n", stats[i], i - 2);
  if (stats[DRAW_RULE + 3] + stats[DRAW_RULE + 4] + stats[DRAW_RULE + 5])
    fprintf(F, "%lu (%lu,%lu) positions win in %d ply.\n",
        stats[DRAW_RULE + 3] + stats[DRAW_RULE + 4] + stats[DRAW_RULE + 5],
        stats[DRAW_RULE + 3], stats[DRAW_RULE + 4], DRAW_RULE + 1);
  for (int i = DRAW_RULE + 6; i < MAX_STATS / 2 + 1; i++)
    if (stats[i])
      fprintf(F, "%lu positions win in %d ply.\n", stats[i], i - 4);
  fprintf(F, "\n");

  uint64_t tot = 0;
  for (int i = 1; i <= DRAW_RULE + 2; i++)
    tot += stats[i];
  fprintf(F, "%lu positions are wins.\n", tot);
  tot = 0;
  for (int i = DRAW_RULE + 3; i < MAX_STATS / 2 + 1; i++)
    tot += stats[i];
  if (tot)
    fprintf(F, "%lu positions are cursed wins.\n", tot);
  if (stats[MAX_STATS / 2 + 1] + stats[MAX_STATS / 2 + 2])
    fprintf(F, "%lu (%lu) positions are draws.\n",
        stats[MAX_STATS / 2 + 1] + stats[MAX_STATS / 2 + 2],
        stats[MAX_STATS / 2 + 1]);
  tot = 0;
  for (int i = DRAW_RULE + 1; i < MAX_STATS / 2 - 3; i++)
    tot += stats[MAX_STATS - 1 - i];
  if (tot)
    fprintf(F, "%lu positions are blessed losses.\n", tot);
  tot = 0;
  for (int i = 0; i <= DRAW_RULE; i++)
    tot += stats[MAX_STATS - 1 -  i];
  fprintf(F, "%lu positions are losses.\n\n", tot);

  for (int i = 0; i < MAX_STATS / 2 - 3; i++)
    if (stats[MAX_STATS - 1 - i])
      fprintf(F, "%lu positions lose in %d ply.\n", stats[MAX_STATS - 1 - i], i);

  tot = 0;
  for (int i = 1; i < MAX_STATS; i++)
    tot += stats[i];
  fprintf(F, "\n%lu legal positions in total.\n", stats[0], tot);
}

void print_max_fens(FILE *F, struct MaxFen *mf)
{
  int w = WHITE ^ flipped, BLACK ^ flipped;

  if (mf->dtz[w][0] > 0)
    fprintf(F, "Longest win for white: %d ply; %s\n", mf->dtz[w][0] / 2,
        mf->fen[w][0]);
  if (mf->dtz[w][1] > 0)
    fprintf(F, "Longest cursed win for white: %d ply; %s\n",
        mf->dtz[w][1] / 2, mf->fen[w][1]);
  if (mf->dtz[b][1] > 0)
    fprintf(F, "Longest cursed win for black: %d ply; %s\n",
        mf->dtz[b][1] / 2, mf->fen[b][1]);
  if (mf->dtz[b][0] > 0)
    fprintf(F, "Longest win for black: %d ply; %s\n", mf->dtz[b][0] / 2,
        mf->fen[b][0]);
  fprintf(F, "\n");
}

// calculate DTZ entropy
static double entropy_helper(uint64_t *stats, uint64_t removed)
{
  uint64_t stats2[4][MAX_VAL] = { 0 };
  for (int i = 0; i <= DRAW_RULE; i++) {
    stats2[0][i] = stats[1 + i];
    stats2[1][i] = stats[MAX_STATS - 1 - i];
  }
  for (int i = 0; i < MAX_VAL; i++) {
    stats2[2][i] = stats[DRAW_RULE + 3 + 2 * i] + stats[DRAW_RULE + 4 + 2 * i];
    stats2[3][i] =  stats[MAX_STATS - DRAW_RULE - 2 - 2 * i]
                  + stats[MAX_STATS - DRAW_RULE - 3 - 2 * i];
  }
  for (int k = 0; k < 4; k++)
    for (int i = 0; i < MAX_VAL; i++)
      for (int j = i + 1; j < MAX_VAL; j++)
        if (stats2[k][i] < stats2[k][j])
          Swap(stats2[k][i], stats2[k][j]);

  for (int i = 0; i < MAX_VAL; i++)
    for (int k = 1; k < 4; k++)
      stats2[0][i] += stats2[k][i];

  stats2[0][0] += removed;

  uint64_t tot = 0;
  for (int i = 0; i < MAX_VAL; i++)
    tot += stats2[0][i];
  double entropy = 0;
  for (int i = 0; i < MAX_VAL && stats2[0][i]; i++) {
    double p = (double)stats2[0][i] / tot;
    entropy += -p * log(p);
  }
  entropy /= log(2.0);
  entropy = entropy * (double)tot / 8.0;

  return entropy;
}

double entropy_one_sided(int stm)
{
  uint64_t stats[MAX_STATS];
  memcpy(stats, g_stats[stm], sizeof stats);

  uint64_t tot = stats[0] + stats[1] + stats[2] + stats[DRAW_RULE + 3] +
    stats[DRAW_RULE + 4] + stats[MAX_STATS / 2 + 1] + stats[MAX_STATS / 2 + 2];
  stats[0] = stats[1] = stats[2] = stats[DRAW_RULE + 3] = stats[DRAW_RULE + 4]
    = stats[MAX_STATS / 2 + 1] = stats[MAX_STATS / 2 + 2] = 0;

  return entropy_helper(stats, tot);
}

double entropy_loss_only(int stm)
{
  uint64_t stats[MAX_STATS];
  memcpy(stats, g_stats[stm], sizeof stats);

  uint64_t tot = 0;
  for (int i = 0; i < MAX_STATS / 2 + 3; i++) {
    tot += stats[i];
    stats[i] = 0;
  }

  return entropy_helper(stats, tot);
}

double entropy_win_only(int stm)
{
  uint64_t stats[MAX_STATS];
  memcpy(stats, g_stats[stm], sizeof stats);

  uint64_t tot = stats[0] + stats[1] + stats[2] + stats[DRAW_RULE + 3] +
    stats[DRAW_RULE + 4] + stats[MAX_STATS / 2 + 1] + stats[MAX_STATS / 2 + 2];
  stats[0] = stats[1] = stats[2] = stats[DRAW_RULE + 3] = stats[DRAW_RULE + 4]
    = stats[MAX_STATS / 2 + 1] = stats[MAX_STATS / 2 + 2] = 0;
  for (int i = MAX_STATS / 2 + 3; i < MAX_STATS; i++) {
    tot += stats[i];
    stats[i] = 0;
  }

  return entropy_helper(stats, tot);
}
