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
#include "stats.h"
#include "types.h"
#include "util.h"
#include "xxhash.h"

uint64_t g_stats[2][MAX_STATS];
XXH128_hash_t wdl_checksum, dtz_checksum[2];

void collect_stats(int stm)
{
  char str[64];
  uint64_t tmp[MAX_STATS];

  uint64_t *stats = g_stats[stm];

  memset(stats, 0, sizeof g_stats[stm]);

#ifndef HAS_PAWNS

  for (int s = 0; s < 462; s++) {
    create_name(str, s, stm, "stats", -1);
    FILE *F = file_open_read(str);
    read_data(F, tmp, sizeof tmp);
    fclose(F);
    for (int i = 0; i < MAX_STATS; i++)
      stats[i] += tmp[i];
  }

#else

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

#endif
}

void print_stats(FILE *F, int stm)
{
  uint64_t *stats = g_stats[stm ^ flipped];

  fprintf(F, "\n%s to move:\n\n", stm == WHITE ? "White" : "Black");

#ifndef HAS_PAWNS
  if (stats[1] + stats[3])
    fprintf(F, "%lu (%lu) positions win in %d ply.\n", stats[1] + stats[3],
        stats[1], 1);
#else
  if (stats[1] + stats[2] + stats[3])
    fprintf(F, "%lu (%lu,%lu) positions win in %d ply.\n",
        stats[1] + stats[2] + stats[3], stats[1], stats[2], 1);
#endif

  for (int i = 4; i <= DRAW_RULE + 2; i++)
    if (stats[i])
      fprintf(F, "%lu positions win in %d ply.\n", stats[i], i - 2);

#ifndef HAS_PAWNS
  if (stats[DRAW_RULE + 3] + stats[DRAW_RULE + 5])
    fprintf(F, "%lu (%lu) positions win in %d ply.\n",
        stats[DRAW_RULE + 3] + stats[DRAW_RULE + 5], stats[DRAW_RULE + 3],
        DRAW_RULE + 1);
#else
  if (stats[DRAW_RULE + 3] + stats[DRAW_RULE + 4] + stats[DRAWRULE + 5])
    fprintf(F, "%lu (%lu,%lu) positions win in %d ply.\n",
        stats[DRAW_RULE + 3] + stats[DRAW_RULE + 4] + stats[DRAW_RULE + 5],
        stats[DRAW_RULE + 3], stats[DRAW_RULE + 4], DRAW_RULE + 1);
#endif

  for (int i = DRAW_RULE + 6; i < MAX_STATS / 2 + 1; i++)
    if (stats[i])
      fprintf(F, "%lu positions win in %d ply.\n", stats[i], i - 4);

  uint64_t tot = 0;
  for (int i = 1; i <= DRAW_RULE + 2; i++)
    tot += stats[i];
  fprintf(F, "\n%lu positions are wins.\n", tot);
  tot = 0;
  for (int i = DRAW_RULE + 3; i < MAX_STATS / 2 + 1; i++)
    tot += stats[i];
  if (tot)
    fprintf(F, "%lu positions are cursed wins.\n", tot);
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
    tot += stats[MAX_STATS - 1 - i];
  fprintf(F, "%lu positions are losses.\n\n", tot);

  for (int i = 0; i < MAX_STATS / 2 - 3; i++)
    if (stats[MAX_STATS - 1 - i])
      fprintf(F, "%lu positions lose in %d ply.\n", stats[MAX_STATS - 1 - i], i);

  tot = 0;
  for (int i = 1; i < MAX_STATS; i++)
    tot += stats[i];
  fprintf(F, "\n%lu legal positions in total.\n", tot);
}

void print_max_fens(FILE *F, struct MaxFen *mf)
{
  int w = WHITE ^ flipped, b = BLACK ^ flipped;

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

static void stats_to_freq(uint64_t stats[MAX_STATS], uint64_t f[2][MAX_VAL])
{
  memset(f, 0, MAX_VAL * 16);

  f[0][1] = stats[1] + stats[2] + stats[3];
  for (int i = 4; i <= DRAW_RULE + 2; i++)
    f[0][i - 2] = stats[i];

  f[0][DRAW_RULE + 1] = stats[DRAW_RULE + 3] + stats[DRAW_RULE + 4];
  for (int i = DRAW_RULE + 5; i <= MAX_STATS / 2; i++)
    f[0][DRAW_RULE + 1 + (i - DRAW_RULE - 5) / 2] += stats[i];

  for (int i = 0; i <= DRAW_RULE; i++)
    f[1][i] = stats[MAX_STATS - 1 - i];

  for (int i = DRAW_RULE + 1; i < MAX_STATS / 2 - 4; i++)
    f[1][DRAW_RULE + 1 + (i - DRAW_RULE - 1) / 2] += stats[MAX_STATS - 1 - i];
}

static XXH128_hash_t freq_to_hash(uint64_t f[2][2][MAX_VAL])
{
  uint64_t v[4 * MAX_VAL];
  int n;
  for (n = MAX_VAL; n >= 1; n--)
    if (f[0][0][n - 1] | f[0][1][n - 1] | f[1][0][n - 1] | f[1][1][n - 1])
      break;
  for (int i = 0; i < n; i++) {
    v[4 * i    ] = f[0][0][i];
    v[4 * i + 1] = f[0][1][i];
    v[4 * i + 2] = f[1][0][i];
    v[4 * i + 3] = f[1][1][i];
  }
  return XXH3_128bits(v, 32 * n);
}

void calc_stats_checksums(void)
{
  uint64_t wdl_counts[2][5] = { 0 };
  for (int stm = 0; stm < 2; stm++) {
    uint64_t *stats = g_stats[stm];
    for (int i = 0; i <= DRAW_RULE; i++)
      wdl_counts[stm][0] += stats[MAX_STATS - 1 - i];
    for (int i = DRAW_RULE + 1; i < MAX_STATS / 2 - 2; i++)
      wdl_counts[stm][1] += stats[MAX_STATS - 1 - i];
    wdl_counts[stm][2] += stats[MAX_STATS / 2] + stats[MAX_STATS / 2 + 1];
    for (int i = DRAW_RULE + 2; i < MAX_STATS / 2; i++)
      wdl_counts[stm][3] += stats[i];
    for (int i = 1; i <= DRAW_RULE + 1; i++)
      wdl_counts[stm][4] += stats[i];
  }
  wdl_checksum = XXH3_128bits(wdl_counts, sizeof wdl_counts);

  uint64_t counts[2][2][MAX_VAL];
  stats_to_freq(g_stats[0], counts[0]);
  stats_to_freq(g_stats[1], counts[1]);
  dtz_checksum[0] = freq_to_hash(counts);

  if (one_sided) {
    memset(counts[one_sided_stm ^ 1], 0, 16 * MAX_VAL);
  }
  else if (wins_only) {
    memset(counts[0][1], 0, 8 * MAX_VAL);
    memset(counts[1][1], 0, 8 * MAX_VAL);
  }
  else {
    memset(counts[0][0], 0, 8 * MAX_VAL);
    memset(counts[1][0], 0, 8 * MAX_VAL);
  }

  dtz_checksum[1] = freq_to_hash(counts);
}

// calculate DTZ entropy
static double entropy_helper(uint64_t *stats, uint64_t removed)
{
  uint64_t freq[4][MAX_VAL] = { 0 };

  for (int i = 1; i <= DRAW_RULE; i++)
    freq[0][i - 1] += stats[2 + i];

  for (int i = DRAW_RULE + 1; i < MAX_STATS / 2 - 3; i++)
    freq[2][(i - DRAW_RULE - 1) / 2] += stats[4 + i];

  freq[1][0] = stats[MAX_STATS - 1];
  for (int i = 1; i <= DRAW_RULE; i++)
    freq[1][i - 1] += stats[MAX_STATS - 1 - i];

  for (int i = DRAW_RULE + 1; i < MAX_STATS / 2 - 3; i++)
    freq[3][(i - DRAW_RULE - 1) / 2] += stats[MAX_STATS - 1 - i];

  for (int k = 0; k < 4; k++)
    for (int i = 0; i < MAX_VAL; i++)
      for (int j = i + 1; j < MAX_VAL; j++)
        if (freq[k][i] < freq[k][j])
          Swap(freq[k][i], freq[k][j]);

  for (int k = 1; k < 4; k++)
    for (int i = 0; i < MAX_VAL; i++)
      freq[0][i] += freq[k][i];

  freq[0][0] += removed;

  uint64_t tot = 0;
  for (int i = 0; i < MAX_VAL; i++)
    tot += freq[0][i];
  double entropy = 0;
  for (int i = 0; i < MAX_VAL && freq[0][i]; i++) {
    double p = (double)freq[0][i] / tot;
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
