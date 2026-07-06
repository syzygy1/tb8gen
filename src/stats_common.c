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
  if (stats[DRAW_RULE + 3] + stats[DRAW_RULE + 4] + stats[DRAW_RULE + 5])
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
