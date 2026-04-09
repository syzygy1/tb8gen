/*
  Copyright (c) 2026 Ronald de Man

  This file is distributed under the terms of the GNU GPL, version 2.
*/

#define NAME(f) EVALUATOR(f,T)

static void NAME(merge_draw_worker)(struct ThreadData *thread)
{
  T *restrict const q = merge_table;
  T n = NAME(mi.v)[MAX_STATS / 2 + 1];
  for (uint64_t idx = thread->begin, end = thread->end; idx < end; idx++)
    q[idx] = n;
}

static void NAME(merge_worker)(struct ThreadData *thread)
{
  T n = NAME(mi.v)[merge_n];

  const uint64_t *restrict p = (uint64_t *)kslice_get_address(-1);
  T *restrict const q = merge_table;
  p += thread->begin >> 6;
  for (uint64_t idx = thread->begin, end = thread->end; idx < end; idx += 64) {
    uint64_t w = *p++;
    while (w)
      q[idx + pop_lsb(&w)] = n;
  }
}

static void NAME(merge_capt_bloss_worker)(struct ThreadData *thread)
{
  const uint64_t *restrict p = (uint64_t *)kslice_get_address(-1);
  T *restrict const q = merge_table;
  p += thread->begin >> 6;
  for (uint64_t idx = thread->begin, end = thread->end; idx < end; idx += 64) {
    uint64_t w = *p++;
    while (w) {
      unsigned bt = pop_lsb(&w);
      if (q[idx + bt] == 1)
        q[idx + bt] = 5;
    }
  }
}

INLINE void NAME(merge_mark_unmoves)(int k, T *restrict const p, Bitboard occ,
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

static void NAME(merge_illegal_worker)(struct ThreadData *thread)
{
  uint32_t sub[MAX_SETS];
  Position pos = g_pos;
  int k = work_set;
  int m = ii.last[k];
  int stm = pos.stm;
  int king_sq = pos.sq[stm ^ 1];

  T *restrict const p = merge_table;

  idx_to_sq_init(thread->begin, sub, &capt_ii[k]);

  for (uint64_t idx = thread->begin, end = thread->end; idx < end;
      idx++, idx_to_sq_inc(sub, &capt_ii[k]))
  {
    Bitboard occ = capt_idx_to_sq(sub, pos.sq, k);
    pos.sq[m] = king_sq;
    NAME(merge_mark_unmoves)(m, p, occ, pos.sq);
  }
}

static void NAME(merge_statistics_worker)(struct ThreadData *thread)
{
  int s = work_slice;
  int t = thread->thread_id;

  T *restrict const p = merge_table;

  if (s < 441) {
    for (uint64_t idx = thread->begin, end = thread->end; idx < end; idx++)
      thread_stats[t][NAME(mi.v_inv)[p[idx]]] += 2;
  } else {
    // Both kings on the diagonal, so count carefully.
    uint32_t sub[MAX_SETS];
    Position pos = g_pos;
    idx_to_sq_init(thread->begin, sub, &ii);

    for (uint64_t idx = thread->begin, end = thread->end; idx < end;
        idx++, idx_to_sq_inc(sub, &ii))
    {
      idx_to_sq(sub, pos.sq);
      mirror_diagonal(pos.sq);
      uint64_t idx2 = sq_to_idx(pos.sq);
      thread_stats[t][NAME(mi.v_inv)[p[idx]]] += idx == idx2 ? 2 : 1;
    }
  }
}

// u16:
// ILLEGAL     = 0
// CAPT_WIN    = 1
// wins          2..101 -> win-in-1 to win-in-100
// CAPT_CWIN   = 102
// wins          103-127  -> win-in-101 to-win-in-125
// CAPT_DRAW   = MAX_STATS/2
// unknown       MAX_STATS/2+1 -> draw
// CAPT_CLOSS  = x
// losses        130..154 -> loss-in-125 to loss-in-101
// losses        MAX_STATS-101..MAX_STATS-2 -> loss-in-100 to loss-in-1
// losses        MAX_STATS-1 -> mate

static void NAME(merge_bitmaps)(int stm, int s)
{
  uint64_t stats[MAX_STATS] = { 0 };

  // DRAW
  run_threaded(NAME(merge_draw_worker), work_g, 0);

  g_pos.stm = stm;
  g_pos.sq[0] = KKSquare[s][0];
  g_pos.sq[1] = KKSquare[s][1];

  // ILLEGAL
  for (int k = 0; k < ii.numsets; k++) {
    if ((g_pos.pt[ii.first[k]] >> 3) != stm)
      continue;
    work_set = k;
    run_threaded(NAME(merge_illegal_worker), work_capt[k], 0);
  }

  // losses
  for (int n = 0; n < max_iteration; n++)
    if (kslice_test(s, stm, "L", n)) {
      kslice_read(-1, s, stm, "L", n);
      merge_n = MAX_STATS - 1 - n;
      run_threaded(NAME(merge_worker), work_g, 0);
      if (!include_losses)
        stats[merge_n] = stat_count(stm, s);
    }

  // wins
  for (int n = 1; n < max_iteration; n++)
    if (kslice_test(s, stm, "W", n)) {
      kslice_read(-1, s, stm, "W", n);
      merge_n = n <= DRAW_RULE ? 1 + n : 2 + n;
      run_threaded(NAME(merge_worker), work_g, 0);
      if (!include_wins)
        stats[merge_n] = stat_count(stm, s);
    }

  // CAPT_WIN
  if (sub_cnt[stm ^ 1][0]) {
    kslice_read(-1, s, stm, "capt/win", -1);
    merge_n = 1;
    run_threaded(NAME(merge_worker), work_g, 0);
  }

  // CAPT_CWIN
  if (sub_cnt[stm ^ 1][1]) {
    kslice_read(-1, s, stm, "capt/cwin", -1);
    merge_n = 2 + DRAW_RULE;
    run_threaded(NAME(merge_worker), work_g, 0);
  }

  // CAPT_DRAW
  if (sub_cnt[stm ^ 1][2]) {
    kslice_read(-1, s, stm, "capt/draw", -1);
    merge_n = MAX_STATS / 2;
    run_threaded(NAME(merge_worker), work_g, 0);
  }

  // CAPT_BLOSS to be added to produce WDL files

  for (int t = 0; t < g_num_threads; t++)
    memset(thread_stats[t], 0, sizeof thread_stats[t]);
  work_slice = s;
  run_threaded(NAME(merge_statistics_worker), work_g, 0);

  uint64_t win_tmp, cwin_tmp, bloss_tmp, loss_tmp;
  win_tmp = stats[2];
  cwin_tmp = stats[DRAW_RULE + 3];
  bloss_tmp = stats[MAX_STATS / 2 + 2];
  loss_tmp = stats[MAX_STATS - 1 - DRAW_RULE];

  for (int t = 0; t < g_num_threads; t++)
    for (int i = 0; i < MAX_STATS; i++)
      stats[i] += thread_stats[t][i];

  // This is a bit hackish.
  if (!include_wins) {
    stats[2] = win_tmp - stats[1];
    stats[DRAW_RULE + 3] = cwin_tmp - stats[DRAW_RULE + 2];
  }

  if (!include_losses) {
    stats[MAX_STATS / 2 + 2] = bloss_tmp;
    stats[MAX_STATS - 1 - DRAW_RULE] = loss_tmp;
  }

  char str[128], tmp[128];
  create_name(str, s, stm, "stats", -1);
  strcat(strcpy(tmp, str), ".tmp");
  FILE *F = fopen(tmp, "wb");
  if (!F) {
    fprintf(stderr, "Could not open %s.\n", tmp);
    exit(EXIT_FAILURE);
  }
  write_data(F, stats, sizeof stats);
  fclose(F);
  rename(tmp, str);

  if (symmetric && stm == BLACK)
    return;

  if (one_sided && stm != one_sided_stm)
    goto start_wdl;

  T z[MAX];
  for (int i = 0; i < MAX; i++)
    z[i] = 0;

  if (one_sided || wins_only) {
    for (int i = 2; i <= DRAW_RULE + 1; i++)
      z[NAME(mi.v)[i]] = NAME(mi.v)[i];
    for (int i = DRAW_RULE + 3; i < MAX_STATS / 2; i++)
      z[NAME(mi.v)[i]] = NAME(mi.v)[i];
  }

  if (one_sided || !wins_only)
    for (int i = 0; i < MAX_STATS / 2 - 2; i++)
      z[NAME(mi.v)[MAX_STATS - 1 - i]] = NAME(mi.v)[MAX_STATS - 1 - i];

  create_name(str, s, stm, "merged/dtz", -1);
  strcat(strcpy(tmp, str), ".tmp");
  F = fopen(tmp, "wb");
  if (!F) {
    fprintf(stderr, "Could not open %s.\n", tmp);
    exit(EXIT_FAILURE);
  }
  NAME(write_data_transform)(F, merge_table, kslice_size * sizeof(T), z);
  fclose(F);
  rename(tmp, str);

start_wdl:
  // 0/1/2/3/4 -> loss/bloss/draw/cwin/win
  // 5/6/7/8 -> capt_bloss/_draw/_cwin/_win+illegal
  uint8_t w[MAX];
  for (int i = 0; i <= DRAW_RULE; i++)
    w[NAME(mi.v)[MAX_STATS - 1 - i]] = 0;
  for (int i = DRAW_RULE + 1; i < MAX_STATS / 2 - 2; i++)
    w[NAME(mi.v)[MAX_STATS - 1 - i]] = 1;
  w[NAME(mi.v)[MAX_STATS / 2 + 1]] = 2;
  w[NAME(mi.v)[MAX_STATS / 2]] = 6;
  for (int i = MAX_STATS / 2 - 1; i >= DRAW_RULE + 3; i--)
    w[NAME(mi.v)[i]] = 3;
  w[NAME(mi.v)[DRAW_RULE + 2]] = 7;
  for (int i = DRAW_RULE + 1; i >= 2; i--)
    w[NAME(mi.v)[i]] = 4;
  w[NAME(mi.v)[1]] = 8;
  w[NAME(mi.v)[0]] = 8;

  // Replace bloss (1) with capt_bloss (5) where appropriate.
  if (sub_cnt[stm ^ 1][3]) {
    kslice_read(-1, s, stm, "capt/bloss", -1);
    run_threaded(NAME(merge_capt_bloss_worker), work_g, 0);
  }

  create_name(str, s, stm, "merged/wdl", -1);
  strcat(strcpy(tmp, str), ".tmp");
  F = fopen(tmp, "wb");
  if (!F) {
    fprintf(stderr, "Could not open %s.\n", tmp);
    exit(EXIT_FAILURE);
  }
  NAME(write_data_transform_to_u8)(F, merge_table, kslice_size, w);
  fclose(F);
  rename(tmp, str);
}
