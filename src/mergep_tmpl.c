/*
  Copyright (c) 2026 Ronald de Man

  This file is distributed under the terms of the GNU GPL, version 2.
*/

#define NAME(f) EVALUATOR(f,T)

static int NAME(init_merge_value_map)(uint64_t *stats, int stm)
{
  int n = 0;
  NAME(mi.v)[0] = n;
  n += (stats[1] != 0);
  NAME(mi.v)[1] = n;
  n += (stats[2] != 0);
  NAME(mi.v)[2] = n;
  if (include_wins) {
    for (int i = 3; i <= DRAW_RULE + 2; i++) {
      n += (stats[i] != 0);
      NAME(mi.v)[i] = n;
    }
    n += (capt_cnt[stm][3] != 0);
    NAME(mi.v)[DRAW_RULE + 3] = n;
    n += (stm == BLACK && pawn_cnt[3] != 0);
    NAME(mi.v)[DRAW_RULE + 4] = n;
    for (int i = DRAW_RULE + 5; i < MAX_STATS / 2 + 1; i++) {
      n += (stats[i] != 0);
      NAME(mi.v)[i] = n;
    }
  } else {
    n += mi.v_wdl[4];
    for (int i = 3; i <= DRAW_RULE + 2; i++)
      NAME(mi.v)[i] = n;
    n += (capt_cnt[stm][3] != 0);
    NAME(mi.v)[DRAW_RULE + 3] = n;
    n += (stm == BLACK && pawn_cnt[3] != 0);
    NAME(mi.v)[DRAW_RULE + 4] = n;
    n += mi.v_wdl[3];
    for (int i = DRAW_RULE + 5; i < MAX_STATS / 2 + 1; i++)
      NAME(mi.v)[i] = n;
  }
  n += (capt_cnt[stm][2] != 0);
  NAME(mi.v)[MAX_STATS / 2 + 1] = n;
  n += (stats[MAX_STATS / 2 + 2] != 0);
  NAME(mi.v)[MAX_STATS / 2 + 2] = n;
  if (include_losses) {
    for (int i = MAX_STATS / 2 - 4; i >= 0; i--) {
      n += (stats[MAX_STATS - 1 - i] != 0);
      NAME(mi.v)[MAX_STATS - 1 - i] = n;
    }
  } else {
    n += mi.v_wdl[1];
    for (int i = MAX_STATS / 2 - 4; i >= DRAW_RULE + 1; i--)
      NAME(mi.v)[MAX_STATS - 1 - i] = n;
    n += mi.v_wdl[0];
    for (int i = DRAW_RULE; i >= 0; i--)
      NAME(mi.v)[MAX_STATS - 1 - i] = n;
  }

  for (int i = 0, j = -1; i < MAX_STATS; i++)
    if (NAME(mi.v)[i] != j)
      NAME(mi.v_inv)[j = NAME(mi.v)[i]] = i;

  return n;
}

static void NAME(merge_transform)(struct ThreadData *thread)
{
  T *restrict const q = merge_table;
  T *w = merge_w;
  for (uint64_t idx = thread->begin, end = thread->end; idx < end; idx++)
    q[idx] = w[q[idx]];
}

static void NAME(merge_draw_worker)(struct ThreadData *thread)
{
  T *restrict const q = merge_table;
  T n = NAME(mi.v)[MAX_STATS / 2 + 2];
  for (uint64_t idx = thread->begin, end = thread->end; idx < end; idx++)
    q[idx] = n;
}

static void NAME(merge_worker)(struct ThreadData *thread)
{
  T n = NAME(mi.v)[merge_n];

  const uint64_t *restrict p = (uint64_t *)k16slice_get_address(-1);
  T *restrict const q = merge_table;
  p += thread->begin >> 6;
  for (uint64_t idx = thread->begin, end = thread->end; idx < end; idx += 64) {
    uint64_t w = *p++;
    while (w)
      q[idx + pop_lsb(&w)] = n;
  }
}

static void NAME(merge_repl_worker)(struct ThreadData *thread)
{
  T r = NAME(mi.v)[merge_r];
  T n = NAME(mi.v)[merge_n];

  const uint64_t *restrict p = (uint64_t *)k16slice_get_address(-1);
  T *restrict const q = merge_table;
  p += thread->begin >> 6;
  for (uint64_t idx = thread->begin, end = thread->end; idx < end; idx += 64) {
    uint64_t w = *p++;
    while (w) {
      unsigned bt = pop_lsb(&w);
      if (q[idx + bt] == r)
        q[idx + bt] = n;
    }
  }
}

static void NAME(merge_capt_bloss_worker)(struct ThreadData *thread)
{
  const uint64_t *restrict p = (uint64_t *)k16slice_get_address(-1);
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
    memcpy(sq2, sq, sizeof sq2);
    sq2[k] = pop_lsb(&b);
    p[sq_to_idx(sq2)] = 0;
  }
}

static void NAME(merge_illegal_worker)(struct ThreadData *thread)
{
  struct IdxState is;
  Position pos = g_pos;
  int k = work_set;
  int m = ri.last[k];
  int stm = pos.stm;
  int king_sq = pos.sq[stm ^ 1];

  T *restrict const p = (T *)merge_table + 8 * k16offset(g_pos.sq);

  idx_state_init(&is, thread->begin, pos.sq, &capt_ri[k]);

  for (uint64_t idx = thread->begin, end = thread->end; idx < end;
      idx++, idx_state_inc(&is, &capt_ri[k]))
  {
    Bitboard occ = idx_state_to_sq(&is, pos.sq, &capt_ri[k]);
    pos.sq[m] = king_sq;
    NAME(merge_mark_unmoves)(m, p, occ, pos.sq);
  }
}

static void NAME(merge_set_illegal_worker)(struct ThreadData *thread)
{
  T *restrict const p = (T *)merge_table + 8 * k16offset(g_pos.sq);

  for (uint64_t idx = thread->begin, end = thread->end; idx < end; idx++)
    p[idx] = 0;
}

static void NAME(merge_statistics_worker)(struct ThreadData *thread)
{
  int t = thread->thread_id;

  T *restrict const p = (T *)merge_table + 8 * k16offset(g_pos.sq);

  for (uint64_t idx = thread->begin, end = thread->end; idx < end; idx++)
    thread_stats[t][NAME(mi.v_inv)[p[idx]]]++;
}

// u16:
// ILLEGAL     = 0
// CAPT_WIN    = 1
// PAWN_WIN    = 2
// wins          3..102 -> win-in-1 to win-in-100
// CAPT_CWIN   = 103
// PAWN_CWIN   = 104
// wins          105-128  -> win-in-101 to-win-in-124
// CAPT_DRAW   = MAX_STATS/2+1
// unknown       MAX_STATS/2+2 -> draw, including PAWN_DRAW
// CAPT_CLOSS  = x
// losses        131..154 -> loss-in-124 to loss-in-101
// losses        MAX_STATS-101..MAX_STATS-2 -> loss-in-100 to loss-in-1
// losses        MAX_STATS-1 -> mate

static void NAME(merge_bitmaps)(int stm, int s)
{
  char str[64];

  create_name_r(str, s, 15, stm, "merged/wdl", -1);
  if (file_exists(str))
    return;

  uint64_t stats[16][MAX_STATS] = { 0 };

  g_pos.stm = stm;

  // DRAW
  run_threaded(NAME(merge_draw_worker), work_g16);

  // losses
  for (int n = 0; n < max_iteration; n++)
    if (   g_stats[stm][MAX_STATS - 1 - n]
        && k16slice_test(s, stm, "L", n)
        && k16slice_read(-1, s, stm, "L", n))
    {
      merge_n = MAX_STATS - 1 - n;
      run_threaded(NAME(merge_worker), work_g16);
      if (!include_losses)
        stat_count(stats, merge_n);
      if (2 * n == mf.dtz[stm ^ 1][0] && !mf.found[stm ^ 1][0])
        find_position(stm ^ 1, s, true, false);
      else if (2 * n == mf.dtz[stm ^ 1][1] && !mf.found[stm ^ 1][1])
        find_position(stm ^ 1, s, true, true);
    }

  // wins
  for (int n = 1; n < max_iteration; n++)
    if (   (n == 1 || n == DRAW_RULE + 1
            || g_stats[stm][2 + n + 2 * (n > DRAW_RULE)])
        && k16slice_test(s, stm, "W", n)
        && k16slice_read(-1, s, stm, "W", n))
    {
      merge_n = n <= DRAW_RULE ? 2 + n : 4 + n;
      run_threaded(NAME(merge_worker), work_g16);
      if (!include_wins)
        stat_count(stats, merge_n);
      if (2 * n + 1 == mf.dtz[stm][0] && !mf.found[stm][0])
        find_position(stm, s, false, false);
      if (2 * n + 1 == mf.dtz[stm][1] && !mf.found[stm][1])
        find_position(stm, s, false, true);
    }

  // CAPT_WIN
  if (g_stats[stm][1] && k16slice_read(-1, s, stm, "capt/win", -1)) {
    merge_r = 3;
    merge_n = 1;
    run_threaded(NAME(merge_repl_worker), work_g16);
  }

  // PAWN_WIN
  if (g_stats[stm][2] && k16slice_read(-1, s, stm, "pawn/win", -1)) {
    merge_r = 3;
    merge_n = 2;
    run_threaded(NAME(merge_repl_worker), work_g16);
  }

  // CAPT_CWIN
  if (capt_cnt[stm][3] && k16slice_read(-1, s, stm, "capt/cwin", -1)) {
    merge_r = DRAW_RULE + 5;
    merge_n = DRAW_RULE + 3;
    run_threaded(NAME(merge_repl_worker), work_g16);
  }

  // PAWN_CWIN
  if (stm == BLACK && pawn_cnt[3] && k16slice_read(-1, s, stm, "pawn/cwin", -1))
  {
    merge_r = DRAW_RULE + 5;
    merge_n = DRAW_RULE + 4;
    run_threaded(NAME(merge_repl_worker), work_g16);
  }

  // CAPT_DRAW
  if (capt_cnt[stm][2] && k16slice_read(-1, s, stm, "capt/draw", -1)) {
    merge_r = MAX_STATS / 2 + 2;
    merge_n = MAX_STATS / 2 + 1;
    run_threaded(NAME(merge_repl_worker), work_g16);
  }

  // ILLEGAL
  for (int r = 0; r < 16; r++) {
    g_pos.sq[0] = KK16Square[s][r][0];
    g_pos.sq[1] = KK16Square[s][r][1];

    if (    is_broken(&g_pos)
        || (    stm == BLACK
            && (pawn_attacks(BLACK, g_pos.sq[2]) & bit(g_pos.sq[0]))))
    {
      run_threaded(NAME(merge_set_illegal_worker), work_g);
      continue;
    }

    for (int k = 0; k < ri.numsets; k++) {
      if ((g_pos.pt[ri.first[k]] >> 3) != stm)
        continue;
      work_set = k;
      run_threaded(NAME(merge_illegal_worker), work_capt[k]);
    }
  }

  // CAPT_BLOSS is added below to produce WDL files

  // Calculate and save stats for each of the 16 K-slices.
  for (int r = 0; r < 16; r++) {
    g_pos.sq[0] = KK16Square[s][r][0];
    g_pos.sq[1] = KK16Square[s][r][1];

    if (is_broken(&g_pos))
      continue;

    for (int t = 0; t < g_num_threads; t++)
      memset(thread_stats[t], 0, sizeof thread_stats[t]);
    run_threaded(NAME(merge_statistics_worker), work_g);

    uint64_t win_tmp, cwin_tmp, bloss_tmp, loss_tmp;
    win_tmp = stats[r][3];
    cwin_tmp = stats[r][DRAW_RULE + 5];
    bloss_tmp = stats[r][MAX_STATS / 2 + 3];
    loss_tmp = stats[r][MAX_STATS - 1 - DRAW_RULE];

    for (int t = 0; t < g_num_threads; t++)
      for (int i = 0; i < MAX_STATS; i++)
        stats[r][i] += thread_stats[t][i];

    // This is a bit hackish.
    if (!include_wins) {
      stats[r][3] = win_tmp - stats[r][1] - stats[r][2];
      stats[r][DRAW_RULE + 5] =
        cwin_tmp - stats[r][DRAW_RULE + 3] - stats[r][DRAW_RULE + 4];
    }

    if (!include_losses) {
      stats[r][MAX_STATS / 2 + 3] = bloss_tmp;
      stats[r][MAX_STATS - 1 - DRAW_RULE] = loss_tmp;
    }

    create_name_r(str, s, r, stm, "stats", -1);
    FILE *F = file_open_write(str);
    write_data(F, stats[r], sizeof stats[r]);
    fclose(F);
    file_rename(str);
  }

  if (symmetric && stm == BLACK)
    return;

  if (!one_sided || stm == one_sided_stm) {
    T z[MAX];
    for (int i = 0; i < MAX; i++)
      z[i] = 0;

    if (one_sided || wins_only) {
      for (int i = 3; i <= DRAW_RULE + 2; i++)
        z[NAME(mi.v)[i]] = NAME(mi.v)[i];
      for (int i = DRAW_RULE + 5; i < MAX_STATS / 2 + 1; i++)
        z[NAME(mi.v)[i]] = NAME(mi.v)[i];
    }

    if (one_sided || !wins_only)
      for (int i = 0; i < MAX_STATS / 2 - 3; i++)
        z[NAME(mi.v)[MAX_STATS - 1 - i]] = NAME(mi.v)[MAX_STATS - 1 - i];

    for (int r = 0; r < 16; r++) {
      g_pos.sq[0] = KK16Square[s][r][0];
      g_pos.sq[1] = KK16Square[s][r][1];

      if (is_broken(&g_pos))
        continue;

      create_name_r(str, s, r, stm, "merged/dtz", -1);
      FILE *F = file_open_write(str);
      NAME(write_data_transform)(F,
          (T *)merge_table + r * 8 * kslice_alloc_size * sizeof(T),
          kslice_size * sizeof(T), z);
      fclose(F);
      file_rename(str);
    }
  }

  // 0/1/2/3/4 -> loss/bloss/draw/cwin/win
  // 5/6/7/8 -> capt_bloss/_draw/_cwin/_win+illegal
  T w[MAX];
  for (int i = 0; i <= DRAW_RULE; i++)
    w[NAME(mi.v)[MAX_STATS - 1 - i]] = 0;
  for (int i = DRAW_RULE + 1; i < MAX_STATS / 2 - 3; i++)
    w[NAME(mi.v)[MAX_STATS - 1 - i]] = 1;
  w[NAME(mi.v)[MAX_STATS / 2 + 2]] = 2; // DRAW including PAWN_DRAW
  w[NAME(mi.v)[MAX_STATS / 2 + 1]] = 6; // CAPT_DRAW
  for (int i = MAX_STATS / 2; i >= DRAW_RULE + 4; i--)
    w[NAME(mi.v)[i]] = 3;
  w[NAME(mi.v)[DRAW_RULE + 3]] = 7;
  for (int i = DRAW_RULE + 2; i >= 2; i--)
    w[NAME(mi.v)[i]] = 4;
  w[NAME(mi.v)[1]] = 8;
  w[NAME(mi.v)[0]] = 8;

  merge_w = w;
  run_threaded(NAME(merge_transform), work_g16);

  // Replace bloss (1) with capt_bloss (5) where appropriate.
  if (capt_cnt[stm][1]) {
    k16slice_read(-1, s, stm, "capt/bloss", -1);
    run_threaded(NAME(merge_capt_bloss_worker), work_g16);
  }

  for (int r = 0; r < 16; r++) {
    g_pos.sq[0] = KK16Square[s][r][0];
    g_pos.sq[1] = KK16Square[s][r][1];

    if (is_broken(&g_pos))
      continue;

    create_name_r(str, s, r, stm, "merged/wdl", -1);
    FILE *F = file_open_write(str);
    NAME(write_data_as_u8)(F, (T *)merge_table + r * 8 * kslice_alloc_size,
        kslice_size);
    fclose(F);
    file_rename(str);
  }
}
