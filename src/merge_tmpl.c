/*
  Copyright (c) 2026 Ronald de Man

  This file is distributed under the terms of the GNU GPL, version 2.
*/

#define NAME(f) EVALUATOR(f,T)

static int NAME(init_merge_value_map)(uint64_t *stats)
{
  int n = 0;
  NAME(mi.v)[0] = n;
  n += (stats[1] != 0);
  NAME(mi.v)[1] = n;
  assert(stats[2] == 0);
  NAME(mi.v)[2] = n;
  if (include_wins) {
    for (int i = 3; i <= DRAW_RULE + 2; i++) {
      n += (stats[i] != 0);
      NAME(mi.v)[i] = n;
    }
    n += (stats[DRAW_RULE + 3] != 0);
    NAME(mi.v)[DRAW_RULE + 3] = n;
    assert(stats[DRAW_RULE + 4] == 0);
    NAME(mi.v)[DRAW_RULE + 4] = n;
    for (int i = DRAW_RULE + 5; i < MAX_STATS / 2 + 1; i++) {
      n += (stats[i] != 0);
      NAME(mi.v)[i] = n;
    }
  } else {
    n += mi.v_wdl[4];
    for (int i = 3; i <= DRAW_RULE + 2; i++)
      NAME(mi.v)[i] = n;
    n += (stats[DRAW_RULE + 3] != 0);
    NAME(mi.v)[DRAW_RULE + 3] = n;
    assert(stats[DRAW_RULE + 4] == 0);
    NAME(mi.v)[DRAW_RULE + 4] = n;
    n += mi.v_wdl[3];
    for (int i = DRAW_RULE + 5; i < MAX_STATS / 2 + 1; i++)
      NAME(mi.v)[i] = n;
  }
  n += (stats[MAX_STATS / 2 + 1] != 0);
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
  bool found = false;
  const uint64_t *restrict p = (uint64_t *)kslice_get_address(-1);
  T *restrict const q = merge_table;
  p += thread->begin >> 6;
  for (uint64_t idx = thread->begin, end = thread->end; idx < end; idx += 64) {
    uint64_t w = *p++;
    while (w) {
      uint64_t cur = idx + pop_lsb(&w);
      if (q[cur] == 1) {
        q[cur] = 5;
        found = true;
      }
    }
  }
  if (found)
    atomic_store_explicit(&found_capt_bloss, true, memory_order_relaxed);
}

INLINE void NAME(merge_mark_unmoves)(int k, T *restrict const p, Bitboard occ,
    const uint8_t *restrict sq)
{
  uint8_t sq2[MAX_PIECES];
  Bitboard b = non_king_piece_moves(g_pos.pt[k], sq[k], occ);
  while (b) {
    memcpy(sq2, sq, sizeof sq2);
    sq2[k] = pop_lsb(&b);
    p[sq_to_idx(sq2)] = 0;
  }
}

INLINE void NAME(merge_mark_ref_unmoves)(int k, T *restrict const p,
    Bitboard occ, const uint8_t *restrict sq)
{
  uint8_t sq2[MAX_PIECES];
  Bitboard b = non_king_piece_moves(g_pos.pt[k], sq[k], occ);
  while (b) {
    memcpy(sq2, sq, sizeof sq2);
    sq2[k] = pop_lsb(&b);
    p[sq_to_idx_ref(sq2)] = 0;
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

  T *restrict const p = merge_table;

  idx_state_init(&is, thread->begin, pos.sq, &capt_ri[k]);

  for (uint64_t idx = thread->begin, end = thread->end; idx < end;
      idx++, idx_state_inc(&is, &capt_ri[k]))
  {
    Bitboard occ = idx_state_to_sq(&is, pos.sq, &capt_ri[k]);
    pos.sq[m] = king_sq;
    NAME(merge_mark_unmoves)(m, p, occ, pos.sq);
  }
}

static void NAME(merge_illegal_ref_worker)(struct ThreadData *thread)
{
  struct IdxState is;
  Position pos = g_pos;
  int k = work_set;
  int m = ri.last[k];
  int stm = pos.stm;
  int king_sq = pos.sq[stm ^ 1];

  T *restrict const p = merge_table;

  idx_state_init(&is, thread->begin, pos.sq, &capt_ri[k]);

  for (uint64_t idx = thread->begin, end = thread->end; idx < end;
      idx++, idx_state_inc(&is, &capt_ri[k]))
  {
    Bitboard occ = idx_state_to_sq(&is, pos.sq, &capt_ri[k]);
    pos.sq[m] = king_sq;
    NAME(merge_mark_ref_unmoves)(m, p, occ, pos.sq);
  }
}

static void NAME(merge_statistics_worker)(struct ThreadData *thread)
{
  int t = thread->thread_id;

  T *restrict const p = merge_table;

  for (uint64_t idx = thread->begin, end = thread->end; idx < end; idx++)
    thread_stats[t][NAME(mi.v_inv)[p[idx]]]++;
}

// u16:
// ILLEGAL     = 0
// CAPT_WIN    = 1
// wins          3..102 -> win-in-1 to win-in-100
// CAPT_CWIN   = 103
// wins          105-127  -> win-in-101 to-win-in-125
// CAPT_DRAW   = MAX_STATS/2+1
// unknown       MAX_STATS/2+2 -> draw
// CAPT_CLOSS  = x
// losses        131..154 -> loss-in-124 to loss-in-101
// losses        MAX_STATS-101..MAX_STATS-2 -> loss-in-100 to loss-in-1
// losses        MAX_STATS-1 -> mate

static void NAME(merge_bitmaps)(int stm, int s)
{
  char str[64];

  create_name(str, s, stm, "merged/wdl", -1);
  if (file_exists(str))
    return;

  uint64_t stats[MAX_STATS] = { 0 };

  // DRAW
  run_threaded(NAME(merge_draw_worker), &work_g_merge_static[s >= 441]);

  g_pos.stm = stm;
  g_pos.sq[0] = KKSquare[s][0];
  g_pos.sq[1] = KKSquare[s][1];

  // ILLEGAL
  for (int k = 0; k < ri.numsets; k++) {
    if ((g_pos.pt[ri.first[k]] >> 3) != stm)
      continue;
    work_set = k;
    if (s < 441)
      run_threaded(NAME(merge_illegal_worker), &work_capt_merge_dynamic[k]);
    else
      run_threaded(NAME(merge_illegal_ref_worker), &work_capt_merge_dynamic[k]);
  }

  // losses
  for (int n = 0; n < max_iteration; n++)
    if (kslice_test(s, stm, "L", n)) {
      kslice_read(-1, s, stm, "L", n);
      merge_n = MAX_STATS - 1 - n;
      run_threaded(NAME(merge_worker), &work_g_merge_dynamic[s >= 441]);
      if (!include_losses)
        stats[merge_n] = kslice_read_count;
      if (2 * n == mf.dtz[stm ^ 1][0] && !mf.found[stm ^ 1][0])
        find_position(s, stm ^ 1, true, false);
      else if (2 * n == mf.dtz[stm ^ 1][1] && !mf.found[stm ^ 1][1])
        find_position(s, stm ^ 1, true, true);
    }

  // wins
  for (int n = 1; n < max_iteration; n++)
    if (kslice_test(s, stm, "W", n)) {
      kslice_read(-1, s, stm, "W", n);
      merge_n = stats_n(n);
      run_threaded(NAME(merge_worker), &work_g_merge_dynamic[s >= 441]);
      if (!include_wins)
        stats[merge_n] = kslice_read_count;
      if (2 * n + 1 == mf.dtz[stm][0] && !mf.found[stm][0])
        find_position(s, stm, false, false);
      if (2 * n + 1 == mf.dtz[stm][1] && !mf.found[stm][1])
        find_position(s, stm, false, true);
    }

  // CAPT_WIN
  if (sub_cnt[stm ^ 1][0]) {
    kslice_read(-1, s, stm, "capt/win", -1);
    merge_n = 1;
    run_threaded(NAME(merge_worker), &work_g_merge_dynamic[s >= 441]);
  }

  // CAPT_CWIN
  if (sub_cnt[stm ^ 1][1]) {
    kslice_read(-1, s, stm, "capt/cwin", -1);
    merge_n = DRAW_RULE + 3;
    run_threaded(NAME(merge_worker), &work_g_merge_dynamic[s >= 441]);
  }

  // CAPT_DRAW
  if (sub_cnt[stm ^ 1][2]) {
    kslice_read(-1, s, stm, "capt/draw", -1);
    merge_n = MAX_STATS / 2 + 1;
    run_threaded(NAME(merge_worker), &work_g_merge_dynamic[s >= 441]);
  }

  // CAPT_BLOSS to be added to produce WDL files

  set_merge_active_threads(&work_g_merge_static[s >= 441]);
  for (int tt = 0; tt < merge_num_active_threads; tt++) {
    int t = merge_active_threads[tt];
    memset(thread_stats[t], 0, sizeof thread_stats[t]);
  }
  run_threaded(NAME(merge_statistics_worker), &work_g_merge_static[s >= 441]);

  uint64_t win_tmp, cwin_tmp, bloss_tmp, loss_tmp;
  win_tmp = stats[3];
  cwin_tmp = stats[DRAW_RULE + 5];
  bloss_tmp = stats[MAX_STATS / 2 + 3];
  loss_tmp = stats[MAX_STATS - 1 - DRAW_RULE];

  for (int tt = 0; tt < merge_num_active_threads; tt++) {
    int t = merge_active_threads[tt];
    for (int i = 0; i < MAX_STATS; i++)
      stats[i] += thread_stats[t][i];
  }

  // This is a bit hackish.
  if (!include_wins) {
    stats[3] = win_tmp - stats[1];
    stats[DRAW_RULE + 5] = cwin_tmp - stats[DRAW_RULE + 3];
  }

  if (!include_losses) {
    stats[MAX_STATS / 2 + 3] = bloss_tmp;
    stats[MAX_STATS - 1 - DRAW_RULE] = loss_tmp;
  }

  create_name(str, s, stm, "stats", -1);
  FILE *F = file_open_write(str);
  write_data(F, stats, sizeof stats);
  fclose(F);
  file_rename(str);

  if (symmetric && stm == BLACK)
    return;

  if (!one_sided || stm == one_sided_stm) {
    T z[MAX] = { 0 };

    if (one_sided || wins_only) {
      for (int i = 3; i <= DRAW_RULE + 2; i++)
        z[NAME(mi.v)[i]] = NAME(mi.v)[i];
      for (int i = DRAW_RULE + 5; i < MAX_STATS / 2 + 1; i++)
        z[NAME(mi.v)[i]] = NAME(mi.v)[i];
    }

    if (one_sided || !wins_only)
      for (int i = 0; i < MAX_STATS / 2 - 3; i++)
        z[NAME(mi.v)[MAX_STATS - 1 - i]] = NAME(mi.v)[MAX_STATS - 1 - i];

    create_name(str, s, stm, "merged/dtz", -1);
    F = file_open_write(str);
    NAME(write_data_transform)(F, merge_table,
        kslice_sizes[s >= 441] * sizeof(T), z);
    fclose(F);
    file_rename(str);
  }

  // 0/1/2/3/4 -> loss/bloss/draw/cwin/win
  // 5/6/7/8 -> capt_bloss/_draw/_cwin/_win+illegal
  T w[MAX];
  for (int i = 0; i <= DRAW_RULE; i++)
    w[NAME(mi.v)[MAX_STATS - 1 - i]] = 0;
  for (int i = DRAW_RULE + 1; i < MAX_STATS / 2 - 3; i++)
    w[NAME(mi.v)[MAX_STATS - 1 - i]] = 1;
  w[NAME(mi.v)[MAX_STATS / 2 + 2]] = 2;
  w[NAME(mi.v)[MAX_STATS / 2 + 1]] = 6;
  for (int i = MAX_STATS / 2; i >= DRAW_RULE + 4; i--)
    w[NAME(mi.v)[i]] = 3;
  w[NAME(mi.v)[DRAW_RULE + 3]] = 7;
  for (int i = DRAW_RULE + 2; i >= 2; i--)
    w[NAME(mi.v)[i]] = 4;
  w[NAME(mi.v)[1]] = 8;
  w[NAME(mi.v)[0]] = 8;

  merge_w = w;
  run_threaded(NAME(merge_transform), &work_g_merge_static[s >= 441]);

  // Replace bloss (1) with capt_bloss (5) where appropriate.
  atomic_store_explicit(&found_capt_bloss, false, memory_order_relaxed);
  if (sub_cnt[stm ^ 1][3]) {
    kslice_read(-1, s, stm, "capt/bloss", -1);
    run_threaded(NAME(merge_capt_bloss_worker), &work_g_merge_dynamic[s >= 441]);
  }

  create_name(str, s, stm, "merged/wdl", -1);
  F = file_open_write(str);
  uint8_t has_capt_bloss =
    atomic_load_explicit(&found_capt_bloss, memory_order_relaxed);
  file_write(&has_capt_bloss, 1, F);
  NAME(write_data_as_u8)(F, (T *)merge_table, kslice_sizes[s >= 441]);
  fclose(F);
  file_rename(str);
}
