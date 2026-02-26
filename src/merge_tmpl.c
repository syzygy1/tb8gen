/*
  Copyright (c) 2026 Ronald de Man

  This file is distributed under the terms of the GNU GPL, version 2.
*/

#define NAME(f) EVALUATOR(f,T)

static void NAME(merge_draw_worker)(struct ThreadData *thread)
{
  T *restrict const q = merge_table;
  for (uint64_t idx = thread->begin, end = thread->end; idx < end; idx++)
    q[idx] = MAX / 2 + 1;
}

static void NAME(merge_worker)(struct ThreadData *thread)
{
  T n = merge_n;

  uint64_t *restrict p = (uint64_t *)kslice_get_address(-1);
  T *restrict const q = merge_table;
  p += thread->begin >> 6;
  for (uint64_t idx = thread->begin, end = thread->end; idx < end; idx += 64) {
    uint64_t w = *p++;
    while (w)
      q[idx + pop_lsb(&w)] = n;
  }
}

static void NAME(merge_capt_worker)(struct ThreadData *thread)
{
  T n = merge_n;

  uint64_t *restrict p = (uint64_t *)kslice_get_address(-1);
  T *restrict const q = merge_table;
  p += thread->begin >> 6;
  for (uint64_t idx = thread->begin, end = thread->end; idx < end; idx += 64) {
    uint64_t w = *p++;
    while (w) {
      unsigned bt = pop_lsb(&w);
      q[idx + bt] = min(q[idx + bt], n);
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

// u8:
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
      merge_n = MAX - 1 - n;
      run_threaded(NAME(merge_worker), work_g, 0);
    }

  // wins
  for (int n = 1; n < max_iteration; n++)
    if (kslice_test(s, stm, "W", n)) {
      kslice_read(-1, s, stm, "W", n);
      merge_n = n <= DRAW_RULE ? 1 + n : 2 + n;
      run_threaded(NAME(merge_worker), work_g, 0);
    }

  // CAPT_WIN
  if (sub_cnt[stm ^ 1][0]) {
    kslice_read(-1, s, stm, "capt/win", -1);
    merge_n = 1;
    run_threaded(NAME(merge_capt_worker), work_g, 0);
  }

  // CAPT_CWIN
  if (sub_cnt[stm ^ 1][1]) {
    kslice_read(-1, s, stm, "capt/cwin", -1);
    merge_n = 2 + DRAW_RULE;
    run_threaded(NAME(merge_capt_worker), work_g, 0);
  }

  // CAPT_DRAW
  if (sub_cnt[stm ^ 1][2]) {
    kslice_read(-1, s, stm, "capt/draw", -1);
    merge_n = MAX / 2;
    run_threaded(NAME(merge_capt_worker), work_g, 0);
  }

  // CAPT_BLOSS to be added to produce WDL files

  for (int t = 0; t < g_num_threads; t++)
    memset(thread_stats[t], 0, sizeof(thread_stats[t]));
  work_slice = s;
  run_threaded(NAME(merge_statistics_worker), work_g, 0);

  uint64_t stats[MAX] = { 0 };
  for (int t = 0; t < g_num_threads; t++)
    for (int i = 0; i < MAX; i++)
      stats[i] += thread_stats[t][i];

  char str[128];
  create_name(str, s, stm, "stats", -1);
  FILE *F = fopen(str, "wb");
  write_data(F, (void *)stats, sizeof(stats));
  fclose(F);
}
