/*
  Copyright (c) 2025, 2026 Ronald de Man

  This file is distributed under the terms of the GNU GPL, version 2.
*/

#define NAME(f) EVALUATOR(f,T)

#define NUM 16

static void NAME(init_source_rank_ri)(struct RankInfo *rank_ri,
    const struct RankInfo *perm_ri)
{
  *rank_ri = ri;
  for (int k = 0; k < ri.numsets; k++) {
    int j = 0;
    for (; j < perm_ri->numsets; j++)
      if (ri.first[k] == perm_ri->first[j])
        break;
    assert(j < perm_ri->numsets);
    rank_ri->perm[k] = j + 1;
  }
}

static void NAME(convert_data_piece)(struct ThreadData *thread)
{
  T *restrict src = convert_data.src;
  T *restrict dst = convert_data.dst;
  struct RankInfo *perm_ri = convert_data.perm_ri;
  struct RankInfo rank_ri;
  struct IdxState2 is;

  uint64_t idx_dec_buf[NUM];

  NAME(init_source_rank_ri)(&rank_ri, perm_ri);
  idx_state2_init(&is, thread->begin, g_pos.sq, perm_ri);

  uint64_t idx = thread->begin, end = thread->end;
  int fill = 0;
  int head = 0; // Next buffered element to consume.

  // Fill pipeline.
  for (; fill < NUM && idx < end;
      fill++, idx++, idx_state2_inc(&is, perm_ri))
  {
    uint64_t idx_dec = perm_rank_bb(is.bb, &rank_ri);
    __builtin_prefetch(src + idx_dec, 0, 3);
    idx_dec_buf[fill] = idx_dec;
  }

  // Steady-state pipeline.
  for (; idx < end; idx++, idx_state2_inc(&is, perm_ri)) {
    uint64_t idx_dec = perm_rank_bb(is.bb, &rank_ri);
    __builtin_prefetch(src + idx_dec, 0, 3);
    dst[idx - NUM] = src[idx_dec_buf[head]];
    idx_dec_buf[head] = idx_dec;
    head = (head + 1) & (NUM - 1);
  }

  // Drain pipeline.
  for (uint64_t out = idx - fill; fill-- > 0; out++) {
    dst[out] = src[idx_dec_buf[head]];
    head = (head + 1) & (NUM - 1);
  }
}

static void NAME(convert_est_data_piece)(struct ThreadData *thread)
{
  T *restrict table = est_data.table;
  int num_cands = est_data.num_cands;
  uint64_t dsize = est_data.dsize;
  T *restrict dst = est_data.dst;
  struct RankInfo rank_ri;
  struct IdxState2 is;

  uint64_t idx_dec_buf[NUM];

  for (int p = 0; p < num_cands; p++) {
    NAME(init_source_rank_ri)(&rank_ri, &try_ri[p]);
    for (int i = thread->begin; i < thread->end; i++) {
      idx_state2_init(&is, segs[i], g_pos.sq, &try_ri[p]);
      int j = 0, fill = 0, head = 0;

      for (; fill < NUM && j < seg_size;
          fill++, j++, idx_state2_inc(&is, &try_ri[p]))
      {
        uint64_t idx_dec = perm_rank_bb(is.bb, &rank_ri);
        __builtin_prefetch(table + idx_dec, 0, 3);
        idx_dec_buf[fill] = idx_dec;
      }

      for (; j < seg_size; j++, idx_state2_inc(&is, &try_ri[p])) {
        uint64_t idx_dec = perm_rank_bb(is.bb, &rank_ri);
        __builtin_prefetch(table + idx_dec, 0, 3);
        dst[p * dsize + i * seg_size + j - NUM] = table[idx_dec_buf[head]];
        idx_dec_buf[head] = idx_dec;
        head = (head + 1) & (NUM - 1);
      }

      for (j -= fill; fill-- > 0; j++) {
        dst[p * dsize + i * seg_size + j] = table[idx_dec_buf[head]];
        head = (head + 1) & (NUM - 1);
      }
    }
  }
}

static void NAME(convert_data_piece_ref)(struct ThreadData *thread)
{
  T *restrict src = convert_data.src;
  T *restrict dst = convert_data.dst;
  struct RankInfo *perm_ri = convert_data.perm_ri;
  struct RankInfo rank_ri;
  struct IdxState2 is;

  uint64_t idx_dec_buf[NUM];

  NAME(init_source_rank_ri)(&rank_ri, perm_ri);
  is.bb[0] = bit(g_pos.sq[0]) | bit(g_pos.sq[1]);

  uint64_t idx = thread->begin, end = thread->end;
  int fill = 0;
  int head = 0; // Next buffered element to consume.

  // Fill pipeline.
  for (; fill < NUM && idx < end; fill++, idx++) {
    unrank_bb_ref(idx, is.bb, perm_ri);
    uint64_t idx_dec = perm_rank_bb_ref(is.bb, &rank_ri);
    __builtin_prefetch(src + idx_dec, 0, 3);
    idx_dec_buf[fill] = idx_dec;
  }

  // Steady-state pipeline.
  for (; idx < end; idx++) {
    unrank_bb_ref(idx, is.bb, perm_ri);
    uint64_t idx_dec = perm_rank_bb_ref(is.bb, &rank_ri);
    __builtin_prefetch(src + idx_dec, 0, 3);
    dst[idx - NUM] = src[idx_dec_buf[head]];
    idx_dec_buf[head] = idx_dec;
    head = (head + 1) & (NUM - 1);
  }

  // Drain pipeline.
  for (uint64_t out = idx - fill; fill-- > 0; out++) {
    dst[out] = src[idx_dec_buf[head]];
    head = (head + 1) & (NUM - 1);
  }
}

static void NAME(convert_est_data_piece_ref)(struct ThreadData *thread)
{
  T *restrict table = est_data.table;
  int num_cands = est_data.num_cands;
  uint64_t dsize = est_data.dsize;
  T *restrict dst = est_data.dst;
  struct RankInfo rank_ri;
  struct IdxState2 is;

  uint64_t idx_dec_buf[NUM];

  is.bb[0] = bit(g_pos.sq[0]) | bit(g_pos.sq[1]);

  for (int p = 0; p < num_cands; p++) {
    NAME(init_source_rank_ri)(&rank_ri, &try_ri[p]);
    for (int i = thread->begin; i < thread->end; i++) {
      int j = 0, fill = 0, head = 0;

      for (; fill < NUM && j < seg_size; fill++, j++) {
        unrank_bb_ref(segs[i] + j, is.bb, &try_ri[p]);
        uint64_t idx_dec = perm_rank_bb_ref(is.bb, &rank_ri);
        __builtin_prefetch(table + idx_dec, 0, 3);
        idx_dec_buf[fill] = idx_dec;
      }

      for (; j < seg_size; j++) {
        unrank_bb_ref(segs[i] + j, is.bb, &try_ri[p]);
        uint64_t idx_dec = perm_rank_bb_ref(is.bb, &rank_ri);
        __builtin_prefetch(table + idx_dec, 0, 3);
        dst[p * dsize + i * seg_size + j - NUM] = table[idx_dec_buf[head]];
        idx_dec_buf[head] = idx_dec;
        head = (head + 1) & (NUM - 1);
      }

      for (j -= fill; fill-- > 0; j++) {
        dst[p * dsize + i * seg_size + j] = table[idx_dec_buf[head]];
        head = (head + 1) & (NUM - 1);
      }
    }
  }
}

#undef NAME
