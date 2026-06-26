/*
  Copyright (c) 2025, 2026 Ronald de Man

  This file is distributed under the terms of the GNU GPL, version 2.
*/

#define NAME(f) EVALUATOR(f,T)

#define NUM 16
static void NAME(convert_data_pawn)(struct ThreadData *thread)
{
  T *restrict src = convert_data.src;
  T *restrict dst = convert_data.dst;
  struct PKIdxInfo *perm_ii = convert_data.perm_ii;
  struct RankInfo rank_ri;
  int8_t king_perm;
  struct PKIdxState is;

  uint64_t idx_dec_buf[NUM];

  int stm = g_pos.stm;

  init_source_rank_ri_pk(&rank_ri, &king_perm, perm_ii, stm);
  pk_idx_state_init(&is, thread->begin, perm_ii, stm);

  uint64_t idx = thread->begin, end = thread->end;
  int fill = 0;
  int head = 0; // Next buffered element to consume.

  // Fill pipeline.
  for (; fill < NUM && idx < end;
      fill++, idx++, pk_idx_state_inc(&is, perm_ii))
  {
    pk_idx_state_to_bb(&is, perm_ii);
    uint64_t idx_dec = pk_bb_to_idx(&is, &rank_ri, king_perm, stm);
    __builtin_prefetch(src + idx_dec, 0, 3);
    idx_dec_buf[fill] = idx_dec;
  }

  // Steady-state pipeline.
  for (; idx < end; idx++, pk_idx_state_inc(&is, perm_ii)) {
    pk_idx_state_to_bb(&is, perm_ii);
    uint64_t idx_dec = pk_bb_to_idx(&is, &rank_ri, king_perm, stm);
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

static void NAME(convert_est_data_pawn)(struct ThreadData *thread)
{
  T *restrict table = est_data.table;
  int num_cands = est_data.num_cands;
  uint32_t dsize = est_data.dsize;
  T *restrict dst = est_data.dst;
  struct RankInfo rank_ri;
  int8_t king_perm;
  struct PKIdxState is;

  uint64_t idx_dec_buf[NUM];

  int stm = g_pos.stm;

  for (int p = 0; p < num_cands; p++) {
    init_source_rank_ri_pk(&rank_ri, &king_perm, &try_ii[p], stm);
    for (int i = thread->begin; i < thread->end; i++) {
      pk_idx_state_init(&is, segs[i], &try_ii[p], stm);
      int j = 0, fill = 0, head = 0;

      for (; fill < NUM && j < seg_size;
          fill++, j++, pk_idx_state_inc(&is, &try_ii[p]))
      {
        pk_idx_state_to_bb(&is, &try_ii[p]);
        uint64_t idx_dec = pk_bb_to_idx(&is, &rank_ri, king_perm, stm);
        __builtin_prefetch(table + idx_dec, 0, 3);
        idx_dec_buf[fill] = idx_dec;
      }

      for (; j < seg_size; j++, pk_idx_state_inc(&is, &try_ii[p])) {
        pk_idx_state_to_bb(&is, &try_ii[p]);
        uint64_t idx_dec = pk_bb_to_idx(&is, &rank_ri, king_perm, stm);
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
