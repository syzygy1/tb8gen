/*
  Copyright (c) 2025, 2026 Ronald de Man

  This file is distributed under the terms of the GNU GPL, version 2.
*/

#define NAME(f) EVALUATOR(f,T)

#define NUM 16
static void NAME(convert_data_piece)(struct ThreadData *thread)
{
  uint8_t sq[MAX_PIECES], sq2[MAX_PIECES];
  T *restrict src = convert_data.src;
  T *restrict dst = convert_data.dst;
  struct P10IdxInfo *perm_ii = convert_data.perm_ii;
  struct P10IdxState is;

  uint64_t idx_dec_buf[NUM];

  int stm = g_pos.stm;
  sq[stm] = g_pos.sq[stm];

  p10_idx_state_init(&is, thread->begin, sq, perm_ii, stm);

  uint64_t idx = thread->begin, end = thread->end;
  int fill = 0;
  int head = 0; // Next buffered element to consume.

  // Fill pipeline.
  for (; fill < NUM && idx < end;
      fill++, idx++, p10_idx_state_inc(&is, perm_ii))
  {
    p10_idx_state_to_sq(&is, sq, perm_ii, stm);
    normalize2(sq, sq2);
    uint64_t idx_dec = p10_sq_to_idx(sq2, sq[stm ^ 1]);
    __builtin_prefetch(src + idx_dec, 0, 3);
    idx_dec_buf[fill] = idx_dec;
  }

  // Steady-state pipeline.
  for (; idx < end; idx++, p10_idx_state_inc(&is, perm_ii)) {
    p10_idx_state_to_sq(&is, sq, perm_ii, stm);
    normalize2(sq, sq2);
    uint64_t idx_dec = p10_sq_to_idx(sq2, sq[stm ^ 1]);
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
  uint32_t dsize = est_data.dsize;
  T *restrict dst = est_data.dst;
  uint8_t sq[MAX_PIECES], sq2[MAX_PIECES];
  struct P10IdxState is;

  uint64_t idx_dec_buf[NUM];

  int stm = g_pos.stm;
  sq[stm] = g_pos.sq[stm];

  for (int p = 0; p < num_cands; p++) {
    for (int i = thread->begin; i < thread->end; i++) {
      p10_idx_state_init(&is, segs[i], sq, &try_ii[p], stm);
      int j = 0, fill = 0, head = 0;

      for (; fill < NUM && j < seg_size;
          fill++, j++, p10_idx_state_inc(&is, &try_ii[p]))
      {
        p10_idx_state_to_sq(&is, sq, &try_ii[p], stm);
        normalize2(sq, sq2);
        uint64_t idx_dec = p10_sq_to_idx(sq2, sq[stm ^ 1]);
        __builtin_prefetch(table + idx_dec, 0, 3);
        idx_dec_buf[fill] = idx_dec;
      }

      for (; j < seg_size; j++, p10_idx_state_inc(&is, &try_ii[p])) {
        p10_idx_state_to_sq(&is, sq, &try_ii[p], stm);
        normalize2(sq, sq2);
        uint64_t idx_dec = p10_sq_to_idx(sq2, sq[stm ^ 1]);
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
