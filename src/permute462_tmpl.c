/*
  Copyright (c) 2025, 2026 Ronald de Man

  This file is distributed under the terms of the GNU GPL, version 2.
*/

#define NAME(f) EVALUATOR(f,T)

#define NUM 16
static void NAME(convert_data_piece)(struct ThreadData *thread)
{
  uint8_t sq[MAX_PIECES];
  T *restrict src = convert_data.src;
  T *restrict dst = convert_data.dst;
  struct IdxInfo *perm_ii = convert_data.perm_ii;
  uint32_t sub[MAX_SETS];

  uint64_t idx_dec_buf[NUM];

  sq[0] = g_pos.sq[0];
  sq[1] = g_pos.sq[1];

  idx_to_sq_init(thread->begin, sub, perm_ii);

  uint64_t idx = thread->begin, end = thread->end;
  int fill = 0;
  int head = 0; // Next buffered element to consume.

  // Fill pipeline.
  for (; fill < NUM && idx < end; fill++, idx++, idx_to_sq_inc(sub, perm_ii)) {
    idx_to_sq_ii(sub, sq, perm_ii);
    uint64_t idx_dec = sq_to_idx(sq);
    __builtin_prefetch(src + idx_dec, 0, 3);
    idx_dec_buf[fill] = idx_dec;
  }

  // Steady-state pipeline.
  for (; idx < end; idx++, idx_to_sq_inc(sub, perm_ii)) {
    idx_to_sq_ii(sub, sq, perm_ii);
    uint64_t idx_dec = sq_to_idx(sq);
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

#if 0
static void NAME(convert_data_pawn)(struct ThreadData *thread)
{
  int n = base_entry.num;
  uint8_t pos[MAX_PIECES];
  T *restrict src = convert_data.src;
  T *restrict dst = convert_data.dst;
  uint8_t *restrict pidx = pidx_list[convert_data.p];
  int rank = convert_data.rank;
  uint64_t idx, end = thread->end;
  struct DecInfo di;
  uint32_t sub[MAX_PIECES];

  uint64_t idx_dec_buf[NUM];

  set_dec_info(&di, &base_entry, convert_data.pcs,
      type_perm_list[convert_data.p], order_list[convert_data.p],
      order2_list[convert_data.p], rank, RANK_ENC);

  decode_init(sub, thread->begin, &di);
  int k = 0, l;
  for (idx = thread->begin; k < NUM && idx < end; idx++, k++) {
    decode_pawn_r_iter(sub, pos, &di, &base_entry, rank);
    uint64_t idx_dec = calc_idx_pawn(pos, pidx, n);
    __builtin_prefetch(&src[idx_dec], 0, 3);
    idx_dec_buf[k] = idx_dec;
  }
  if (idx >= end) {
    l = k;
    goto finish;
  }

  for (; idx < end; idx++, k++) {
    k &= NUM - 1;
    decode_pawn_r_iter(sub, pos, &di, &base_entry, rank);
    uint64_t idx_dec = calc_idx_pawn(pos, pidx, n);
    __builtin_prefetch(&src[idx_dec], 0, 3);
    dst[idx - NUM] = src[idx_dec_buf[k]];
    idx_dec_buf[k] = idx_dec;
  }
  l = NUM;

finish:
  for (idx -= l; l > 0; l--, idx++, k++) {
    k &= NUM - 1;
    dst[idx] = src[idx_dec_buf[k]];
  }
}
#endif

static void NAME(convert_est_data_piece)(struct ThreadData *thread)
{
  T *restrict table = est_data.table;
  int num_cands = est_data.num_cands;
  uint32_t dsize = est_data.dsize;
  T *restrict dst = est_data.dst;
  uint8_t sq[MAX_PIECES];
  uint32_t sub[MAX_SETS];

  uint64_t idx_dec_buf[NUM];

  sq[0] = g_pos.sq[0];
  sq[1] = g_pos.sq[1];

  for (int p = 0; p < num_cands; p++) {
    for (int i = thread->begin; i < thread->end; i++) {
      idx_to_sq_init(segs[i], sub, &try_ii[p]);
      int j = 0, fill = 0, head = 0;

      for (; fill < NUM && j < seg_size;
          fill++, j++, idx_to_sq_inc(sub, &try_ii[p]))
      {
        idx_to_sq_ii(sub, sq, &try_ii[p]);
        uint64_t idx_dec = sq_to_idx(sq);
        __builtin_prefetch(table + idx_dec, 0, 3);
        idx_dec_buf[fill] = idx_dec;
      }

      for (; j < seg_size; j++, idx_to_sq_inc(sub, &try_ii[p])) {
        idx_to_sq_ii(sub, sq, &try_ii[p]);
        uint64_t idx_dec = sq_to_idx(sq);
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

#if 0
static void NAME(convert_est_data_pawn)(struct ThreadData *thread)
{
  int i, j, m, p, q, r;
  T *restrict table = est_data.table;
  int num_cands = est_data.num_cands;
  uint8_t *restrict pcs = est_data.pcs;
  uint32_t dsize = est_data.dsize;
  T *restrict dst = est_data.dst;
  int rank = est_data.rank;
  uint64_t idx;
  int n = base_entry.num;
  uint8_t pos[MAX_PIECES];
  struct DecInfo di;
  uint32_t sub[MAX_PIECES];

  uint64_t idx_cache[MAX_CANDS];
  
  for (i = thread->begin; i < thread->end; i++)
  {
    for (p = 0; p < num_cands;) {
      for (q = p + 1; q < num_cands; q++) {
	for (m = 0; m < num_types; m++)
	  if (cmp[type_perm_list[trylist[p]][m]] != cmp[type_perm_list[trylist[q]][m]]) break;
	if (m < num_types) break;
      }
      int l = trylist[p];
      set_dec_info(&di, &base_entry, pcs, type_perm_list[l],
          order_list[l], order2_list[l], rank, RANK_ENC);
      // prefetch for j = 0
      decode_init(sub, segs[i], &di);
      decode_pawn_r_iter(sub, pos, &di, &base_entry, rank);
      for (r = p; r < q; r++) {
	l = trylist[r];
        idx = calc_idx_pawn(pos, pidx_list[l], n);
	__builtin_prefetch(&table[idx], 0, 3);
	idx_cache[r] = idx;
      }
      for (j = 1; j < seg_size; j++) {
        decode_pawn_r_iter(sub, pos, &di, &base_entry, rank);
	for (r = p; r < q; r++) {
	  l = trylist[r];
          idx = calc_idx_pawn(pos, pidx_list[l], n);
	  __builtin_prefetch(&table[idx], 0, 3);
	  dst[r * dsize + i * seg_size + j - 1] = table[idx_cache[r]];
	  idx_cache[r] = idx;
	}
      }
      for (r = p; r < q; r++)
	dst[r * dsize + i * seg_size + j - 1] = table[idx_cache[r]];
      p = q;
    }
  }
}
#endif
#undef NAME
