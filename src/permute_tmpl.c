/*
  Copyright (c) 2025, 2026 Ronald de Man

  This file is distributed under the terms of the GNU GPL, version 2.
*/

#define NAME(f) EVALUATOR(f,T)

#define NUM 16
static void NAME(convert_data_piece)(struct ThreadData *thread)
{
  uint8_t pos[MAX_PIECES];
  T *restrict src = convert_data.src;
  T *restrict dst = convert_data.dst;
  uint8_t *restrict pidx = pidx_list[convert_data.p];
  uint64_t idx, end = thread->end;
  struct DecInfo di;
  uint32_t sub[MAX_PIECES];
#ifdef RAMGEN
  T *restrict v = permute_v;
#endif

  uint64_t idx_dec_buf[NUM];

  set_dec_info(&di, &tb_entry, convert_data.pcs,
      type_perm_list[convert_data.p], order_list[convert_data.p], 0xf, -1,
      LT_PIECE);

  decode_init(sub, thread->begin, &di);
  int k = 0, l;
  for (idx = thread->begin; k < NUM && idx < end; idx++, k++) {
    decode_piece_iter(sub, pos, &di, &tb_entry);
    uint64_t idx_dec = calc_idx_piece(pos, pidx);
    __builtin_prefetch(&src[idx_dec], 0, 3);
    idx_dec_buf[k] = idx_dec;
  }
  if (idx >= end) {
    l = k;
    goto finish;
  }

  for (; idx < end; idx++, k++) {
    k &= NUM - 1;
    decode_piece_iter(sub, pos, &di, &tb_entry);
    uint64_t idx_dec = calc_idx_piece(pos, pidx);
    __builtin_prefetch(&src[idx_dec], 0, 3);
#ifndef RAMGEN
    dst[idx - NUM] = src[idx_dec_buf[k]];
#else
    dst[idx - NUM] = v[src[idx_dec_buf[k]]];
#endif
    idx_dec_buf[k] = idx_dec;
  }
  l = NUM;

finish:
  for (idx -= l; l > 0; l--, idx++, k++) {
    k &= NUM - 1;
#ifndef RAMGEN
    dst[idx] = src[idx_dec_buf[k]];
#else
    dst[idx] = v[src[idx_dec_buf[k]]];
#endif
  }
}

static void NAME(convert_est_data_piece)(struct ThreadData *thread)
{
  int i, j, m, p, q, r;
  T *restrict table = est_data.table;
  int num_cands = est_data.num_cands;
  uint8_t *restrict pcs = est_data.pcs;
  uint32_t dsize = est_data.dsize;
  T *restrict dst = est_data.dst;
  uint64_t idx;
  uint8_t pos[MAX_PIECES];
  struct DecInfo di;
  uint32_t sub[MAX_PIECES];
#ifdef RAMGEN
  T *restrict v = permute_v;
#endif

  uint64_t idx_cache[MAX_CANDS];
 
  for (i = thread->begin; i < thread->end; i++) {
    for (p = 0; p < num_cands;) {
      // Look for permutations that can be treated together in terms of
      // decode_piece().
      for (q = p + 1; q < num_cands; q++) {
	for (m = 0; m < num_types; m++)
	  if (pcs[type_perm_list[trylist[p]][m]] != pcs[type_perm_list[trylist[q]][m]]) break;
	if (m < num_types) break;
      }
      int l = trylist[p];
      set_dec_info(&di, &tb_entry, pcs, type_perm_list[l],
          order_list[l], 0xf, -1, LT_PIECE);
      // prefetch for j = 0
      decode_init(sub, segs[i], &di);
      decode_piece_iter(sub, pos, &di, &tb_entry);
      for (r = p; r < q; r++) {
	l = trylist[r];
        idx = calc_idx_piece(pos, pidx_list[l]);
	__builtin_prefetch(&table[idx], 0, 3);
	idx_cache[r] = idx;
      }
      for (j = 1; j < seg_size; j++) {
	// prefetch for j, copy for j - 1
        decode_piece_iter(sub, pos, &di, &tb_entry);
	for (r = p; r < q; r++) {
	  l = trylist[r];
          idx = calc_idx_piece(pos, pidx_list[l]);
	  __builtin_prefetch(&table[idx], 0, 3);
#ifndef RAMGEN
	  dst[r * dsize + i * seg_size + j - 1] = table[idx_cache[r]];
#else
	  dst[r * dsize + i * seg_size + j - 1] = v[table[idx_cache[r]]];
#endif
	  idx_cache[r] = idx;
	}
      }
      for (r = p; r < q; r++)
#ifndef RAMGEN
	dst[r * dsize + i * seg_size + j - 1] = table[idx_cache[r]];
#else
	dst[r * dsize + i * seg_size + j - 1] = v[table[idx_cache[r]]];
#endif
      p = q;
    }
  }
}
#undef NAME
