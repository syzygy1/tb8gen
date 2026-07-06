/*
  Copyright (c) 2026 Ronald de Man

  This file is distributed under the terms of the GNU GPL, version 2.
*/

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <immintrin.h>

#include "compress.h"
#include "defs.h"
#include "stats.h"
#include "rgenerate.h"
#include "tbrgen.h"
#include "threads.h"
#include "types.h"
#include "util.h"
#include "hash/xxhash.h"

uint64_t g_stats[2][MAX_STATS];
XXH128_hash_t wdl_checksum, dtz_checksum[2];

#include "stats_common.c"

static uint64_t (*per_thread_stats)[MAX_STATS];

INLINE bool is_symmetric(Bitboard b)
{
  return flip_main_diag(b) == b;
}

struct RIdxStateS {
  Bitboard occ[8];
  uint32_t sub[MAX_SETS + 1];
  bool sym[MAX_SETS + 1];
};

static bool ridx_state_sym_init(struct RIdxStateS *is, uint64_t idx,
    const struct RankInfo *ri)
{
  memset(is, 0, sizeof *is);
  for (int k = ri->numsets - 1; k >= 0; k--)
    idx = divmod_recip(idx, ri->factor[k], ri->recip[k], &is->sub[k + 1]);
  is->sub[0] = idx;
  is->occ[0] = bit(KKSquare[idx][0]) | bit(KKSquare[idx][1]);
  is->sym[0] = is_symmetric(is->occ[0]);
  for (int i = 0; i < ri->numsets; i++) {
    Bitboard occ = unrank_binomial(is->sub[i + 1], ri->mult[i], is->occ[i]);
    is->occ[i + 1] = is->occ[i] | occ;
    is->sym[i + 1] = is->sym[i] && is_symmetric(occ);
  }
  return is->sym[ri->numsets];
}

INLINE bool ridx_state_sym_inc(struct RIdxStateS *is,
    const struct RankInfo *ri)
{
  uint32_t *const restrict sub = is->sub;
  int i = ri->numsets;

  for (;;) {
    i--;
    if (++sub[i + 1] < ri->factor[i]) break;
    sub[i + 1] = 0;
    if (i == 0) {
      sub[0]++;
      is->occ[0] = bit(KKSquare[sub[0]][0]) | bit(KKSquare[sub[0]][1]);
      break;
    }
  }

  Bitboard occ = is->occ[i];
  for (; i < ri->numsets; i++) {
    occ |= unrank_binomial(is->sub[i + 1], ri->mult[i], occ);
    is->occ[i + 1] = occ;
    is->sym[i + 1] = is->sym[i] && is_symmetric(occ);
  }
  return is->sym[ri->numsets];
}

static void count_stats_worker(struct ThreadData *thread)
{
  uint64_t *restrict stats = per_thread_stats[thread->thread_id];
  uint8_t *restrict table = g_table[g_pos.stm];

  uint64_t idx = thread->begin, end = thread->end;

  if (end <= table_diagonal) {
    for (; idx < end; idx++)
      stats[table[idx]] += 2;
    return;
  }

  for (; idx < table_diagonal; idx++)
    stats[table[idx]] += 2;

  struct RIdxStateS is;
  bool symmetric = ridx_state_sym_init(&is, idx, &ri);
  for (; idx < end; idx++, symmetric = ridx_state_sym_inc(&is, &ri)) {
#if 0
    __m512i x = _mm512_load_si512((__m512i *)is.bb);
    __m512i y = flip_main_8xbb(x);
    if (_mm512_cmpeq_epi64_mask(x, y) == 0xff)
      stats[table[idx]] += 2;
    else
      stats[table[idx]]++;
#else
    stats[table[idx]] += 1 + symmetric;
#endif
  }
}

void collect_stats(int stm)
{
  per_thread_stats = aligned_alloc(64, g_num_threads * MAX_STATS * 8);
  if (!per_thread_stats)
    out_of_mem();
  memset(per_thread_stats, 0, g_num_threads * MAX_STATS * 8);
  memset(g_stats[stm], 0, sizeof g_stats[stm]);

  g_pos.stm = stm;
  run_threaded(count_stats_worker, &work_g_static);

  for (int t = 0; t < g_num_threads; t++) {
    uint64_t *restrict stats = per_thread_stats[t];
    g_stats[stm][0] += stats[255]; // illegal
    g_stats[stm][1] += stats[254]; // capt_win
    g_stats[stm][2] += stats[253]; // pawn_win
    for (int i = 1; i <= DRAW_RULE; i++)
      g_stats[stm][i + 2] += stats[253 - i];
    for (int i = 0; i <= DRAW_RULE; i++)
      g_stats[stm][MAX_STATS - 1 - i] += stats[i + 1];
    g_stats[stm][MAX_STATS / 2 + 1] += stats[127]; // capt_draw
    g_stats[stm][MAX_STATS / 2 + 2] += stats[0]; // unresolved -> draw
  }
  for (int i = 0; i < MAX_STATS; i++)
    g_stats[stm][i] >>= 1;

  free(per_thread_stats);
}
