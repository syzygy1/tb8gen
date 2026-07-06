/*
  Copyright (c) 2026 Ronald de Man

  This file is distributed under the terms of the GNU GPL, version 2.
*/

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "compress.h"
#include "defs.h"
#include "reduce.h"
#include "rstats.h"
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

static int stats_first_win;
static int stats_win_limit;
static int stats_first_loss;
static int stats_loss_limit;
static bool stats_include_epoch0_special;
static bool stats_include_unresolved;

static void collect_stats_range(int stm);

void reset_stats(void)
{
  memset(g_stats, 0, sizeof g_stats);
}

static bool stats_want_win(int n)
{
  return n >= stats_first_win && n < stats_win_limit;
}

static bool stats_want_loss(int n)
{
  return n >= stats_first_loss && n < stats_loss_limit;
}

static int byte_to_stat(uint8_t b)
{
  switch (b) {
  case RAM_UNRESOLVED:
    return stats_include_unresolved ? MAX_STATS / 2 + 2 : -1;
  case RAM_ILLEGAL:
    return stats_include_epoch0_special ? 0 : -1;
  case RAM_CAPT_WIN:
    return stats_include_epoch0_special ? 1 : -1;
  case RAM_PAWN_WIN:
    return epoch == 0 && stats_include_epoch0_special ? 2 : -1;
  case RAM_CAPT_CWIN:
    return epoch == 0 && stats_include_epoch0_special ? DRAW_RULE + 3 : -1;
  case RAM_PAWN_CWIN:
    return epoch == 0 && stats_include_epoch0_special ? DRAW_RULE + 4 : -1;
  case RAM_CAPT_DRAW:
    return stats_include_epoch0_special ? MAX_STATS / 2 + 1 : -1;
  case RAM_PAWN_DRAW:
    return stats_include_epoch0_special ? MAX_STATS / 2 + 2 : -1;
  default:
    break;
  }

  if (epoch == 0) {
    if (b >= RAM_LOSS_IN_0 && b < RAM_PAWN_DRAW) {
      int n = b - RAM_LOSS_IN_0;
      return stats_want_loss(n) ? loss_to_stat(n) : -1;
    }
    if (b > RAM_CAPT_DRAW && b < RAM_CAPT_CWIN) {
      int n = 253 - b - 2;
      return stats_want_win(n) ? win_to_stat(n) : -1;
    }
    if (b > RAM_CAPT_CWIN && b < RAM_PAWN_WIN) {
      int n = 253 - b;
      return stats_want_win(n) ? win_to_stat(n) : -1;
    }
    return -1;
  }

  if (b > RAM_REDUCED_BLOSS && b < RAM_CAPT_DRAW) {
    int n = b + reduce_cnt_loss[epoch - 1];
    return stats_want_loss(n) ? loss_to_stat(n) : -1;
  }
  if (b > RAM_CAPT_DRAW && b < RAM_REDUCED_CWIN) {
    int n = reduce_cnt_win[epoch - 1] - b;
    return stats_want_win(n) ? win_to_stat(n) : -1;
  }

  return -1;
}

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
  if (epoch == 0) {
    stats_first_win = 1;
    stats_first_loss = 0;
  } else {
    stats_first_win = reduce_cnt_win[epoch - 1] - 250;
    stats_first_loss = reduce_cnt_loss[epoch - 1] + 3;
  }
  stats_win_limit = MAX_STATS / 2 - 3;
  stats_loss_limit = MAX_STATS / 2 - 3;
  stats_include_epoch0_special = epoch == 0;
  stats_include_unresolved = true;
  collect_stats_range(stm);
}

void collect_stats_before_reduce(int stm, int n)
{
  if (epoch == 0) {
    stats_first_win = 1;
    stats_first_loss = 0;
  } else {
    stats_first_win = reduce_cnt_win[epoch - 1] - 250;
    stats_first_loss = reduce_cnt_loss[epoch - 1] + 3;
  }
  stats_win_limit = n - 1;
  stats_loss_limit = n;
  stats_include_epoch0_special = epoch == 0;
  stats_include_unresolved = false;
  collect_stats_range(stm);
}

static void collect_stats_range(int stm)
{
  per_thread_stats = aligned_alloc(64, g_num_threads * MAX_STATS * 8);
  if (!per_thread_stats)
    out_of_mem();
  memset(per_thread_stats, 0, g_num_threads * MAX_STATS * 8);

  g_pos.stm = stm;
  run_threaded(count_stats_worker, &work_g_static);

  uint64_t add[MAX_STATS] = { 0 };
  for (int t = 0; t < g_num_threads; t++) {
    uint64_t *restrict stats = per_thread_stats[t];
    for (int b = 0; b < 256; b++) {
      int s = byte_to_stat(b);
      if (s >= 0)
        add[s] += stats[b];
    }
  }
  for (int i = 0; i < MAX_STATS; i++)
    g_stats[stm][i] += add[i] >> 1;

  free(per_thread_stats);
}
