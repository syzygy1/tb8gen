/*
  Copyright (c) 2011-2017, 2025, 2026 Ronald de Man

  This file is distributed under the terms of the GNU GPL, version 2.
*/

#include <assert.h>
#include <fcntl.h>
#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "defs.h"
#include "types.h"
#include "compress.h"
#include "huffman.h"
#include "index.h"
#include "kslicep.h"
#include "permutep.h"
#include "probe.h"
#include "tb8gen.h"
#include "threads.h"
#include "util.h"

struct PIdxInfo {
  int numsets;
  uint32_t factor[MAX_SETS + 2];
  int first[MAX_SETS + 2];
  int mult[MAX_SETS + 2];
};

struct PIdxState {
  uint32_t sub[MAX_SETS + 2];
  Bitboard occ[MAX_SETS + 2];
  alignas(64) Bitboard bb[8];
  int n;
};

struct PIdxInfo p_ii;

uint64_t tb_size;

static constexpr int MAX_PERMS = 1*2*3*4*5*6*7;
static constexpr int MAX_CANDS = 6*7;

static struct Work *work_convert;
static struct Work *work_est;

static int num_set_perms;
static uint8_t set_perm_list[MAX_PERMS][MAX_SETS];

static uint64_t compest[MAX_PERMS];

static uint8_t set_pt[MAX_SETS];

static int trylist[MAX_CANDS];
static struct PIdxInfo try_ii[MAX_CANDS];
static struct PIdxInfo best_ii;

static uint8_t perm_tmp[MAX_SETS];

void p_idx_state_init(struct PIdxState *is, uint64_t idx,
    const struct PIdxInfo *ii)
{
  for (int k = ii->numsets - 1; k > 0; k--) {
    is->sub[k] = idx % ii->factor[k];
    idx /= ii->factor[k];
  }
  is->sub[0] = idx;
  is->n = 0;
  is->bb[0] = bit(g_slice.sq[2]);
  is->occ[0] = is->bb[0];
}

INLINE void p_idx_state_inc(struct PIdxState *is, const struct PIdxInfo *ii)
{
  uint32_t *restrict sub = is->sub;
  int i = ii->numsets - 1;
  for (; ++sub[i] >= ii->factor[i] && i > 0; i--)
    sub[i] = 0;
  is->n = i;
}

void p_idx_state_to_bb(struct PIdxState *is, const struct PIdxInfo *ii)
{
  int i = is->n;
  Bitboard occ = is->occ[i];
  for (; i < ii->numsets; i++) {
    is->occ[i] = occ;
    is->bb[i + 1] = unrank_binomial2(is->sub[i], ii->mult[i], occ);
    occ |= is->bb[i + 1];
  }
}

static void init_source_rank_ri_p(struct RankInfo *rank_ri,
    int8_t king_perm[2], const struct PIdxInfo *perm_ii)
{
  *rank_ri = ri;
  king_perm[0] = king_perm[1] = -1;
  for (int j = 0; j < perm_ii->numsets; j++) {
    if (perm_ii->first[j] == 0)
      king_perm[0] = j + 1;
    else if (perm_ii->first[j] == 1)
      king_perm[1] = j + 1;
  }
  assert(king_perm[0] >= 0 && king_perm[1] >= 0);

  for (int k = 0; k < ri.numsets; k++) {
    int j = 0;
    for (; j < perm_ii->numsets; j++)
      if (ri.first[k] == perm_ii->first[j])
        break;
    assert(j < perm_ii->numsets);
    rank_ri->perm[k] = j + 1;
  }
}

uint64_t p_bb_to_idx(const struct PIdxState *is,
    const struct RankInfo *rank_ri, const int8_t king_perm[2])
{
  int sq0 = lsb(is->bb[king_perm[0]]);
  int sq1 = lsb(is->bb[king_perm[1]]);
  int psq = g_slice.sq[2];
  Bitboard occ = bit(sq0) | bit(sq1) | bit(psq);

  int s0 = (sq0 > psq);
  int s1 = (sq1 > psq) + (sq1 > sq0);
  uint64_t idx = (sq0 - s0) * 62 + (sq1 - s1);
  return perm_rank_bb_from(is->bb, idx, 0, occ, rank_ri);
}

static void generate_set_perms_helper(int n, int k)
{
  if (k == n) {
    for (int i = 0; i < n; i++)
      set_perm_list[num_set_perms][i] = perm_tmp[i];
    num_set_perms++;
    return;
  }

  for (int i = 0; i < n; i++)
    if (perm_tmp[i] == 0xff) {
      perm_tmp[i] = k;
      generate_set_perms_helper(n, k + 1);
      perm_tmp[i] = 0xff;
    }
}

static void generate_set_perms(int n)
{
  num_set_perms = 0;
  for (int i = 0; i < n; i++)
    perm_tmp[i] = 0xff;
  generate_set_perms_helper(n, 0);
}

static int num_segs, seg_size;
static uint64_t *restrict segs = nullptr;

static constexpr int NUM_SEGS = 1000;
static constexpr int SEG_SIZE = 64 * 256;

static void generate_test_list(uint64_t size, int n)
{
  if (segs) free(segs);

  if (n <= 3 || size <= 100000) {
    // 1 entry covering whole table
    num_segs = 1;
    segs = malloc(sizeof(uint64_t));
    seg_size = size;
    segs[0] = 0;
    return;
  }

  if (n <= 3 || (NUM_SEGS + 1) * SEG_SIZE >= size) {
    num_segs = 100;
    seg_size = size / num_segs;
  } else {
    num_segs = n > 5 ? NUM_SEGS : NUM_SEGS / 2;
    seg_size = n > 5 ? SEG_SIZE : SEG_SIZE / 1;
  }

  uint64_t max = size - num_segs * (uint64_t)seg_size + 1;
  segs = malloc(sizeof(uint64_t) * num_segs);
  for (int i = 0; i < num_segs; i++)
    segs[i] = llrand() % max;
  for (int i = 0; i < num_segs; i++)
    for (int j = i + 1; j < num_segs; j++)
      if (segs[i] > segs[j])
        Swap(segs[i], segs[j]);
  for (int i = 0; i < num_segs; i++)
    segs[i] += i * seg_size;
}

void init_permute_pawn_p(void)
{
  assume(ri.numsets + 2 <= MAX_SETS);
  //int stm = g_slice.stm;
  p_ii.numsets = ri.numsets + 2;
  p_ii.first[0] = 0; // stm ?
  p_ii.mult[0] = 1;
  p_ii.first[1] = 1; // stm^1 ?
  p_ii.mult[1] = 1;
  for (int i = 0; i < ri.numsets; i++) {
    p_ii.first[i + 2] = ri.first[i];
    p_ii.mult[i + 2] = ri.mult[i];
  }
  for (int i = 0, n = 63; i < p_ii.numsets; i++) {
    p_ii.factor[i] = Binomial[p_ii.mult[i]][n];
    n -= p_ii.mult[i];
  }

  generate_set_perms(p_ii.numsets);

  for (int i = 0; i < p_ii.numsets; i++)
    set_pt[i] = g_pos.pt[p_ii.first[i]];

  tb_size = 63*62 * kslice_size;
  generate_test_list(tb_size, g_pos.num - 1);
  if (!work_convert)
    work_convert = create_work(g_total_work, tb_size, 0);
}

static struct {
  void *src;
  void *dst;
  struct PIdxInfo *perm_ii;
  int rank;
} convert_data;

static struct {
  void *table;
  void *dst;
  int num_cands;
  uint64_t dsize;
  int rank;
} est_data;

#define T u8
#include "permutep_tmpl.c"
#undef T

#define T u16
#include "permutep_tmpl.c"
#undef T

static void estimate_compression_pawn(void *table, int num_cands, bool wide,
    bool wdl)
{
  uint64_t dsize = num_segs * seg_size;
  void *dst = malloc((num_cands * dsize + 1) * (1 + wide));
  est_data.table = table;
  est_data.dst = dst;
  est_data.num_cands = num_cands;
  est_data.dsize = dsize;
  uint8_t *dst0 = dst;

  if (!wide)
    run_threaded(convert_est_data_pawn_u8, work_est);
  else
    run_threaded(convert_est_data_pawn_u16, work_est);

  uint64_t csize;

  for (int p = 0; p < num_cands; p++) {
    int num_syms;
    int64_t *freq = construct_pairs(dst0, dsize, 20, 100, &num_syms, wide, wdl);
    struct HuffCode *c = create_code(freq, num_syms);
    csize = calc_size(c);
    free_code(c);
    printf("[%2d]", p);
    printf(";");
    for (int i = 0; i < p_ii.numsets; i++) {
      int k = set_perm_list[trylist[p]][i];
      printf(" %c", PieceChar[set_pt[k]]);
    }
    printf("; %"PRIu64"\n", csize);
    compest[trylist[p]] = csize;
    dst0 += (wide ? 2 : 1) * dsize;
  }

  free(dst);
}

static int64_t estimate_compression(void *table, int *bestp, bool wide,
    bool wdl)
{
  int i, j, k, p, q;
  int num_cands, bp = 0;
  uint64_t best = 0;
  uint8_t bestperm[MAX_SETS];

  if (g_compress_type == 1) {
    *bestp = 0;
    return 0;
  }

  if (!work_est)
    work_est = create_work(g_total_work, num_segs, 0);
  else
    work_refill_units(work_est, g_total_work, num_segs, 0);

  for (i = 0; i < num_set_perms; i++)
    compest[i] = 0;

  for (k = 0; k < p_ii.numsets - 1; k++) {
    int pos0 = p_ii.numsets - k - 2;
    int pos1 = p_ii.numsets - k - 1;
    best = UINT64_MAX;
    num_cands = 0;
    for (p = 0; p < p_ii.numsets; p++) {
      for (i = pos1 + 1; i < p_ii.numsets; i++)
	if (p == bestperm[i]) break;
      if (i < p_ii.numsets) continue;
      for (q = 0; q < p_ii.numsets; q++) {
	if (q == p) continue;
	for (i = pos1 + 1; i < p_ii.numsets; i++)
	  if (q == bestperm[i]) break;
	if (i < p_ii.numsets) continue;
	// look for permutation ending with p,q,bestperm[pos1+1..p_ii.numsets-1]
	for (i = 0; i < num_set_perms; i++) {
	  for (j = pos1 + 1; j < p_ii.numsets; j++)
	    if (set_perm_list[i][j] != bestperm[j]) break;
	  if (j < p_ii.numsets) continue;
	  if (   set_perm_list[i][pos0] == p
              && set_perm_list[i][pos1] == q) break;
	}
	if (i < num_set_perms) {
	  if (compest[i]) {
	    if (compest[i] < best) {
	      best = compest[i];
	      bp = i;
	    }
	  } else
	    trylist[num_cands++] = i;
	}
      }
    }
    for (i = 0; i < num_cands; i++)
      for (j = i + 1; j < num_cands; j++)
	if (trylist[i] > trylist[j])
          Swap(trylist[i], trylist[j]);
    for (i = 0; i < num_cands; i++) {
      try_ii[i].numsets = p_ii.numsets;
      for (j = 0; j < p_ii.numsets; j++) {
        p = set_perm_list[trylist[i]][j];
        try_ii[i].mult[j] = p_ii.mult[p];
        try_ii[i].first[j] = p_ii.first[p];
      }
      for (int l = 0, n = 63; l < try_ii[i].numsets; l++) {
        try_ii[i].factor[l] = Binomial[try_ii[i].mult[l]][n];
        n -= try_ii[i].mult[l];
      }
    }
    estimate_compression_pawn(table, num_cands, wide, wdl);
    for (i = 0; i < num_cands; i++) {
      if (compest[trylist[i]] < best) {
	best = compest[trylist[i]];
	bp = trylist[i];
      }
    }
    bestperm[pos1] = set_perm_list[bp][pos1];
  }
  *bestp = bp;

  return best;
}

void permute_pawn_p(void *tb_table, void *table, uint8_t *best, int type, 
    bool wide)
{
  int bestp;

  estimate_compression(table, &bestp, wide, type == WDL);

  for (int i = 0; i < p_ii.numsets; i++)
    best[i] = set_perm_list[bestp][i];

  best_ii.numsets = p_ii.numsets;
  for (int i = 0; i < p_ii.numsets; i++) {
    int k = best[i];
    best_ii.mult[i] = p_ii.mult[k];
    best_ii.first[i] = p_ii.first[k];
  }
  for (int i = 0, n = 63; i < best_ii.numsets; i++) {
    best_ii.factor[i] = Binomial[best_ii.mult[i]][n];
    n -= best_ii.mult[i];
  }

  printf("\nbest permutation: ");
  for (int i = 0; i < p_ii.numsets; i++) {
    for (int j = 0; j < best_ii.mult[i]; j++)
      printf("%c", PieceChar[set_pt[best[i]]]);
  }
  printf("\n");

  convert_data.src = table;
  convert_data.dst = tb_table;
  convert_data.perm_ii = &best_ii;

  if (g_compress_type == 1)
    return;

  if (!wide)
    run_threaded(convert_data_pawn_u8, work_convert);
  else
    run_threaded(convert_data_pawn_u16, work_convert);
}
