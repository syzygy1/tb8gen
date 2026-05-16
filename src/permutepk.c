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
#include "indexp.h"
#include "kslicep.h"
#include "permutepk.h"
#include "probe.h"
#include "tb8genp.h"
#include "threads.h"
#include "util.h"

struct PKIdxInfo {
  int numsets;
  uint32_t factor[MAX_SETS + 1];
  int first[MAX_SETS + 1];
  int mult[MAX_SETS + 1];
};

struct IdxInfo pk_ii;

extern uint64_t tb_size;

static constexpr int MAX_PERMS = 1*2*3*4*5*6*7;
static constexpr int MAX_CANDS = 6*7;

static uint64_t *restrict work_convert;
static uint64_t *restrict work_est;

static int num_set_perms;
static uint8_t set_perm_list[MAX_PERMS][MAX_SETS];

static uint64_t compest[MAX_PERMS];

static uint8_t set_pt[MAX_SETS];

static int trylist[MAX_CANDS];
static struct PKIdxInfo try_ii[MAX_CANDS];
static struct PKIdxInfo best_ii;

static uint8_t perm_tmp[MAX_SETS];

void pk_idx_to_sq_init(uint64_t idx, uint32_t *restrict sub,
    const struct PKIdxInfo *ii)
{
  for (int k = ii->numsets - 1; k > 0; k--) {
    sub[k] = idx % ii->factor[k];
    idx /= ii->factor[k];
  }
  sub[0] = idx;
}

INLINE void pk_idx_to_sq_inc(uint32_t *sub, const struct PKIdxInfo *ii)
{
  for (int i = ii->numsets - 1; i >= 0 && ++sub[i] >= ii->factor[i]; i--)
    sub[i] = 0;
}

void pk_idx_to_sq(uint32_t *sub, uint8_t *restrict sq,
    const struct PKIdxInfo *ii, int stm)
{
  Bitboard occ = bit(sq[stm]) | bit(sq[2]);
  for (int i = 0; i < ii->numsets; i++)
    occ = unrank_binomial(sub[i], ii->mult[i], sq + ii->first[i], occ);
}

uint64_t pk_sq_to_idx(uint8_t *restrict sq, int stm)
{
  Bitboard occ = bit(sq[stm]) | bit(sq[2]);

  uint64_t idx = 0;
  for (int k = 0; k < pk_ii.numsets; k++) {
    int i = pk_ii.first[k];
    sort_squares(pk_ii.mult[k], &sq[i]);
    size_t s = 0;
    Bitboard occ2 = occ;
    for (int j = 0; j < pk_ii.mult[k]; i++, j++) {
      int rank = rank_among_free(sq[i], occ);
      occ2 |= bit(sq[i]);
      s += Binomial[j + 1][rank];
    }
    idx = idx * pk_ii.factor[k] + s;
    occ = occ2;
  }

  return idx;
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

void init_permute_pawn_pk(int stm)
{
  assume(ii.numsets + 1 <= MAX_SETS);
  g_pos.stm = stm;
  pk_ii.numsets = ii.numsets + 1;
  pk_ii.first[0] = stm ^ 1;
  pk_ii.mult[0] = 1;
  for (int i = 0; i < ii.numsets; i++) {
    pk_ii.first[i + 1] = ii.first[i];
    pk_ii.mult[i + 1] = ii.mult[i];
  }
  for (int i = 0, n = 62; i < pk_ii.numsets; i++) {
    pk_ii.factor[i] = Binomial[pk_ii.mult[i]][n];
    n -= pk_ii.mult[i];
  }

  generate_set_perms(pk_ii.numsets);

  for (int i = 0; i < pk_ii.numsets; i++)
    set_pt[i] = g_pos.pt[pk_ii.first[i]];

  tb_size = 62 * kslice_size;
  generate_test_list(tb_size, g_pos.num - 2);
  if (!work_convert)
    work_convert = create_work(g_total_work, tb_size, 0);
}

static struct {
  void *src;
  void *dst;
  struct PKIdxInfo *perm_ii;
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
#include "permutepk_tmpl.c"
#undef T

#define T u16
#include "permutepk_tmpl.c"
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

  if (num_segs > 1) {
    if (!wide)
      run_threaded(convert_est_data_pawn_u8, work_est, 0);
    else
      run_threaded(convert_est_data_pawn_u16, work_est, 0);
  }
  else {
    if (!wide)
      run_single(convert_est_data_pawn_u8, work_est, 0);
    else
      run_single(convert_est_data_pawn_u16, work_est, 0);
  }

  uint64_t csize;

  for (int p = 0; p < num_cands; p++) {
    int num_syms;
    int64_t *freq = construct_pairs(dst0, dsize, 20, 100, &num_syms, wide, wdl);
    struct HuffCode *c = create_code(freq, num_syms);
    csize = calc_size(c);
    free_code(c);
    printf("[%2d]", p);
    printf("; perm:");
    for (int i = 0; i < pk_ii.numsets; i++)
      printf(" %2d", set_perm_list[trylist[p]][i]);
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
    work_est = alloc_work(g_total_work);
  fill_work(g_total_work, num_segs, 0, work_est);

  for (i = 0; i < num_set_perms; i++)
    compest[i] = 0;

  for (k = 0; k < pk_ii.numsets - 1; k++) {
    best = UINT64_MAX;
    num_cands = 0;
    for (p = 0; p < pk_ii.numsets; p++) {
      for (i = 0; i < k; i++)
	if (p == bestperm[i]) break;
      if (i < k) continue;
      for (q = 0; q < pk_ii.numsets; q++) {
	if (q == p) continue;
	for (i = 0; i < k; i++)
	  if (q == bestperm[i]) break;
	if (i < k) continue;
	// look for permutation starting with bestperm[0..k-1],p,q
	for (i = 0; i < num_set_perms; i++) {
	  for (j = 0; j < k; j++)
	    if (set_perm_list[i][j] != bestperm[j]) break;
	  if (j < k) continue;
	  if (   set_perm_list[i][k]   == p
              && set_perm_list[i][k+1] == q) break;
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
      try_ii[i].numsets = pk_ii.numsets;
      for (j = 0; j < pk_ii.numsets; j++) {
        p = set_perm_list[trylist[i]][j];
        try_ii[i].mult[j] = pk_ii.mult[p];
        try_ii[i].first[j] = pk_ii.first[p];
      }
      for (int l = 0, n = 62; l < try_ii[i].numsets; l++) {
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
    bestperm[k] = set_perm_list[bp][k];
  }
  *bestp = bp;

  return best;
}

void permute_pawn_pk(void *tb_table, void *table, uint8_t *best, int type, 
    bool wide)
{
  int bestp;

  estimate_compression(table, &bestp, wide, type == WDL);

  for (int i = 0; i < pk_ii.numsets; i++)
    best[i] = set_perm_list[bestp][i];

  best_ii.numsets = pk_ii.numsets;
  for (int i = 0; i < pk_ii.numsets; i++) {
    int k = best[i];
    best_ii.mult[i] = pk_ii.mult[k];
    best_ii.first[i] = pk_ii.first[k];
  }
  for (int i = 0, n = 62; i < best_ii.numsets; i++) {
    best_ii.factor[i] = Binomial[best_ii.mult[i]][n];
    n -= best_ii.mult[i];
  }

  printf("\nbest permutation: ");
  for (int i = 0; i < pk_ii.numsets; i++) {
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
    run_threaded(convert_data_pawn_u8, work_convert, 1);
  else
    run_threaded(convert_data_pawn_u16, work_convert, 1);
}
