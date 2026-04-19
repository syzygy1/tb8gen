/*
  Copyright (c) 2011-2017, 2025, 2026 Ronald de Man

  This file is distributed under the terms of the GNU GPL, version 2.
*/

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
#include "kslice.h"
#include "permute462.h"
#include "probe.h"
#include "tb8gen.h"
#include "threads.h"
#include "util.h"

static constexpr int MAX_PERMS = 1*2*3*4*5*6*7;
static constexpr int MAX_CANDS = 6*7;

static uint64_t *restrict work_convert = nullptr;
static uint64_t *restrict work_est = nullptr;

static int num_sets, num_set_perms;
static uint8_t set_perm_list[MAX_PERMS][MAX_SETS];

static uint64_t compest[MAX_PERMS];

static uint8_t set_pt[MAX_SETS];

static int trylist[MAX_CANDS];
static struct IdxInfo try_ii[MAX_CANDS];
static struct IdxInfo best_ii;

static uint8_t perm_tmp[MAX_SETS];

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

  if (n <= 4 || (NUM_SEGS + 1) * SEG_SIZE >= size) {
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

void init_permute_piece_462(void)
{
  num_sets = ii.numsets;

  generate_set_perms(num_sets);

  for (int i = 0; i < num_sets; i++)
    set_pt[i] = g_pos.pt[ii.first[i]];

  generate_test_list(kslice_size, g_pos.num - 2);

  work_convert = create_work(g_total_work, kslice_size, 0);
}

static struct {
  void *src;
  void *dst;
  struct IdxInfo *perm_ii;
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
#include "permute462_tmpl.c"
#undef T

#define T u16
#include "permute462_tmpl.c"
#undef T

static void estimate_compression_piece(void *table, int num_cands, bool wide,
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
      run_threaded(convert_est_data_piece_u8, work_est, 0);
    else
      run_threaded(convert_est_data_piece_u16, work_est, 0);
  }
  else {
    if (!wide)
      run_single(convert_est_data_piece_u8, work_est, 0);
    else
      run_single(convert_est_data_piece_u16, work_est, 0);
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
    for (int i = 0; i < num_sets; i++)
      printf(" %2d", set_perm_list[trylist[p]][i]);
    printf("; %"PRIu64"\n", csize);
    compest[trylist[p]] = csize;
    dst0 += (wide ? 2 : 1) * dsize;
  }

  free(dst);
}

static int64_t estimate_compression(void *table, int *bestp, int rank,
    bool wide, bool wdl)
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

  for (k = 0; k < num_sets - 1; k++) {
    best = UINT64_MAX;
    num_cands = 0;
    for (p = 0; p < num_sets; p++) {
      for (i = 0; i < k; i++)
	if (p == bestperm[i]) break;
      if (i < k) continue;
      for (q = 0; q < num_sets; q++) {
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
      try_ii[i].numsets = num_sets;
      for (j = 0; j < num_sets; j++) {
        p = set_perm_list[trylist[i]][j];
        try_ii[i].mult[j] = ii.mult[p];
        try_ii[i].first[j] = ii.first[p];
        try_ii[i].last[j] = ii.last[p];
      }
      calc_factors(&try_ii[i]);
    }
    if (rank < 0)
      estimate_compression_piece(table, num_cands, wide, wdl);
//    else
//      estimate_compression_pawn(table, pcs, rank, num_cands, wide);
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

void permute_piece_462(void *tb_table, void *table, uint8_t *best, int type, 
    bool wide)
{
  int bestp;

  estimate_compression(table, &bestp, -1, wide, type == WDL);

  for (int i = 0; i < num_sets; i++)
    best[i] = set_perm_list[bestp][i];

  best_ii.numsets = num_sets;
  for (int i = 0; i < num_sets; i++) {
    int k = best[i];
    best_ii.mult[i] = ii.mult[k];
    best_ii.first[i] = ii.first[k];
  }
  calc_factors(&best_ii);

  printf("\nbest permutation: ");
  for (int i = 0; i < num_sets; i++) {
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
    run_threaded(convert_data_piece_u8, work_convert, 1);
  else
    run_threaded(convert_data_piece_u16, work_convert, 1);
}
