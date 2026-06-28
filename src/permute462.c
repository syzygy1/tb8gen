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

static struct Work *work_convert = nullptr;
static struct Work *work_est = nullptr;
static uint64_t slice_size = 0;

static int num_sets, num_set_perms;
static uint8_t set_perm_list[MAX_PERMS][MAX_SETS];

static uint64_t compest[MAX_PERMS];

static int trylist[MAX_CANDS];
static struct RankInfo try_ri[MAX_CANDS];
static struct RankInfo best_ri;

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
  num_sets = ri.numsets;

  generate_set_perms(num_sets);

  generate_test_list(kslice_size, g_pos.num - 2);
  work_convert = create_work(g_total_work, kslice_size, 0);
  slice_size = kslice_size;
}

static struct {
  void *src;
  void *dst;
  struct RankInfo *perm_ri;
} convert_data;

static struct {
  void *table;
  void *dst;
  int num_cands;
  uint64_t dsize;
} est_data;

#define T u8
#include "permute462_tmpl.c"
#undef T

#define T u16
#include "permute462_tmpl.c"
#undef T

static void estimate_compression_piece(void *table, int num_cands, bool wide,
    bool wdl, bool reflection)
{
  uint64_t dsize = num_segs * seg_size;
  void *dst = malloc((num_cands * dsize + 1) * (1 + wide));
  est_data.table = table;
  est_data.dst = dst;
  est_data.num_cands = num_cands;
  est_data.dsize = dsize;
  uint8_t *dst0 = dst;

  if (reflection) {
    if (!wide)
      run_threaded(convert_est_data_piece_ref_u8, work_est);
    else
      run_threaded(convert_est_data_piece_ref_u16, work_est);
  } else {
    if (!wide)
      run_threaded(convert_est_data_piece_u8, work_est);
    else
      run_threaded(convert_est_data_piece_u16, work_est);
  }

  uint64_t csize;

  for (int p = 0; p < num_cands; p++) {
    int num_syms;
    int64_t *freq = construct_pairs(dst0, dsize, 20, 100, &num_syms, wide, wdl);
    struct HuffCode *c = create_code(freq, num_syms);
    csize = calc_size(c);
    free_code(c);
    printf("[%2d]", p);
    printf(";");
    for (int i = 0; i < num_sets; i++) {
      int k = set_perm_list[trylist[p]][i];
      printf(" %c", PieceChar[g_set_pt[k]]);
    }
    printf("; %"PRIu64"\n", csize);
    compest[trylist[p]] = csize;
    dst0 += (wide ? 2 : 1) * dsize;
  }

  free(dst);
}

static int64_t estimate_compression(void *table, int *bestp, bool wide,
    bool wdl, bool reflection)
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

  for (k = 0; k < num_sets - 1; k++) {
    int pos0 = num_sets - k - 2;
    int pos1 = num_sets - k - 1;
    best = UINT64_MAX;
    num_cands = 0;
    for (p = 0; p < num_sets; p++) {
      for (i = pos1 + 1; i < num_sets; i++)
	if (p == bestperm[i]) break;
      if (i < num_sets) continue;
      for (q = 0; q < num_sets; q++) {
	if (q == p) continue;
	for (i = pos1 + 1; i < num_sets; i++)
	  if (q == bestperm[i]) break;
	if (i < num_sets) continue;
	// look for permutation ending with p,q,bestperm[pos1+1..num_sets-1]
	for (i = 0; i < num_set_perms; i++) {
	  for (j = pos1 + 1; j < num_sets; j++)
	    if (set_perm_list[i][j] != bestperm[j]) break;
	  if (j < num_sets) continue;
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
      uint8_t m[MAX_SETS] = { 0 };
      for (j = 0; j < num_sets; j++) {
        p = set_perm_list[trylist[i]][j];
        m[j] = ri.mult[p];
      }
      try_ri[i] = rank_info_62[rank_mult(m)];
      for (j = 0; j < num_sets; j++) {
        p = set_perm_list[trylist[i]][j];
        try_ri[i].first[j] = ri.first[p];
      }
    }
    estimate_compression_piece(table, num_cands, wide, wdl, reflection);
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

void permute_piece_462(void *tb_table, void *table, uint8_t *best, int type, 
    bool wide)
{
  int s = KKMap[g_slice.sq[0]][g_slice.sq[1]];
  bool reflection = s >= 441;
  uint64_t size = kslice_sizes[reflection];

  if (size != slice_size) {
    generate_test_list(size, g_pos.num - 2);
    if (work_convert)
      work_free(work_convert);
    work_convert = create_work(g_total_work, size, 0);
    slice_size = size;
  }

  int bestp;

  estimate_compression(table, &bestp, wide, type == WDL, reflection);

  for (int i = 0; i < num_sets; i++)
    best[i] = set_perm_list[bestp][i];

  uint8_t m[MAX_SETS] = { 0 };
  for (int i = 0; i < num_sets; i++)
    m[i] = ri.mult[best[i]];
  best_ri = rank_info_62[rank_mult(m)];
  for (int i = 0; i < num_sets; i++)
    best_ri.first[i] = ri.first[best[i]];

  printf("\nbest permutation: ");
  for (int i = 0; i < num_sets; i++) {
    for (int j = 0; j < best_ri.mult[i]; j++)
      printf("%c", PieceChar[g_set_pt[best[i]]]);
  }
  printf("\n");

  convert_data.src = table;
  convert_data.dst = tb_table;
  convert_data.perm_ri = &best_ri;

  if (g_compress_type == 1)
    return;

  if (reflection) {
    if (!wide)
      run_threaded(convert_data_piece_ref_u8, work_convert);
    else
      run_threaded(convert_data_piece_ref_u16, work_convert);
  } else {
    if (!wide)
      run_threaded(convert_data_piece_u8, work_convert);
    else
      run_threaded(convert_data_piece_u16, work_convert);
  }
}
