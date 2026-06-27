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
#include "kslice.h"
#include "permute10.h"
#include "probe.h"
#include "tb8gen.h"
#include "threads.h"
#include "util.h"

struct P10IdxInfo {
  int numsets;
  int k2;
  uint32_t factor[MAX_SETS + 1];
  int first[MAX_SETS];
  int mult[MAX_SETS];
};

struct P10IdxState {
  uint32_t sub[MAX_SETS + 1];
  Bitboard occ[MAX_SETS + 1];
  alignas(64) Bitboard bb[8];
  uint8_t sq[2];
  int n;
};

extern uint64_t tb_size;

static constexpr int MAX_PERMS = 1*2*3*4*5*6*7;
static constexpr int MAX_CANDS = 6*7;

static struct Work *work_convert;
static struct Work *work_est;

static int num_set_perms;
static uint8_t set_perm_list[MAX_PERMS][MAX_SETS];

static uint64_t compest[MAX_PERMS];

static uint8_t set_pt[MAX_SETS + 1];

static int trylist[MAX_CANDS];
static struct P10IdxInfo try_ii[MAX_CANDS];
static struct P10IdxInfo best_ii;

static uint8_t perm_tmp[MAX_SETS];

static uint8_t P10Square[64];
static bool P10Ref[64];
static uint8_t InvSquare[64];

static int32_t current_size = -1;

void p10_idx_state_init(struct P10IdxState *is, uint64_t idx,
    const struct P10IdxInfo *ii, int stm)
{
  for (int k = ii->numsets - 1; k > 0; k--) {
    is->sub[k] = idx % ii->factor[k];
    idx /= ii->factor[k];
  }
  is->sub[0] = idx;
  is->sq[stm] = g_slice.sq[stm];
  is->sq[stm ^ 1] = InvSquare[is->sub[ii->k2]];
  is->n = 0;
  is->bb[0] = bit(is->sq[0]) | bit(is->sq[1]);
  is->occ[0] = is->bb[0];
}

INLINE void p10_idx_state_inc(struct P10IdxState *is,
    const struct P10IdxInfo *ii)
{
  uint32_t *restrict sub = is->sub;
  int i = ii->numsets - 1;
  for (; i >= 0 && ++sub[i] >= ii->factor[i]; i--)
    sub[i] = 0;
  is->n = i <= ii->k2 ? 0 : i;
}

void p10_idx_state_to_bb(struct P10IdxState *is, const struct P10IdxInfo *ii,
    int stm)
{
  int i = is->n;
  Bitboard occ;

  if (i == 0) {
    is->sq[stm ^ 1] = InvSquare[is->sub[ii->k2]];
    occ = bit(is->sq[0]) | bit(is->sq[1]);
    is->bb[0] = occ;
    is->occ[0] = occ;
  } else
    occ = is->occ[i];

  for (; i < ii->numsets; i++) {
    is->occ[i] = occ;
    if (ii->mult[i] == 0) {
      is->bb[i + 1] = 0;
      continue;
    }
    is->bb[i + 1] = unrank_binomial2(is->sub[i], ii->mult[i], occ);
    occ |= is->bb[i + 1];
  }
}

static void init_source_rank_ri_10(struct RankInfo *rank_ri,
    const struct P10IdxInfo *perm_ii)
{
  *rank_ri = ri;
  for (int k = 0; k < ri.numsets; k++) {
    int j = 0;
    for (; j < perm_ii->numsets; j++)
      if (perm_ii->mult[j] && ri.first[k] == perm_ii->first[j])
        break;
    assert(j < perm_ii->numsets);
    rank_ri->perm[k] = j + 1;
  }
}

uint64_t p10_bb_to_idx(struct P10IdxState *is,
    const struct RankInfo *rank_ri, int k2sq)
{
  assert(P10Square[k2sq] != 0xff);

  alignas(64) Bitboard bb[8];
  const Bitboard *set_bb = is->bb;
  int t = KK_transform[is->sq[0]][is->sq[1]];
  if (t) {
    transform_set_bb(t, is->bb, bb);
    set_bb = bb;
  }

  uint64_t idx = P10Square[k2sq];
  return  !P10Ref[k2sq]
        ? perm_rank_bb_from(set_bb, idx, 0, set_bb[0], rank_ri)
        : idx * kslice_size + perm_rank_bb_ref(set_bb, rank_ri);
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

void init_permute_piece_10(int k)
{
  static const int num[] = { 58, 58, 58, 55, 55, 55, 33, 30, 30, 30 };

  if (current_size < 0) {
    generate_set_perms(ri.numsets + 1);

    for (int i = 0; i < ri.numsets; i++)
      set_pt[i + 1] = g_set_pt[i];
    set_pt[0] = g_pos.pt[g_slice.stm ^ 1];
  }

  if (num[k] != current_size) {
    current_size = num[k];
    tb_size = current_size * kslice_size;
    generate_test_list(tb_size, g_pos.num - 2);
    work_convert = create_work(g_total_work, tb_size, 0);
  }

  int n = 0;
  for (int l = 0; l < 64; l ++)
    if (KKIdx[k][l] >= 0) {
      P10Square[l] = n;
      P10Ref[l] = KKIdx[k][l] >= 441;
      InvSquare[n++] = l;
    } else
      P10Square[l] = 0xff;
  assert(n == current_size);
}

static struct {
  void *src;
  void *dst;
  struct P10IdxInfo *perm_ii;
} convert_data;

static struct {
  void *table;
  void *dst;
  int num_cands;
  uint64_t dsize;
} est_data;

#define T u8
#include "permute10_tmpl.c"
#undef T

#define T u16
#include "permute10_tmpl.c"
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

  if (!wide)
    run_threaded(convert_est_data_piece_u8, work_est);
  else
    run_threaded(convert_est_data_piece_u16, work_est);

  uint64_t csize;

  for (int p = 0; p < num_cands; p++) {
    int num_syms;
    int64_t *freq = construct_pairs(dst0, dsize, 20, 100, &num_syms, wide, wdl);
    struct HuffCode *c = create_code(freq, num_syms);
    csize = calc_size(c);
    free_code(c);
    printf("[%2d]", p);
    printf(";");
    for (int i = 0; i < ri.numsets + 1; i++) {
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

  for (k = 0; k < ri.numsets; k++) {
    int numsets = ri.numsets + 1;
    int pos0 = numsets - k - 2;
    int pos1 = numsets - k - 1;
    best = UINT64_MAX;
    num_cands = 0;
    for (p = 0; p < numsets; p++) {
      for (i = pos1 + 1; i < numsets; i++)
	if (p == bestperm[i]) break;
      if (i < numsets) continue;
      for (q = 0; q < numsets; q++) {
	if (q == p) continue;
	for (i = pos1 + 1; i < numsets; i++)
	  if (q == bestperm[i]) break;
	if (i < numsets) continue;
	// look for permutation ending with p,q,bestperm[pos1+1..numsets-1]
	for (i = 0; i < num_set_perms; i++) {
	  for (j = pos1 + 1; j < numsets; j++)
	    if (set_perm_list[i][j] != bestperm[j]) break;
	  if (j < numsets) continue;
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
      try_ii[i].numsets = ri.numsets + 1;
      for (j = 0; j < ri.numsets + 1; j++) {
        p = set_perm_list[trylist[i]][j];
        try_ii[i].mult[j] = p == 0 ? 0 : ri.mult[p - 1];
        try_ii[i].first[j] = p == 0 ? 0 : ri.first[p - 1];
        if (p == 0)
          try_ii[i].k2 = j;
      }
      for (int l = 0, n = 62; l < try_ii[i].numsets; l++) {
        int m = try_ii[i].mult[l];
        try_ii[i].factor[l] = m == 0 ? current_size : Binomial[m][n];
        n -= m;
      }
    }
    estimate_compression_piece(table, num_cands, wide, wdl);
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

void permute_piece_10(void *tb_table, void *table, uint8_t *best, int type, 
    bool wide)
{
  int bestp;

  estimate_compression(table, &bestp, wide, type == WDL);

  for (int i = 0; i < ri.numsets + 1; i++)
    best[i] = set_perm_list[bestp][i];

  best_ii.numsets = ri.numsets + 1;
  for (int i = 0; i < best_ii.numsets; i++) {
    int k = best[i];
    best_ii.mult[i] = k == 0 ? 0 : ri.mult[k - 1];
    best_ii.first[i] = k == 0 ? 0 : ri.first[k - 1];
    if (k == 0)
      best_ii.k2 = i;
  }
  for (int i = 0, n = 62; i < best_ii.numsets; i++) {
    int m = best_ii.mult[i];
    best_ii.factor[i] = m == 0 ? current_size : Binomial[m][n];
    n -= best_ii.mult[i];
  }

  printf("\nbest permutation: ");
  for (int i = 0; i < ri.numsets + 1; i++) {
    for (int j = 0; j < max(1, best_ii.mult[i]); j++)
      printf("%c", PieceChar[set_pt[best[i]]]);
  }
  printf("\n");

  convert_data.src = table;
  convert_data.dst = tb_table;
  convert_data.perm_ii = &best_ii;

  if (g_compress_type == 1)
    return;

  if (!wide)
    run_threaded(convert_data_piece_u8, work_convert);
  else
    run_threaded(convert_data_piece_u16, work_convert);
}
