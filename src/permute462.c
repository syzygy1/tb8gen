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
#include "permute462.h"
#include "probe.h"
#include "tb8gen.h"
#include "threads.h"
#include "util.h"

static constexpr int MAX_PERMS = 1*2*3*4*5*6*7;
static constexpr int MAX_CANDS = 6*7;

static uint64_t *restrict work_convert = NULL;
static uint64_t *restrict work_est = NULL;

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
static uint64_t *restrict segs = NULL;

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

#if 0
// sq_to_idx()
INLINE uint64_t calc_idx_piece(uint8_t *restrict sq, uint8_t *restrict pidx,
    int n)
{
  // need to think about this
  //
  uint8_t sq2[MAX_PIECES];
  for (int i = 0; i < n; i++)
    sq2[i] = sq[pidx[i]];
  int s = KKMap[sq2[0]][sq2[1]];
  if (s < 0)
    return 462 * kslice_size;
    // table[462 * kslice_size] is set to the don't care value.
  else {
    // FIXME: improve efficiency
    uint8_t sq3[MAX_PIECES];
    normalize(sq2, sq3);
    return (unsigned)s * kslice_size + sq_to_idx(sq3);
  }
}
#endif

#if 0
INLINE uint64_t calc_idx_pawn(uint8_t *pos, uint8_t *pidx, int n)
{
  uint64_t idx = pos[pidx[1]];
  for (int m = 2; m < n; m++)
    idx = (idx << 6) | pos[pidx[m]];
  idx ^= sq_mask_pawn[pos[pidx[0]]];
  return idx;
}
#endif

void init_permute_piece_462(void)
{
  num_sets = ii.numsets;

  generate_set_perms(num_sets);

  for (int i = 0; i < num_sets; i++)
    set_pt[i] = g_pos.pt[ii.first[i]];

  generate_test_list(kslice_size, g_pos.num - 2);

  work_convert = create_work(g_total_work, kslice_size, 0);
}

#if 0
void init_permute_pawn(uint8_t *pcs, uint8_t *pt)
{
  int i, j, j0, j1, k, m, l;
  int tidx[16];
  int pivtype;
  int numpcs = g_pos.num;

  for (i = 0; i < numpcs; i++)
    pw[i] = (pt[i] == WPAWN) ? 0x38 : 0x00;
  uint64_t pw_mask = 0;
  for (i = 1; i < numpcs; i++)
    pw_mask |= (uint64_t)pw[i] << (6 * (numpcs - i - 1));

  if (pw[0])
    for (i = 1; i < 4; i++)
      Swap(flip0x40[i], flip0x40[7 - i]);

  for (i = 0; i < 64; i++) {
    sq_mask_pawn[i] = (i & 0x04) ? mask_a1h1 : 0ULL;
//    sq_mask_pawn[i] |= flip0x40[i >> 3];
    sq_mask_pawn[i] |= (uint64_t)(i & 0x03) << (6 * (numpcs - 1));
    sq_mask_pawn[i] ^= pw_mask;
  }

  for (i = 0, k = 0; i < 16; i++)
    if (pcs[i]) type[k++] = i;
  num_types = k;

  pivtype = pt[0];
  base_entry.pawns[0] = pcs[pivtype == WPAWN ? WPAWN : BPAWN];
  base_entry.pawns[1] = pcs[pivtype == WPAWN ? BPAWN : WPAWN];

  base_entry.num = 0;
  for (i = 0; i < 16; i++)
    base_entry.num += pcs[i];

  generate_type_perms(num_types);

  for (i = 0; i < num_types; i++) {
    for (j = 0; j < base_entry.num; j++)
      if (pt[j] == type[i]) break;
    tidx[type[i]] = j;
  }

  for (i = 0; i < num_type_perms; i++) {
    for (j0 = 0; j0 < num_types; j0++)
      if (type_perm_list[i][j0] == pivtype) break;
    for (j1 = 0; j1 < num_types; j1++)
      if (type_perm_list[i][j1] == (pivtype ^ 0x08)) break;
    if (j1 == num_types) j1 = 0x0f;
    order_list[i] = j0;
    order2_list[i] = j1;
    for (m = 0; m < base_entry.pawns[0]; m++)
      piece_perm_list[i][m] = tidx[pivtype] + m;
    for (; m < base_entry.pawns[0] + base_entry.pawns[1]; m++)
      piece_perm_list[i][m] = tidx[pivtype ^ 0x08] + (m - base_entry.pawns[0]);
    for (k = 0; k < num_types; k++)
      if (k != j0 && k != j1)
	for (l = 0; l < pcs[type_perm_list[i][k]]; l++)
	  piece_perm_list[i][m++] = tidx[type_perm_list[i][k]] + l;
  }

  for (i = 0; i < 16; i++)
    cmp[i] = pcs[i];
  cmp[WPAWN] = 13;
  cmp[BPAWN] = 14;

  for (i = 0; i < num_type_perms; i++)
    for (j = i + 1; j < num_type_perms; j++) {
      for (k = 0; k < num_types; k++)
	if (cmp[type_perm_list[i][k]] != cmp[type_perm_list[j][k]]) break;
      if (   k < num_types
          && cmp[type_perm_list[i][k]] > cmp[type_perm_list[j][k]])
      {
	for (k = 0; k < num_types; k++)
          Swap(type_perm_list[i][k], type_perm_list[j][k]);
	for (k = 0; k < base_entry.num; k++)
          Swap(piece_perm_list[i][k], piece_perm_list[j][k]);
        Swap(order_list[i], order_list[j]);
        Swap(order2_list[i], order2_list[j]);
      }
    }

  for (i = 0; i < num_type_perms; i++)
    for (j = 0; j < base_entry.num; j++)
      pidx_list[i][piece_perm_list[i][j]] = j;

  work_convert = alloc_work(g_total_work);
}

void *init_permute_rank(uint8_t *pcs, int rank, void *tb_table, bool wide)
{
  tb_size = set_dec_info(&dec_info, &base_entry, pcs, type_perm_list[0],
      order_list[0], order2_list[0], rank, RANK_ENC);
  printf("tb_size = %"PRIu64"\n", tb_size);

  generate_test_list(tb_size, base_entry.num);

  if (!tb_table && !(tb_table = malloc((tb_size + 1) * (wide ? 2 : 1))))
    out_of_mem();

  fill_work(g_total_work, tb_size, 0, work_convert);

  return tb_table;
}
#endif

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

#if 0
void estimate_compression_pawn(void *table, uint8_t *pcs, int rank,
    int num_cands, bool wide)
{
  uint64_t dsize = num_segs * seg_size;
  void *dst = malloc((num_cands * dsize + 1) * (wide ? 2 : 1));
  est_data.table = table;
  est_data.pcs = pcs;
  est_data.dst = dst;
  est_data.num_cands = num_cands;
  est_data.dsize = dsize;
  est_data.rank = rank;
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
    int64_t *freq = construct_pairs(dst0, dsize, 20, 100, wide, &num_syms);
    struct HuffCode *c = create_code(freq, num_syms);
    csize = calc_size(c);
    free_code(c);
    printf("[%2d] order: %d", p, order_list[trylist[p]]);
    printf("; perm:");
    for (int i = 0; i < num_types; i++)
      printf(" %2d", type_perm_list[trylist[p]][i]);
    printf("; %"PRIu64"\n", csize);
    compest[trylist[p]] = csize;
    dst0 += (wide ? 2 : 1) * dsize;
  }

  free(dst);
}
#endif

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
    best_ii.last[i] = ii.last[k];
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

#if 0
void permute_pawn_dtm(void *tb_table, void *table, uint8_t *best, int rank,
    void *v, bool wide)
{
  int i;

//  permute_v = v;

  int bestp;

  estimate_compression(table, &bestp, pcs, rank, wide);

  for (i = 0; i < base_entry.num; i++)
    best[i] = pt[piece_perm_list[bestp][i]];
  best[0] |= order_list[bestp] << 4;
  best[1] |= order2_list[bestp] << 4;

  printf("best order: %d", best[0] >> 4);
  if ((best[1] >> 4) < 0x0f)
    printf(" %d", best[1] >> 4);
  printf("\nbest permutation:");
  for (i = 0; i < base_entry.num; i++)
    printf(" %d", best[i] & 0x0f);
  printf("\n");

  if (g_compress_type == 1) return;

  convert_data.src = table;
  convert_data.dst = tb_table;
  convert_data.pcs = pcs;
  convert_data.p = bestp;
  convert_data.rank = rank;

  if (!wide)
    run_threaded(convert_data_pawn_u8, work_convert, 1);
  else
    run_threaded(convert_data_pawn_u16, work_convert, 1);
}
#endif
