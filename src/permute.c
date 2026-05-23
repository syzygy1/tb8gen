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
#include "permute.h"
#include "probe.h"
#include "tb8gen.h"
#include "threads.h"
#include "util.h"

uint64_t tb_size;

//extern uint64_t sq_mask[64];

static constexpr int MAX_PERMS = 1*2*3*4*5*6*7;
static constexpr int MAX_CANDS = 6*7;

static struct Work *work_convert = nullptr;
static struct Work *work_est = nullptr;

static struct TbEntry tb_entry;
static struct DecInfo dec_info;

static uint8_t order_list[MAX_PERMS];
static uint8_t order2_list[MAX_PERMS];

static uint8_t piece_perm_list[MAX_PERMS][MAX_PIECES];
static uint8_t pidx_list[MAX_PERMS][MAX_PIECES];

static int num_types, num_type_perms;
static uint8_t type[MAX_PIECES];
static uint8_t type_perm_list[MAX_PERMS][MAX_PIECES];

static uint64_t compest[MAX_PERMS];

static int trylist[MAX_CANDS];

#if 0
static int pw[MAX_PIECES];
static int cmp[16];
#endif

static uint8_t perm_tmp[MAX_PIECES];

static void generate_type_perms2(int n, int k)
{
  if (k == n) {
    for (int i = 0; i < n; i++)
      type_perm_list[num_type_perms][i] = perm_tmp[i];
    num_type_perms++;
    return;
  }

  for (int i = 0; i < n; i++)
    if (perm_tmp[i] == 0xff) {
      perm_tmp[i] = type[k];
      generate_type_perms2(n, k + 1);
      perm_tmp[i] = 0xff;
    }
}

static void generate_type_perms(int n)
{
  num_type_perms = 0;
  for (int i = 0; i < n; i++)
    perm_tmp[i] = 0xff;
  generate_type_perms2(n, 0);
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

INLINE uint64_t calc_idx_piece(uint8_t *restrict sq, uint8_t *restrict pidx,
    int n)
{
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

void init_permute_piece(uint8_t *pcs, uint8_t *pt)
{
  int i, j, k, m, l;
  int tidx[16];

  k = 0;
  for (int i = 0; i < 16; i++)
    if (pcs[i]) type[k++] = i;
  num_types = k;

  k = 0;
  for (int i = 0; i < 16; i++)
    if (pcs[i] == 1) k++;
  tb_entry.kk_enc = (k == 2);
#if 0
  if (k >= 3) tb_entry.kk_enc = false;
  else if (k == 2) tb_entry.kk_enc = true;
  else { /* only possible for suicide chess */
    k = 16;
    for (int i = 0; i < 16; i++)
      if (pcs[i] < k && pcs[i] > 1)
        k = pcs[i];
    tb_entry.kk_enc = 1 + k;
  }
#endif

  tb_entry.num = 0;
  for (int i = 0; i < 16; i++)
    tb_entry.num += pcs[i];

  generate_type_perms(num_types);

  for (i = 0; i < num_types; i++) {
    for (j = 0; j < tb_entry.num; j++)
      if (pt[j] == type[i]) break;
    tidx[type[i]] = j;
  }

  if (!tb_entry.kk_enc) { /* 111 */
    for (i = 0; i < num_type_perms;) {
      for (j = num_types - 3; j >= 0; j--)
        if (   pcs[type_perm_list[i][j    ]] == 1
            && pcs[type_perm_list[i][j + 1]] == 1
            && pcs[type_perm_list[i][j + 2]] == 1) break;
      if (j < 0) {
        num_type_perms--;
        for (k = 0; k < num_types; k++)
          type_perm_list[i][k] = type_perm_list[num_type_perms][k];
      } else {
        piece_perm_list[i][2] = tidx[type_perm_list[i][j]];
        piece_perm_list[i][1] = tidx[type_perm_list[i][j + 1]];
        piece_perm_list[i][0] = tidx[type_perm_list[i][j + 2]];
        for (k = 0, m = 3; k < num_types;)
          if (k != j) {
            for (l = 0; l < pcs[type_perm_list[i][k]]; l++)
              piece_perm_list[i][m++] = tidx[type_perm_list[i][k]] + l;
            k++;
          } else
            k += 3;
        order_list[i] = j;
        i++;
      }
    }
  }
  else {
    for (i = 0; i < num_type_perms;) {
      for (j = num_types - 2; j >= 0; j--)
        if (   pcs[type_perm_list[i][j    ]] == 1
            && pcs[type_perm_list[i][j + 1]] == 1) break;
      if (j < 0) {
        num_type_perms--;
        for (k = 0; k < num_types; k++)
          type_perm_list[i][k] = type_perm_list[num_type_perms][k];
      } else {
        piece_perm_list[i][1] = tidx[type_perm_list[i][j]];
        piece_perm_list[i][0] = tidx[type_perm_list[i][j + 1]];
        for (k = 0, m = 2; k < num_types;)
          if (k != j) {
            for (l = 0; l < pcs[type_perm_list[i][k]]; l++)
              piece_perm_list[i][m++] = tidx[type_perm_list[i][k]] + l;
            k++;
          } else
            k += 2;
        order_list[i] = j;
        i++;
      }
    }
  }
#if 0
 else { /* 2, or 3, or 4, or higher; only possible for suicide chess */
    int p = tb_entry.kk_enc - 1;
    for (i = 0; i < num_type_perms;) {
      for (j = num_types - 1; j >= 0; j--)
        if (pcs[type_perm_list[i][j]] == p) break;
      for (k = 0; k < p; k++)
        piece_perm_list[i][k] = tidx[type_perm_list[i][j]] + k;
      for (k = 0, m = p; k < num_types; k++)
        if (k != j) {
          for (l = 0; l < pcs[type_perm_list[i][k]]; l++)
            piece_perm_list[i][m++] = tidx[type_perm_list[i][k]] + l;
        }
      order_list[i] = j;
      i++;
    }
  }
#endif

  for (i = 0; i < num_type_perms; i++)
    for (j = i + 1; j < num_type_perms; j++) {
      for (k = 0; k < num_types; k++)
        if (pcs[type_perm_list[i][k]] != pcs[type_perm_list[j][k]]) break;
      if (   k < num_types
          && pcs[type_perm_list[i][k]] > pcs[type_perm_list[j][k]])
      {
        for (k = 0; k < num_types; k++)
          Swap(type_perm_list[i][k], type_perm_list[j][k]);
        for (k = 0; k < tb_entry.num; k++)
          Swap(piece_perm_list[i][k], piece_perm_list[j][k]);
        Swap(order_list[i], order_list[j]);
      }
    }

  for (i = 0; i < num_type_perms; i++)
    order2_list[i] = 0xf;

  for (i = 0; i < num_type_perms; i++)
    for (j = 0; j < tb_entry.num; j++)
      pidx_list[i][piece_perm_list[i][j]] = j;

  tb_size = set_dec_info(&dec_info, &tb_entry, pcs, type_perm_list[0],
      order_list[0], 0xf, -1, LT_PIECE);
  printf("tb_size = %"PRIu64"\n", tb_size);

  generate_test_list(tb_size, tb_entry.num);

  work_convert = create_work(g_total_work, tb_size, 0);
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
  tb_entry.pawns[0] = pcs[pivtype == WPAWN ? WPAWN : BPAWN];
  tb_entry.pawns[1] = pcs[pivtype == WPAWN ? BPAWN : WPAWN];

  tb_entry.num = 0;
  for (i = 0; i < 16; i++)
    tb_entry.num += pcs[i];

  generate_type_perms(num_types);

  for (i = 0; i < num_types; i++) {
    for (j = 0; j < tb_entry.num; j++)
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
    for (m = 0; m < tb_entry.pawns[0]; m++)
      piece_perm_list[i][m] = tidx[pivtype] + m;
    for (; m < tb_entry.pawns[0] + tb_entry.pawns[1]; m++)
      piece_perm_list[i][m] = tidx[pivtype ^ 0x08] + (m - tb_entry.pawns[0]);
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
	for (k = 0; k < tb_entry.num; k++)
          Swap(piece_perm_list[i][k], piece_perm_list[j][k]);
        Swap(order_list[i], order_list[j]);
        Swap(order2_list[i], order2_list[j]);
      }
    }

  for (i = 0; i < num_type_perms; i++)
    for (j = 0; j < tb_entry.num; j++)
      pidx_list[i][piece_perm_list[i][j]] = j;

  work_convert = alloc_work(g_total_work);
}

void *init_permute_rank(uint8_t *pcs, int rank, void *tb_table, bool wide)
{
  tb_size = set_dec_info(&dec_info, &tb_entry, pcs, type_perm_list[0],
      order_list[0], order2_list[0], rank, RANK_ENC);
  printf("tb_size = %"PRIu64"\n", tb_size);

  generate_test_list(tb_size, tb_entry.num);

  if (!tb_table && !(tb_table = malloc((tb_size + 1) * (wide ? 2 : 1))))
    out_of_mem();

  fill_work(g_total_work, tb_size, 0, work_convert);

  return tb_table;
}
#endif

static struct {
  void *src;
  void *dst;
  uint8_t *pcs;
  int p;
  int rank;
} convert_data;

static struct {
  void *table;
  uint8_t *pcs;
  void *dst;
  int num_cands;
  uint32_t dsize;
  int rank;
} est_data;

#define T u8
#include "permute_tmpl.c"
#undef T

#define T u16
#include "permute_tmpl.c"
#undef T

static void estimate_compression_piece(void *table, uint8_t *pcs, int num_cands,
    bool wide, bool wdl)
{
  uint32_t dsize = num_segs * seg_size;
  void *dst = malloc((num_cands * dsize + 1) * (wide ? 2 : 1));
  est_data.table = table;
  est_data.pcs = pcs;
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

#if 0
static void estimate_compression_pawn(void *table, uint8_t *pcs, int rank,
    int num_cands, bool wide)
{
  uint32_t dsize = num_segs * seg_size;
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

static uint64_t estimate_compression(void *table, int *bestp, uint8_t *pcs,
    int rank, bool wide, bool wdl)
{
  int i, j, k, p, q;
  int num_cands, bp = 0;
  uint64_t best = 0;
  uint8_t bestperm[MAX_PIECES];

  if (g_compress_type == 1) {
    *bestp = 0;
    return 0;
  }

  if (!work_est)
    work_est = create_work(g_total_work, num_segs, 0);
  else
    work_refill_units(work_est, g_total_work, num_segs, 0);

  for (i = 0; i < num_type_perms; i++)
    compest[i] = 0;

  for (k = 0; k < num_types - 1; k++) {
    best = UINT64_MAX;
    num_cands = 0;
    for (p = 0; p < num_types; p++) {
      for (i = 0; i < k; i++)
	if (type[p] == bestperm[i]) break;
      if (i < k) continue;
      for (q = 0; q < num_types; q++) {
	if (q == p) continue;
	for (i = 0; i < k; i++)
	  if (type[q] == bestperm[i]) break;
	if (i < k) continue;
	// look for permutation starting with bestperm[0..k-1],p,q
	for (i = 0; i < num_type_perms; i++) {
	  for (j = 0; j < k; j++)
	    if (type_perm_list[i][j] != bestperm[j]) break;
	  if (j < k) continue;
	  if (   type_perm_list[i][k]   == type[p]
              && type_perm_list[i][k+1] == type[q]) break;
	}
	if (i < num_type_perms) {
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
	if (trylist[i] > trylist[j]) {
	  int tmp = trylist[i];
	  trylist[i] = trylist[j];
	  trylist[j] = tmp;
	}
    if (rank < 0)
      estimate_compression_piece(table, pcs, num_cands, wide, wdl);
//    else
//      estimate_compression_pawn(table, pcs, rank, num_cands, wide);
    for (i = 0; i < num_cands; i++) {
      if (compest[trylist[i]] < best) {
	best = compest[trylist[i]];
	bp = trylist[i];
      }
    }
    bestperm[k] = type_perm_list[bp][k];
  }
  *bestp = bp;

  return best;
}

void permute_piece_wdl(void *tb_table, uint8_t *pcs, uint8_t *pt,
    void *table, uint8_t *best)
{
  int bestp;

  estimate_compression(table, &bestp, pcs, -1, false, true);

  for (int i = 0; i < tb_entry.num; i++)
    best[i] = pt[piece_perm_list[bestp][i]];
  best[0] |= order_list[bestp] << 4;

  printf("best order: %d", best[0] >> 4);
  printf("\nbest permutation:");
  for (int i = 0 ;i < tb_entry.num; i++)
    printf(" %d", best[i] & 0x0f);
  printf("\n");

  convert_data.src = table;
  convert_data.dst = tb_table;
  convert_data.pcs = pcs;
  convert_data.p = bestp;

  if (g_compress_type == 1)
    return;

  run_threaded(convert_data_piece_u8, work_convert, 1);
}

void permute_piece_dtz(void *tb_table, uint8_t *pcs, uint8_t *pt, void *table,
    uint8_t *best, bool wide)
{
  int bestp;
  estimate_compression(table, &bestp, pcs, -1, wide, false);

  for (int i = 0; i < tb_entry.num; i++)
    best[i] = pt[piece_perm_list[bestp][i]];
  best[0] |= order_list[bestp] << 4;

  printf("best order: %d", best[0] >> 4);
  printf("\nbest permutation:");
  for (int i = 0; i < tb_entry.num; i++)
    printf(" %d", best[i] & 0x0f);
  printf("\n");

  convert_data.src = table;
  convert_data.dst = tb_table;
  convert_data.pcs = pcs;
  convert_data.p = bestp;

  if (g_compress_type == 1)
    return;

  if (!wide)
    run_threaded(convert_data_piece_u8, work_convert, 1);
  else
    run_threaded(convert_data_piece_u16, work_convert, 1);
}

#if 0
void permute_pawn_dtm(void *tb_table, uint8_t *pcs, uint8_t *pt, void *table,
    uint8_t *best, int rank, void *v, bool wide)
{
  int i;

//  permute_v = v;

  int bestp;

  estimate_compression(table, &bestp, pcs, rank, wide);

  for (i = 0; i < tb_entry.num; i++)
    best[i] = pt[piece_perm_list[bestp][i]];
  best[0] |= order_list[bestp] << 4;
  best[1] |= order2_list[bestp] << 4;

  printf("best order: %d", best[0] >> 4);
  if ((best[1] >> 4) < 0x0f)
    printf(" %d", best[1] >> 4);
  printf("\nbest permutation:");
  for (i = 0; i < tb_entry.num; i++)
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
