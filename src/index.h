/*
  Copyright (c) 2026 Ronald de Man

  This file is distributed under the terms of the GNU GPL, version 2.
*/

#ifndef INDEX_H
#define INDEX_H

#include <assert.h>
#include <inttypes.h>
#include <x86intrin.h>

#include "defs.h"
#include "probe.h"
#include "tb8gen.h"
#include "types.h"

struct IdxInfo {
  int numsets;   // number of sets of like pieces, excluding kings.
  uint64_t size;
  uint32_t factor[MAX_SETS]; // total number of placements for a set
  int first[MAX_SETS];          // index of first piece of each set
  int mult[MAX_SETS];           // number of like pieces in each set
  int last[MAX_SETS];
};

extern struct IdxInfo ii, capt_ii[MAX_SETS];

extern int pc_to_set[MAX_PIECES];

INLINE void sort_squares(int n, uint8_t *restrict sq)
{
  for (int i = 0; i < n; i++)
    for (int j = i + 1; j < n; j++)
      if (sq[i] > sq[j])
        Swap(sq[i], sq[j]);
}

INLINE int rank_among_free(uint8_t sq, Bitboard occ)
{
  return sq - popcnt(occ & ((1ULL << sq) - 1));
}

INLINE Bitboard unrank_binomial(uint64_t idx, int n, uint8_t *restrict sq,
    Bitboard occ)
{
  if (n == 0)
    return occ;

  Bitboard b = ~occ;
  for (int i = n - 1; i > 0; i--) {
    int r = i;
    while (idx >= Binomial[i + 1][r + 1])
      r++;
    idx -= Binomial[i + 1][r];
    Bitboard b1 = _pdep_u64(1ULL << r, b);
    sq[i] = lsb(b1);
    occ |= b1;
  }
  Bitboard b1 = _pdep_u64(1ULL << idx, b);
  sq[0] = lsb(b1);
  occ |= b1;

  return occ;
}

INLINE void idx_to_sq_inc(uint32_t *sub, const struct IdxInfo *ii)
{
  for (int i = ii->numsets - 1; ++sub[i] >= ii->factor[i] && i > 0; i--)
    sub[i] = 0;
}

INLINE void idx_to_sq_add(uint64_t v, uint32_t *restrict sub,
    const struct IdxInfo *ii)
{
  for (int i = ii->numsets - 1; v && i >= 0; i--) {
    uint64_t sum = sub[i] + v;
    uint32_t factor = ii->factor[i];
    if (sum < factor) {
      sub[i] = sum;
      break;
    }
    if (sum < 2 * (uint64_t)factor) {
      sub[i] = sum - factor;
      v = 1;
    } else {
      sub[i] = sum % factor;
      v = sum / factor;
    }
  }
}

// Mirror wK to A1-D1-D4 and, if wK on A1-D4, then bK to A1-H1-H8.
INLINE void normalize(const uint8_t *restrict sq, uint8_t *restrict sq2)
{
  assume(g_pos.num <= MAX_PIECES);

  for (int i = 0; i < g_pos.num; i++)
    sq2[i] = sq[i] ^ MirrorMask[sq[0]];

  if (FlipTest[sq[0]][sq[1]])
    for (int i = 0; i < g_pos.num; i++)
      sq2[i] = FlipDiag[sq2[i]];
}

INLINE void mirror_diagonal(uint8_t *restrict sq)
{
  assume(g_pos.num <= MAX_PIECES);

  for (int i = 2; i < g_pos.num; i++)
    sq[i] = FlipDiag[sq[i]];
}

void init_tables(void);
void calc_factors(struct IdxInfo *ii);
uint64_t sq_to_idx(uint8_t *sq);
uint64_t capt_sq_to_idx(uint8_t *sq, int k);
void idx_to_sq_init(uint64_t idx, uint32_t *sub, const struct IdxInfo *ii);
Bitboard idx_to_sq(uint32_t *sub, uint8_t *sq);
void idx_to_sq_ii(uint32_t *sub, uint8_t *sq, const struct IdxInfo *ii);
Bitboard capt_idx_to_sq(uint32_t *sub, uint8_t *sq, int k);

#endif
