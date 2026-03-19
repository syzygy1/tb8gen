/*
  Copyright (c) 2026 Ronald de Man

  This file is distributed under the terms of the GNU GPL, version 2.
*/

#ifndef INDEX_H
#define INDEX_H

#include <inttypes.h>
#include <x86intrin.h>

#include "defs.h"
#include "probe.h"
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

extern int16_t KKMap[64][64];
extern uint8_t MirrorMask[64];
extern bool FlipTest[64][64];

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
  for (int i = ii->numsets - 1; ++sub[i] >= ii->factor[i]; i--)
    sub[i] = 0;
}

// FIXME: make sure that v and sub[] never overflow
// probably just insert a check: if v too big, then do as in init().
INLINE void idx_to_sq_add(uint32_t v, uint32_t *restrict sub,
    const struct IdxInfo *ii)
{
  int i = ii->numsets;
  while (v) {
    sub[--i] += v;
    v = 0;
    while (sub[i] >= ii->factor[i]) {
      sub[i] -= ii->factor[i];
      v++;
    }
  }
}

// Mirror wK to A1-D1-D4 and, if wK on A1-D4, then bK to A1-H1-H8.
INLINE void normalize(const uint8_t *restrict sq, uint8_t *restrict sq2)
{
  for (int i = 0; i < MAX_PIECES; i++)
    sq2[i] = sq[i] ^ MirrorMask[sq[0]];

  if (FlipTest[sq[0]][sq[1]])
    for (int i = 0; i < MAX_PIECES; i++)
      sq2[i] = FlipDiag[sq2[i]];
}

INLINE void mirror_diagonal(uint8_t *restrict sq)
{
  for (int i = 2; i < MAX_PIECES; i++)
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
