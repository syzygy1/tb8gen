/*
  Copyright (c) 2026 Ronald de Man

  This file is distributed under the terms of the GNU GPL, version 2.
*/

#include <stdint.h>
#include <x86intrin.h>

#include "defs.h"
#include "index.h"
#include "probe.h"

int16_t KKMap[64][64];
uint8_t MirrorMask[64];
bool FlipTest[64][64];

struct IdxInfo ii, capt_ii[MAX_SETS];
int pc_to_set[MAX_PIECES];

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

// We expect a normalized position.
INLINE uint64_t sq_to_idx_helper(uint8_t *restrict sq, struct IdxInfo *ii)
{
  uint64_t idx = 0;
  Bitboard occ = bit(sq[0]) | bit(sq[1]);

  for (int k = 0; k < ii->numsets; k++) {
    int i = ii->first[k];
    sort_squares(ii->mult[k], &sq[i]);
    size_t s = 0;
    Bitboard occ2 = occ;
    for (int j = 0; j < ii->mult[k]; i++, j++) {
      int rank = rank_among_free(sq[i], occ);
      occ2 |= bit(sq[i]);
      s += Binomial[j + 1][rank];
    }
    idx += s * ii->factor[k];
    occ = occ2;
  }

  return idx;
}

uint64_t sq_to_idx(uint8_t *restrict sq)
{
  return sq_to_idx_helper(sq, &ii);
}

uint64_t capt_sq_to_idx(uint8_t *restrict sq, int k)
{
  return sq_to_idx_helper(sq, &capt_ii[k]);
}

INLINE Bitboard idx_to_sq_unpack(uint32_t *sub, uint8_t *sq, struct IdxInfo *ii)
{
  Bitboard occ = bit(sq[0]) | bit(sq[1]);
  for (int i = 0; i < ii->numsets; i++)
    occ = unrank_binomial(sub[i], ii->mult[i], sq + ii->first[i], occ);
  return occ;
}

void idx_to_sq_init(uint64_t idx, uint32_t *sub, struct IdxInfo *ii)
{
  for (int k = 0; k < ii->numsets; k++) {
    sub[k] = idx / ii->factor[k];
    idx -= sub[k] * ii->factor[k];
  }
}

Bitboard idx_to_sq(uint32_t *sub, uint8_t *restrict sq)
{
  return idx_to_sq_unpack(sub, sq, &ii);
}

Bitboard capt_idx_to_sq(uint32_t *sub, uint8_t *restrict sq, const int k)
{
  return idx_to_sq_unpack(sub, sq, &capt_ii[k]);
}

void calc_factors(struct IdxInfo *ii)
{
  for (int i = 0; i < MAX_SETS; i++)
    ii->factor[i] = 0;

  for (int i = 0, n = 62; i < ii->numsets; i++) {
    ii->subfactor[i] = subfactor(ii->mult[i], n);
    n -= ii->mult[i];
  }

  uint64_t f = 1;
  for (int i = ii->numsets - 1; i >= 0; i--) {
    ii->factor[i] = f;
    f *= ii->subfactor[i];
  }
  ii->size = f;
  // Increase subfactor[0] to ensure we can go slightly beyond the end
  // without hanging in sq_to_idx_add().
  ii->subfactor[0] += 64;
}

void init_tables(void)
{
  static constexpr Bitboard A1D1D4 = 0x080c0e0full;
  static constexpr Bitboard A1D4   = 0x08040201ull;
  static constexpr Bitboard LOWER  = 0x80c0e0f0f8fcfeffull;

  for (int s = 0; s < 64; s++)
    MirrorMask[s] = ((s & 0x04) ? 0x07 : 0x00) | ((s & 0x20) ? 0x38 : 0x00);

  for (int i = 0; i < 64; i++)
    for (int j = 0; j < 64; j++) {
      int s1 = i ^ MirrorMask[i];
      int s2 = j ^ MirrorMask[i];
      if (!(bit(s1) & A1D1D4) || ((bit(s1) & A1D4) && !(bit(s2) & LOWER))) {
        FlipTest[i][j] = true;
        s1 = FlipDiag[s1];
        s2 = FlipDiag[s2];
      }
      KKMap[i][j] = KKIdx[Triangle[s1]][s2];
    }
}
