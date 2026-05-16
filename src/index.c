/*
  Copyright (c) 2026 Ronald de Man

  This file is distributed under the terms of the GNU GPL, version 2.
*/

#include <stdint.h>

#include "defs.h"
#include "index.h"
#include "probe.h"
#include "types.h"

struct IdxInfo ii, capt_ii[MAX_SETS];
int pc_to_set[MAX_PIECES];
Bitboard Unrank2[62 * 61 / 2], Unrank3[62 * 61 * 60 / 6];

#ifndef HAS_PAWNS
static constexpr bool has_pawns = false;
#else
static constexpr bool has_pawns = true;
#endif

void init_unrank(void)
{
  int idx = 0;
  for (int s1 = 1; s1 < 62; s1++)
    for (int s0 = 0; s0 < s1; s0++)
      Unrank2[idx++] = bit(s0) | bit(s1);

  idx = 0;
  for (int s2 = 2; s2 < 62; s2++)
    for (int s1 = 1; s1 < s2; s1++)
      for (int s0 = 0; s0 < s1; s0++)
        Unrank3[idx++] = bit(s0) | bit(s1) | bit(s2);
}

// We expect a normalized position.
INLINE uint64_t sq_to_idx_helper(uint8_t *restrict sq, const struct IdxInfo *ii)
{
  uint64_t idx = 0;
  Bitboard occ = has_pawns ? bit(sq[0]) | bit(sq[1]) | bit(sq[2])
                           : bit(sq[0]) | bit(sq[1]);

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
    idx = idx * ii->factor[k] + s;
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

INLINE Bitboard idx_to_sq_unpack(uint32_t *sub, uint8_t *restrict sq,
    const struct IdxInfo *ii)
{
  Bitboard occ = has_pawns ? bit(sq[0]) | bit(sq[1]) | bit(sq[2])
                           : bit(sq[0]) | bit(sq[1]);
  for (int i = 0; i < ii->numsets; i++)
    occ = unrank_binomial(sub[i], ii->mult[i], sq + ii->first[i], occ);
  return occ;
}

void idx_to_sq_init(uint64_t idx, uint32_t *restrict sub,
    const struct IdxInfo *ii)
{
  for (int k = ii->numsets - 1; k >= 0; k--)
    idx = divmod_u64_u32_recip(idx, ii->factor[k], ii->recip[k], &sub[k]);
}

Bitboard idx_to_sq(uint32_t *sub, uint8_t *restrict sq)
{
  return idx_to_sq_unpack(sub, sq, &ii);
}

void idx_to_sq_ii(uint32_t *sub, uint8_t *restrict sq, const struct IdxInfo *ii)
{
  idx_to_sq_unpack(sub, sq, ii);
}

Bitboard capt_idx_to_sq(uint32_t *sub, uint8_t *restrict sq, const int k)
{
  return idx_to_sq_unpack(sub, sq, &capt_ii[k]);
}

void calc_factors(struct IdxInfo *ii)
{
  for (int i = 0, n = has_pawns ? 61 : 62; i < ii->numsets; i++) {
    ii->factor[i] = Binomial[ii->mult[i]][n];
    ii->recip[i] = (((__uint128_t)1 << 64) + ii->factor[i] - 1) / ii->factor[i];
    n -= ii->mult[i];
  }

  uint64_t f = 1;
  for (int i = ii->numsets - 1; i >= 0; i--)
    f *= ii->factor[i];
  ii->size = f;
}
