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
#ifndef HAS_PAWNS
#include "tb8gen.h"
#else
#include "tb8genp.h"
#endif
#include "types.h"

#ifndef HAS_PAWNS
static constexpr bool has_pawns = false;
#else
static constexpr bool has_pawns = true;
#endif

struct IdxInfo {
  int numsets;   // number of sets of like pieces, excluding kings.
  int first[MAX_SETS];          // index of first piece of each set
  int mult[MAX_SETS];           // number of like pieces in each set
  uint32_t factor[MAX_SETS]; // total number of placements for a set
  uint64_t recip[MAX_SETS];
  uint64_t size;
  int last[MAX_SETS];
};

struct IdxState {
  uint32_t sub[MAX_SETS];
  Bitboard occ[MAX_SETS];
  int n;
};

extern struct IdxInfo ii, capt_ii[MAX_SETS];
extern int pc_to_set[MAX_PIECES];
extern Bitboard Unrank2[62 * 61 / 2], Unrank3[62 * 61 * 60 / 6];

#define SORT2(a, b) do { \
  if ((b) < (a))         \
    Swap(a, b);          \
} while (0)

INLINE void sort3(uint8_t *x)
{
  SORT2(x[0], x[1]);
  SORT2(x[1], x[2]);
  SORT2(x[0], x[1]);
}

INLINE void sort4(uint8_t *x)
{
  SORT2(x[0], x[1]);
  SORT2(x[2], x[3]);
  SORT2(x[0], x[2]);
  SORT2(x[1], x[3]);
  SORT2(x[1], x[2]);
}

INLINE void sort_squares(int n, uint8_t *restrict x)
{
  assume(n <= MAX_PIECES - 2);
  switch (n) {
  case 0:
  case 1:
    break;

  case 2:
    SORT2(x[0], x[1]);
    break;

  case 3:
    sort3(x);
    break;

  case 4:
    sort4(x);
    break;

  default:
    // insertion sort
    for (int i = 1; i < n; i++) {
      int v = x[i];
      int j = i;
      while (j > 0 && v < x[j - 1]) {
        x[j] = x[j - 1];
        j--;
      }
      x[j] = v;
    }
  }
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

  if (n == 1) {
    Bitboard b1 = _pdep_u64(1ULL << idx, b);
    occ |= b1;
    sq[0] = lsb(b1);
  }
  else if (n == 2) {
    Bitboard b1 = _pdep_u64(Unrank2[idx], b);
    occ |= b1;
    sq[0] = pop_lsb(&b1);
    sq[1] = lsb(b1);
  }
  else if (n == 3) {
    Bitboard b1 = _pdep_u64(Unrank3[idx], b);
    occ |= b1;
    sq[0] = pop_lsb(&b1);
    sq[1] = pop_lsb(&b1);
    sq[2] = lsb(b1);
  }
  else {
    Bitboard b1 = 0;
    for (int i = n - 1; i > 0; i--) {
      int r = i;
      while (idx >= Binomial[i + 1][r + 1])
        r++;
      idx -= Binomial[i + 1][r];
      b1 |= bit(r);
    }
    b1 = _pdep_u64(b1 | bit(idx), b);
    occ |= b1;
    while (b1)
      *sq++ = pop_lsb(&b1);
  }

  return occ;
}

INLINE void idx_to_sq_inc(uint32_t *sub, const struct IdxInfo *ii)
{
  for (int i = ii->numsets - 1; ++sub[i] >= ii->factor[i] && i > 0; i--)
    sub[i] = 0;
}

// Valid if x <= 2^N, d-1 <= 2^l and N + l <= 64.
// This should not be a problem even for 9-piece tables.
// See https://gmplib.org/~tege/divcnst-pldi94.pdf
INLINE uint64_t divmod_recip(uint64_t x, uint32_t d, uint64_t recip,
    uint32_t *rem)
{
  uint64_t q = ((__uint128_t)x * recip) >> 64;
  uint64_t r = x - q * d;

  *rem = (uint32_t)r;
  return q;
}

INLINE uint64_t recip(uint64_t f)
{
  assert(f != 1);
  return (((__uint128_t)1 << 64) + f - 1) / f;
}

INLINE void idx_to_sq_add(uint64_t v, uint32_t *restrict sub,
    const struct IdxInfo *restrict ii)
{
  int i = ii->numsets;

  while (i > 0) {
    uint64_t s = (uint64_t)sub[--i] + v;
    uint32_t f = ii->factor[i];

    if (s < f) {
      sub[i] = s;
      return;
    }

    v = divmod_u64_u32_recip(s, f, ii->recip[i], &sub[i]);
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

INLINE void mirror_diagonal2(uint8_t *restrict sq, uint8_t *restrict sq2)
{
  assume(g_pos.num <= MAX_PIECES);

  sq2[0] = sq[0];
  sq2[1] = sq[1];
  for (int i = 2; i < g_pos.num; i++)
    sq2[i] = FlipDiag[sq[i]];
}

INLINE void idx_state_inc(struct IdxState *is, const struct IdxInfo *ii)
{
  uint32_t *restrict sub = is->sub;
  int i = ii->numsets - 1;
  for (; ++sub[i] >= ii->factor[i] && i > 0; i--)
    sub[i] = 0;
  is->n = i;
}

INLINE void idx_state_add(struct IdxState *is, uint64_t v,
    const struct IdxInfo *restrict ii)
{
  uint32_t *restrict sub = is->sub;
  int i = ii->numsets;

  while (i > 0) {
    uint64_t s = (uint64_t)sub[--i] + v;
    uint32_t f = ii->factor[i];

    if (s < f) {
      sub[i] = s;
      is->n = i;
      return;
    }

    v = divmod_u64_u32_recip(s, f, ii->recip[i], &sub[i]);
  }
}

INLINE Bitboard idx_state_to_sq(struct IdxState *is, uint8_t *restrict sq,
    const struct IdxInfo *ii)
{
  int i = is->n;
  Bitboard occ = is->occ[i];
  for (; i < ii->numsets; i++) {
    is->occ[i] = occ;
    occ = unrank_binomial(is->sub[i], ii->mult[i], sq + ii->first[i], occ);
  }
  return occ;
}

void init_unrank(void);
void calc_factors(struct IdxInfo *ii);
uint64_t sq_to_idx(uint8_t *sq);
uint64_t capt_sq_to_idx(uint8_t *sq, int k);
void idx_to_sq_init(uint64_t idx, uint32_t *sub, const struct IdxInfo *ii);
Bitboard idx_to_sq(uint32_t *sub, uint8_t *sq);
void idx_to_sq_ii(uint32_t *sub, uint8_t *sq, const struct IdxInfo *ii);
Bitboard capt_idx_to_sq(uint32_t *sub, uint8_t *sq, int k);
void idx_state_init(struct IdxState *is, uint64_t idx, uint8_t *restrict sq,
    const struct IdxInfo *ii);

#endif
