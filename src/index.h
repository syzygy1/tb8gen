/*
  Copyright (c) 2026 Ronald de Man

  This file is distributed under the terms of the GNU GPL, version 2.
*/

#ifndef INDEX_H
#define INDEX_H

#include <assert.h>
#include <inttypes.h>
#include <string.h>
#include <x86intrin.h>

#include "defs.h"
#include "tb8gen.h"
#include "types.h"

static constexpr int MAX_MULT = MAX_PIECES - (has_pawns ? 3 : 2);

struct RankInfo {
  uint8_t numsets;
  uint8_t first[MAX_SETS];
  uint8_t mult[MAX_SETS];
  uint8_t last[MAX_SETS];
  uint8_t transition_id[MAX_SETS];
  uint32_t factor[MAX_SETS];
  uint64_t recip[MAX_SETS];
  uint64_t sizes[2];
};

// Indexing used by permute10.c deviates a bit.
// Now the non-leading king is permuted together with the non-pawn pieces.
struct RankInfo10 {
  uint8_t numsets;
  uint8_t k2;                 // Set index of the non-leading king.
  uint8_t first[MAX_SETS];
  uint8_t mult[MAX_SETS];
  uint32_t factor[MAX_SETS];
};

extern struct RankInfo rank_info_61[32];
extern struct RankInfo rank_info_62[64];
extern struct RankInfo rank_info_63[64];

struct IdxState {
  uint32_t sub[MAX_SETS];
  Bitboard occ[MAX_SETS];
  int n;
};

extern struct RankInfo ri, capt_ri[MAX_SETS];
extern int pc_to_set[MAX_PIECES];
extern Bitboard Unrank2[62 * 61 / 2], Unrank3[62 * 61 * 60 / 6];
extern uint32_t Binomial[8][64];
extern uint64_t MirrorMask[64];
extern bool FlipTest[64][64];
extern const int16_t KKIdx[10][64];
extern uint8_t KKSquare[462][2];
extern int16_t KKMap[64][64];

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

#undef SORT2

INLINE Bitboard unrank_binomial(uint64_t idx, int n, uint8_t *restrict sq,
    Bitboard occ)
{
  if (n == 0)
    return occ;

  assume(n > 0 && n <= MAX_MULT);

  Bitboard b = ~occ;

  if (n == 1) {
    Bitboard b1 = _pdep_u64(bit(idx), b);
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
    int r = popcnt(b) - 1;
    for (int i = n - 1; i > 0; i--) {
      while (idx < Binomial[i + 1][r])
        r--;
      idx -= Binomial[i + 1][r];
      b1 |= bit(r);
      r--;
    }
    b1 = _pdep_u64(b1 | bit(idx), b);
    occ |= b1;
    while (b1)
      *sq++ = pop_lsb(&b1);
  }

  return occ;
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

INLINE int rank_among_free(uint8_t sq, Bitboard occ)
{
  return sq - popcnt(occ & ((1ULL << sq) - 1));
}

// We expect a normalized position.
INLINE uint64_t sq_to_idx_helper(const uint8_t *restrict sq, uint64_t idx,
    Bitboard occ, const struct RankInfo *ii)
{
  for (int k = 0; k < ii->numsets; k++) {
    size_t s;
    int i = ii->first[k];
    int m = ii->mult[k];
    assume(m > 0 && m <= MAX_MULT);
    if (m == 1) {
      s = rank_among_free(sq[i], occ);
      occ |= bit(sq[i]);
    } else if (m == 2) {
      Bitboard b = ~occ;
      Bitboard b1 = bit(sq[i]) | bit(sq[i + 1]);
      occ |= b1;
      b1 = _pext_u64(b1, b);
      s = pop_lsb(&b1);
      s += Binomial[2][lsb(b1)];
    } else {
      Bitboard b = ~occ, b1 = 0;
      for (int j = 0; j < m; j++)
        b1 |= bit(sq[i + j]);
      occ |= b1;
      b1 = _pext_u64(b1, b);
      s = 0;
      for (int j = 1; b1; j++)
        s += Binomial[j][pop_lsb(&b1)];
    }
    idx = idx * ii->factor[k] + s;
  }

  return idx;
}

INLINE uint64_t mirror_diagonal_u64(uint64_t v)
{
  return  ((v & 0x0707070707070707ULL) << 3)
        | ((v & 0x3838383838383838ULL) >> 3);
}

INLINE void mirror_diagonal(uint8_t *sq)
{
  uint64_t v;
  memcpy(&v, sq, 8);
  v = mirror_diagonal_u64(v);
  memcpy(sq, &v, 8);
}

INLINE void normalize(uint8_t *sq)
{
  uint64_t v;
  memcpy(&v, sq, 8);
  v ^= MirrorMask[sq[0]];
  if (FlipTest[sq[0]][sq[1]])
    v = mirror_diagonal_u64(v);
  memcpy(sq, &v, 8);
}

INLINE void normalize_quadrant(uint8_t *sq)
{
  uint64_t v;
  memcpy(&v, sq, 8);
  v ^= MirrorMask[sq[0]];
  memcpy(sq, &v, 8);
}

INLINE void mirror_diagonal2(const uint8_t *restrict sq, uint8_t *restrict sq2)
{
  uint64_t v;
  memcpy(&v, sq, 8);
  v = mirror_diagonal_u64(v);
  memcpy(sq2, &v, 8);
}

INLINE void normalize2(const uint8_t *restrict sq, uint8_t *restrict sq2)
{
  uint64_t v;
  memcpy(&v, sq, 8);
  v ^= MirrorMask[sq[0]];
  if (FlipTest[sq[0]][sq[1]])
    v = mirror_diagonal_u64(v);
  memcpy(sq2, &v, 8);
}

INLINE void idx_state_inc(struct IdxState *is, const struct RankInfo *ii)
{
  uint32_t *restrict sub = is->sub;
  int i = ii->numsets - 1;
  for (; ++sub[i] >= ii->factor[i] && i > 0; i--)
    sub[i] = 0;
  is->n = i;
}

INLINE void idx_state_add(struct IdxState *is, uint64_t v,
    const struct RankInfo *restrict ii)
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

    v = divmod_recip(s, f, ii->recip[i], &sub[i]);
  }
}

INLINE Bitboard idx_state_to_sq(struct IdxState *is, uint8_t *restrict sq,
    const struct RankInfo *ii)
{
  int i = is->n;
  Bitboard occ = is->occ[i];
  for (; i < ii->numsets; i++) {
    is->occ[i] = occ;
    occ = unrank_binomial(is->sub[i], ii->mult[i], sq + ii->first[i], occ);
  }
  return occ;
}

void init_ranking(void);
int rank_mult(uint8_t mult[MAX_SETS]);

void calc_factors(struct RankInfo *ri, int n);
uint64_t sq_to_idx(const uint8_t *sq);
uint64_t sq_to_idx_ref(const uint8_t *sq);
uint64_t capt_sq_to_idx(const uint8_t *sq, int k);
void idx_state_init(struct IdxState *is, uint64_t idx, uint8_t *restrict sq,
    const struct RankInfo *ri);

uint64_t rank_trivial_from(const uint8_t *restrict sq, int k, Bitboard occ,
    const uint8_t *restrict first, const struct RankInfo *ri);
uint64_t rank_reflection(const uint8_t *restrict sq, Bitboard occ,
    const uint8_t *restrict first, const struct RankInfo *ri);
Bitboard unrank_reflection(uint64_t idx, uint8_t *restrict sq, Bitboard occ,
    const struct RankInfo *ri);

#endif
