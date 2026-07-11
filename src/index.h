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
#include "movegen.h"
#include "types.h"

static constexpr int MAX_MULT = MAX_PIECES - (has_pawns ? 3 : 2);

struct RankInfo {
  uint8_t numsets;
  uint8_t first[MAX_SETS];
  uint8_t mult[MAX_SETS];
  uint8_t last[MAX_SETS];
  uint8_t perm[MAX_SETS];
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
  alignas(64) Bitboard occ[8];
  Bitboard bb[8];
  uint32_t sub[MAX_SETS + 1];
  uint8_t sq[3];
};

struct Slice {
  uint8_t sq[3];
  int stm;
};

#ifdef HAS_PAWNS
INLINE bool is_broken(struct Slice *s)
{
  return   (king_mask(s->sq[0]) & bit(s->sq[1]))
        || s->sq[0] == s->sq[2] || s->sq[1] == s->sq[2];
}
#endif

extern struct Slice g_slice;

extern struct RankInfo ri, capt_ri[MAX_SETS];
extern int8_t g_sets[2][8];
extern uint8_t g_set_pt[8];
extern int8_t g_piece_set[2][8];
extern Bitboard Unrank2[63 * 62 / 2], Unrank3[63 * 62 * 61 / 6];
extern uint32_t Binomial[8][64];
extern uint64_t MirrorMask[64];
extern bool FlipTest[64][64];
extern const int16_t KKIdx[10][64];
extern uint8_t KKSquare[463][2];
extern int16_t KKMap[64][64];
extern uint8_t KK_transform[64][64];

INLINE __m512i flip_vertical_8xbb(__m512i x)
{
  const __m128i mask128 =
    _mm_setr_epi8( 7, 6, 5, 4, 3, 2, 1, 0, 15,14,13,12,11,10, 9, 8);
  return _mm512_shuffle_epi8(x, _mm512_broadcast_i64x2(mask128));
}

INLINE __m512i flip_horizontal_8xbb(__m512i x)
{
  const __m512i v = _mm512_set1_epi64(0x8040201008040201ULL);
  return _mm512_gf2p8affine_epi64_epi8(x, v, 0);
}

INLINE __m512i rotate270_8xbb(__m512i x)
{
  const __m512i v = _mm512_set1_epi64(0x8040201008040201ULL);
  return _mm512_gf2p8affine_epi64_epi8(v, x, 0);
}

INLINE __m512i flip_main_8xbb(__m512i x)
{
  x = flip_vertical_8xbb(x);
  x = rotate270_8xbb(x);
  return x;
}

INLINE __m512i flip_anti_8xbb(__m512i x)
{
  x = rotate270_8xbb(x);
  x = flip_vertical_8xbb(x);
  return x;
}

INLINE __m512i rotate90_8xbb(__m512i x)
{
  x = flip_vertical_8xbb(x);
  x = rotate270_8xbb(x);
  x = flip_vertical_8xbb(x);
  return x;
}

INLINE __m512i rotate180_8xbb(__m512i x)
{
  x = flip_vertical_8xbb(x);
  x = flip_horizontal_8xbb(x);
  return x;
}

INLINE Bitboard unrank_binomial(uint64_t idx, int n, Bitboard occ)
{
  if (n == 0)
    return 0;

  Bitboard b = ~occ;

  if (n == 1)
    return _pdep_u64(bit(idx), b);
  else if (n == 2)
    return _pdep_u64(Unrank2[idx], b);
  else if (n == 3)
    return _pdep_u64(Unrank3[idx], b);

  Bitboard b1 = 0;
  int r = popcnt(b) - 1;
  for (int i = n - 1; i > 0; i--) {
    while (idx < Binomial[i + 1][r])
      r--;
    idx -= Binomial[i + 1][r];
    b1 |= bit(r);
    r--;
  }
  return _pdep_u64(b1 | bit(idx), b);
}

INLINE uint64_t rank_combination(Bitboard subset, Bitboard universe)
{
  uint64_t dense = _pext_u64(subset, universe);

  uint64_t r = 0;
  for (int j = 1; dense; j++)
    r += Binomial[j][pop_lsb(&dense)];
  return r;
}

// Valid if x <= 2^N, d-1 <= 2^L and N + L <= 64.
// This should not be a problem even for 9-piece tables.
// See https://gmplib.org/~tege/divcnst-pldi94.pdf
INLINE uint64_t divmod_recip(uint64_t x, uint32_t d, uint64_t recip,
    uint32_t *rem)
{
  // FIXME: look into avoiding this special case (caller-side).
  if (d == 1) {
    *rem = 0;
    return x;
  }

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

INLINE Bitboard flip_main_diag(Bitboard b)
{
  Bitboard t;

  t  = (b ^ (b << 28)) & UINT64_C(0x0f0f0f0f00000000);
  b ^= t ^ (t >> 28);

  t  = (b ^ (b << 14)) & UINT64_C(0x3333000033330000);
  b ^= t ^ (t >> 14);

  t  = (b ^ (b << 7)) & UINT64_C(0x5500550055005500);
  b ^= t ^ (t >> 7);

  return b;
}

// Still used in rank_reflection() and in probe.c.
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

INLINE Bitboard idx_state_inc(struct IdxState *is, const struct RankInfo *ri)
{
  uint32_t *const restrict sub = is->sub;
  int i = ri->numsets;

  for (;;) {
    i--;
    if (++sub[i] < ri->factor[i]) break;
    if (i == 0)
      return 0;
    sub[i] = 0;
  }

  Bitboard occ = is->occ[i];
  for (; i < ri->numsets; i++) {
    occ |= is->bb[i + 1] = unrank_binomial(is->sub[i], ri->mult[i], occ);
    is->occ[i + 1] = occ;
  }
  return occ;
}

INLINE Bitboard idx_state_add(struct IdxState *is, uint64_t v,
    const struct RankInfo *restrict ri)
{
  uint32_t *const restrict sub = is->sub;
  int i = ri->numsets;

  for (;;) {
    uint64_t s = (uint64_t)sub[--i] + v;
    uint32_t f = ri->factor[i];

    if (s < f) {
      sub[i] = s;
      break;
    }
    if (i == 0)
      return 0;

    v = divmod_recip(s, f, ri->recip[i], &sub[i]);
  }

  Bitboard occ = is->occ[i];
  for (; i < ri->numsets; i++) {
    occ |= is->bb[i + 1] = unrank_binomial(is->sub[i], ri->mult[i], occ);
    is->occ[i + 1] = occ;
  }
  return occ;
}

#ifdef HAS_PAWNS
struct PawnIdxState {
  uint32_t sub[MAX_SETS];
  Bitboard inv_diff_bb;
};

INLINE Bitboard idx_state_pawn_add(struct IdxState *is, uint64_t v,
    struct PawnIdxState *restrict pis8, struct PawnIdxState *restrict pis16,
    const struct RankInfo *ri)
{
  uint32_t *const restrict sub = is->sub;
  int i = ri->numsets;

  for (;;) {
    uint64_t s = (uint64_t)sub[--i] + v;
    uint32_t f = ri->factor[i];

    if (s < f) {
      sub[i] = s;
      break;
    }
    if (i == 0)
      return 0;

    v = divmod_recip(s, f, ri->recip[i], &sub[i]);
  }

  Bitboard occ = is->occ[i];
  for (; i < ri->numsets; i++) {
    occ |= is->bb[i + 1] = unrank_binomial(is->sub[i], ri->mult[i], occ);
    is->occ[i + 1] = occ;
    pis8->sub[i] =
      rank_combination(is->bb[i + 1], is->occ[i] ^ pis8->inv_diff_bb);
    if (pis16)
      pis16->sub[i] =
        rank_combination(is->bb[i + 1], is->occ[i] ^ pis16->inv_diff_bb);
  }
  return occ;
}

INLINE uint64_t pawn_idx(struct PawnIdxState *pis, const struct RankInfo *ri)
{
  uint64_t idx = 0;
  for (int i = 0; i < ri->numsets; i++)
    idx = idx * ri->factor[i] + pis->sub[i];
  return idx;
}
#endif

INLINE bool idx_state_legal(const struct IdxState *is, int stm, Bitboard occ)
{
  int ksq = is->sq[stm ^ 1];
  if (has_pawns && stm == BLACK && (pawn_attacks(BLACK, is->sq[2]) & bit(ksq)))
    return false;
  for (int i = 0; g_sets[stm][i] >= 0; i++) {
    int k = g_sets[stm][i];
    Bitboard b = non_king_piece_attacks(g_set_pt[k], ksq, occ);
    if (b & is->bb[k + 1])
      return false;
  }
  return true;
}

INLINE void idx_state_to_sq(const struct IdxState *is, uint8_t *restrict sq,
    const struct RankInfo *ri)
{
  for (int k = 0; k < ri->numsets; k++) {
    Bitboard b = is->bb[k + 1];
    int i = ri->first[k];
    while (b)
      sq[i++] = pop_lsb(&b);
  }
}

void init_ranking(void);
int rank_mult(uint8_t mult[MAX_SETS]);

void transform_set_bb(int t, Bitboard *set_bb, Bitboard *set_bb2);
uint64_t rank_bb(const Bitboard *set_bb, const struct RankInfo *ri);
uint64_t rank_bb_from(const Bitboard *set_bb, uint64_t idx, int k, Bitboard occ,
    const struct RankInfo *ri);
uint64_t perm_rank_bb(const Bitboard *set_bb, const struct RankInfo *ri);
uint64_t perm_rank_bb_from(const Bitboard *set_bb, uint64_t idx, int k,
    Bitboard occ, const struct RankInfo *ri);
uint64_t rank_bb_ref(const Bitboard *set_bb, const struct RankInfo *ri);
uint64_t perm_rank_bb_ref(const Bitboard *set_bb, const struct RankInfo *ri);
uint64_t unrank_bb_ref(uint64_t idx, Bitboard *set_bb,
    const struct RankInfo *ri);
Bitboard idx_state_init(struct IdxState *is, uint64_t idx,
    const  uint8_t *restrict sq, const struct RankInfo *ri, const bool ref);
bool idx_state_mate(struct IdxState *is, int stm, Bitboard occ);
bool idx_state_has_legal_moves(struct IdxState *is, int stm, Bitboard occ);
#ifdef HAS_PAWNS
void pawn_idx_init(struct PawnIdxState *pis, Bitboard diff_bb,
    const struct IdxState *is, const struct RankInfo *ri);
#endif

void calc_factors(struct RankInfo *ri, int n);

uint64_t rank_trivial_from(const uint8_t *restrict sq, int k, Bitboard occ,
    const uint8_t *restrict first, const struct RankInfo *ri);
uint64_t rank_reflection(const uint8_t *restrict sq, Bitboard occ,
    const uint8_t *restrict first, const struct RankInfo *ri);

#endif
