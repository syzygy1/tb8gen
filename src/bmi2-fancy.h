#ifndef BMI2_FANCY_H
#define BMI2_FANCY_H

#ifdef BMI2_FANCY

#include <immintrin.h>

extern Bitboard RookMasks[64], RookMasks2[64];
extern Bitboard BishopMasks[64], BishopMasks2[64];
extern uint16_t *RookAttacks[64];
extern uint16_t *BishopAttacks[64];

INLINE Bitboard bishop_attacks(int sq, Bitboard occ)
{
  return _pdep_u64(BishopAttacks[sq][_pext_u64(occ, BishopMasks[sq])],
      BishopMasks2[sq]);
}

INLINE Bitboard rook_attacks(int sq, Bitboard occ)
{
  return _pdep_u64(RookAttacks[sq][_pext_u64(occ, RookMasks[sq])],
      RookMasks2[sq]);
}

INLINE Bitboard queen_attacks(int sq, Bitboard occ)
{
  return bishop_attacks(sq, occ) | rook_attacks(sq, occ);
}

#endif

#endif
