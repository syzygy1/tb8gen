/*
  Copyright (c) 2011-2013, 2025 Ronald de Man

  This file is distributed under the terms of the GNU GPL, version 2.
*/

#ifndef BMI2_PLAIN_H
#define BMI2_PLAIN_H

#ifdef BMI2_PLAIN

#include <immintrin.h>

extern Bitboard RookMasks[64];
extern Bitboard BishopMasks[64];
extern Bitboard *RookAttacks[64];
extern Bitboard *BishopAttacks[64];

INLINE Bitboard bishop_attacks(int sq, Bitboard occ)
{
  return BishopAttacks[sq][_pext_u64(occ, BishopMasks[sq])];
}

INLINE Bitboard rook_attacks(int sq, Bitboard occ)
{
  return RookAttacks[sq][_pext_u64(occ, RookMasks[sq])];
}

INLINE Bitboard queen_attacks(int sq, Bitboard occ)
{
  return bishop_attacks(sq, occ) | rook_attacks(sq, occ);
}

#endif

#endif
