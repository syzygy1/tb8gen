/*
  Copyright (c) 2011-2013, 2025 Ronald de Man

  This file is distributed under the terms of the GNU GPL, version 2.
*/

#ifndef MAGIC_H
#define MAGIC_H

#ifdef MAGIC

extern Bitboard AttackTable[97264];

struct Magic {
  Bitboard *data;
  Bitboard mask;
  uint64_t magic;
};

extern struct Magic BishopMagic[64];
extern struct Magic RookMagic[64];

INLINE Bitboard bishop_attacks(int sq, Bitboard occ)
{
  struct Magic *mag = &BishopMagic[sq];
  return mag->data[((occ & mag->mask) * mag->magic) >> (64-9)];
}

INLINE Bitboard rook_attacks(int sq, Bitboard occ)
{
  struct Magic *mag = &RookMagic[sq];
  return mag->data[((occ & mag->mask) * mag->magic) >> (64-12)];
}

INLINE Bitboard queen_attacks(int sq, Bitboard occ)
{
  return bishop_attacks(sq, occ) | rook_attacks(sq, occ);
}

#endif

#endif
