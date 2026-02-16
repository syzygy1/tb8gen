/*
  Copyright (c) 2011-2013, 2025, 2026 Ronald de Man

  This file is distributed under the terms of the GNU GPL, version 2.
*/

#include "movegen.h"

Bitboard Bit[64];
Bitboard KnightAttacks[64], KingAttacks[64];
Bitboard PawnAttacks[2][64];
Bitboard SidesMask[64];

static signed char PawnDelta[2][2][2] = {
  { {  7,  15 }, {  9,  17 } },
  { { -9, -17 }, { -7, -15 } }
};

static signed char KnightDelta[8][2] = {
  { -17, -33 }, { -15, -31 }, { -10, -18 }, { -6, -14 },
  {   6,  14 }, {  10,  18 }, {  15,  31 }, { 17,  33 }
};

static signed char BishopDelta[4][2] = {
  { -9, -17 }, { -7, -15 }, { 7, 15 }, { 9, 17 }
};

static signed char RookDelta[4][2] = {
  { -8, -16 }, { -1, -1 }, { 1, 1 }, { 8, 16 }
};

static signed char KingDelta[8][2] = {
  { -9, -17 }, { -8, -16 }, { -7, -15 }, { -1, -1 },
  {  1,   1 }, {  7,  15 }, {  8,  16 }, {  9, 17 }
};

INLINE bool valid(int sq, signed char delta[2])
{
  int sq88 = sq + (sq & ~7);
  return !((sq88 + delta[1]) & 0x88);
}

#include "magic.c"
//#include "hyper.c"
//#include "bmi2.c"

const char PieceChar[] = " PNBRQK  pnbrqk ";

static Bitboard calc_attacks(int sq, signed char delta[][2], int num)
{
  Bitboard bb = 0;

  for (int d = 0; d < num; d++)
    if (valid(sq, delta[d]))
      bb |= bit(sq + delta[d][0]);

  return bb;
}

void init_movegen(void)
{
  for (int sq = 0; sq < 64; sq++)
    Bit[sq] = 1ULL << sq;

  for (int sq = 0; sq < 64; sq++) {
    PawnAttacks[WHITE][sq] = calc_attacks(sq, PawnDelta[WHITE], 2);
    PawnAttacks[BLACK][sq] = calc_attacks(sq, PawnDelta[BLACK], 2);
    KnightAttacks[sq]      = calc_attacks(sq, KnightDelta,      8);
    KingAttacks[sq]        = calc_attacks(sq, KingDelta,        8);
  }

  init_sliding_attacks();
}

bool has_legal_moves(Position *pos)
{
  for (int i = 0; i < pos->num; i++) {
    if ((pos->pt[i] >> 3) != pos->stm)
      continue;
    int from = pos->sq[i];
    Bitboard b = piece_moves(pos->pt[i], from, pos->occ);
    while (b) {
      int to = pop_lsb(&b);
      if (do_move(pos, from, to, i)) {
        undo_move(pos, from, to, i);
        return true;
      }
    }
  }
  return false;
}

bool has_legal_caps(Position *pos)
{
  for (int i = 0; i < pos->num; i++) {
    if ((pos->pt[i] >> 3) != pos->stm)
      continue;
    int from = pos->sq[i];
    Bitboard b = piece_attacks(pos->pt[i], from, pos->occ) & pos->occ;
    while (b) {
      int to = pop_lsb(&b);
      int j = piece_idx(pos, to);
      if (!((pos->pt[i] ^ pos->pt[j]) & 8)) continue;
      if (do_capture(pos, from, to, i, j)) {
        undo_capture(pos, from, to, i, j);
        return true;
      }
    }
  }
  return false;
}
