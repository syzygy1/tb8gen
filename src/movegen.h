/*
  Copyright (c) 2011-2013, 2025 Ronald de Man

  This file is distributed under the terms of the GNU GPL, version 2.
*/

#ifndef MOVEGEN_H
#define MOVEGEN_H

#include <stdbit.h>
#include <stddef.h>

#include "defs.h"
#include "types.h"

extern const char PieceChar[];

struct Position {
  Bitboard occ;
  uint8_t sq[8], pt[8];
  int num, stm;
  int8_t pcs[2][8];
#ifdef HAS_PAWNS
  int8_t pawns[2][8];
#endif
};

typedef struct Position Position;

extern Bitboard Bit[64], KnightAttacks[64], KingAttacks[64];
extern Bitboard PawnAttacks[2][64], SidesMask[64];

#ifdef ATOMIC
extern Bitboard AtomMask[64];
#endif

INLINE Bitboard bit(int sq)
{
  return 1ULL << sq;
//  return Bit[sq];
}

INLINE int lsb(Bitboard b)
{
  return stdc_trailing_zeros(b);
}

INLINE int pop_lsb(Bitboard *b)
{
  int s = lsb(*b);
  *b &= *b - 1;
  return s;
}

INLINE int popcnt(Bitboard b)
{
  return stdc_count_ones(b);
}

#include "magic.h"
//#include "hyper.h"
//#include "bmi2.h"

INLINE Bitboard pawn_attacks(int c, int sq)
{
  return PawnAttacks[c][sq];
}

INLINE Bitboard knight_attacks(int sq)
{
  return KnightAttacks[sq];
}

INLINE Bitboard king_attacks(int sq)
{
  return KingAttacks[sq];
}

INLINE Bitboard piece_attacks(int pt, int sq, Bitboard occ)
{
  switch (pt & 7) {
  case PAWN:
#ifndef HAS_PAWNS
    unreachable();
#endif
    return pawn_attacks(pt >> 3, sq);
  case KNIGHT:
    return knight_attacks(sq);
  case BISHOP:
    return bishop_attacks(sq, occ);
  case ROOK:
    return rook_attacks(sq, occ);
  case QUEEN:
    return queen_attacks(sq, occ);
  case KING:
    return king_attacks(sq);
  default:
    unreachable();
  }
}

// only used in rtbgen / rtbver
INLINE Bitboard non_king_piece_attacks(int pt, int sq, Bitboard occ)
{
  switch (pt & 0x07) {
  case KNIGHT:
    return knight_attacks(sq);
  case BISHOP:
    return bishop_attacks(sq, occ);
  case ROOK:
    return rook_attacks(sq, occ);
  case QUEEN:
    return queen_attacks(sq, occ);
  default:
    unreachable();
  }
}

INLINE Bitboard non_king_piece_moves(int pt, int sq, Bitboard occ)
{
  return non_king_piece_attacks(pt, sq, occ) & ~occ;
}

INLINE bool rank18(int sq)
{
  return (unsigned)(sq - 8) >= 48;
}

INLINE Bitboard piece_moves(int pt, int sq, Bitboard occ)
{
  switch (pt & 7) {
  case PAWN:
#ifndef HAS_PAWNS
    unreachable();
#endif
    Bitboard b = 0;
    int fwd = (pt & 8) ? -8 : 8;
    if (!(bit(sq + fwd) & occ)) {
      b = bit(sq + fwd);
      if (rank18(sq - fwd) && !(bit(sq + 2 * fwd) & occ))
	b |= bit(sq + 2 * fwd);
    }
    return b;
  case KNIGHT:
    return knight_attacks(sq) & ~occ;
  case BISHOP:
    return bishop_attacks(sq, occ) & ~occ;
  case ROOK:
    return rook_attacks(sq, occ) & ~occ;
  case QUEEN:
    return queen_attacks(sq, occ) & ~occ;
  default:
    return king_attacks(sq) & ~occ;
  }
}

#ifdef HAS_PAWNS
extern uint8_t KingIndex[2];
#if 0
INLINE uint8_t king_index(int stm)
{
  return KingIndex[stm];
}
#else
INLINE uint8_t king_index(Position *pos, int stm)
{
  for (int i = 0; ; i++)
    if (pos->pt[i] == (stm << 3 | KING))
      return i;
}
#endif
#else
INLINE uint8_t king_index(Position *pos, int stm)
{
  (void)pos;
  return stm;
}
#endif

INLINE bool is_attacked_by(Position *pos, int sq, int stm)
{
  Bitboard b = bit(sq);
  for (int i = 0; i < pos->num; i++)
    if (   (pos->pt[i] >> 3) == stm
        && (piece_attacks(pos->pt[i], pos->sq[i], pos->occ) & b))
      return true;
  return false;
}

INLINE bool my_king_attacked(Position *pos)
{
  return is_attacked_by(pos, pos->sq[king_index(pos, pos->stm)], pos->stm ^ 1);
}

INLINE bool opp_king_attacked(Position *pos)
{
  return is_attacked_by(pos, pos->sq[king_index(pos, pos->stm ^ 1)], pos->stm);
}

INLINE int piece_idx(Position *pos, int sq)
{
  for (int i = 0; i < pos->num; i++)
    if (pos->sq[i] == sq) return i;
  assume(false); // signal to the compiler that this path is unreachable
  return 0;
}

INLINE void undo_capture(Position *pos, int from, int to, int i, int j)
{
  pos->pt[pos->num] = pos->pt[j];   // move back piece in j to last position
  pos->sq[pos->num] = pos->sq[j];
  pos->pt[j] = pos->pt[++pos->num]; // restore captured piece type
  pos->sq[j] = to;
  pos->sq[i] = from;
  pos->occ ^= bit(from);
  pos->stm ^= 1;
}

INLINE bool do_capture(Position *pos, int from, int to, int i, int j)
{
  pos->occ ^= bit(from);
  pos->sq[i] = to;
  pos->pt[pos->num--] = pos->pt[j]; // save captured piece type
  pos->pt[j] = pos->pt[pos->num];   // move last piece in list to position j
  pos->sq[j] = pos->sq[pos->num];
  pos->stm ^= 1;

  if (opp_king_attacked(pos)) {
    undo_capture(pos, from, to, i, j);
    return false;
  }

  return true;
}

INLINE void undo_move(Position *pos, int from, int to, int i)
{
  pos->occ ^= bit(from) ^ bit(to);
  pos->sq[i] = from;
  pos->stm ^= 1;
}

INLINE bool do_move(Position *pos, int from, int to, int i)
{
  pos->occ ^= bit(from) ^ bit(to);
  pos->sq[i] = to;
  pos->stm ^= 1;
  if (opp_king_attacked(pos)) {
    undo_move(pos, from, to, i);
    return false;
  }
  return true;
}

#ifdef HAS_PAWNS
// to is square of captured pawn
INLINE void undo_ep_capture(Position *pos, int from, int to, int i, int j)
{
  int ep = to ^ 8;
  pos->pt[pos->num] = pos->pt[j];   // move back piece in j to last position
  pos->sq[pos->num] = pos->sq[j];
  pos->pt[j] = pos->pt[++pos->num]; // restore captured piece type
  pos->sq[j] = to;
  pos->sq[i] = from;
  pos->occ ^= bit(from) ^ bit(to) ^ bit(ep);
  pos->stm ^= 1;
}

// to is square of captured pawn
INLINE bool do_ep_capture(Position *pos, int from, int to, int i, int j)
{
  int ep = to ^ 8;
  pos->occ ^= bit(from) ^ bit(to) ^ bit(ep);
  pos->sq[i] = ep;
  pos->pt[pos->num--] = pos->pt[j]; // save captured piece type
  pos->pt[j] = pos->pt[pos->num];   // move last piece in list to position j
  pos->sq[j] = pos->sq[pos->num];
  pos->stm ^= 1;

  if (opp_king_attacked(pos)) {
    undo_ep_capture(pos, from, to, i, j);
    return false;
  }

  return true;
}
#endif

bool has_legal_moves(Position *pos);
bool has_legal_caps(Position *pos);

void init_movegen(void);
#endif
