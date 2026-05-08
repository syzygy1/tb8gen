/*
  Copyright (c) 2026 Ronald de Man

  This file is distributed under the terms of the GNU GPL, version 2.
*/

#ifndef GENERATEP_H
#define GENERATEP_H

#include <inttypes.h>
#include <stdint.h>

#include "movegen.h"
#include "types.h"

extern uint64_t capt_cnt[2][5], pawn_cnt[5];
extern int max_iteration;

INLINE bool is_broken(struct Position *pos)
{
  return   (king_mask(pos->sq[0]) & bit(pos->sq[1]))
        || pos->sq[0] == pos->sq[2] || pos->sq[1] == pos->sq[2];
}

INLINE int stats_n(int n)
{
  return 2 + n + 2 * (n > DRAW_RULE);
}

void generate(void);

#endif
