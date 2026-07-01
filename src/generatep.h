/*
  Copyright (c) 2026 Ronald de Man

  This file is distributed under the terms of the GNU GPL, version 2.
*/

#ifndef GENERATEP_H
#define GENERATEP_H

#include <inttypes.h>
#include <stdint.h>

#include "index.h"
#include "movegen.h"
#include "types.h"

extern uint64_t capt_cnt[2][5], pawn_cnt[5];
extern int max_iteration;

INLINE int stats_n(int n)
{
  return 2 + n + 2 * (n > DRAW_RULE);
}

void generate(void);
void init_generation_work(void);
void delete_intermediate_slices(void);
void cleanup_generation(void);

#endif
