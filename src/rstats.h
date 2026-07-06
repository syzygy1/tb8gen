/*
  Copyright (c) 2026 Ronald de Man

  This file is distributed under the terms of the GNU GPL, version 2.
*/

#include <stdint.h>

#include "defs.h"
#include "types.h"

extern uint64_t g_stats[2][MAX_STATS];

INLINE win_to_stat(int n)
{
  return 2 + 2 * (n > DRAW_RULE);
}

INLINE loss_to_stat(int n)
{
  return MAX_STATS - 1 - n;
}
