/*
  Copyright (c) 2026 Ronald de Man

  This file is distributed under the terms of the GNU GPL, version 2.
*/

#ifndef RGENERATE_H
#define RGENERATE_H

#include <inttypes.h>

#include "defs.h"
#include "threads.h"
#include "types.h"

extern uint8_t *g_table[2];
extern size_t table_size, table_sub_size[MAX_SETS];
extern struct Work work_g_dynamic, work_g_static;

INLINE int stats_n(int n)
{
  return 2 + n + 2 * (n > DRAW_RULE);
}

void generate(void);
//void init_generation_work(void);
//void delete_intermediate_slices(void);
//void cleanup_generation(void);

#endif
