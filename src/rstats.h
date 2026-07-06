/*
  Copyright (c) 2026 Ronald de Man

  This file is distributed under the terms of the GNU GPL, version 2.
*/

#ifndef RSTATS_H
#define RSTATS_H

#include <stdint.h>

#include "defs.h"

void reset_stats(void);
void collect_stats(int stm);
void collect_stats_before_reduce(int stm, int n);
extern uint64_t g_stats[2][MAX_STATS];

#endif
