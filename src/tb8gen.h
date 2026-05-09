/*
  Copyright (c) 2026 Ronald de Man

  This file is distributed under the terms of the GNU GPL, version 2.
*/

#ifndef TB8GEN_H
#define TB8GEN_H

#include <stdint.h>

#include "stats.h"

extern bool one_sided, wins_only;
extern int one_sided_stm;

extern struct MaxFen mf;

extern char *g_tablename;
extern bool symmetric;
extern Position g_pos;
extern uint64_t *work_g, *work_capt[MAX_SETS];
extern uint64_t g_stats[2][MAX_STATS];
extern bool g_use_rans;
extern bool used_rans;
extern const char *typename[3];
extern bool g_cleanup;

static constexpr int g_num_pawns = 0;

#endif
