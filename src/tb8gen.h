/*
  Copyright (c) 2026 Ronald de Man

  This file is distributed under the terms of the GNU GPL, version 2.
*/

#ifndef TB8GEN_H
#define TB8GEN_H

#include <stdint.h>

#include "stats.h"

extern struct MaxFen mf;

extern int one_sided_stm;
extern bool one_sided, wins_only;

extern bool g_use_rans;
extern bool used_rans;
extern bool symmetric;
extern bool g_cleanup;

extern char *g_tablename;
extern char *g_output_dir;
extern Position g_pos;
extern uint64_t *work_g, *work_capt[MAX_SETS];
extern uint64_t g_stats[2][MAX_STATS];
extern const char *typename[3];

static constexpr int g_num_pawns = 0;

#endif
