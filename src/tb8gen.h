/*
  Copyright (c) 2026 Ronald de Man

  This file is distributed under the terms of the GNU GPL, version 2.
*/

#ifndef TB8GEN_H
#define TB8GEN_H

#include <stdint.h>

#include "movegen.h"
#include "stats.h"

#ifdef HAS_PAWNS
struct DtzFormat {
  bool one_sided, wins_only;
  int one_sided_stm;
};

extern struct DtzFormat dtz_format[24];
extern bool flipped;
extern int g_num_pawns;
extern char pawnstr[24][3];
extern struct Work *work_g16;
static constexpr bool has_pawns = true;
#else
static constexpr bool flipped = false;
static constexpr int g_num_pawns = 0;
static constexpr bool has_pawns = false;
#endif

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
extern int8_t g_sets[2][8];
extern uint8_t g_set_type[8];
extern struct Work *work_g, *work_capt[MAX_SETS];
extern uint64_t g_stats[2][MAX_STATS];
extern const char *typename[3];

#endif
