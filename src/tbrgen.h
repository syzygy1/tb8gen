/*
  Copyright (c) 2026 Ronald de Man

  This file is distributed under the terms of the GNU GPL, version 2.
*/

#ifndef TBRGEN_H
#define TBRGEN_H

#include <stdint.h>

#include "movegen.h"
#include "rstats.h"

#ifdef HAS_PAWNS
struct DtzFormat {
  bool one_sided, wins_only;
  int one_sided_stm;
};

extern struct DtzFormat dtz_format[24];
extern bool flipped;
extern int g_num_pawns;
extern char pawnstr[24][3];
#else
static constexpr bool flipped = false;
static constexpr int g_num_pawns = 0;
#endif

extern struct MaxFen mf;

extern int one_sided_stm;
extern bool one_sided, wins_only;

extern bool g_use_rans;
extern bool used_rans;
extern bool symmetric;

extern char *g_tablename;
extern char *g_output_dir;
extern size_t table_size;
extern size_t table_diagonal;

#endif
