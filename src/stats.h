/*
  Copyright (c) 2026 Ronald de Man

  This file is distributed under the terms of the GNU GPL, version 2.
*/

#ifndef STATS_H
#define STATS_H

#include <stdio.h>

#include "tb8gen.h"
#include "hash/xxhash.h"

struct MaxFen {
  int dtz[2][2];
  char fen[2][2][48];
  bool found[2][2];
};

extern uint64_t g_stats[2][MAX_STATS];
extern XXH128_hash_t wdl_checksum, dtz_checksum[2];

void collect_stats(int stm);
void print_stats(FILE *F, int stm);
void print_max_fens(FILE *F, struct MaxFen *mf);
void calc_stats_checksums(void);
double entropy_one_sided(int stm);
double entropy_loss_only(int stm);
double entropy_win_only(int stm);

#endif
