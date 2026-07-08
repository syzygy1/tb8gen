/*
  Copyright (c) 2026 Ronald de Man

  This file is distributed under the terms of the GNU GPL, version 2.
*/

#ifndef RSTATS_H
#define RSTATS_H

#include <stdint.h>

#include "defs.h"
#include "hash/xxhash.h"

struct MaxFen {
  int dtz[2][2];
  char fen[2][2][48];
  bool found[2][2];
};

extern uint64_t g_stats[2][MAX_STATS];

void reset_stats(void);
void collect_stats(int stm);
void collect_stats_before_reduce(int stm, int n);
void print_stats(FILE *F, int stm);
void print_max_fens(FILE *F, struct MaxFen *mf);

XXH128_hash_t freq_to_hash(int n, uint64_t *f0, uint64_t *f1, uint64_t *f2,
    uint64_t *f3);
void calc_stats_checksums(void);
#ifdef HAS_PAWNS
void calc_partial_stats_checksum(int q, uint64_t stats[2][MAX_STATS]);
#endif

double entropy_one_sided(int stm);
double entropy_loss_only(int stm);
double entropy_win_only(int stm);

#endif
