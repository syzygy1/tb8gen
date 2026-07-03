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
#ifdef HAS_PAWNS
extern XXH128_hash_t dtz_partial_checksum[24];
#endif

void collect_stats(int stm);
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
