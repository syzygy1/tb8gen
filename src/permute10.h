/*
  Copyright (c) 2026 Ronald de Man

  This file is distributed under the terms of the GNU GPL, version 2.
*/

#ifndef PERMUTE10_H
#define PERMUTE10_H

#include <inttypes.h>

void init_permute_piece_10(int k);
void permute_piece_10(void *tb_table, void *table, uint8_t *best, int type,
    bool wide);

void init_permute_pawn_10(void);
void *init_permute_rank_10(int rank, void *tb_table, bool wide);
void permute_pawn_dtz_10(void *tb_table, void *table, uint8_t *best, int rank,
    bool wide);

#endif
