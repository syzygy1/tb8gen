/*
  Copyright (c) 2026 Ronald de Man

  This file is distributed under the terms of the GNU GPL, version 2.
*/

#ifndef PERMUTE462_H
#define PERMUTE462_H

#include <inttypes.h>

void init_permute_piece_462(void);
void permute_piece_462(void *tb_table, void *table, uint8_t *best, int type,
    bool wide);

void init_permute_pawn_462(void);
void *init_permute_rank_462(int rank, void *tb_table, bool wide);
void permute_pawn_dtz_462(void *tb_table, void *table, uint8_t *best, int rank,
    bool wide);

#endif
