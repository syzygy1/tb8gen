/*
  Copyright (c) 2026 Ronald de Man

  This file is distributed under the terms of the GNU GPL, version 2.
*/

#ifndef PERMUTEP_H
#define PERMUTEP_H

#include <inttypes.h>

void init_permute_pawn_p(void);
void permute_pawn_p(void *tb_table, void *table, uint8_t *best, int type,
    bool wide);

#endif
