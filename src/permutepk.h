/*
  Copyright (c) 2026 Ronald de Man

  This file is distributed under the terms of the GNU GPL, version 2.
*/

#ifndef PERMUTEPK_H
#define PERMUTEPK_H

#include <inttypes.h>

void init_permute_pawn_pk(void);
void permute_pawn_pk(void *tb_table, void *table, uint8_t *best, int type,
    bool wide);

#endif
