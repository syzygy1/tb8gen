/*
  Copyright (c) 2026 Ronald de Man

  This file is distributed under the terms of the GNU GPL, version 2.
*/

#ifndef PERMUTE10_H
#define PERMUTE10_H

#include <inttypes.h>

void init_permute_piece_10(int k);
void permute_piece_10(void *tb_table, void *table, uint8_t *best, int type,
    bool wide, void *v);

#endif
