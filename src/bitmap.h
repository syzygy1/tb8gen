/*
  Copyright (c) 2026 Ronald de Man

  This file is distributed under the terms of the GNU GPL, version 2.
*/

#ifndef BITMAP_H
#define BITMAP_H

#include <stdio.h>

#include "types.h"

constexpr int TABLE_READ = 0;
constexpr int TABLE_WRITE = 1;

struct Table {
  FILE *F;
  int type;
};

struct Bitmap {
  FILE *F;
};

void table_write(struct Table *table, int stm, int kslice);
void table_read(struct Table *table, int stm, int kslice);
int table_get(struct Table *table);
void table_put(struct Table *table, int v);
void table_close(struct Table *table);

#endif
