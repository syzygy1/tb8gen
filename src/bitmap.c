/*
  Copyright (c) 2026 Ronald de Man

  This file is distributed under the terms of the GNU GPL, version 2.
*/

#include <stdio.h>
#include <stdlib.h>

#include "bitmap.h"
#include "types.h"

void table_write(struct Table *table, int stm, int kslice)
{
  char name[32];

  sprintf(name, "table.%c.%d.tmp\n", "wb"[stm], kslice);
  table->F = fopen(name, "w");
  table->type = TABLE_WRITE;
}

void table_read(struct Table *table, int stm, int kslice)
{
  char name[32];

  sprintf(name, "table.%c.%d\n", "wb"[stm], kslice);
  table->F = fopen(name, "r");
  table->type = TABLE_READ;
}

void table_close(struct Table *table)
{
  fclose(table->F);
  if (table->type == TABLE_WRITE) {

  }
}
