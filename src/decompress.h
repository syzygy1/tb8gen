/*
  Copyright (c) 2011-2016, 2018, 2025 Ronald de Man

  This file is distributed under the terms of the GNU GPL, version 2.
*/

#ifndef DECOMPRESS_H
#define DECOMPRESS_H

#include "probe.h"

void decompress_table(struct TbTable2 *table, uint8_t *decomp_table,
    size_t size);
//void print_decompression_info(struct EncInfo *ei);
//void print_freq_table(struct EncInfo *ei);

#endif
