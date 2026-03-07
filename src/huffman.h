/*
  Copyright (c) 2026 Ronald de Man

  This file is distributed under the terms of the GNU GPL, version 2.
*/

#ifndef HUFFMAN_H
#define HUFFMAN_H

#include <inttypes.h>

#include "defs.h"

struct HuffCode {
  int64_t *freq;
  int64_t nfreq[MAXSYMB];
  int map[MAXSYMB];
  int inv[MAXSYMB];
  int length[MAXSYMB];
  int num_syms, num, max_len, min_len;
  int base[33];
  int offset[33];
};

struct HuffCode *create_code(int64_t *freq, int num_syms);
uint64_t calc_size(struct HuffCode *c);
void free_code(struct HuffCode *c);

#endif
