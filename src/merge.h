/*
  Copyright (c) 2026 Ronald de Man

  This file is distributed under the terms of the GNU GPL, version 2.
*/

#ifndef MERGE_H
#define MERGE_H

struct MergeInfo {
  bool v_wdl[5];
  bool wide;
  union {
    uint8_t v_u8[MAX_STATS];
    uint16_t v_u16[MAX_STATS];
  };
  union {
    uint16_t v_inv_u8[256];
    uint16_t v_inv_u16[MAX_STATS];
  };
};

extern struct MergeInfo mi;

void merge(int merge_type);
void collect_stats(int stm);
void print_stats(int stm);
double entropy_one_sided(int stm), entropy_loss_only(int stm),
       entropy_win_only(int stm);

#endif
