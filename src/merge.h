/*
  Copyright (c) 2026 Ronald de Man

  This file is distributed under the terms of the GNU GPL, version 2.
*/

#ifndef MERGE_H
#define MERGE_H

enum { MERGE_SAVE, MERGE_COMPRESS };

void merge(int merge_type);
void collect_stats(int stm);
double entropy_one_sided(int stm), entropy_loss_only(int stm),
       entropy_win_only(int stm);

#endif
