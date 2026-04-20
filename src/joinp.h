/*
  Copyright (c) 2026 Ronald de Man

  This file is distributed under the terms of the GNU GPL, version 2.
*/

#ifndef JOINP_H
#define JOINP_H

extern uint8_t g_dist_format;

void compress_slice_pk(void);
void compress_slice_p(void);
void join_slices_pk(void);
void join_slices_p(void);

#endif
