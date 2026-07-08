/*
  Copyright (c) 2011-2018, 2024, 2026 Ronald de Man

  This file is distributed under the terms of the GNU GPL, version 2.
*/

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "defs.h"
#include "reduce.h"
#include "rgenerate.h"
#include "rstats.h"
#include "tbrgen.h"
#include "types.h"
#include "util.h"

int epoch;
int reduce_cnt_win[64];
int reduce_cnt_loss[64];

void *transform_v;

void transform(struct ThreadData *thread)
{
  uint8_t *restrict p = g_table[WHITE];
  uint8_t *restrict q = g_table[BLACK];
  uint8_t *v = transform_v;

  for (uint64_t idx = thread->begin, end = thread->end; idx < end; idx++) {
    p[idx] = v[p[idx]];
    q[idx] = v[q[idx]];
  }
}

#if 0
void transform_table_u8(struct ThreadData *thread)
{
  uint64_t idx;
  uint64_t end = thread->end;
  uint8_t *v = transform;
  uint8_t *table = transform_tbl_u8;

  for (idx = thread->begin; idx < end; idx++)
    table[idx] = v[table[idx]];
}

void transform_table_u16(struct ThreadData *thread)
{
  uint64_t idx;
  uint64_t end = thread->end;
  uint16_t *v = transform_v_u16;
  uint16_t *table = transform_tbl_u16;

  for (idx = thread->begin; idx < end; idx++)
    table[idx] = v[table[idx]];
}
#endif

static void init_save_map(uint8_t v[256], int n)
{
  memset(v, 0, 256);

  int first_win = 1;
  int first_loss = 0;

  if (epoch == 0) {
    v[RAM_ILLEGAL] = RAM_ILLEGAL;
    v[RAM_CAPT_WIN] = RAM_CAPT_WIN;
    v[RAM_PAWN_WIN] = RAM_PAWN_WIN;
    v[RAM_CAPT_CWIN] = RAM_CAPT_CWIN;
    v[RAM_PAWN_CWIN] = RAM_PAWN_CWIN;
    v[RAM_CAPT_DRAW] = RAM_CAPT_DRAW;
    v[RAM_CAPT_BLOSS] = RAM_CAPT_BLOSS;
    v[RAM_PAWN_DRAW] = RAM_PAWN_DRAW;
  } else {
    first_win = reduce_cnt_win[epoch - 1] - 250;
    first_loss = reduce_cnt_loss[epoch - 1] + 4;
  }

  for (int i = first_win; i < n - 1; i++) {
    uint8_t b = win_to_byte(i);
    if (b != RAM_REDUCED_CWIN && b != RAM_REDUCED_WIN)
      v[b] = b;
  }
  for (int i = first_loss; i < n; i++) {
    uint8_t b = loss_to_byte(i);
    if (   b != RAM_REDUCED_LOSS
        && b != RAM_REDUCED_CAPT_BLOSS
        && b != RAM_REDUCED_BLOSS)
      v[b] = b;
  }
}

static void save_table(uint8_t *table, int stm, int n)
{
  char name[64];
  uint8_t v[256];

  sprintf(name, "%s.%c.%d", g_tablename, "wb"[stm], epoch);

  FILE *F;
  if (!(F = fopen(name, "wb"))) {
    fprintf(stderr, "Could not open %s for writing.\n", name);
    exit(EXIT_FAILURE);
  }

  init_save_map(v, n);
  write_data_transform_u8(F, table, table_size, v);
  fclose(F);
}

static void init_reduce_map(uint8_t v[256], int n)
{
  memset(v, 0, 256);

  int next_reduce_cnt_win = n + 249;
  int next_reduce_cnt_loss = n - 4;

  if (epoch == 0) {
    v[RAM_UNRESOLVED] = RAM_UNRESOLVED;
    v[RAM_ILLEGAL] = RAM_ILLEGAL;
    v[RAM_CAPT_WIN] = RAM_CAPT_WIN;
    v[RAM_CAPT_DRAW] = RAM_CAPT_DRAW;
    v[RAM_CAPT_CWIN] = RAM_REDUCED_CAPT_CWIN;
    v[RAM_CAPT_BLOSS] = RAM_REDUCED_CAPT_BLOSS;
    v[RAM_PAWN_WIN] = RAM_REDUCED_WIN;
    v[RAM_PAWN_CWIN] = RAM_REDUCED_CWIN;
    v[RAM_PAWN_DRAW] = RAM_UNRESOLVED;

    for (int i = 1; i <= DRAW_RULE; i++)
      v[win_to_byte(i)] = RAM_REDUCED_WIN;
    for (int i = DRAW_RULE + 1; i < n - 1; i++)
      v[win_to_byte(i)] = RAM_REDUCED_CWIN;

    for (int i = 0; i <= DRAW_RULE; i++)
      v[loss_to_byte(i)] = RAM_REDUCED_LOSS;
    for (int i = DRAW_RULE + 1; i < n; i++)
      v[loss_to_byte(i)] = RAM_REDUCED_BLOSS;
  } else {
    int first_win = reduce_cnt_win[epoch - 1] - 250;
    int first_loss = reduce_cnt_loss[epoch - 1] + 4;

    v[RAM_UNRESOLVED] = RAM_UNRESOLVED;
    v[RAM_REDUCED_LOSS] = RAM_REDUCED_LOSS;
    v[RAM_REDUCED_CAPT_BLOSS] = RAM_REDUCED_CAPT_BLOSS;
    v[RAM_REDUCED_BLOSS] = RAM_REDUCED_BLOSS;
    v[RAM_CAPT_DRAW] = RAM_CAPT_DRAW;
    v[RAM_REDUCED_CWIN] = RAM_REDUCED_CWIN;
    v[RAM_REDUCED_CAPT_CWIN] = RAM_REDUCED_CAPT_CWIN;
    v[RAM_REDUCED_WIN] = RAM_REDUCED_WIN;
    v[RAM_CAPT_WIN] = RAM_CAPT_WIN;
    v[RAM_ILLEGAL] = RAM_ILLEGAL;

    for (int i = first_loss; i <= DRAW_RULE; i++)
      v[loss_to_byte(i)] = RAM_REDUCED_LOSS;
    for (int i = max(first_loss, DRAW_RULE + 1); i < n; i++)
      v[loss_to_byte(i)] = RAM_REDUCED_BLOSS;
    for (int i = first_win; i <= DRAW_RULE && i < n - 1; i++)
      v[win_to_byte(i)] = RAM_REDUCED_WIN;
    for (int i = max(first_win, DRAW_RULE + 1); i < n - 1; i++)
      v[win_to_byte(i)] = RAM_REDUCED_CWIN;
  }

  v[win_to_byte(n - 1)] = next_reduce_cnt_win - (n - 1);
  v[win_to_byte(n)] = next_reduce_cnt_win - n;
  v[loss_to_byte(n)] = n - next_reduce_cnt_loss;
  v[loss_to_byte(n + 1)] = n + 1 - next_reduce_cnt_loss;
}

void reduce_tables(int n)
{
  collect_stats_before_reduce(WHITE, n);
  collect_stats_before_reduce(BLACK, n);

  save_table(g_table[WHITE], WHITE, n);
  if (!symmetric)
    save_table(g_table[BLACK], BLACK, n);

  uint8_t v[256];
  init_reduce_map(v, n);
  transform_v = v;
  run_threaded(transform, &work_g_static);

  reduce_cnt_win[epoch] = n + 249;
  reduce_cnt_loss[epoch] = n - 4;
  epoch++;
}

#if 0
static void store_table(uint8_t *table, int stm)
{
  FILE *F;
  char name[64];

  sprintf(name, "%s.%c", tablename, "wb"[stm]);

  if (!(F = fopen(name, "wb"))) {
    fprintf(stderr, "Could not open %s for writing.\n", name);
    exit(EXIT_FAILURE);
  }

  write_data(F, table, 0, size);

  fclose(F);
}
#endif

void unlink_table(int stm)
{
  char name[64];

  sprintf(name, "%s.%c", g_tablename, "wb"[stm]);
  unlink(name);
}

void unlink_saves(int stm)
{
  char name[64];

  for (int k = 0; k < epoch; k++) {
    sprintf(name, "%s.%c.%d", g_tablename, "wb"[stm], k);
    unlink(name);
  }
}

#define T u8
#define MAX 256
#include "reduce_tmpl.c"
#undef MAX
#undef T

#define T u16
#define MAX MAX_STATS
#include "reduce_tmpl.c"
#undef MAX
#undef T
