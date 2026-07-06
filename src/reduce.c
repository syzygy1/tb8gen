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

static void save_table(uint8_t *table, int stm, int n)
{
  FILE *F;
  char name[64];
  uint8_t v[256];

  sprintf(name, "%s.%c.%d", g_tablename, "wb"[stm], epoch);

  if (!(F = fopen(name, "wb"))) {
    fprintf(stderr, "Could not open %s for writing.\n", name);
    exit(EXIT_FAILURE);
  }

  memset(&v, 0, sizeof v);

  if (epoch == 0) {
    v[1] = 0;
    for (int i = 1; i < n; i++)
      v[loss_to_byte(i)] = 256 - i;
    for (int i = 1; i < n; i++)
      v[win_to_byte(i)] = i;
  }

  write_data_transform_u8(F, table, table_size, v);

  fclose(F);
}

void reduce_tables(int n, int epoch)
{
  collect_stats(WHITE);
  collect_stats(BLACK);

  save_table(g_table[WHITE], WHITE, n);
  if (!symmetric)
    save_table(g_table[BLACK], BLACK, n);

  uint8_t v[256] = { 0 };

  if (epoch == 0) {
    v[0] = 0; // unresolved
    v[255] = 255; // illegal
    v[254] = 254; // capt_win
    v[253] = 253; // pawn_win -> reduced win
    for (int i = 1; i <= DRAW_RULE; i++)
      v[win_to_byte(i)] = 253;
    v[253 - DRAW_RULE - 1] = 252; // capt_cwin -> reduced capt_cwin
    v[253 - DRAW_RULE - 2] = 251; // pawn_cwin -> redcuded
    for (int i = DRAW_RULE + 1; i < n; i++)
      v[win_to_byte(i)] = 251;

    for(int i = 0; i <= DRAW_RULE; i++)
      v[loss_to_byte(i)] = 1; // loss -> reduced loss
    for (int i = DRAW_RULE + 1; i < n; i++)
      v[loss_to_byte(i)] = 2; // bloss -> reduced bloss
  } else {
    // to be done
  }


  transform_v = v;
  run_threaded(transform, &work_g_static);

  epoch++;
  // update reduce_cnt_win[], reduce_cnt_loss[]
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
