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

struct MergeInfo mi[2];

int epoch;
int reduce_cnt_win[64];
int reduce_cnt_loss[64];

static uint8_t *transform_table;
static uint8_t *transform_v8;

void reduce_worker(struct ThreadData *thread)
{
  uint8_t *restrict table = transform_table;
  uint8_t *v = transform_v8;

  for (uint64_t idx = thread->begin, end = thread->end; idx < end; idx++)
    table[idx] = v[table[idx]];
}

static void init_reduce_map(uint8_t v[256], int n);

static void init_save_map(uint8_t v[256], int n)
{
  memset(v, 0, 256);
  if (epoch == 0) {
    // L0 - L100
    for (int i = 0; i <= DRAW_RULE; i++)
      v[RAM_LOSS_IN_0 + i] = 255 - i;
    // L101 - L(n-1)
    for (int i = DRAW_RULE + 1; i <= n - 1; i++)
      v[RAM_LOSS_IN_0 + i + 1] = 255 - i;
    // W1 - W100
    for (int i = 1; i <= DRAW_RULE; i++)
      v[RAM_PAWN_WIN - i] = i;
    // W101 - W(n-2)
    for (int i = DRAW_RULE + 1; i <= n - 2; i++)
      v[RAM_PAWN_CWIN - (i - DRAW_RULE)] = i;
  } else {
    int first_loss = reduce_cnt_loss[epoch - 1] + 4;
    // L(first_loss) - L(n-1)
    for (int i = 0; i <= n - 1 - first_loss; i++)
      v[4 + i] = 255 - i;
    int first_win = reduce_cnt_win[epoch - 1] - 250;
    // W(first_win) - W(n-2)
    for (int i = 0; i <= n - 2 - first_win; i++)
      v[250 - i] = 1 + i;
  }
}

static void save_table(uint8_t *table, int stm, int n)
{
  char name[64];
  uint8_t v[256];

  sprintf(name, "%s.%c.%d", g_tablename, "wb"[stm], epoch);

  FILE *F = file_open_write(name);
  init_save_map(v, n);
  write_data_transform_u8(F, table, table_size, v);
  fclose(F);
  file_rename(name);
}

static int ram_byte_to_stat(uint8_t b)
{
  if (epoch != 0) {
    switch (b) {
    case RAM_REDUCED_LOSS:
      return loss_to_stat(0);
    case RAM_REDUCED_CAPT_BLOSS:
    case RAM_REDUCED_BLOSS:
      return loss_to_stat(DRAW_RULE + 1);
    case RAM_REDUCED_CWIN:
      return DRAW_RULE + 4;
    case RAM_REDUCED_CAPT_CWIN:
      return DRAW_RULE + 3;
    case RAM_REDUCED_WIN:
      return 2;
    default:
      break;
    }
  }

  switch (b) {
  case RAM_UNRESOLVED:
  case RAM_PAWN_DRAW:
    return MAX_STATS / 2 + 2;
  case RAM_ILLEGAL:
    return epoch == 0 ? 0 : -1;
  case RAM_CAPT_WIN:
    return epoch == 0 ? 1 : -1;
  case RAM_PAWN_WIN:
    return epoch == 0 ? 2 : -1;
  case RAM_CAPT_CWIN:
    return epoch == 0 ? DRAW_RULE + 3 : -1;
  case RAM_PAWN_CWIN:
    return epoch == 0 ? DRAW_RULE + 4 : -1;
  case RAM_CAPT_BLOSS:
    return epoch == 0 ? loss_to_stat(DRAW_RULE + 1) : -1;
  case RAM_CAPT_DRAW:
    return MAX_STATS / 2 + 1;
  default:
    break;
  }

  if (epoch == 0) {
    if (b >= RAM_LOSS_IN_0 && b < RAM_PAWN_DRAW) {
      int n = b - RAM_LOSS_IN_0 - (b > RAM_CAPT_BLOSS);
      return n < MAX_STATS / 2 - 3 ? loss_to_stat(n) : -1;
    }
    if (b > RAM_CAPT_DRAW && b < RAM_CAPT_CWIN) {
      int n = 253 - b - 2;
      return n > 0 && n < MAX_STATS / 2 - 3 ? win_to_stat(n) : -1;
    }
    if (b > RAM_CAPT_CWIN && b < RAM_PAWN_WIN) {
      int n = 253 - b;
      return n > 0 && n < MAX_STATS / 2 - 3 ? win_to_stat(n) : -1;
    }
    return -1;
  }

  if (b > RAM_REDUCED_BLOSS && b < RAM_CAPT_DRAW) {
    int n = b + reduce_cnt_loss[epoch - 1];
    return n < MAX_STATS / 2 - 3 ? loss_to_stat(n) : -1;
  }
  if (b > RAM_CAPT_DRAW && b < RAM_REDUCED_CWIN) {
    int n = reduce_cnt_win[epoch - 1] - b;
    return n > 0 && n < MAX_STATS / 2 - 3 ? win_to_stat(n) : -1;
  }

  return -1;
}

// For epoch == 0:
// - Losses L0 - L100 collapse to RAM_REDUCED_LOSS (byte 1).
// - Losses L101 - L(n-1) collapse to RAM_REDUCED_BLOSS (byte 3).
// - Wins W1 - W100 collapse to RAM_REDUCED_WIN (byte 253).
// - Wins W101 - W(n-2) collapse to RAM_REDUCED_CWIN (byte 251).
//
// For epoch > 0, earlier reductions have already compacted the old ranges.
// The code only collapses the newly active range, starting at:
// - first_loss = reduce_cnt_loss[epoch - 1] + 4
// - first_win = reduce_cnt_win[epoch - 1] - 250
//
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

  if (!g_only_generate) {
    save_table(g_table[WHITE], WHITE, n);
    if (!symmetric)
      save_table(g_table[BLACK], BLACK, n);
  }

  uint8_t v[256];
  init_reduce_map(v, n);
  transform_v8 = v;
  transform_table = g_table[WHITE];
  run_threaded(reduce_worker, &work_g_static);
  transform_table = g_table[BLACK];
  run_threaded(reduce_worker, &work_g_static);

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

// mi.v[]    :stats -> ram_byte / ram_u16
// mi.inv_v[]:ram_byte -> stats

static uint16_t stats_to_val[2][MAX_STATS];
static uint16_t val_to_stats[2][MAX_STATS];
bool dtz_wide[2];

void create_stats_to_val(int stm)
{
  // Determine whether we need to store wins, losses or both in the merged
  // in-RAM table.
  bool wins, losses;
  if (one_sided)
    wins = losses = stm == one_sided_stm;
  else {
    wins = wins_only;
    losses = !wins_only;
  }

  uint64_t *stats = g_stats[stm];
  uint16_t *s_to_v = stats_to_val[stm];

  memset(s_to_v, 0, sizeof stats_to_val[stm]);
  int n = 0;
  if (wins) {
    for (int i = 3; i <= DRAW_RULE + 2; i++) {
      n += (stats[i] != 0);
      s_to_v[i] = n;
    }
    for (int i = DRAW_RULE + 5; i < MAX_STATS / 2 + 1; i++) {
      n += (stats[i] != 0);
      s_to_v[i] = n;
    }
  }
  if (losses) {
    for (int i = MAX_STATS / 2 - 4; i >= 0; i--) {
      n += (stats[MAX_STATS - 1 - i] != 0);
      s_to_v[MAX_STATS - 1 - i] = n;
    }
  }

  memset(val_to_stats[stm], 0, sizeof val_to_stats[stm]);

  for (int i = 0, j = -1; i < MAX_STATS; i++)
    if (s_to_v[i] != j)
      val_to_stats[stm][j = s_to_v[i]] = i;
  dtz_wide[stm] = n >= 256;
}

static void *transform_src, *transform_dst;
static uint16_t *transform_v;

static void transform_u8_worker(struct ThreadData *thread)
{
  uint8_t *restrict table = transform_src;
  uint16_t *restrict v = transform_v;

  for (uint64_t idx = thread->begin, end = thread->end; idx < end; idx++)
    table[idx] = v[table[idx]];
}

static void transform_u16_worker(struct ThreadData *thread)
{
  uint8_t *restrict src = transform_src;
  uint16_t *restrict dst = transform_dst;
  uint16_t *restrict v = transform_v;

  for (uint64_t idx = thread->begin, end = thread->end; idx < end; idx++)
    dst[idx] = v[src[idx]];
}

static void transform_to_val(int stm, uint8_t *src, void *dst, uint16_t *v)
{
  transform_src = src;
  transform_dst = dst;
  transform_v = v;

  if (dtz_wide[stm]) {
    run_threaded(transform_u16_worker, &work_g_static);
  } else {
    assert(src == dst);
    run_threaded(transform_u8_worker, &work_g_static);
  }
}

void reconstruct_pass(int stm, void *table, uint16_t v[256], int k)
{
  char name[64];
  sprintf(name, "%s.%c.%d", g_tablename, "wb"[stm], k);

  FILE *F = file_open_read(name);

  if (dtz_wide[stm]) {
    read_data_transform_or_u8_to_u16(F, table, table_size, v);
  } else {
    uint8_t v8[256];
    for (int i = 0; i < 256; i++)
      v8[i] = v[i];
    read_data_transform_or_u8(F, table, table_size, v8);
  }
}

// If dtz_wide[stm], then dst must point to an u16 buffer.
// If !dtz_wide[stm], then we must have dst == src.
void reconstruct_table(int stm, void *src, void *dst)
{
  uint16_t v[256] = { 0 };  // byte -> stats[]

  if (epoch == 0) {
    // byte -> stats[]
    for (int i = 0; i <= DRAW_RULE; i++)
      v[RAM_LOSS_IN_0 + i] = MAX_STATS - 1 - i;
    for (int i = DRAW_RULE + 1; RAM_LOSS_IN_0 + i + 1 < RAM_PAWN_DRAW; i++)
      v[RAM_LOSS_IN_0 + i + 1] = MAX_STATS - 1 - i;
    for (int i = 1; i <= DRAW_RULE; i++)
      v[RAM_PAWN_WIN - i] = i + 2;
    for (int i = DRAW_RULE + 1; RAM_PAWN_CWIN - (i - DRAW_RULE) > RAM_CAPT_DRAW;
        i++)
      v[RAM_PAWN_CWIN - (i - DRAW_RULE)] = i + 4;

    // stats -> reduced values
    for (int i = 0; i < 256; i++)
      v[i] = stats_to_val[stm][v[i]];

    transform_to_val(stm, src, dst, v);
    return;
  }

  // epoch > 0
  int first_loss = reduce_cnt_loss[epoch - 1] + 4;
  // L(first_loss) - L(n-1)
  for (int i = 0; 4 + i < RAM_PAWN_DRAW; i++)
    v[4 + i] = MAX_STATS - 1 - (first_loss + i);
  int first_win = reduce_cnt_win[epoch - 1] - 250;
  // W(first_win) - W(n-2)
  for (int i = 0; 250 - i > RAM_CAPT_DRAW; i++)
    v[250 - i] = 4 + (first_win + i);

  for (int i = 0; i < 256; i++)
    v[i] = stats_to_val[stm][v[i]];

  transform_to_val(stm, src, dst, v);

  // k == 0
  memset(v, 0, sizeof v);
  for (int i = 0; 255 - i >= 128; i++)
    v[255 - i] = MAX_STATS - 1 - i;
  for (int i = 1; i <= DRAW_RULE; i++)
    v[1 + i] = 2 + i;
  for (int i = DRAW_RULE + 1; 1 + i < 128; i++)
    v[1 + i] = 4 + i;
  reconstruct_pass(stm, dst, v, 0);

  for (int k = 1 ; k < epoch; k++) {
    memset(v, 0, sizeof v);
    int first_loss = reduce_cnt_loss[k - 1] + 4;
    for (int i = 0; 255 - i >= 128; i++)
      v[255 - i] = MAX_STATS - 1 - (first_loss + i);
    int first_win = reduce_cnt_win[k - 1] - 250;
    for (int i = 0; 1 + i < 128; i++)
      v[1 + i] = 4 + (first_win + i);
    reconstruct_pass(stm, dst, v, k);
  }
}



// Instead of doing this, first look at g_stats[] to determine all
// the wins and/or losses that need a "slot", and create a mapping
// between the stats[] values that need to survive (win values if
// wins == true, loss values if losses == true) and a "minimal" range,
// similar to how this is done in
// sort_values() is called on stats[], so this connects.
//
// I guess we could already round values > DRAW_RULE, but maybe don't
// do that immediately.
//
// Depending on the minimal range, we need u8 or u16.
// If u16, we need to persist the other side to disk if the other side
// still must be compressed (unless we know that RAM is plenty).
//
// Probably not here but in reconstruct_table():
// Then create a mapping to transform the RAM table, via their corresponding
// stats value, into this minimal range.
// We do the same for each of the saved epochs.
// This should be easy since we know exactly to which stats[] they map.

#if 0
  mi.wide = tot_vals > 256;

  if (!mi.wide) {
    // One byte suffices.

    int n = init_merge_value_map_u8(stats, stm);
    if (n > 255) {
      fprintf(stderr, "Internal error.\n");
      exit(EXIT_FAILURE);
    }

    merge_table = alloc_huge(sizeof(u8) * table_size);
    if (!merge_table)
      out_of_mem();

    reconstruct_table_u8(stm);

    free(merge_table);

  } else {

    init_merge_value_map_u16(stats, stm);

    merge_table = alloc_huge(sizeof(u16) * table_size);
    if (!merge_table)
      out_of_mem();

    reconstruct_table_u16(stm);

    free(merge_table);

#endif
