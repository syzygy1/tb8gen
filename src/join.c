/*
  Copyright (c) 2026 Ronald de Man

  This file is distributed under the terms of the GNU GPL, version 2.
*/

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "compress.h"
#include "defs.h"
#include "merge.h"
#include "movegen.h"
#include "permute.h"
#include "permute462.h"
#include "tb8gen.h"
#include "threads.h"
#include "types.h"
#include "util.h"

static constexpr int minfreq = 8;

static void *join_table, *tb_table;
static uint8_t *join_pt, *join_pcs;
bool join_wide, tb_wide;
struct DtzMap dtzmap;

static int sort_list(uint64_t *freq, uint16_t *map, uint16_t *inv_map)
{
  int num;

  for (int i = 0; i < MAX_VAL; i++)
    map[i] = inv_map[i] = 0;

  num = 0;
  for (int i = 0; i < MAX_VAL; i++)
    if (freq[i])
      map[num++] = i;

  for (int i = 0; i < num; i++)
    for (int j = i + 1; j < num; j++)
      if (freq[map[i]] < freq[map[j]])
        Swap(map[i], map[j]);

  for (int i = 0; i < num; i++)
    inv_map[map[i]] = i;

  return num;
}

static void sort_values(int stm, uint64_t *stats, struct DtzMap *dtzmap)
{
  uint64_t freq[4][MAX_VAL];
  uint16_t (*map)[MAX_VAL] = dtzmap->map;
  uint16_t (*inv_map)[MAX_VAL] = dtzmap->inv_map;

  dtzmap->stm = stm;

  for (int j = 0; j < 4; j++)
    for (int i = 0; i < MAX_VAL; i++)
      freq[j][i] = 0;

  for (int i = 0; i < 4; i++)
    dtzmap->num[i] = 0;

  if (one_sided || wins_only) {
    for (int i = 1; i <= DRAW_RULE; i++)
      freq[0][i - 1] += stats[1 + i];
    dtzmap->num[0] = sort_list(freq[0], map[0], inv_map[0]);

    for (int i = DRAW_RULE + 1; i < MAX_STATS / 2 - 2; i++)
      freq[2][(i - DRAW_RULE - 1) / 2] += stats[2 + i];
    dtzmap->num[2] = sort_list(freq[2], map[2], inv_map[2]);
  }

  if (one_sided || !wins_only) {
    freq[1][0] = stats[MAX_STATS - 1];
    for (int i = 1; i <= DRAW_RULE; i++)
      freq[1][i - 1] += stats[MAX_STATS - 1 - i];
    dtzmap->num[1] = sort_list(freq[1], map[1], inv_map[1]);

    for (int i = DRAW_RULE + 1; i < MAX_STATS / 2 - 2; i++)
      freq[3][(i - DRAW_RULE - 1) / 2] += stats[MAX_STATS - 1 - i];
    dtzmap->num[3] = sort_list(freq[3], map[3], inv_map[3]);
  }

  int num = 1;
  for (int i = 0; i < 4; i++)
    num = max(num, dtzmap->num[i]);
  dtzmap->max_num = num;
  dtzmap->wide = num >= 240;

  static uint64_t tot_pos[] = {
    518,
    31332,
    1911252,
    114675120,
    6765832080,
    392418260640
  };
  uint64_t tot = tot_pos[g_pos.num - 2] / 7000ULL;

  int i;
  for (i = 0; i < num; i++) {
    uint64_t f = 0;
    for (int j = 0; j < 4; j++)
      if (i < dtzmap->num[j])
        f += freq[j][map[j][i]];
    if (f < tot) break;
  }
  dtzmap->high_freq_max = i;
}

static void prepare_dtz_map(uint16_t *v, struct DtzMap *map)
{
  uint16_t (*inv_map)[MAX_VAL] = map->inv_map;
  int num = map->max_num;

  // The value 'num' will be compressed as don't care.
  for (int i = 0; i < MAX_STATS; i++)
    v[i] = num;

  if (one_sided || wins_only) {
    for (int i = 1; i <= DRAW_RULE; i++)
      v[1 + i] = inv_map[0][i - 1];
    for (int i = DRAW_RULE + 1; i < MAX_STATS / 2 - 2; i++)
      v[2 + i] = inv_map[2][(i - DRAW_RULE - 1) / 2];
  }
  if (one_sided || !wins_only) {
    v[MAX_STATS - 1] = inv_map[1][0];
    for (int i = 1; i <= DRAW_RULE; i++)
      v[MAX_STATS - 1 - i] = inv_map[1][i - 1];
    for (int i = DRAW_RULE + 1; i < MAX_STATS / 2 - 2; i++)
      v[MAX_STATS - 1 - i] = inv_map[3][(i - DRAW_RULE - 1) / 2];
  }
}

static void read_merge_info(int stm)
{
  char str[128];
  sprintf(str, "merge_info.%c", "wb"[stm]);
  FILE *F = fopen(str, "rb");
  if (!F) {
    fprintf(stderr, "Could not open %s.\n", str);
    exit(EXIT_FAILURE);
  }
  file_read(&mi, sizeof(mi), F);
  fclose(F);
}

// Join everything in one table per side.
static void join_wdl(int stm, struct tb_handle *G)
{
  char str[128];
  uint8_t *table = join_table;

  // Broken positions with adjacent kings are mapped to the index one
  // beyond the end of the table. We need to set this value to don't care.
  table[462 * kslice_size] = 8;

  for (int s = 0; s < 462; s++) {
    create_name(str, s, stm, "merged/wdl", -1);
    FILE *F = fopen(str, "rb");
    if (!F) {
      fprintf(stderr, "Could not open %s.\n", str);
      exit(EXIT_FAILURE);
    }
    read_data(F, table + s * kslice_size, kslice_size);
    fclose(F);
  }

  read_merge_info(stm);
  compress_init_wdl(mi.v_wdl);

  uint8_t best[MAX_PIECES];
  printf("Find optimal premutation for %ctm/wdl.\n", "wb"[stm]);
  permute_piece_wdl(tb_table, join_pcs, join_pt, table, best);
  printf("Compressing date for %ctm/wdl.\n", "wb"[stm]);
  compress_tb(G, -1, tb_table, tb_size, best, minfreq, false);
}

static void join_dtz(int stm, struct tb_handle *G)
{
  uint16_t v[MAX_STATS];
  char str[128];

  read_merge_info(stm);

  sort_values(stm, g_stats[stm], &dtzmap);

  if (dtzmap.wide != join_wide) {
    join_wide = dtzmap.wide;
    free(join_table);
    join_table = alloc_huge((462 * kslice_size + 1) * (1 + join_wide));
    if (!join_table)
      out_of_mem();
  }

  if (dtzmap.wide != tb_wide) {
    tb_wide = dtzmap.wide;
    free(tb_table);
    tb_table = alloc_huge((tb_size + 1) * (1 + tb_wide));
    if (!tb_table)
      out_of_mem();
  }

  prepare_dtz_map(v, &dtzmap);

  compress_alloc_dtz(tb_wide);

  if (!mi.wide) {

    if (tb_wide) {
      fprintf(stderr, "Internal error.\n");
      exit(EXIT_FAILURE);
    }

    uint8_t w[256];
    for (int i = 255; i >= 0; i--)
      w[i] = v[mi.v_inv_u8[i]];

    uint8_t *table = join_table;
    table[462 * kslice_size] = dtzmap.max_num;

    for (int s = 0; s < 462; s++) {
      create_name(str, s, stm, "merged/dtz", -1);
      FILE *F = fopen(str, "rb");
      if (!F) {
        fprintf(stderr, "Could not open %s.\n", str);
        exit(EXIT_FAILURE);
      }
      read_data_transform_u8(F, table + s * kslice_size, kslice_size, w);
      fclose(F);
    }

  } else if (!tb_wide) {

    uint8_t w[MAX_STATS];
    for (int i = 0; i < MAX_STATS; i++)
      w[i] = v[i];

    uint8_t *table = join_table;
    table[462 * kslice_size] = dtzmap.max_num;

    for (int s = 0; s < 462; s++) {
      create_name(str, s, stm, "merged/dtz", -1);
      FILE *F = fopen(str, "rb");
      if (!F) {
        fprintf(stderr, "Could not open %s.\n", str);
        exit(EXIT_FAILURE);
      }
      read_data_transform_to_u8_u16(F, table + s * kslice_size,
          kslice_size * 2, w);
      fclose(F);
    }

  } else {

    uint16_t *table = join_table;
    table[462 * kslice_size] = dtzmap.max_num;

    for (int s = 0; s < 462; s++) {
      create_name(str, s, stm, "merged/dtz", -1);
      FILE *F = fopen(str, "rb");
      if (!F) {
        fprintf(stderr, "could not open %s.\n", str);
        exit(EXIT_FAILURE);
      }
      read_data_transform_u16(F, table + s * kslice_size, kslice_size * 2, v);
      fclose(F);
    }

  }

  compress_init_dtz(&dtzmap, tb_wide);

  uint8_t best[MAX_PIECES];
  printf("Find optimal premutation for %ctm/dtz.\n", "wb"[stm]);
  permute_piece_dtz(tb_table, join_pcs, join_pt, join_table, best, tb_wide);
  printf("Compressing date for %ctm/dtz.\n", "wb"[stm]);
  compress_tb(G, -1, tb_table, tb_size, best, minfreq, tb_wide);

  compress_free_dtz();
}

void join_slices(uint8_t *pcs, uint8_t *pt)
{
  join_pcs = pcs;
  join_pt = pt;

  init_permute_piece(pcs, pt);

  join_table = alloc_huge(462 * kslice_size + 1);
  tb_table = alloc_huge(tb_size + 1);
  if (!join_table || !tb_table)
    out_of_mem();
  join_wide = tb_wide = false;

  compress_alloc_wdl();
  struct tb_handle *G = create_tb(g_tablename, WDL, 6);

  join_wdl(WHITE, G);
  if (!symmetric)
    join_wdl(BLACK, G);

  merge_tb(G);
  compress_free_wdl();

  G = create_tb(g_tablename, DTZ, 10);

  if (!one_sided || one_sided_stm == WHITE)
    join_dtz(WHITE, G);

  if (   (one_sided && one_sided_stm == BLACK)
      || (!one_sided && !symmetric))
    join_dtz(BLACK, G);

  free(join_table);
  free(tb_table);

  merge_tb(G);
}

// Join everything in one table per side.
static void join_wdl_462(int stm)
{
  uint64_t stats[MAX_STATS];
  char str[128];
  uint8_t *table = join_table;

  create_dir(-1, stm, "wdl");

  for (int s = 0; s < 462; s++) {
    g_pos.sq[0] = KKSquare[s][0];
    g_pos.sq[1] = KKSquare[s][1];
    g_pos.stm = stm;

    create_name(str, s, stm, "merged/wdl", -1);
    FILE *F = fopen(str, "rb");
    if (!F) {
      fprintf(stderr, "Could not open %s.\n", str);
      exit(EXIT_FAILURE);
    }
    read_data(F, table, kslice_size);
    fclose(F);

    create_name(str, s, stm, "stats", -1);
    F = fopen(str, "rb");
    if (!F) {
      fprintf(stderr, "Could not open %s.\n", str);
      exit(EXIT_FAILURE);
    }
    read_data(F, stats, sizeof(stats));
    fclose(F);

    bool v_wdl[5] = { 0 };
    for (int i = 2; i <= DRAW_RULE + 1; i++)
      if (stats[i]) {
        v_wdl[4] = true;
        break;
      }
    for (int i = DRAW_RULE + 3; i < MAX_STATS / 2; i++)
      if (stats[i]) {
        v_wdl[3] = true;
        break;
      }
    v_wdl[2] = (bool)stats[MAX_STATS / 2 + 1];
    for (int i = DRAW_RULE + 1; i < MAX_STATS / 2 - 2; i++)
      if (stats[i]) {
        v_wdl[1] = true;
        break;
      }
    for (int i = 0; i <= DRAW_RULE; i++)
      if (stats[i]) {
        v_wdl[0] = true;
        break;
      }

    compress_init_wdl(v_wdl);

    uint8_t best[MAX_SETS];
    printf("Find optimal premutation for %ctm/wdl, slice = %d.\n", "wb"[stm],
        s);
    permute_piece_462(tb_table, table, best, WDL, false);
    printf("Compressing date for %ctm/wdl, slice = %d.\n", "wb"[stm], s);
    compress_data_462(s, stm, WDL, tb_table, kslice_size, best, minfreq, false);
  }
}

static void join_dtz_462(int stm)
{
  uint16_t v[MAX_STATS];
  char str[128];

  read_merge_info(stm);

  sort_values(stm, g_stats[stm], &dtzmap);

  if (dtzmap.wide != join_wide) {
    join_wide = dtzmap.wide;
    free(join_table);
    join_table = alloc_huge(kslice_size * (1 + join_wide));
    if (!join_table)
      out_of_mem();
  }

  if (dtzmap.wide != tb_wide) {
    tb_wide = dtzmap.wide;
    free(tb_table);
    tb_table = alloc_huge((tb_size + 1) * (1 + tb_wide));
    if (!tb_table)
      out_of_mem();
  }

  prepare_dtz_map(v, &dtzmap);

  compress_alloc_dtz(tb_wide);

  if (!mi.wide) {

    if (tb_wide) {
      fprintf(stderr, "Internal error.\n");
      exit(EXIT_FAILURE);
    }

    uint8_t w[256];
    for (int i = 255; i >= 0; i--)
      w[i] = v[mi.v_inv_u8[i]];

    uint8_t *table = join_table;

    for (int s = 0; s < 462; s++) {
      create_name(str, s, stm, "merged/dtz", -1);
      FILE *F = fopen(str, "rb");
      if (!F) {
        fprintf(stderr, "Could not open %s.\n", str);
        exit(EXIT_FAILURE);
      }
      read_data_transform_u8(F, table, kslice_size, w);
      fclose(F);
    }

  } else if (!tb_wide) {

    uint8_t w[MAX_STATS];
    for (int i = 0; i < MAX_STATS; i++)
      w[i] = v[i];

    uint8_t *table = join_table;
    table[kslice_size] = dtzmap.max_num;

    for (int s = 0; s < 462; s++) {
      create_name(str, s, stm, "merged/dtz", -1);
      FILE *F = fopen(str, "rb");
      if (!F) {
        fprintf(stderr, "Could not open %s.\n", str);
        exit(EXIT_FAILURE);
      }
      read_data_transform_to_u8_u16(F, table, kslice_size * 2, w);
      fclose(F);
    }

  } else {

    uint16_t *table = join_table;
    table[kslice_size] = dtzmap.max_num;

    for (int s = 0; s < 462; s++) {
      create_name(str, s, stm, "merged/dtz", -1);
      FILE *F = fopen(str, "rb");
      if (!F) {
        fprintf(stderr, "could not open %s.\n", str);
        exit(EXIT_FAILURE);
      }
      read_data_transform_u16(F, table, kslice_size * 2, v);
      fclose(F);
    }

  }

  compress_init_dtz(&dtzmap, tb_wide); // we must somehow relay dtzmap!

  uint8_t best[MAX_PIECES];
  printf("Find optimal premutation for %ctm/dtz.\n", "wb"[stm]);
  permute_piece_462(tb_table, join_table, best, DTZ, tb_wide);
  printf("Compressing date for %ctm/dtz.\n", "wb"[stm]);
//  compress_tb_462(G, -1, tb_table, tb_size, best, minfreq, tb_wide);

  compress_free_dtz();
}

void join_slices_462(void)
{
  init_permute_piece_462();

  join_table = alloc_huge(kslice_size);
  tb_table = alloc_huge(tb_size + 1);
  if (!join_table || !tb_table)
    out_of_mem();
  join_wide = tb_wide = false;

  compress_alloc_wdl();

  join_wdl_462(WHITE);
  if (!symmetric)
    join_wdl_462(BLACK);

//  merge_tb_462(G);
  compress_free_wdl();

  if (!one_sided || one_sided_stm == WHITE)
    join_dtz_462(WHITE);

  if (   (one_sided && one_sided_stm == BLACK)
      || (!one_sided && !symmetric))
    join_dtz_462(BLACK);

  free(join_table);
  free(tb_table);

//  merge_tb_462(G);
}
