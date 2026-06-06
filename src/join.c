/*
  Copyright (c) 2026 Ronald de Man

  This file is distributed under the terms of the GNU GPL, version 2.
*/

#include <fcntl.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#include "checksum.h"
#include "compress.h"
#include "defs.h"
#include "index.h"
#include "kslice.h"
#include "merge.h"
#include "movegen.h"
#include "permute.h"
#include "permute10.h"
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

static char *create_final_name(int type)
{
  size_t len = strlen(g_output_dir) + strlen(g_tablename)
      + strlen(suffix[type]) + 2;
  char *name = malloc(len);
  if (!name)
    out_of_mem();
  snprintf(name, len, "%s/%s%s", g_output_dir, g_tablename, suffix[type]);
  return name;
}

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

static void sort_values(int stm, uint64_t *stats, struct DtzMap *dtzmap,
    uint64_t tot_pos)
{
  uint64_t freq[4][MAX_VAL] = { 0 };
  uint16_t (*map)[MAX_VAL] = dtzmap->map;
  uint16_t (*inv_map)[MAX_VAL] = dtzmap->inv_map;

  dtzmap->stm = stm;

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

  uint64_t tot = tot_pos / 7000ULL;

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

static void prepare_wdl_map(uint64_t *stats, bool *v, bool has_capt_bloss)
{
  for (int i = 0; i < 5; i++)
    v[i] = false;

  for (int i = 2; i <= DRAW_RULE + 1; i++)
    if (stats[i]) {
      v[4] = true;
      break;
    }
  for (int i = DRAW_RULE + 3; i < MAX_STATS / 2; i++)
    if (stats[i]) {
      v[3] = true;
      break;
    }
  v[2] = (bool)stats[MAX_STATS / 2 + 1];
  for (int i = DRAW_RULE + 1; i < MAX_STATS / 2 - 2; i++)
    if (stats[MAX_STATS - 1 - i]) {
      v[1] = true;
      break;
    }
  for (int i = 0; i <= DRAW_RULE; i++)
    if (stats[MAX_STATS - 1 - i]) {
      v[0] = true;
      break;
    }

  bool dc[4] = {
    has_capt_bloss, stats[MAX_STATS / 2], stats[DRAW_RULE + 2], true
  };

  int i, j;
  for (i = 0; i < 4; i++)
    if (dc[i]) break;
  for (j = 0; j < 5; j++)
    if (v[j]) break;
  if (j > i + 1)
    v[0] = true;
}

static bool read_wdl_slice(const char *name, uint8_t *table, uint64_t size)
{
  uint8_t has_capt_bloss;
  FILE *F = file_open_read(name);
  file_read(&has_capt_bloss, 1, F);
  read_data(F, table, size);
  fclose(F);
  return has_capt_bloss != 0;
}

static void read_merge_info(int stm)
{
  char str[64];
  sprintf(str, "merge_info.%c", "wb"[stm]);
  FILE *F = file_open_read(str);
  read_data(F, &mi, sizeof mi);
  fclose(F);
}

void delete_slices_dir(int stm, const char *name)
{
  if (!g_cleanup) return;

  for (int s = 0; s < 462; s++)
    kslice_delete(s, stm, name, -1);
  char str[64];
  sprintf(str, "%s/%c", name, "wb"[stm]);
  rmdir(str);
}

void final_cleanup(void)
{
  if (!g_cleanup) return;

  delete_slices_dir(WHITE, "stats");
  delete_slices_dir(BLACK, "stats");
  rmdir("stats");
  unlink("done_wdl");
  unlink("generate_info");
  unlink("merge_info.w");
  unlink("merge_info.b");
  unlink("maxfens");
}

// Join everything in one table per side.
static void join_wdl(int stm, struct tb_handle *G)
{
  char str[64];
  uint8_t *table = join_table;
  bool has_capt_bloss = false;

  // Broken positions with adjacent kings are mapped to the index one
  // beyond the end of the table. We need to set this value to don't care.
  table[462 * kslice_size] = 8;

  for (int s = 0; s < 462; s++) {
    create_name(str, s, stm, "merged/wdl", -1);
    if (read_wdl_slice(str, table + s * kslice_size, kslice_sizes[s >= 441]))
      has_capt_bloss = true;
  }

  read_merge_info(stm);

  bool v_wdl[5];
  prepare_wdl_map(g_stats[stm], v_wdl, has_capt_bloss);
  compress_init_wdl(mi.v_wdl);

  uint8_t best[MAX_PIECES];
  printf("Find optimal permutation for %ctm/wdl.\n", "wb"[stm]);
  permute_piece_wdl(tb_table, join_pcs, join_pt, table, best);
  printf("Compressing data for %ctm/wdl.\n", "wb"[stm]);
  compress_tb(G, -1, tb_table, tb_size, best, minfreq, false);
}

static void join_dtz(int stm, struct tb_handle *G)
{
  uint16_t v[MAX_STATS];
  char str[64];

  read_merge_info(stm);

  sort_values(stm, g_stats[stm], &dtzmap,
      441 * kslice_sizes[0] + 21 * kslice_sizes[1]);

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
    for (int i = 0; i < 256; i++)
      w[i] = v[mi.v_inv_u8[i]];

    uint8_t *table = join_table;
    table[462 * kslice_size] = dtzmap.max_num;

    for (int s = 0; s < 462; s++) {
      create_name(str, s, stm, "merged/dtz", -1);
      FILE *F = file_open_read(str);
      read_data_transform_u8(F, table + s * kslice_size,
          kslice_sizes[s >= 441], w);
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
      FILE *F = file_open_read(str);
      read_data_transform_to_u8_u16(F, table + s * kslice_size,
          kslice_sizes[s >= 441] * 2, w);
      fclose(F);
    }

  } else {

    uint16_t *table = join_table;
    table[462 * kslice_size] = dtzmap.max_num;

    for (int s = 0; s < 462; s++) {
      create_name(str, s, stm, "merged/dtz", -1);
      FILE *F = file_open_read(str);
      read_data_transform_u16(F, table + s * kslice_size,
          kslice_sizes[s >= 441] * 2, v);
      fclose(F);
    }

  }

  compress_init_dtz(&dtzmap, tb_wide);

  uint8_t best[MAX_PIECES];
  printf("Find optimal permutation for %ctm/dtz.\n", "wb"[stm]);
  permute_piece_dtz(tb_table, join_pcs, join_pt, join_table, best, tb_wide);
  printf("Compressing data for %ctm/dtz.\n", "wb"[stm]);
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

  if (!file_exists("done_wdl")) {
    compress_alloc_wdl();
    struct tb_handle *G = create_tb(g_tablename, WDL, 6);

    join_wdl(WHITE, G);
    if (!symmetric)
      join_wdl(BLACK, G);

    merge_tb(G);

    compress_free_wdl();

    create_empty("done_wdl");

    if (g_cleanup) {
      delete_slices_dir(WHITE, "merged/wdl");
      if (!symmetric)
        delete_slices_dir(BLACK, "merged/wdl");
      rmdir("merged/wdl");
    }
  }

  struct tb_handle *G = create_tb(g_tablename, DTZ, 10);

  if (!one_sided || one_sided_stm == WHITE)
    join_dtz(WHITE, G);

  if (   (one_sided && one_sided_stm == BLACK)
      || (!one_sided && !symmetric))
    join_dtz(BLACK, G);

  free(join_table);
  free(tb_table);

  merge_tb(G);

  if (!g_cleanup) return;

  if (!one_sided || one_sided_stm == WHITE)
    delete_slices_dir(WHITE, "merged/dtz");
  if (   (one_sided && one_sided_stm == BLACK)
      || (!one_sided && !symmetric))
    delete_slices_dir(BLACK, "merged/dtz");
  rmdir("merged/dtz");
  rmdir("merged");

  final_cleanup();
}

// Compress each of the 462 K-slices per side individually.
static void join_wdl_462(int stm)
{
  uint64_t stats[MAX_STATS];
  char str[64];
  uint8_t *table = join_table;

  create_dir(-1, stm, "wdl");

  for (int s = 0; s < 462; s++) {
    char name[64];
    create_name(name, s, stm, "wdl", -1);
    if (file_exists(name))
      continue;

    g_pos.sq[0] = KKSquare[s][0];
    g_pos.sq[1] = KKSquare[s][1];
    g_pos.stm = stm;

    create_name(str, s, stm, "merged/wdl", -1);
    bool has_capt_bloss = read_wdl_slice(str, table, kslice_sizes[s >= 441]);

    create_name(str, s, stm, "stats", -1);
    FILE *F = file_open_read(str);
    read_data(F, stats, sizeof stats);
    fclose(F);

    bool v_wdl[5];
    prepare_wdl_map(stats, v_wdl, has_capt_bloss);
    compress_init_wdl(v_wdl);

    uint8_t best[MAX_SETS];
    printf("Find optimal permutation for %ctm/wdl, slice = %d.\n", "wb"[stm],
        s);
    permute_piece_462(tb_table, table, best, WDL, false);
    printf("Compressing data for %ctm/wdl, slice = %d.\n", "wb"[stm], s);
    compress_data_slice(name, stm, WDL, tb_table, kslice_sizes[s >= 441], best,
        minfreq, false, false);

    if (!g_cleanup) continue;

    kslice_delete(s, stm, "merged/wdl", -1);
  }
}

static void join_dtz_462(int stm)
{
  uint64_t stats[MAX_STATS];
  uint16_t v[MAX_STATS];
  char str[64];

  create_dir(-1, stm, "dtz");

  read_merge_info(stm);

  sort_values(stm, g_stats[stm], &dtzmap,
      441 * kslice_sizes[0] + 21 * kslice_sizes[1]);

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
    tb_table = alloc_huge((kslice_size + 1) * (1 + tb_wide));
    if (!tb_table)
      out_of_mem();
  }

  if (!mi.wide && tb_wide) {
    fprintf(stderr, "Internal error.\n");
    exit(EXIT_FAILURE);
  }

  compress_alloc_dtz(tb_wide);

  for (int s = 0; s < 462; s++) {
    char name[64];
    create_name(name, s, stm, "dtz", -1);
    if (file_exists(name))
      continue;

    g_pos.sq[0] = KKSquare[s][0];
    g_pos.sq[1] = KKSquare[s][1];
    g_pos.stm = stm;

    create_name(str, s, stm, "stats", -1);
    FILE *F = file_open_read(str);
    read_data(F, stats, sizeof stats);
    fclose(F);

    sort_values(stm, stats, &dtzmap, kslice_sizes[s >= 441]);
    prepare_dtz_map(v, &dtzmap);

    if (!mi.wide) {

      uint8_t w[256];
      for (int i = 0; i < 256; i++)
        w[i] = v[mi.v_inv_u8[i]];

      uint8_t *table = join_table;

      create_name(str, s, stm, "merged/dtz", -1);
      FILE *F = file_open_read(str);
      read_data_transform_u8(F, table, kslice_sizes[s >= 441], w);
      fclose(F);

    } else if (!tb_wide) {

      uint8_t w[MAX_STATS];
      for (int i = 0; i < MAX_STATS; i++)
        w[i] = v[i];

      uint8_t *table = join_table;

      create_name(str, s, stm, "merged/dtz", -1);
      FILE *F = file_open_read(str);
      read_data_transform_to_u8_u16(F, table, kslice_sizes[s >= 441] * 2, w);
      fclose(F);

    } else {

      uint16_t *table = join_table;

      create_name(str, s, stm, "merged/dtz", -1);
      FILE *F = file_open_read(str);
      read_data_transform_u16(F, table, kslice_sizes[s >= 441] * 2, v);
      fclose(F);

    }

    compress_init_dtz(&dtzmap, tb_wide);

    uint8_t best[MAX_SETS];
    printf("Find optimal permutation for %ctm/dtz, slice = %d.\n", "wb"[stm],
        s);
    permute_piece_462(tb_table, join_table, best, DTZ, tb_wide);
    printf("Compressing data for %ctm/dtz, slice = %d.\n", "wb"[stm], s);
    compress_data_slice(name, stm, DTZ, tb_table, kslice_sizes[s >= 441], best,
        minfreq, tb_wide, false);

    if (!g_cleanup) continue;

    kslice_delete(s, stm, "merged/dtz", -1);
  }

  compress_free_dtz();
}

void join_final_462(int type)
{
  char str[64];
  struct stat st;
  uint64_t slice_size[924];
  bool has_stm[2] = {
    type == WDL || !one_sided || one_sided_stm == WHITE,
    !symmetric && (type == WDL || !one_sided || one_sided_stm == BLACK)
  };

  uint8_t *buf = calloc(1, 2 * 462 * sizeof(uint64_t) + 1000);
  if (!buf)
    out_of_mem();

  uint8_t *p = buf;
  write_le_u32(p, magic2[type]);
  p += 4;
  *p++ = 1; // version number
  *p++ = g_pos.num;
  for (int i = 2; i < g_pos.num; i++)
    *p++ = (g_pos.pt[i] & 7) | ((g_pos.pt[i] & 8) << 4);
  *p++ = LT_PIECE_KK;
  if (type != WDL) {
    uint8_t dist_format = (has_stm[WHITE] && has_stm[BLACK]) ? TWO_SIDED : 0;
    if (one_sided)
      dist_format |= WTM_OR_BTM | (one_sided_stm == WHITE ? WTM_ONLY : 0);
    else
      dist_format |= WIN_OR_LOSS | (wins_only ? WIN_ONLY : 0);
    *p++ = dist_format;
  }
  p = buf + (((p - buf) + 7) & ~7);

  int num = 0;
  for (int s = 0; s < 462; s++)
    for (int stm = 0; stm < 2; stm++) {
      if (!has_stm[stm]) continue;
      create_name(str, s, stm, typename[type], -1);
      if (stat(str, &st) < 0) {
        fprintf(stderr, "Could not access %s.\n", str);
        exit(EXIT_FAILURE);
      }
      slice_size[num++] = st.st_size;
    }

  uint64_t offset = (p - buf) + num * sizeof(uint64_t);
  uint64_t *p64 = (uint64_t *)p;
  for (int i = 0; i < num; i++)
    if (slice_size[i] < 64) {
      p64[i] = offset;
      offset += slice_size[i];
    }
  offset = (offset + 0x3f) & ~0x3fULL;
  for (int i = 0; i < num; i++)
    if (slice_size[i] >= 64) {
      p64[i] = offset;
      offset += slice_size[i];
    }

  char *fname = create_final_name(type);
  char *tmp = malloc(strlen(fname) + 5);
  if (!tmp)
    out_of_mem();
  strcat(strcpy(tmp, fname), ".tmp");
  int fd = creat(tmp, 0666);
  if (fd < 0) {
    fprintf(stderr, "Could not open %s for writing.\n", tmp);
    exit(EXIT_FAILURE);
  }

  write_data_fd(fd, buf, (p - buf) + num * sizeof(uint64_t));

  free(buf);

  int n = 0;
  for (int s = 0; s < 462; s++) {
    for (int stm = 0; stm < 2; stm++) {
      if (!has_stm[stm]) continue;
      if (slice_size[n] < 64) {
        create_name(str, s, stm, typename[type], -1);
        int slice_fd = open(str, O_RDONLY);
        if (slice_fd < 0) {
          fprintf(stderr, "Could not open %s.\n", str);
          exit(EXIT_FAILURE);
        }
        copy_data_fd(slice_fd, fd, slice_size[n]);
        close(slice_fd);
      }
      n++;
    }
  }
  off_t size = lseek(fd, 0, SEEK_CUR);
  if (size == (off_t)-1) {
    fprintf(stderr, "Could not lseek().\n");
    exit(EXIT_FAILURE);
  }
  size_t pad = (0x40 - (size & 0x3f)) & 0x3f;
  if (pad) {
    char zeros[64] = { 0 };
    write_data_fd(fd, zeros, pad);
  }
  n = 0;
  for (int s = 0; s < 462; s++) {
    for (int stm = 0; stm < 2; stm++) {
      if (!has_stm[stm]) continue;
      if (slice_size[n] >= 64) {
        create_name(str, s, stm, typename[type], -1);
        int slice_fd = open(str, O_RDONLY);
        if (slice_fd < 0) {
          fprintf(stderr, "Could not open %s.\n", str);
          exit(EXIT_FAILURE);
        }
        copy_data_fd(slice_fd, fd, slice_size[n]);
        close(slice_fd);
      }
      n++;
    }
  }
  close(fd);
  add_checksum(tmp);
  file_rename(fname);
  free(tmp);
  free(fname);

  if (type == WDL)
    create_empty("done_wdl");

  if (!g_cleanup) return;

  for (int s = 0; s < 462; s++)
    for (int stm = 0; stm < 2; stm++) {
      if (!has_stm[stm]) continue;
        kslice_delete(s, stm, typename[type], -1);
    }

  delete_dir(-1, typename[type]);
}

void join_slices_462(void)
{
  init_permute_piece_462();

  join_table = alloc_huge(kslice_size);
  tb_table = alloc_huge(kslice_size + 1);
  if (!join_table || !tb_table)
    out_of_mem();
  join_wide = tb_wide = false;

  if (!file_exists("done_wdl")) {
    compress_alloc_wdl();

    join_wdl_462(WHITE);
    if (!symmetric)
      join_wdl_462(BLACK);

    compress_free_wdl();

    join_final_462(WDL);

    if (g_cleanup) {
      rmdir("merged/wdl/w");
      rmdir("merged/wdl/b");
      rmdir("merged/wdl");
    }
  }

  if (!one_sided || one_sided_stm == WHITE)
    join_dtz_462(WHITE);

  if ((one_sided && one_sided_stm == BLACK) || (!one_sided && !symmetric))
    join_dtz_462(BLACK);

  free(join_table);
  free(tb_table);

  join_final_462(DTZ);

  if (!g_cleanup) return;

  rmdir("merged/dtz/w");
  rmdir("merged/dtz/b");
  rmdir("merged/dtz");
  rmdir("merged");

  final_cleanup();
}

// Join everything in 10 tables per side.
static void join_wdl_10(int stm)
{
  char str[64];
  uint8_t *table = join_table;

  create_dir(-1, stm, "wdl");
  g_pos.stm = stm;

  for (int k = 0; k < 10; k++) {
    char name[64];
    create_name_10(name, k, stm, "wdl");
    if (file_exists(name))
      continue;

    init_permute_piece_10(k);
    int num = 0;
    uint64_t stats[MAX_STATS] = { 0 };
    bool has_capt_bloss = false;

    g_pos.sq[stm] = InvTriangle[k];

    for (int l = 0; l < 64; l++) {
      if (KKIdx[k][l] < 0)
        continue;

      g_pos.sq[stm ^ 1] = l;

      int s = KKMap[g_pos.sq[0]][g_pos.sq[1]];
      create_name(str, s, stm, "merged/wdl", -1);
      if (read_wdl_slice(str, table + num * kslice_size, kslice_sizes[s >= 441]))
        has_capt_bloss = true;

      create_name(str, s, stm, "stats", -1);
      FILE *F = file_open_read(str);
      uint64_t tmp[MAX_STATS];
      read_data(F, tmp, sizeof tmp);
      fclose(F);
      for (int i = 0; i < MAX_STATS; i++)
        stats[i] += tmp[i];

      num++;
    }

    bool v_wdl[5];
    prepare_wdl_map(stats, v_wdl, has_capt_bloss);
    compress_init_wdl(v_wdl);

    uint8_t best[MAX_SETS];
    printf("Find optimal permutation for %ctm/wdl, slice = %d.\n", "wb"[stm],
        k);
    permute_piece_10(tb_table, table, best, WDL, false);
    printf("Compressing data for %ctm/wdl, slice = %d.\n", "wb"[stm], k);
    compress_data_slice(name, stm, WDL, tb_table, num * kslice_size, best,
        minfreq, false, true);

    if (!g_cleanup) continue;

    // Delete merged slices which are no longer needed.
    for (int l = 0; l < 64; l++) {
      if (KKIdx[k][l] < 0)
        continue;

      g_pos.sq[stm ^ 1] = l;
      int s = KKMap[g_pos.sq[0]][g_pos.sq[1]];

      kslice_delete(s, stm, "merged/wdl", -1);
    }
  }
}

static void join_dtz_10(int stm)
{
  uint16_t v[MAX_STATS];
  char str[64];

  create_dir(-1, stm, "dtz");

  read_merge_info(stm);

  // Call sort_values() to set dtzmap.wide.
  sort_values(stm, g_stats[stm], &dtzmap,
      441 * kslice_sizes[0] + 21 * kslice_sizes[1]);

  if (dtzmap.wide != join_wide) {
    join_wide = dtzmap.wide;
    free(join_table);
    join_table = alloc_huge(58 * kslice_size * (1 + join_wide));
    if (!join_table)
      out_of_mem();
  }

  if (dtzmap.wide != tb_wide) {
    tb_wide = dtzmap.wide;
    free(tb_table);
    tb_table = alloc_huge((58 * kslice_size + 1) * (1 + tb_wide));
    if (!tb_table)
      out_of_mem();
  }

  if (!mi.wide && tb_wide) {
    fprintf(stderr, "Internal error.\n");
    exit(EXIT_FAILURE);
  }

  compress_alloc_dtz(tb_wide);
  g_pos.stm = stm;

  for (int k = 0; k < 10; k++) {
    char name[64];
    create_name_10(name, k, stm, "dtz");
    if (file_exists(name))
      continue;

    init_permute_piece_10(k);
    int num = 0;
    uint64_t stats[MAX_STATS] = { 0 };
    g_pos.sq[stm] = InvTriangle[k];

    for (int l = 0; l < 64; l++) {
      if (KKIdx[k][l] < 0)
        continue;

      g_pos.sq[stm ^ 1] = l;
      int s = KKMap[g_pos.sq[0]][g_pos.sq[1]];

      create_name(str, s, stm, "stats", -1);
      FILE *F = file_open_read(str);
      uint64_t tmp[MAX_STATS];
      read_data(F, tmp, sizeof tmp);
      fclose(F);
      for (int i = 0; i < MAX_STATS; i++)
        stats[i] += tmp[i];

      num++;
    }

    // FIXME: should perhaps double contribution of s >= 441 slices to stats.
    // - or properly count them and reduce num * kslice_size accordingly
    sort_values(stm, stats, &dtzmap, num * kslice_size);
    prepare_dtz_map(v, &dtzmap);

    int n = 0;

    if (!mi.wide) {

      uint8_t w[256];
      for (int i = 0; i < 256; i++)
        w[i] = v[mi.v_inv_u8[i]];

      uint8_t *table = join_table;

      for (int l = 0; l < 64; l++) {
        if (KKIdx[k][l] < 0)
          continue;

        g_pos.sq[stm ^ 1] = l;
        int s = KKMap[g_pos.sq[0]][g_pos.sq[1]];

        create_name(str, s, stm, "merged/dtz", -1);
        FILE *F = file_open_read(str);
        read_data_transform_u8(F, table + n * kslice_size,
            kslice_sizes[s >= 441], w);
        fclose(F);

        n++;
      }

    } else if (!tb_wide) {

      uint8_t w[MAX_STATS];
      for (int i = 0; i < MAX_STATS; i++)
        w[i] = v[i];

      uint8_t *table = join_table;

      for (int l = 0; l < 64; l++) {
        if (KKIdx[k][l] < 0)
          continue;

        g_pos.sq[stm ^ 1] = l;
        int s = KKMap[g_pos.sq[0]][g_pos.sq[1]];

        create_name(str, s, stm, "merged/dtz", -1);
        FILE *F = file_open_read(str);
        read_data_transform_to_u8_u16(F, table + n * kslice_size,
            kslice_sizes[s >= 441] * 2, w);
        fclose(F);

        n++;
      }

    } else {

      uint16_t *table = join_table;

      for (int l = 0; l < 64; l++) {
        if (KKIdx[k][l] < 0)
          continue;

        g_pos.sq[stm ^ 1] = l;
        int s = KKMap[g_pos.sq[0]][g_pos.sq[1]];

        create_name(str, s, stm, "merged/dtz", -1);
        FILE *F = file_open_read(str);
        read_data_transform_u16(F, table + n * kslice_size,
            kslice_sizes[s >= 441] * 2, v);
        fclose(F);

        n++;
      }

    }

    compress_init_dtz(&dtzmap, tb_wide);

    uint8_t best[MAX_SETS];
    printf("Find optimal permutation for %ctm/dtz, slice = %d.\n", "wb"[stm],
        k);
    permute_piece_10(tb_table, join_table, best, DTZ, tb_wide);
    printf("Compressing data for %ctm/dtz, slice = %d.\n", "wb"[stm], k);
    compress_data_slice(name, stm, DTZ, tb_table, num * kslice_size, best,
        minfreq, tb_wide, true);

    if (!g_cleanup) continue;

    // Delete slices which are no longer needed.
    for (int l = 0; l < 64; l++) {
      if (KKIdx[k][l] < 0)
        continue;

      g_pos.sq[stm ^ 1] = l;
      int s = KKMap[g_pos.sq[0]][g_pos.sq[1]];

      kslice_delete(s, stm, "merged/dtz", -1);
    }
  }

  compress_free_dtz();
}

static void join_final_10(int type)
{
  char str[64];
  struct stat st;
  uint64_t slice_size[20];
  bool has_stm[2] = {
    type == WDL || !one_sided || one_sided_stm == WHITE,
    !symmetric && (type == WDL || !one_sided || one_sided_stm == BLACK)
  };

  uint8_t *buf = calloc(1, 2 * 10 * sizeof(uint64_t) + 1000);
  if (!buf)
    out_of_mem();

  uint8_t *p = buf;
  write_le_u32(p, magic2[type]);
  p += 4;
  *p++ = 1; // version number
  *p++ = g_pos.num;
  for (int i = 2; i < g_pos.num; i++)
    *p++ = (g_pos.pt[i] & 7) | ((g_pos.pt[i] & 8) << 4);
  *p++ = LT_PIECE_K;
  if (type != WDL) {
    uint8_t dist_format = (has_stm[WHITE] && has_stm[BLACK]) ? TWO_SIDED : 0;
    if (one_sided)
      dist_format |= WTM_OR_BTM | (one_sided_stm == WHITE ? WTM_ONLY : 0);
    else
      dist_format |= WIN_OR_LOSS | (wins_only ? WIN_ONLY : 0);
    *p++ = dist_format;
  }
  p = buf + (((p - buf) + 7) & ~7);

  int num = 0;
  for (int k = 0; k < 10; k++)
    for (int stm = 0; stm < 2; stm++) {
      if (!has_stm[stm]) continue;
      create_name_10(str, k, stm, typename[type]);
      if (stat(str, &st) < 0) {
        fprintf(stderr, "Could not access %s.\n", str);
        exit(EXIT_FAILURE);
      }
      slice_size[num++] = st.st_size;
    }

  uint64_t offset = (p - buf) + num * sizeof(uint64_t);
  uint64_t *p64 = (uint64_t *)p;
  for (int i = 0; i < num; i++)
    if (slice_size[i] < 64) {
      p64[i] = offset;
      offset += slice_size[i];
    }
  offset = (offset + 0x3f) & ~0x3fULL;
  for (int i = 0; i < num; i++)
    if (slice_size[i] >= 64) {
      p64[i] = offset;
      offset += slice_size[i];
    }

  char *fname = create_final_name(type);
  char *tmp = malloc(strlen(fname) + 5);
  if (!tmp)
    out_of_mem();
  strcat(strcpy(tmp, fname), ".tmp");
  int fd = creat(tmp, 0666);
  if (fd < 0) {
    fprintf(stderr, "Could not open %s for writing.\n", tmp);
    exit(EXIT_FAILURE);
  }

  write_data_fd(fd, buf, (p - buf) + num * sizeof(uint64_t));

  free(buf);

  int n = 0;
  for (int k = 0; k < 10; k++) {
    for (int stm = 0; stm < 2; stm++) {
      if (!has_stm[stm]) continue;
      if (slice_size[n] < 64) {
        create_name_10(str, k, stm, typename[type]);
        int slice_fd = open(str, O_RDONLY);
        if (slice_fd < 0) {
          fprintf(stderr, "Could not open %s.\n", str);
          exit(EXIT_FAILURE);
        }
        copy_data_fd(slice_fd, fd, slice_size[n]);
        close(slice_fd);
      }
      n++;
    }
  }
  off_t size = lseek(fd, 0, SEEK_CUR);
  if (size == (off_t)-1) {
    fprintf(stderr, "Could not lseek().\n");
    exit(EXIT_FAILURE);
  }
  size_t pad = (0x40 - (size & 0x3f)) & 0x3f;
  if (pad) {
    char zeros[64] = { 0 };
    write_data_fd(fd, zeros, pad);
  }
  n = 0;
  for (int k = 0; k < 10; k++) {
    for (int stm = 0; stm < 2; stm++) {
      if (!has_stm[stm]) continue;
      if (slice_size[n] >= 64) {
        create_name_10(str, k, stm, typename[type]);
        int slice_fd = open(str, O_RDONLY);
        if (slice_fd < 0) {
          fprintf(stderr, "Could not open %s.\n", str);
          exit(EXIT_FAILURE);
        }
        copy_data_fd(slice_fd, fd, slice_size[n]);
        close(slice_fd);
      }
      n++;
    }
  }
  close(fd);
  add_checksum(tmp);
  file_rename(fname);
  free(tmp);
  free(fname);

  create_empty("done_wdl");

  if (!g_cleanup) return;

  for (int k = 0; k < 10; k++)
    for (int stm = 0; stm < 2; stm++) {
      if (!has_stm[stm]) continue;
      create_name_10(str, k, stm, typename[type]);
      unlink(str);
    }

  delete_dir(-1, typename[type]);
}

void join_slices_10(void)
{
  join_table = alloc_huge(58 * kslice_size);
  tb_table = alloc_huge(58 * kslice_size + 1);
  if (!join_table || !tb_table)
    out_of_mem();
  join_wide = tb_wide = false;

  if (!file_exists("done_wdl")) {
    compress_alloc_wdl();

    join_wdl_10(WHITE);
    if (!symmetric)
      join_wdl_10(BLACK);

    compress_free_wdl();

    join_final_10(WDL);

    if (g_cleanup) {
      rmdir("merged/wdl/w");
      rmdir("merged/wdl/b");
      rmdir("merged/wdl");
    }
  }

  if (!one_sided || one_sided_stm == WHITE)
    join_dtz_10(WHITE);

  if ((one_sided && one_sided_stm == BLACK) || (!one_sided && !symmetric))
    join_dtz_10(BLACK);

  free(join_table);
  free(tb_table);

  join_final_10(DTZ);

  if (!g_cleanup) return;

  rmdir("merged/dtz/w");
  rmdir("merged/dtz/b");
  rmdir("merged/dtz");
  rmdir("merged");

  final_cleanup();
}
