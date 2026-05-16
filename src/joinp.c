/*
  Copyright (c) 2026 Ronald de Man

  This file is distributed under the terms of the GNU GPL, version 2.
*/

#include <fcntl.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>

#include "checksum.h"
#include "compress.h"
#include "defs.h"
#include "generatep.h"
#include "index.h"
#include "kslicep.h"
#include "merge.h"
#include "movegen.h"
#include "permute.h"
#include "permutepk.h"
#include "permutep.h"
#include "tb8genp.h"
#include "threads.h"
#include "types.h"
#include "util.h"

static constexpr int minfreq = 8;

static void *join_table, *tb_table;
bool join_wide, tb_wide;
uint8_t g_dist_format;
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

static void sort_values(int stm, uint64_t *stats, struct DtzMap *dtzmap,
    uint64_t tot_pos)
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
      freq[0][i - 1] += stats[2 + i];
    dtzmap->num[0] = sort_list(freq[0], map[0], inv_map[0]);

    for (int i = DRAW_RULE + 1; i < MAX_STATS / 2 - 3; i++)
      freq[2][(i - DRAW_RULE - 1) / 2] += stats[4 + i];
    dtzmap->num[2] = sort_list(freq[2], map[2], inv_map[2]);
  }

  if (one_sided || !wins_only) {
    freq[1][0] = stats[MAX_STATS - 1];
    for (int i = 1; i <= DRAW_RULE; i++)
      freq[1][i - 1] += stats[MAX_STATS - 1 - i];
    dtzmap->num[1] = sort_list(freq[1], map[1], inv_map[1]);

    for (int i = DRAW_RULE + 1; i < MAX_STATS / 2 - 3; i++)
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
      v[2 + i] = inv_map[0][i - 1];
    for (int i = DRAW_RULE + 1; i < MAX_STATS / 2 - 3; i++)
      v[4 + i] = inv_map[2][(i - DRAW_RULE - 1) / 2];
  }
  if (one_sided || !wins_only) {
    v[MAX_STATS - 1] = inv_map[1][0];
    for (int i = 1; i <= DRAW_RULE; i++)
      v[MAX_STATS - 1 - i] = inv_map[1][i - 1];
    for (int i = DRAW_RULE + 1; i < MAX_STATS / 2 - 3; i++)
      v[MAX_STATS - 1 - i] = inv_map[3][(i - DRAW_RULE - 1) / 2];
  }
}

static void read_merge_info(int stm)
{
  char str[64];
  sprintf(str, "merge_info.%c", "wb"[stm]);
  FILE *F = file_open_read(str);
  read_data(F, &mi, sizeof mi);
  fclose(F);
}

static void delete_slices_dir(int stm, const char *name, int psq)
{
  if (!g_cleanup) return;

  char str[64];
  int q = PawnFlip[0][psq];
  for (int k1 = 0; k1 < 64; k1++) {
    if (k1 == psq) continue;
    for (int k2 = 0; k2 < 64; k2++) {
      if (k2 == psq || (king_mask(k1) & bit(k2))) continue;
      sprintf(str, "../%s/", pawnstr[q]);
      create_name_sq(str + strlen(str), k1, k2, stm, name, -1);
      unlink(str);
    }
  }
}

static void delete_merged(int psq)
{
  char str[64];
  int q = PawnFlip[0][psq];
  sprintf(str, "../%s/merged/wdl", pawnstr[q]);
  delete_dir(-1, str);
  sprintf(str, "../%s/merged", pawnstr[q]);
  rmdir(str);
}

static void pawn_slice_cleanup(void)
{
  if (!g_cleanup) return;

  delete_dir(-1, "merged/dtz");
  if (g_pos.sq[2] >= 48) {
    delete_dir(-1, "merged/wdl");
    rmdir("merged");
  }
}

// Compress one pawn-slice as 63 tables of 62 K-slices.
static void join_wdl_pk(int stm)
{
  init_permute_pawn_pk(stm);

  char str[64];
  uint8_t *table = join_table;

  create_dir(-1, stm, "../wdl");

  for (int k1 = 0; k1 < 64; k1++) {
    if (k1 == g_pos.sq[2]) continue;

    char name[64], tmp[64];
    create_name_sq(tmp, g_pos.sq[2], k1, stm, "wdl", -1);
    strcat(strcpy(name, "../"), tmp);
    if (file_exists(name))
      continue;

    uint64_t stats[MAX_STATS] = { 0 };
    g_pos.sq[stm] = k1;

    int num = 0;

    for (int k2 = 0; k2 < 64; k2++) {
      if (k2 == k1 || k2 == g_pos.sq[2]) continue;
      g_pos.sq[stm ^ 1] = k2;

      if (king_attacks(k1) & bit(k2)) {
        memset(table + num * kslice_size, 8, kslice_size);
      } else {
        create_name_sq(str, g_pos.sq[0], g_pos.sq[1], stm, "merged/wdl", -1);
        FILE *F = file_open_read(str);
        read_data(F, table + num * kslice_size, kslice_size);
        fclose(F);

        create_name_sq(str, g_pos.sq[0], g_pos.sq[1], stm, "stats", -1);
        F = file_open_read(str);
        uint64_t tmp[MAX_STATS];
        read_data(F, tmp, sizeof tmp);
        fclose(F);
        for (int i = 0; i < MAX_STATS; i++)
          stats[i] += tmp[i];
      }
      num++;
    }

    bool v_wdl[5] = { 0 };
    for (int i = 2; i <= DRAW_RULE + 2; i++)
      if (stats[i]) {
        v_wdl[4] = true;
        break;
      }
    for (int i = DRAW_RULE + 4; i < MAX_STATS / 2 + 1; i++)
      if (stats[i]) {
        v_wdl[3] = true;
        break;
      }
    v_wdl[2] = (bool)stats[MAX_STATS / 2 + 2];
    for (int i = DRAW_RULE + 1; i < MAX_STATS / 2 - 3; i++)
      if (stats[MAX_STATS - 1 - i]) {
        v_wdl[1] = true;
        break;
      }
    for (int i = 0; i <= DRAW_RULE; i++)
      if (stats[MAX_STATS - 1 - i]) {
        v_wdl[0] = true;
        break;
      }

    compress_init_wdl(v_wdl);

    uint8_t best[MAX_SETS];
    printf("Find optimal permutation for %ctm/wdl, slice = (%d,%d).\n",
        "wb"[stm], g_pos.sq[2], g_pos.sq[stm]);
    permute_pawn_pk(tb_table, table, best, WDL, false);
    printf("Compressing data for %ctm/wdl, slice = (%d,%d).\n", "wb"[stm],
        g_pos.sq[2], g_pos.sq[stm]);
    compress_data_slice(name, stm, WDL, tb_table, 62 * kslice_size, best,
        minfreq, false, false);

    if (!g_cleanup) continue;

    // FIXME: do not delete merged/wdl too soon; needed for pawn moves.
    for (int k2 = 0; k2 < 64; k2++) {
      g_pos.sq[stm ^ 1] = k2;
      if (is_broken(&g_pos)) continue;

      create_name_sq(str, g_pos.sq[0], g_pos.sq[1], stm, "merged/wdl", -1);
//      unlink(str);
    }
  }
}

static void join_dtz_pk(int stm)
{
  uint16_t v[MAX_STATS];
  char str[64];

  create_dir(-1, stm, "../dtz");

  read_merge_info(stm);

  g_dist_format = one_sided ? 0 : WIN_OR_LOSS | (wins_only ? WIN_ONLY : 0);

  // FIXME: 63*62 ?
  sort_values(stm, g_stats[stm], &dtzmap, 63*62 * kslice_size);

  if (dtzmap.wide != join_wide) {
    join_wide = dtzmap.wide;
    free(join_table);
    join_table = alloc_huge(62 * kslice_size * (1 + join_wide));
    if (!join_table)
      out_of_mem();
  }

  if (dtzmap.wide != tb_wide) {
    tb_wide = dtzmap.wide;
    free(tb_table);
    tb_table = alloc_huge((62 * kslice_size + 1) * (1 + tb_wide));
    if (!tb_table)
      out_of_mem();
  }

  if (!mi.wide && tb_wide) {
    fprintf(stderr, "Internal error.\n");
    exit(EXIT_FAILURE);
  }

  compress_alloc_dtz(tb_wide);

  init_permute_pawn_pk(stm);
  for (int k1 = 0; k1 < 64; k1++) {
    if (k1 == g_pos.sq[2]) continue;

    char name[64], tmp[64];
    create_name_sq(tmp, g_pos.sq[2], k1, stm, "dtz", -1);
    strcat(strcpy(name, "../"), tmp);
    if (file_exists(name))
      continue;

    uint64_t stats[MAX_STATS] = { 0 };
    g_pos.sq[stm] = k1;

    for (int k2 = 0; k2 < 64; k2++) {
      if (k2 == k1 || k2 == g_pos.sq[2]) continue;

      g_pos.sq[stm ^ 1] = k2;

      if (!(king_attacks(k1) & bit(k2))) {
        create_name_sq(str, g_pos.sq[0], g_pos.sq[1], stm, "stats", -1);
        FILE *F = file_open_read(str);
        uint64_t tmp[MAX_STATS];
        read_data(F, tmp, sizeof tmp);
        fclose(F);
        for (int i = 0; i < MAX_STATS; i++)
          stats[i] += tmp[i];
      }
    }

    sort_values(stm, stats, &dtzmap, 62 * kslice_size);
    prepare_dtz_map(v, &dtzmap);

    int n = 0;

    if (!mi.wide) {

      uint8_t w[256];
      for (int i = 0; i < 256; i++)
        w[i] = v[mi.v_inv_u8[i]];

      uint8_t *table = join_table;

      for (int k2 = 0; k2 < 64; k2++) {
        if (k2 == k1 || k2 == g_pos.sq[2]) continue;
        
        g_pos.sq[stm ^ 1] = k2;

        if (king_attacks(k1) & bit(k2)) {
          memset(table + n * kslice_size, v[0], kslice_size);
        } else {
          create_name_sq(str, g_pos.sq[0], g_pos.sq[1], stm, "merged/dtz", -1);
          FILE *F = file_open_read(str);
          read_data_transform_u8(F, table + n * kslice_size, kslice_size, w);
          fclose(F);
        }
        n++;
      }

    } else if (!tb_wide) {

      uint8_t w[MAX_STATS];
      for (int i = 0; i < MAX_STATS; i++)
        w[i] = v[i];

      uint8_t *table = join_table;

      for (int k2 = 0; k2 < 64; k2++) {
        if (k2 == k1 || k2 == g_pos.sq[2]) continue;

        g_pos.sq[stm ^ 1] = k2;

        if (king_attacks(k1) & bit(k2)) {
          memset(table + n * kslice_size, v[0], kslice_size);
        } else {
          create_name_sq(str, g_pos.sq[0], g_pos.sq[1], stm, "merged/dtz", -1);
          FILE *F = file_open_read(str);
          read_data_transform_to_u8_u16(F, table + n * kslice_size,
              kslice_size * 2, w);
          fclose(F);
        }
        n++;
      }

    } else {

      uint16_t *table = join_table;

      for (int k2 = 0; k2 < 64; k2++) {
        if (k2 == k1 || k2 == g_pos.sq[2]) continue;

        g_pos.sq[stm ^ 1] = k2;

        if (king_attacks(k1) & bit(k2)) {
          for (size_t i = 0; i < kslice_size; i++)
            table[n * kslice_size + i] = v[0];
        } else {
          create_name_sq(str, g_pos.sq[0], g_pos.sq[1], stm, "merged/dtz", -1);
          FILE *F = file_open_read(str);
          read_data_transform_u16(F, table + n * kslice_size, kslice_size * 2,
              v);
          fclose(F);
        }
        n++;
      }

    }

    compress_init_dtz(&dtzmap, tb_wide);

    uint8_t best[MAX_SETS];
    printf("Find optimal permutation for %ctm/dtz, slice = %d.\n", "wb"[stm],
        k1);
    permute_pawn_pk(tb_table, join_table, best, DTZ, tb_wide);
    printf("Compressing data for %ctm/dtz, slice = %d.\n", "wb"[stm], k1);
    compress_data_slice(name, stm, DTZ, tb_table, 62 * kslice_size, best,
        minfreq, tb_wide, false);

    if (!g_cleanup) continue;

    for (int k2 = 0; k2 < 64; k2++) {
      g_pos.sq[stm ^ 1] = k2;
      if (is_broken(&g_pos)) continue;

      create_name_sq(str, g_pos.sq[0], g_pos.sq[1], stm, "merged/dtz", -1);
      unlink(str);
      create_name_sq(str, g_pos.sq[0], g_pos.sq[1], stm, "stats", -1);
      unlink(str);
    }
  }

  compress_free_dtz();
}

void join_final_pk(int type)
{
  char str[64];
  struct stat st;
  uint64_t slice_size[2 * 1512];
  bool has_stm[24][2];

  for (int i = 0; i < 24; i++) {
    has_stm[i][WHITE] = type == WDL || !dtz_format[i].one_sided
        || dtz_format[i].one_sided_stm == WHITE;
    has_stm[i][BLACK] = !symmetric && (type == WDL || !dtz_format[i].one_sided
        || dtz_format[i].one_sided_stm == BLACK);
  }

  uint8_t *buf = calloc(1, 2 * 1512 * sizeof(uint64_t) + 1000);
  if (!buf)
    out_of_mem();

  uint8_t *p = buf;
  write_le_u32(p, magic2[type]);
  p += 4;
  *p++ = 1; // version number
  *p++ = g_pos.num;
  for (int i = 2; i < g_pos.num; i++)
    *p++ = (g_pos.pt[i] & 7) | ((g_pos.pt[i] & 8) << 4);
  *p++ = LT_PAWN_PK;
  p = buf + (((p - buf) + 7) & ~7);

  int sides = symmetric ? 1 : 2;
  int num = 0;
  for (int q = 0; q < 24; q++) {
    int psq = InvPawnFlip[0][q];
    for (int k1 = 0; k1 < 64; k1++) {
      if (k1 == psq) continue;
      for (int stm = 0; stm < sides; stm++) {
        if (!has_stm[q][stm]) {
          slice_size[num++] = 0;
          continue;
        }
        create_name_sq(str, psq, k1, stm, typename[type], -1);
        if (stat(str, &st) < 0) {
          fprintf(stderr, "Could not access %s.\n", str);
          exit(EXIT_FAILURE);
        }
        slice_size[num++] = st.st_size;
      }
    }
  }
  assert(num == 1512 * sides);

  uint64_t offset = (p - buf) + num * sizeof(uint64_t);
  uint64_t *p64 = (uint64_t *)p;
  for (int i = 0; i < num; i++)
    if (slice_size[i] < 64) {
      p64[i] = slice_size[i] > 0 ? offset : 0;
      offset += slice_size[i];
    }
  offset = (offset + 0x3f) & ~0x3fULL;
  for (int i = 0; i < num; i++)
    if (slice_size[i] >= 64) {
      p64[i] = offset;
      offset += slice_size[i];
    }

  char fname[64], tmp[64];
  sprintf(fname, "../%s%s", g_tablename, suffix[type]);
  strcat(strcpy(tmp, fname), ".tmp");
  int fd = creat(tmp, 0666);
  if (fd < 0) {
    fprintf(stderr, "Could not open %s for writing.\n", tmp);
    exit(EXIT_FAILURE);
  }

  write_data_fd(fd, buf, (p - buf) + num * sizeof(uint64_t));

  free(buf);

  int n = 0;
  for (int q = 0; q < 24; q++) {
    int psq = InvPawnFlip[0][q];
    for (int k1 = 0; k1 < 64; k1++) {
      if (k1 == psq) continue;
      for (int stm = 0; stm < sides; stm++) {
        if (slice_size[n] > 0 && slice_size[n] < 64) {
          create_name_sq(str, psq, k1, stm, typename[type], -1);
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
  }
  off_t size = lseek(fd, 0, SEEK_CUR);
  if (size == (off_t)-1) {
    fprintf(stderr, "Could not lseek().\n");
    exit(EXIT_FAILURE);
  }
  size_t pad = (-size) & 0x3f;
  if (pad) {
    char zeros[64] = { 0 };
    write_data_fd(fd, zeros, pad);
  }
  n = 0;
  for (int q = 0; q < 24; q++) {
    int psq = InvPawnFlip[0][q];
    for (int k1 = 0; k1 < 64; k1++) {
      if (k1 == psq) continue;
      for (int stm = 0; stm < sides; stm++) {
        if (slice_size[n] >= 64) {
          create_name_sq(str, psq, k1, stm, typename[type], -1);
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
  }
  close(fd);
  add_checksum(tmp);
  rename(tmp, fname);

  if (!g_cleanup) return;

  for (int q = 0; q < 24; q++) {
    int psq = InvPawnFlip[0][q];
    for (int k1 = 0; k1 < 64; k1++) {
      if (k1 == psq) continue;
      for (int stm = 0; stm < sides; stm++) {
        create_name_sq(str, psq, k1, stm, typename[type], -1);
        unlink(str);
      }
    }
  }

  delete_dir(-1, typename[type]);
}

void compress_slice_pk(void)
{
  join_table = alloc_huge(62 * kslice_size);
  tb_table = alloc_huge(62 * kslice_size + 1);
  if (!join_table || !tb_table)
    out_of_mem();
  join_wide = tb_wide = false;

  compress_alloc_wdl();

  join_wdl_pk(WHITE);
  if (!symmetric)
    join_wdl_pk(BLACK);

  compress_free_wdl();

  if (!one_sided || one_sided_stm == WHITE)
    join_dtz_pk(WHITE);

  if (   (one_sided && one_sided_stm == BLACK)
      || (!one_sided && !symmetric))
    join_dtz_pk(BLACK);

  free(join_table);
  free(tb_table);

  if (!g_cleanup) return;

  pawn_slice_cleanup();
}

void join_slices_pk(void)
{
  join_final_pk(WDL);
  join_final_pk(DTZ);
}

// Join everything in 10 tables per side.
static void join_wdl_p(int stm)
{
  char name[64], tmp[64];
  create_name_p(tmp, g_pos.sq[2], stm, "wdl");
  strcat(strcpy(name, "../"), tmp);
  if (file_exists(name))
    return;

  uint64_t stats[MAX_STATS] = { 0 };
  uint8_t *table = join_table;
  char str[64];

  create_dir(-1, stm, "../wdl");

  g_pos.stm = stm;
  int num = 0;
  for (int k1 = 0; k1 < 64; k1++) {
    if (k1 == g_pos.sq[2])
      continue;

    g_pos.sq[0] = k1;    // FIXME: sq[stm] ?

    for (int k2 = 0; k2 < 64; k2++) {
      if (k2 == k1 || k2 == g_pos.sq[2])
        continue;

      g_pos.sq[1] = k2;

      if (king_attacks(k1) & bit(k2)) {
        memset(table + num * kslice_size, 8, kslice_size);
// FIXME: increase stats[] ?
      } else {

        create_name_sq(str, k1, k2, stm, "merged/wdl", -1);
        FILE *F = file_open_read(str);
        read_data(F, table + num * kslice_size, kslice_size);
        fclose(F);

        create_name_sq(str, k1, k2, stm, "stats", -1);
        F = file_open_read(str);
        uint64_t tmp[MAX_STATS];
        read_data(F, tmp, sizeof tmp);
        fclose(F);
        for (int i = 0; i < MAX_STATS; i++)
          stats[i] += tmp[i];
      }

      num++;
    }
  }

  bool v_wdl[5] = { 0 };
  for (int i = 2; i <= DRAW_RULE + 2; i++)
    if (stats[i]) {
      v_wdl[4] = true;
      break;
    }
  for (int i = DRAW_RULE + 4; i < MAX_STATS / 2 + 1; i++)
    if (stats[i]) {
      v_wdl[3] = true;
      break;
    }
  v_wdl[2] = (bool)stats[MAX_STATS / 2 + 2];
  for (int i = DRAW_RULE + 1; i < MAX_STATS / 2 - 3; i++)
    if (stats[MAX_STATS - 1 - i]) {
      v_wdl[1] = true;
      break;
    }
  for (int i = 0; i <= DRAW_RULE; i++)
    if (stats[MAX_STATS - 1 - i]) {
      v_wdl[0] = true;
      break;
    }

  compress_init_wdl(v_wdl);

  uint8_t best[MAX_SETS];
  printf("Find optimal permutation for %ctm/wdl.\n", "wb"[stm]);
  permute_pawn_p(tb_table, table, best, WDL, false);
  printf("Compressing data for %ctm/wdl.\n", "wb"[stm]);
  compress_data_slice(name, stm, WDL, tb_table, tb_size, best, minfreq, false,
      true);

  if (!g_cleanup) return;

  if (g_pos.sq[2] >= 16 && g_pos.sq[2] < 40) {
    delete_slices_dir(stm, "merged/wdl", g_pos.sq[2] - 8);
    if (stm == BLACK)
      delete_merged(g_pos.sq[2] - 8);
  }
  else if (g_pos.sq[2] >= 48) {
    delete_slices_dir(stm, "merged/wdl", g_pos.sq[2] - 16);
    delete_slices_dir(stm, "merged/wdl", g_pos.sq[2] - 8);
    delete_slices_dir(stm, "merged/wdl", g_pos.sq[2]);
    if (stm == BLACK) {
      delete_merged(g_pos.sq[2] - 16);
      delete_merged(g_pos.sq[2] - 8);
    }
  }
}

static void join_dtz_p(int stm)
{
  char name[64], tmp[64];
  create_name_p(tmp, g_pos.sq[2], stm, "dtz");
  strcat(strcpy(name, "../"), tmp);
  if (file_exists(name))
    return;

  uint64_t stats[MAX_STATS] = { 0 };
  uint16_t v[MAX_STATS];
  char str[64];

  create_dir(-1, stm, "../dtz");

  read_merge_info(stm);

  g_dist_format = one_sided ? 0 : WIN_OR_LOSS | (wins_only ? WIN_ONLY : 0);

  // FIXME: lower 63*62 here?
  sort_values(stm, g_stats[stm], &dtzmap, 63*62 * kslice_size);

  if (dtzmap.wide != join_wide) {
    join_wide = dtzmap.wide;
    free(join_table);
    join_table = alloc_huge(63*62 * kslice_size * (1 + join_wide));
    if (!join_table)
      out_of_mem();
  }

  if (dtzmap.wide != tb_wide) {
    tb_wide = dtzmap.wide;
    free(tb_table);
    tb_table = alloc_huge((63*62 * kslice_size + 1) * (1 + tb_wide));
    if (!tb_table)
      out_of_mem();
  }

  if (!mi.wide && tb_wide) {
    fprintf(stderr, "Internal error.\n");
    exit(EXIT_FAILURE);
  }

  compress_alloc_dtz(tb_wide);

  g_pos.stm = stm;
  for (int k1 = 0; k1 < 64; k1++) {
    if (k1 == g_pos.sq[2])
      continue;

    g_pos.sq[0] = k1;

    for (int k2 = 0; k2 < 64; k2++) {
      if (k2 == k1 || k2 == g_pos.sq[2])
        continue;

      g_pos.sq[1] = k2;

      if (king_attacks(k1) & bit(k2)) {
        // increase stats[] ?
      } else {

        create_name_sq(str, k1, k2, stm, "stats", -1);
        FILE *F = file_open_read(str);
        uint64_t tmp[MAX_STATS];
        read_data(F, tmp, sizeof tmp);
        fclose(F);
        for (int i = 0; i < MAX_STATS; i++)
          stats[i] += tmp[i];
      }
    }
  }

  sort_values(stm, stats, &dtzmap, 62*63 * kslice_size);
  prepare_dtz_map(v, &dtzmap);

  int n = 0;

  if (!mi.wide) {

    uint8_t w[256];
    for (int i = 0; i < 256; i++)
      w[i] = v[mi.v_inv_u8[i]];

    uint8_t *table = join_table;

    for (int k1 = 0; k1 < 64; k1++) {
      if (k1 == g_pos.sq[2])
        continue;
      g_pos.sq[0] = k1;

      for (int k2 = 0; k2 < 64; k2++) {
        if (k2 == k1 || k2 == g_pos.sq[2])
          continue;
        g_pos.sq[1] = k2;

        if (king_attacks(k1) & bit(k2)) {
          memset(table + n * kslice_size, v[0], kslice_size);
        } else {
          create_name_sq(str, k1, k2, stm, "merged/dtz", -1);
          FILE *F = file_open_read(str);
          read_data_transform_u8(F, table + n * kslice_size, kslice_size, w);
          fclose(F);
        }
        n++;
      }
    }

  } else if (!tb_wide) {

    uint8_t w[MAX_STATS];
    for (int i = 0; i < MAX_STATS; i++)
      w[i] = v[i];

    uint8_t *table = join_table;

    for (int k1 = 0; k1 < 64; k1++) {
      if (k1 == g_pos.sq[2])
        continue;
      g_pos.sq[0] = k1;

      for (int k2 = 0; k2 < 64; k2++) {
        if (k2 == k1 || k2 == g_pos.sq[2])
          continue;
        g_pos.sq[1] = k2;

        if (king_attacks(k1) & bit(k2)) {
          memset(table + n * kslice_size, v[0], kslice_size);
        } else {
          create_name_sq(str, k1, k2, stm, "merged/dtz", -1);
          FILE *F = file_open_read(str);
          read_data_transform_to_u8_u16(F, table + n * kslice_size,
              kslice_size * 2, w);
          fclose(F);
        }
        n++;
      }
    }

  } else {

    uint16_t *table = join_table;

    for (int k1 = 0; k1 < 64; k1++) {
      if (k1 == g_pos.sq[2])
        continue;
      g_pos.sq[0] = k1;

      for (int k2 = 0; k2 < 64; k2++) {
        if (k2 == k1 || k2 == g_pos.sq[2])
          continue;
        g_pos.sq[1] = k2;

        if (king_attacks(k1) & bit(k2)) {
          for (size_t i = 0; i < kslice_size; i++)
            table[n * kslice_size + i] = v[0];
        } else {
          create_name_sq(str, k1, k2, stm, "merged/dtz", -1);
          FILE *F = file_open_read(str);
          read_data_transform_u16(F, table + n * kslice_size, kslice_size * 2,
              v);
          fclose(F);
        }
        n++;
      }
    }

  }

  compress_init_dtz(&dtzmap, tb_wide);

  uint8_t best[MAX_SETS];
  printf("Find optimal permutation for %ctm/dtz.\n", "wb"[stm]);
  permute_pawn_p(tb_table, join_table, best, DTZ, tb_wide);
  printf("Compressing data for %ctm/dtz.\n", "wb"[stm]);
  compress_data_slice(name, stm, DTZ, tb_table, tb_size, best, minfreq, tb_wide,
      true);

  compress_free_dtz();

  if (!g_cleanup) return;

  delete_slices_dir(stm, "merged/dtz", g_pos.sq[2]);
}

static void join_final_p(int type)
{
  char str[64];
  struct stat st;
  uint64_t slice_size[2 * 24];
  bool has_stm[24][2];

  for (int i = 0; i < 24; i++) {
    has_stm[i][WHITE] = type == WDL || !dtz_format[i].one_sided
        || dtz_format[i].one_sided_stm == WHITE;
    has_stm[i][BLACK] = !symmetric && (type == WDL || !dtz_format[i].one_sided
        || dtz_format[i].one_sided_stm == BLACK);
  }

  uint8_t *buf = calloc(1, 2 * 24 * sizeof(uint64_t) + 1000);
  if (!buf)
    out_of_mem();

  uint8_t *p = buf;
  write_le_u32(p, magic2[type]);
  p += 4;
  *p++ = 1; // version number
  *p++ = g_pos.num;
  for (int i = 2; i < g_pos.num; i++)
    *p++ = (g_pos.pt[i] & 7) | ((g_pos.pt[i] & 8) << 4);
  *p++ = LT_PAWN_P;
  p = buf + (((p - buf) + 7) & ~7);

  int sides = symmetric ? 1 : 2;
  int num = 0;
  for (int q = 0; q < 24; q++) {
    int psq = InvPawnFlip[0][q];
    for (int stm = 0; stm < sides; stm++) {
      if (!has_stm[q][stm]) {
        slice_size[num++] = 0;
        continue;
      }
      create_name_p(str, psq, stm, typename[type]);
      if (stat(str, &st) < 0) {
        fprintf(stderr, "Could not access %s.\n", str);
        exit(EXIT_FAILURE);
      }
      slice_size[num++] = st.st_size;
    }
  }
  assert(num == 24 * sides);

  uint64_t offset = (p - buf) + num * sizeof(uint64_t);
  uint64_t *p64 = (uint64_t *)p;
  for (int i = 0; i < num; i++)
    if (slice_size[i] < 64) {
      p64[i] = slice_size[i] > 0 ? offset : 0;
      offset += slice_size[i];
    }
  offset = (offset + 0x3f) & ~0x3fULL;
  for (int i = 0; i < num; i++)
    if (slice_size[i] >= 64) {
      p64[i] = offset;
      offset += slice_size[i];
    }

  char fname[64], tmp[64];
  sprintf(fname, "../%s%s", g_tablename, suffix[type]);
  strcat(strcpy(tmp, fname), ".tmp");
  int fd = creat(tmp, 0666);
  if (fd < 0) {
    fprintf(stderr, "Could not open %s for writing.\n", tmp);
    exit(EXIT_FAILURE);
  }

  write_data_fd(fd, buf, (p - buf) + num * sizeof(uint64_t));

  free(buf);

  int n = 0;
  for (int q = 0; q < 24; q++) {
    int psq = InvPawnFlip[0][q];
    for (int stm = 0; stm < sides; stm++) {
      if (slice_size[n] > 0 && slice_size[n] < 64) {
        create_name_p(str, psq, stm, typename[type]);
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
  size_t pad = (-size) & 0x3f;
  if (pad) {
    char zeros[64] = { 0 };
    write_data_fd(fd, zeros, pad);
  }
  n = 0;
  for (int q = 0; q < 24; q++) {
    int psq = InvPawnFlip[0][q];
    for (int stm = 0; stm < sides; stm++) {
      if (slice_size[n] >= 64) {
        create_name_p(str, psq, stm, typename[type]);
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
  rename(tmp, fname);

  if (!g_cleanup) return;

  for (int q = 0; q < 24; q++) {
    int psq = InvPawnFlip[0][q];
    for (int stm = 0; stm < sides; stm++) {
      create_name_p(str, psq, stm, typename[type]);
      unlink(str);
    }
  }

  delete_dir(-1, typename[type]);
}

void compress_slice_p(void)
{
  join_table = alloc_huge(63*62 * kslice_size);
  tb_table = alloc_huge(63*62 * kslice_size + 1);
  if (!join_table || !tb_table)
    out_of_mem();
  join_wide = tb_wide = false;

  init_permute_pawn_p();

  compress_alloc_wdl();

  join_wdl_p(WHITE);
  if (!symmetric)
    join_wdl_p(BLACK);

  compress_free_wdl();

  if (!one_sided || one_sided_stm == WHITE)
    join_dtz_p(WHITE);

  if (   (one_sided && one_sided_stm == BLACK)
      || (!one_sided && !symmetric))
    join_dtz_p(BLACK);

  if (g_pos.sq[2] >= 48)
    unlink("merged");

  free(join_table);
  free(tb_table);

  if (!g_cleanup) return;

  pawn_slice_cleanup();
}

void join_slices_p(void)
{
  join_final_p(WDL);
  join_final_p(DTZ);
}
