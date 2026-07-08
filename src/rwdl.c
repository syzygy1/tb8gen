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
#include <unistd.h>

#include "checksum.h"
#include "compress.h"
#include "defs.h"
#include "index.h"
#include "kslice.h"
#include "movegen.h"
#include "permute.h"
#include "permute10.h"
#include "permute462.h"
#include "probe.h"
#include "reduce.h"
#include "rgenerate.h"
#include "rstats.h"
#include "rwdl.h"
#include "tb8gen.h"
#include "threads.h"
#include "types.h"
#include "util.h"
#include "hash/xxhash.h"

static constexpr int minfreq = 8;

static void *join_table, *tb_table;
static uint8_t ram_wdl_map[256];

extern XXH128_hash_t wdl_checksum;

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

static void init_ram_wdl_map(void)
{
  for (int i = 0; i < 256; i++)
    ram_wdl_map[i] = 8;

  ram_wdl_map[RAM_UNRESOLVED] = 2;
  ram_wdl_map[RAM_ILLEGAL] = 8;
  ram_wdl_map[RAM_CAPT_WIN] = 8;
  ram_wdl_map[RAM_CAPT_CWIN] = 7;
  ram_wdl_map[RAM_CAPT_DRAW] = 6;
  ram_wdl_map[RAM_CAPT_BLOSS] = 5;
  ram_wdl_map[RAM_PAWN_WIN] = 4;
  ram_wdl_map[RAM_PAWN_CWIN] = 3;
  ram_wdl_map[RAM_PAWN_DRAW] = 2;

  if (epoch == 0) {
    for (int n = 0; n <= DRAW_RULE + 23; n++)
      ram_wdl_map[loss_to_byte(n)] = n <= DRAW_RULE ? 0 : 1;
    for (int n = 1; n <= DRAW_RULE + 23; n++)
      ram_wdl_map[win_to_byte(n)] = n <= DRAW_RULE ? 4 : 3;
  } else {
    ram_wdl_map[RAM_REDUCED_LOSS] = 0;
    ram_wdl_map[RAM_REDUCED_CAPT_BLOSS] = 6;
    ram_wdl_map[RAM_REDUCED_BLOSS] = 1;
    ram_wdl_map[RAM_REDUCED_CWIN] = 3;
    ram_wdl_map[RAM_REDUCED_CAPT_CWIN] = 7;
    ram_wdl_map[RAM_REDUCED_WIN] = 4;

    int first_loss = reduce_cnt_loss[epoch - 1] + 4;
    int first_win = reduce_cnt_win[epoch - 1] - 250;
    int loss_limit = min(MAX_STATS / 2 - 3, reduce_cnt_loss[epoch - 1] + 127);
    int win_limit = min(MAX_STATS / 2 - 3, reduce_cnt_win[epoch - 1] - 127);

    for (int n = first_loss; n < loss_limit; n++) {
      uint8_t b = loss_to_byte(n);
      if (b > RAM_REDUCED_BLOSS && b < RAM_CAPT_DRAW)
        ram_wdl_map[b] = n <= DRAW_RULE ? 0 : 1;
    }
    for (int n = first_win; n < win_limit; n++) {
      uint8_t b = win_to_byte(n);
      if (b > RAM_CAPT_DRAW && b < RAM_REDUCED_CWIN)
        ram_wdl_map[b] = n <= DRAW_RULE ? 4 : 3;
    }
  }
}

static void update_wdl_vals(bool vals[5], bool dc[4], uint8_t v)
{
  if (v < 5)
    vals[v] = true;
  else if (v < 9)
    dc[v - 5] = true;
}

static void finish_wdl_vals(bool vals[5], bool dc[4])
{
  dc[3] = true;

  int i, j;
  for (i = 0; i < 4; i++)
    if (dc[i]) break;
  for (j = 0; j < 5; j++)
    if (vals[j]) break;
  if (j > i + 1)
    vals[0] = true;
}

static void prepare_wdl_map(uint64_t *stats, bool *vals, bool has_capt_bloss)
{
  for (int i = 0; i < 5; i++)
    vals[i] = false;

  for (int i = 2; i <= DRAW_RULE + 2; i++)
    if (stats[i]) {
      vals[4] = true;
      break;
    }
  if (stats[1])
    vals[4] = true;
  for (int i = DRAW_RULE + 4; i < MAX_STATS / 2 + 1; i++)
    if (stats[i]) {
      vals[3] = true;
      break;
    }
  if (stats[DRAW_RULE + 3])
    vals[3] = true;
  vals[2] = stats[MAX_STATS / 2 + 1] || stats[MAX_STATS / 2 + 2];
  for (int i = DRAW_RULE + 1; i < MAX_STATS / 2 - 3; i++)
    if (stats[MAX_STATS - 1 - i]) {
      vals[1] = true;
      break;
    }
  for (int i = 0; i <= DRAW_RULE; i++)
    if (stats[MAX_STATS - 1 - i]) {
      vals[0] = true;
      break;
    }

  bool dc[4] = { false, false, false, true };
  finish_wdl_vals(vals, dc);
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

static void transform_table(uint8_t *dst, const uint8_t *src, uint64_t size,
    bool vals[5], bool dc[4])
{
  for (uint64_t idx = 0; idx < size; idx++) {
    uint8_t v = ram_wdl_map[src[idx]];
    dst[idx] = v;
    update_wdl_vals(vals, dc, v);
  }
}

static bool transform_slice_stats(uint8_t *dst, const uint8_t *src, int s,
    uint64_t stats[MAX_STATS])
{
  bool vals[5] = { 0 };
  bool dc[4] = { 0 };
  const uint8_t *slice = src + (uint64_t)s * kslice_size;

  if (s < 441) {
    transform_table(dst, slice, kslice_size, vals, dc);
    for (uint64_t idx = 0; idx < kslice_size; idx++) {
      int stat = ram_byte_to_stat(slice[idx]);
      if (stat >= 0)
        stats[stat]++;
    }
    return dc[0];
  }

  Bitboard bb[8];
  bb[0] = bit(KKSquare[s][0]) | bit(KKSquare[s][1]);
  for (uint64_t idx = 0; idx < kslice_sizes[1]; idx++) {
    unrank_bb_ref(idx, bb, &ri);
    uint64_t full_idx = rank_bb(bb, &ri);
    uint8_t b = slice[full_idx];
    uint8_t v = ram_wdl_map[b];
    dst[idx] = v;
    update_wdl_vals(vals, dc, v);
    int stat = ram_byte_to_stat(b);
    if (stat >= 0)
      stats[stat]++;
  }

  return dc[0];
}

static void compress_full_wdl(int stm, struct tb_handle *G, uint8_t *pcs,
    uint8_t *pt)
{
  bool vals[5] = { 0 };
  bool dc[4] = { 0 };
  uint8_t *table = join_table;

  transform_table(table, g_table[stm], table_size, vals, dc);
  table[table_size] = 8;
  prepare_wdl_map(g_stats[stm], vals, dc[0]);
  compress_init_wdl(vals);

  uint8_t best[MAX_PIECES];
  printf("Find optimal permutation for %ctm/wdl.\n", "wb"[stm]);
  permute_piece_wdl(tb_table, pcs, pt, table, best);
  printf("Compressing data for %ctm/wdl.\n", "wb"[stm]);
  compress_tb(G, -1, tb_table, tb_size, best, minfreq, false);
}

static void join_full_wdl(uint8_t *pcs, uint8_t *pt)
{
  init_permute_piece(pcs, pt);

  join_table = alloc_huge(table_size + 1);
  tb_table = alloc_huge(tb_size + 1);
  if (!join_table || !tb_table)
    out_of_mem();

  compress_alloc_wdl();
  struct tb_handle *G = create_tb(g_tablename, WDL, 6);
  compress_full_wdl(WHITE, G, pcs, pt);
  if (!symmetric)
    compress_full_wdl(BLACK, G, pcs, pt);
  merge_tb(G);
  compress_free_wdl();

  free(join_table);
  free(tb_table);
}

static void join_final_462_wdl(void)
{
  char str[64];
  struct stat st;
  uint64_t slice_size[924];
  bool has_stm[2] = { true, !symmetric };

  uint8_t *buf = calloc(1, 2 * 462 * sizeof(uint64_t) + 1000);
  if (!buf)
    out_of_mem();

  uint8_t *p = buf;
  write_le_u32(p, magic2[WDL]);
  p += 4;
  memcpy(p, &wdl_checksum, sizeof wdl_checksum);
  p += sizeof wdl_checksum;
  *p++ = 1;
  *p++ = g_pos.num;
  for (int i = 2; i < g_pos.num; i++)
    *p++ = (g_pos.pt[i] & 7) | ((g_pos.pt[i] & 8) << 4);
  *p++ = LT_PIECE_KK;
  p = buf + (((p - buf) + 7) & ~7);

  int num = 0;
  for (int s = 0; s < 462; s++)
    for (int stm = 0; stm < 2; stm++) {
      if (!has_stm[stm]) continue;
      create_name(str, s, stm, "wdl", -1);
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

  char *fname = create_final_name(WDL);
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

  for (int pass = 0; pass < 2; pass++) {
    if (pass == 1) {
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
    }
    int n = 0;
    for (int s = 0; s < 462; s++)
      for (int stm = 0; stm < 2; stm++) {
        if (!has_stm[stm]) continue;
        bool small = slice_size[n] < 64;
        if ((pass == 0) == small) {
          create_name(str, s, stm, "wdl", -1);
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
  add_xxhash(tmp);
  file_rename(fname);
  free(tmp);
  free(fname);
}

static void compress_462_wdl_side(int stm)
{
  create_dir(-1, stm, "wdl");
  for (int s = 0; s < 462; s++) {
    bool vals[5] = { 0 };
    uint64_t stats[MAX_STATS] = { 0 };
    uint8_t *table = join_table;

    g_slice.sq[0] = KKSquare[s][0];
    g_slice.sq[1] = KKSquare[s][1];
    g_slice.stm = stm;
    bool has_capt_bloss = transform_slice_stats(table, g_table[stm], s, stats);
    prepare_wdl_map(stats, vals, has_capt_bloss);
    compress_init_wdl(vals);

    uint8_t best[MAX_SETS];
    printf("Find optimal permutation for %ctm/wdl, slice = %d.\n", "wb"[stm],
        s);
    permute_piece_462(tb_table, table, best, WDL, false);
    printf("Compressing data for %ctm/wdl, slice = %d.\n", "wb"[stm], s);

    char name[64];
    create_name(name, s, stm, "wdl", -1);
    compress_data_slice(name, stm, WDL, tb_table, kslice_sizes[s >= 441], best,
        minfreq, false, false);
  }
}

static void join_462_wdl(void)
{
  init_permute_piece_462();

  join_table = alloc_huge(kslice_size);
  tb_table = alloc_huge(kslice_size + 1);
  if (!join_table || !tb_table)
    out_of_mem();

  compress_alloc_wdl();
  compress_462_wdl_side(WHITE);
  if (!symmetric)
    compress_462_wdl_side(BLACK);
  compress_free_wdl();
  join_final_462_wdl();

  free(join_table);
  free(tb_table);
}

static void join_final_10_wdl(void)
{
  char str[64];
  struct stat st;
  uint64_t slice_size[20];
  bool has_stm[2] = { true, !symmetric };

  uint8_t *buf = calloc(1, 2 * 10 * sizeof(uint64_t) + 1000);
  if (!buf)
    out_of_mem();

  uint8_t *p = buf;
  write_le_u32(p, magic2[WDL]);
  p += 4;
  memcpy(p, &wdl_checksum, sizeof wdl_checksum);
  p += sizeof wdl_checksum;
  *p++ = 1;
  *p++ = g_pos.num;
  for (int i = 2; i < g_pos.num; i++)
    *p++ = (g_pos.pt[i] & 7) | ((g_pos.pt[i] & 8) << 4);
  *p++ = LT_PIECE_K;
  p = buf + (((p - buf) + 7) & ~7);

  int num = 0;
  for (int k = 0; k < 10; k++)
    for (int stm = 0; stm < 2; stm++) {
      if (!has_stm[stm]) continue;
      create_name_10(str, k, stm, "wdl");
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

  char *fname = create_final_name(WDL);
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

  for (int pass = 0; pass < 2; pass++) {
    if (pass == 1) {
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
    }
    int n = 0;
    for (int k = 0; k < 10; k++)
      for (int stm = 0; stm < 2; stm++) {
        if (!has_stm[stm]) continue;
        bool small = slice_size[n] < 64;
        if ((pass == 0) == small) {
          create_name_10(str, k, stm, "wdl");
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
  add_xxhash(tmp);
  file_rename(fname);
  free(tmp);
  free(fname);
}

static void compress_10_wdl_side(int stm)
{
  create_dir(-1, stm, "wdl");
  for (int k = 0; k < 10; k++) {
    init_permute_piece_10(k);
    memset(join_table, 8, 58 * kslice_size);

    int num = 0;
    uint64_t stats[MAX_STATS] = { 0 };
    bool has_capt_bloss = false;

    g_slice.stm = stm;
    g_slice.sq[stm] = InvTriangle[k];
    for (int l = 0; l < 64; l++) {
      if (KKIdx[k][l] < 0)
        continue;
      g_slice.sq[stm ^ 1] = l;
      int s = KKMap[g_slice.sq[0]][g_slice.sq[1]];
      if (transform_slice_stats((uint8_t *)join_table
            + (uint64_t)num * kslice_size, g_table[stm], s, stats))
        has_capt_bloss = true;
      num++;
    }

    bool vals[5];
    prepare_wdl_map(stats, vals, has_capt_bloss);
    compress_init_wdl(vals);

    uint8_t best[MAX_SETS];
    printf("Find optimal permutation for %ctm/wdl, slice = %d.\n", "wb"[stm],
        k);
    permute_piece_10(tb_table, join_table, best, WDL, false);
    printf("Compressing data for %ctm/wdl, slice = %d.\n", "wb"[stm], k);

    char name[64];
    create_name_10(name, k, stm, "wdl");
    compress_data_slice(name, stm, WDL, tb_table, num * kslice_size, best,
        minfreq, false, true);
  }
}

static void join_10_wdl(void)
{
  join_table = alloc_huge(58 * kslice_size);
  tb_table = alloc_huge(58 * kslice_size + 1);
  if (!join_table || !tb_table)
    out_of_mem();

  compress_alloc_wdl();
  compress_10_wdl_side(WHITE);
  if (!symmetric)
    compress_10_wdl_side(BLACK);
  compress_free_wdl();
  join_final_10_wdl();

  free(join_table);
  free(tb_table);
}

void rjoin_wdl(uint8_t *pcs, uint8_t *pt, int layout)
{
  init_ram_wdl_map();

  switch (layout) {
  case 0:
    join_full_wdl(pcs, pt);
    break;
  case 1:
    join_10_wdl();
    break;
  case 2:
    join_462_wdl();
    break;
  default:
    unreachable();
  }

  create_empty("done_wdl");
}
