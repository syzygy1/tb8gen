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
#include "rpermute.h"
#include "rpermute10.h"
#include "rpermute462.h"
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

static void *dtz_table, *tb_table;
static uint8_t ram_wdl_map[256];
struct DtzMap dtzmap;

extern XXH128_hash_t wdl_checksum;

static struct Work work_stats;

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
    ram_wdl_map[RAM_REDUCED_CAPT_BLOSS] = 5;
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

static _Atomic bool found_capt_bloss;
static uint8_t capt_bloss_val;
static void *work_table;
static uint64_t (*per_thread_stats)[MAX_STATS];

static void collect_stats_u8_worker(struct ThreadData *thread)
{
  uint64_t *restrict stats = per_thread_stats[thread->thread_id];
  uint8_t *restrict table = work_table;

  for (uint64_t idx = thread->begin, end = thread->end; idx < end; idx++)
    stats[table[idx]]++;
}

static void collect_stats_u16_worker(struct ThreadData *thread)
{
  uint64_t *restrict stats = per_thread_stats[thread->thread_id];
  uint16_t *restrict table = work_table;

  for (uint64_t idx = thread->begin, end = thread->end; idx < end; idx++)
    stats[table[idx]]++;
}

static int byte_to_wdl_stat(uint8_t b)
{
  // Capture blessed losses are handled as WDL don't-cares through
  // find_capt_bloss(), not as an ordinary blessed-loss value in a slice.
  if (epoch == 0 && b == RAM_CAPT_BLOSS)
    return -1;

  int s = byte_to_stat(b);
  if (s >= 0 || epoch == 0)
    return s;

  switch (b) {
  case RAM_REDUCED_LOSS:
    return MAX_STATS - 1;
  case RAM_REDUCED_BLOSS:
    return loss_to_stat(DRAW_RULE + 1);
  case RAM_REDUCED_CWIN:
    return DRAW_RULE + 4;
  case RAM_REDUCED_CAPT_CWIN:
    return DRAW_RULE + 3;
  case RAM_REDUCED_WIN:
    return 2;
  default:
    return -1;
  }
}

static void wdl_collect_stats(uint8_t *restrict table, uint64_t *restrict stats)
{
  if (!per_thread_stats)
    per_thread_stats = aligned_alloc(64, g_num_threads * MAX_STATS * 8);
  if (!per_thread_stats)
    out_of_mem();
  memset(per_thread_stats, 0, g_num_threads * MAX_STATS * 8);
  work_table = table;
  run_threaded(collect_stats_u8_worker, &work_stats);

  for (int t = 0; t < g_num_threads; t++)
    for (int b = 0; b < 256; b++) {
      int s = byte_to_wdl_stat(b);
      if (s >= 0)
        stats[s] += per_thread_stats[t][b];
    }
}

static void find_capt_bloss_worker(struct ThreadData *thread)
{
  if (atomic_load_explicit(&found_capt_bloss, memory_order_relaxed))
    return;

  uint8_t v = capt_bloss_val;
  uint8_t *restrict table = work_table;

  for (uint64_t idx = thread->begin, end = thread->end; idx < end; idx++)
    if (table[idx] == v) {
      atomic_store_explicit(&found_capt_bloss, true, memory_order_relaxed);
      break;
    }
}

static bool find_capt_bloss(uint8_t *table)
{
  capt_bloss_val = epoch == 0 ? RAM_CAPT_BLOSS : RAM_REDUCED_CAPT_BLOSS;
  work_table = table;
  atomic_store_explicit(&found_capt_bloss, false, memory_order_relaxed);
  run_threaded(find_capt_bloss_worker, &work_stats);
  return atomic_load_explicit(&found_capt_bloss, memory_order_relaxed);
}

static void dtz_collect_stats_u8(int stm, uint8_t *restrict table,
    uint64_t *restrict stats)
{
  if (!per_thread_stats)
    per_thread_stats = aligned_alloc(64, g_num_threads * MAX_STATS * 8);
  if (!per_thread_stats)
    out_of_mem();
  memset(per_thread_stats, 0, g_num_threads * MAX_STATS * 8);
  work_table = table;

  run_threaded(collect_stats_u8_worker, &work_stats);

  for (int t = 0; t < g_num_threads; t++)
    for (int i = 1; i < 256; i++)
      stats[val_to_stats[stm][i]] += per_thread_stats[t][i];
  stats[0] = 0;
}

static void dtz_collect_stats_u16(int stm, uint16_t *restrict table,
    uint64_t *restrict stats)
{
  if (!per_thread_stats)
    per_thread_stats = aligned_alloc(64, g_num_threads * MAX_STATS * 8);
  if (!per_thread_stats)
    out_of_mem();
  memset(per_thread_stats, 0, g_num_threads * MAX_STATS * 8);
  work_table = table;

  run_threaded(collect_stats_u16_worker, &work_stats);

  for (int t = 0; t < g_num_threads; t++)
    for (int i = 1; i < MAX_STATS; i++)
      stats[val_to_stats[stm][i]] += per_thread_stats[t][i];
  stats[0] = 0;
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

static void prepare_wdl_map(uint64_t *stats, bool *v, bool has_capt_bloss)
{
  for (int i = 0; i < 5; i++)
    v[i] = false;

  for (int i = 2; i <= DRAW_RULE + 2; i++)
    if (stats[i]) {
      v[4] = true;
      break;
    }
  for (int i = DRAW_RULE + 4; i < MAX_STATS / 2 + 1; i++)
    if (stats[i]) {
      v[3] = true;
      break;
    }
  v[2] = (bool)stats[MAX_STATS / 2 + 2];
  for (int i = DRAW_RULE + 1; i < MAX_STATS / 2 - 3; i++)
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
    has_capt_bloss, stats[MAX_STATS / 2 + 1], stats[DRAW_RULE + 3], true
  };

  int i, j;
  for (i = 0; i < 4; i++)
    if (dc[i]) break;
  for (j = 0; j < 5; j++)
    if (v[j]) break;
  if (j > i + 1)
    v[0] = true;
}

static void compress_full_wdl(int stm, struct tb_handle *G, uint8_t *pcs,
    uint8_t *pt)
{
  bool vals[5];
  uint8_t *table = g_table[stm];

  work_init(&work_stats, table_size, 0, WORK_DYNAMIC, 16, 512);
  bool has_capt_bloss = find_capt_bloss(g_table[stm]);
  table[table_size] = RAM_ILLEGAL;
  prepare_wdl_map(g_stats[stm], vals, has_capt_bloss);
  compress_init_wdl(vals);

  uint8_t best[MAX_PIECES];
  printf("Find optimal permutation for %ctm/wdl.\n", "wb"[stm]);
  permute_piece_wdl(tb_table, pcs, pt, table, best, ram_wdl_map);
  printf("Compressing data for %ctm/wdl.\n", "wb"[stm]);
  compress_tb(G, -1, tb_table, tb_size, best, minfreq, false);
}

static void join_full_wdl(uint8_t *pcs, uint8_t *pt)
{
  init_permute_piece(pcs, pt);

  tb_table = alloc_huge(tb_size + 1);
  if (!tb_table)
    out_of_mem();

  compress_alloc_wdl();
  struct tb_handle *G = create_tb(g_tablename, WDL, 6);
  compress_full_wdl(WHITE, G, pcs, pt);
  if (!symmetric)
    compress_full_wdl(BLACK, G, pcs, pt);
  merge_tb(G);
  compress_free_wdl();

  free(tb_table);
}

static void compress_full_dtz(int stm, struct tb_handle *G, uint8_t *pcs,
    uint8_t *pt)
{
  uint16_t v[MAX_STATS];
  uint8_t w8[256];
  uint16_t w16[MAX_STATS];
  void *vv;

  sort_values(stm, g_stats[stm], &dtzmap, 441 * ri.sizes[0] + 21 * ri.sizes[1]);
  prepare_dtz_map(v, &dtzmap);
  if (!dtz_wide[stm]) {
    w8[0] = dtzmap.max_num;
    for (int i = 1; i < 256; i++)
      w8[i] = v[val_to_stats[stm][i]];
    vv = w8;
  } else {
    w16[0] = dtzmap.max_num;
    for (int i = 1; i < MAX_STATS; i++)
      w16[i] = v[val_to_stats[stm][i]];
    vv = w16;
  }
  compress_init_dtz(&dtzmap, dtz_wide[stm]);

  uint8_t best[MAX_PIECES];
  printf("Find optimal permutation for %ctm/dtz.\n", "wb"[stm]);
  permute_piece_dtz(tb_table, pcs, pt, dtz_table, best, dtz_wide[stm], vv);
  printf("Compressing data for %ctm/dtz.\n", "wb"[stm]);
  compress_tb(G, -1, tb_table, tb_size, best, minfreq, dtz_wide[stm]);
}

#if 0
static void join_full_dtz(int stm, uint8_t *pcs, uint8_t *pt)
{
  init_permute_piece(pcs, pt);

  compress_alloc_wdl();

  struct tb_handle *G = create_tb(g_tablename, WDL, 6);

  if (!one_sided || one_sided_stm == WHITE)
    compress_full_dtz(WHITE, G, pcs, pt);

  if (   (one_sided && one_sided_stm == BLACK)
      || (!one_sided && !symmetric))
    compress_full_dtz(BLACK, G, pcs, pt);

  merge_tb(G);
  compress_free_wdl();

  free(tb_table);
}
#endif

static void join_final_462(int type)
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

  if (type == WDL) {
    memcpy(p, &wdl_checksum, sizeof wdl_checksum);
    p += sizeof wdl_checksum;
  } else if (type == DTZ) {
    memcpy(p, &dtz_checksum, sizeof dtz_checksum);
    p += sizeof dtz_checksum;
  }

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
  add_xxhash(tmp);
  file_rename(fname);
  free(tmp);
  free(fname);
}

static void compress_462_wdl(int stm)
{
  work_init(&work_stats, ri.sizes[0], 0, WORK_DYNAMIC, 16, 512);

  create_dir(-1, stm, "wdl");
  for (int s = 0; s < 462; s++) {
    g_slice.sq[0] = KKSquare[s][0];
    g_slice.sq[1] = KKSquare[s][1];
    g_slice.stm = stm;

    uint8_t *restrict table = g_table[stm] + s * ri.sizes[0];

    uint64_t stats[MAX_STATS] = { 0 };
    wdl_collect_stats(table, stats);
    bool has_capt_bloss = find_capt_bloss(table);
    bool vals[5];
    prepare_wdl_map(stats, vals, has_capt_bloss);

    compress_init_wdl(vals);

    uint8_t best[MAX_SETS];
    printf("Find optimal permutation for %ctm/wdl, slice = %d.\n", "wb"[stm],
        s);
    permute_piece_462(tb_table, table, best, WDL, false, ram_wdl_map);
    printf("Compressing data for %ctm/wdl, slice = %d.\n", "wb"[stm], s);

    char name[64];
    create_name(name, s, stm, "wdl", -1);
    compress_data_slice(name, stm, WDL, tb_table, kslice_sizes[s >= 441], best,
        minfreq, false, false);
  }
}

// This should be invoked with dtz_table pointing to an u8 or u16
// reconstructed table. The table's values should be stats[] values when
// mapped through mi[stm].v_inv_u8/u16[].
static void compress_462_dtz(int stm)
{
  uint16_t v[MAX_STATS];
  uint8_t w8[256];
  uint16_t w16[MAX_STATS];

  work_init(&work_stats, ri.sizes[0], 0, WORK_DYNAMIC, 16, 512);

  create_dir(-1, stm, "dtz");

  for (int s = 0; s < 462; s++) {
    g_slice.sq[0] = KKSquare[s][0];
    g_slice.sq[1] = KKSquare[s][1];
    g_slice.stm = stm;

    void *tbl;
    void *vv;

    if (!dtz_wide[stm]) {

      uint8_t *table = (uint8_t *)dtz_table + s * ri.sizes[0];

      uint64_t stats[MAX_STATS] = { 0 };
      dtz_collect_stats_u8(stm, table, stats);
      sort_values(stm, stats, &dtzmap, ri.sizes[0]);
      prepare_dtz_map(v, &dtzmap);
      compress_init_dtz(&dtzmap, false);

      w8[0] = dtzmap.max_num;
      for (int i = 1; i < 256; i++)
        w8[i] = v[val_to_stats[stm][i]];
      tbl = table;
      vv = w8;

    } else {

      uint16_t *table = (uint16_t *)dtz_table + s * ri.sizes[0];

      uint64_t stats[MAX_STATS] = { 0 };
      dtz_collect_stats_u16(stm, table, stats);
      sort_values(stm, stats, &dtzmap, ri.sizes[0]);
      prepare_dtz_map(v, &dtzmap);
      compress_init_dtz(&dtzmap, true);

      w16[0] = dtzmap.max_num;
      for (int i = 1; i < MAX_STATS; i++)
        w16[i] = v[val_to_stats[stm][i]];
      tbl = table;
      vv = w16;

    }

    char name[64];
    create_name(name, s, stm, "dtz", -1);

    uint8_t best[MAX_SETS];
    printf("Find optimal permutation for %ctm/dtz, slice = %d.\n", "wb"[stm],
        s);
    permute_piece_462(tb_table, tbl, best, DTZ, dtz_wide[stm], vv);
    printf("Compressing data for %ctm/dtz, slice = %d.\n", "wb"[stm], s);
    compress_data_slice(name, stm, DTZ, tb_table, ri.sizes[s >= 441], best,
        minfreq, dtz_wide[stm], false);
  }
}

static void join_462_wdl(void)
{
  init_permute_piece_462();

  tb_table = alloc_huge(kslice_size + 1);
  if (!tb_table)
    out_of_mem();

  compress_alloc_wdl();
  compress_462_wdl(WHITE);
  if (!symmetric)
    compress_462_wdl(BLACK);
  compress_free_wdl();
  join_final_462(WDL);

  free(tb_table);
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

  if (type == WDL) {
    memcpy(p, &wdl_checksum, sizeof wdl_checksum);
    p += sizeof wdl_checksum;
  } else if (type == DTZ) {
    memcpy(p, &dtz_checksum, sizeof dtz_checksum);
    p += sizeof dtz_checksum;
  }

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
  add_xxhash(tmp);
  file_rename(fname);
  free(tmp);
  free(fname);
}

static void compress_10_wdl(int stm)
{
  work_init(&work_stats, ri.sizes[0], 0, WORK_DYNAMIC, 16, 512);

  create_dir(-1, stm, "wdl");
  g_slice.stm = stm;
  for (int k = 0; k < 10; k++) {
    init_permute_piece_10(k);

    int num = 0;
    uint64_t stats[MAX_STATS] = { 0 };
    bool has_capt_bloss = false;

    g_slice.sq[stm] = InvTriangle[k];
    for (int l = 0; l < 64; l++) {
      if (KKIdx[k][l] < 0)
        continue;
      g_slice.sq[stm ^ 1] = l;
      int s = KKMap[g_slice.sq[0]][g_slice.sq[1]];
      uint8_t *restrict table = g_table[stm] + s * ri.sizes[0];
      wdl_collect_stats(table, stats);
      has_capt_bloss = has_capt_bloss || find_capt_bloss(table);
      num++;
    }

    bool vals[5];
    prepare_wdl_map(stats, vals, has_capt_bloss);
    compress_init_wdl(vals);

    uint8_t best[MAX_SETS];
    printf("Find optimal permutation for %ctm/wdl, slice = %d.\n", "wb"[stm],
        k);
    permute_piece_10(tb_table, g_table[stm], best, WDL, false, ram_wdl_map);
    printf("Compressing data for %ctm/wdl, slice = %d.\n", "wb"[stm], k);

    char name[64];
    create_name_10(name, k, stm, "wdl");
    compress_data_slice(name, stm, WDL, tb_table, num * ri.sizes[0], best,
        minfreq, false, true);
  }
}

static void compress_10_dtz(int stm)
{
  uint16_t v[MAX_STATS];
  uint8_t w8[256];
  uint16_t w16[MAX_STATS];
  void *vv;

  work_init(&work_stats, ri.sizes[0], 0, WORK_DYNAMIC, 16, 512);

  create_dir(-1, stm, "dtz");

  g_slice.stm = stm;
  for (int k = 0; k < 10; k++) {
    init_permute_piece_10(k);

    int num = 0;
    uint64_t stats[MAX_STATS] = { 0 };
    g_slice.sq[stm] = InvTriangle[k];

    if (!dtz_wide[stm]) {

      for (int l = 0; l < 64; l++) {
        if (KKIdx[k][l] < 0)
          continue;
        g_slice.sq[stm ^ 1] = l;
        int s = KKMap[g_slice.sq[0]][g_slice.sq[1]];
        uint8_t *restrict table = (uint8_t *)dtz_table + s * ri.sizes[0];
        uint64_t slice_stats[MAX_STATS] = { 0 };
        dtz_collect_stats_u8(stm, table, slice_stats);
        for (int i = 0; i < MAX_STATS; i++)
          stats[i] += slice_stats[i] * (1 + (s >= 441));
        num++;
      }

      sort_values(stm, stats, &dtzmap, num * ri.sizes[0]);
      prepare_dtz_map(v, &dtzmap);
      compress_init_dtz(&dtzmap, false);

      w8[0] = dtzmap.max_num;
      for (int i = 1; i < 256; i++)
        w8[i] = v[val_to_stats[stm][i]];
      vv = w8;

    } else {

      for (int l = 0; l < 64; l++) {
        if (KKIdx[k][l] < 0)
          continue;
        g_slice.sq[stm ^ 1] = l;
        int s = KKMap[g_slice.sq[0]][g_slice.sq[1]];
        uint16_t *restrict table = (uint16_t *)dtz_table + s * ri.sizes[0];
        uint64_t slice_stats[MAX_STATS] = { 0 };
        dtz_collect_stats_u16(stm, table, slice_stats);
        for (int i = 0; i < MAX_STATS; i++)
          stats[i] += slice_stats[i] * (1 + (s >= 441));
        num++;
      }

      sort_values(stm, stats, &dtzmap, num * ri.sizes[0]);
      prepare_dtz_map(v, &dtzmap);
      compress_init_dtz(&dtzmap, true);

      w16[0] = dtzmap.max_num;
      for (int i = 1; i < MAX_STATS; i++)
        w16[i] = v[val_to_stats[stm][i]];
      vv = w16;

    }

    char name[64];
    create_name_10(name, k, stm, "dtz");

    uint8_t best[MAX_SETS];
    printf("Find optimal permutation for %ctm/dtz, slice = %d.\n", "wb"[stm],
        k);
    permute_piece_10(tb_table, dtz_table, best, DTZ, dtz_wide[stm], vv);
    printf("Compressing data for %ctm/dtz, slice = %d.\n", "wb"[stm], k);
    compress_data_slice(name, stm, DTZ, tb_table, num * ri.sizes[0], best,
        minfreq, dtz_wide[stm], true);
  }
}


static void join_10_wdl(void)
{
  tb_table = alloc_huge(58 * kslice_size + 1);
  if (!tb_table)
    out_of_mem();

  compress_alloc_wdl();
  compress_10_wdl(WHITE);
  if (!symmetric)
    compress_10_wdl(BLACK);
  compress_free_wdl();
  join_final_10(WDL);

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
}

static struct tb_handle *dtz_G;
static uint8_t *dtz_pcs;
static uint8_t *dtz_pt;

static void join_dtz_side(int stm, int layout)
{
  if (dtz_wide[stm]) {
    dtz_table = alloc_huge(2 * (table_size + 1));
    if (!dtz_table)
      out_of_mem();
    reconstruct_table(stm, g_table[stm], dtz_table);
    free(g_table[stm]);
    g_table[stm] = nullptr;
  } else {
    reconstruct_table(stm, g_table[stm], g_table[stm]);
    dtz_table = g_table[stm];
  }
  compress_alloc_dtz(dtz_wide[stm]);
  switch (layout) {
  case 0:
    work_init(&work_stats, table_size, 0, WORK_DYNAMIC, 16, 512);
    compress_full_dtz(stm, dtz_G, dtz_pcs, dtz_pt);
    break;
  case 1:
    compress_10_dtz(stm);
    break;
  case 2:
    compress_462_dtz(stm);
    break;
  }
  compress_free_dtz();
  if (dtz_wide[stm]) {
    free(dtz_table);
    dtz_table = nullptr;
  } else {
    free(g_table[stm]);
    g_table[stm] = nullptr;
  }
}

void rjoin_dtz(uint8_t *pcs, uint8_t *pt, int layout)
{
  dtz_pcs = pcs;
  dtz_pt = pt;
  switch (layout) {
  case 0:
    init_permute_piece(pcs, pt);
    tb_table = alloc_huge(2 * (tb_size + 1));
    break;
  case 1:
    tb_table = alloc_huge(2 * (58 * kslice_size + 1));
    break;
  case 2:
    init_permute_piece_462();
    tb_table = alloc_huge(2 * (kslice_size + 1));
    break;
  default:
    unreachable();
  }
  if (!tb_table)
    out_of_mem();
  if (layout == 0)
    dtz_G = create_tb(g_tablename, DTZ, 10);

  bool dtz_w = !one_sided || one_sided_stm == WHITE;
  bool dtz_b = one_sided ? one_sided_stm == BLACK : !symmetric;

  if (dtz_w)
    create_stats_to_val(WHITE);
  else
    free(g_table[WHITE]);
  if (dtz_b)
    create_stats_to_val(BLACK);
  else
    free(g_table[BLACK]);
  if (dtz_w && dtz_b && dtz_wide[WHITE] && !dtz_wide[BLACK]) {
    join_dtz_side(BLACK, layout);
    join_dtz_side(WHITE, layout);
  } else {
    if (dtz_w) {
      if (dtz_b && dtz_wide[BLACK]) {
        // save_table(BLACK);
        join_dtz_side(WHITE, layout);
        // load_table(BLACK);
      } else {
        join_dtz_side(WHITE, layout);
      }
    }
    if (dtz_b)
      join_dtz_side(BLACK, layout);
  }

  switch (layout) {
  case 0:
    merge_tb(dtz_G);
    break;
  case 1:
    join_final_10(DTZ);
    break;
  case 2:
    join_final_462(DTZ);
    break;
  }
  free(tb_table);
}
