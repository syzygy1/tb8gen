/*
  Copyright (c) 2026 Ronald de Man

  This file is distributed under the terms of the GNU GPL, version 2.
*/

#include <assert.h>
#include <inttypes.h>
#include <string.h>
#include <stdlib.h>

#include "defs.h"
#include "movegen.h"
#include "kslice.h"
#include "probe.h"
#include "tb8gen.h"
#include "threads.h"
#include "types.h"
#include "util.h"

struct KSliceManager manager[2][462];

uint8_t *kslice_buf[20];
uint8_t *kslice_sub_buf[19];
size_t sub_offset[MAX_SETS];
static bool kslice_in_use[19];
//static bool kslice_dirty[19];
int8_t kslice_slot[463];
uint64_t kslice_cache_lines;
size_t sub_size[2];
static uint64_t *work_cl, *work_clc;
static uint64_t *work_sub_cl[2];

static int flip(int s)
{
  int wk = KKSquare[s][0];
  int bk = KKSquare[s][1];
  return KKMap[bk][wk];
}

static void init_kslice_manager(void)
{
  int j, s = 0;
  Bitboard loaded = 0;
  for (int i = 0; i < 10; i++) {
    for (int bk = 0; bk < 64; bk++) {
      if (KKIdx[i][bk] < 0) continue;
      manager[BLACK][s].kslice = KKIdx[i][bk];
      Bitboard b = (bit(bk) | king_attacks(bk)) & ~loaded;
      j = 0;
      while (b) {
        int sq = pop_lsb(&b);
        if (KKIdx[i][sq] >= 0) {
          manager[BLACK][s].in[j++] = KKIdx[i][sq];
          loaded |= bit(sq);
          assert(j < 11);
        }
      }
      manager[BLACK][s].in[j] = -1;
      b = bk == 63 ? loaded : bk >= 9 ? (loaded & ((1ULL << (bk - 8)) - 1)) : 0;
      j = 0;
      while (b) {
        int sq = pop_lsb(&b);
        manager[BLACK][s].out[j++] = KKIdx[i][sq];
        loaded ^= bit(sq);
        assert(j < 11);
      }
      manager[BLACK][s].out[j] = -1;
      s++;
    }
    assert(!loaded);
  }
  assert(s == 462);

  for (int s = 0; s < 462; s++) {
    manager[WHITE][s].kslice = flip(manager[BLACK][s].kslice);
    for (j = 0; manager[BLACK][s].in[j] >= 0; j++)
      manager[WHITE][s].in[j] = flip(manager[BLACK][s].in[j]);
    manager[WHITE][s].in[j] = -1;
    for (j = 0; manager[BLACK][s].out[j] >= 0; j++)
      manager[WHITE][s].out[j] = flip(manager[BLACK][s].out[j]);
    manager[WHITE][s].out[j] = -1;
  }
}

struct KSliceManager *kslice_get_manager(int stm, int i)
{
  return &manager[stm][i];
}

// Convert number of bits to number of bytes rounded up to cache lines.
INLINE size_t bits_to_aligned(size_t size)
{
  size = (size + 7) >> 3;
  return (size + 0x3f) & ~0x3f;
}

void kslice_setup(void)
{
  init_kslice_manager();
  size_t size = bits_to_aligned(kslice_size);
  for (int i = 0; i < 20; i++) {
    kslice_buf[i] = alloc_huge(size);
    if (!kslice_buf[i])
      out_of_mem();
  }
  kslice_slot[0] = 19;
  for (int i = 0; i < 462; i++)
    kslice_slot[i + 1] = -1;
  kslice_cache_lines = size >> 6;
  work_cl = create_work(g_total_work, kslice_cache_lines, 0);
  work_clc = create_work(g_total_work, kslice_cache_lines - 1, 0);
  sub_size[0] = sub_size[1] = 0;
  for (int i = 0; i < ii.numsets; i++) {
    int stm = g_pos.pt[ii.first[i]] >> 3;
    sub_offset[i] = sub_size[stm];
    sub_size[stm] += bits_to_aligned(kslice_sub_size[i]);
  }
  size = max(sub_size[WHITE], sub_size[BLACK]);
  for (int i = 0; i < 19; i++) {
    kslice_sub_buf[i] = alloc_huge(size);
    if (!kslice_sub_buf[i])
      out_of_mem();
  }
  work_sub_cl[WHITE] = create_work(g_total_work, sub_size[WHITE] >> 6, 0);
  work_sub_cl[BLACK] = create_work(g_total_work, sub_size[BLACK] >> 6, 0);
}

void kslice_cleanup(void)
{
  for (int i = 0; i < 20; i++)
    if (kslice_buf[i])
      free(kslice_buf[i]);
  for (int i = 0; i < 19; i++)
    if (kslice_sub_buf[i])
      free(kslice_sub_buf[i]);
}

void kslice_reserve(int s)
{
  assert(kslice_slot[s + 1] < 0);
  for (int i = 0; i < 19; i++) {
    if (kslice_in_use[i]) continue;
    kslice_in_use[i] = true;
    kslice_slot[s + 1] = i;
    return;
  }
  assert(false);
}

void kslice_release(int s)
{
  assert(kslice_slot[s + 1] >= 0);
  kslice_in_use[kslice_slot[s + 1]] = false;
  kslice_slot[s + 1] = -1;
}

static void *work_p, *work_q;

static void set_worker(struct ThreadData *thread)
{
  uint8_t *restrict p = work_p;

  memset(p + (thread->begin << 6), 0xff, (thread->end - thread->begin) << 6);
}

void kslice_set_addr(void *p)
{
  work_p = p;
  run_threaded(set_worker, work_cl, 0);
}

void kslice_set(int s)
{
  work_p = kslice_get_address(s);
  run_threaded(set_worker, work_cl, 0);
}

static void clear_worker(struct ThreadData *thread)
{
  uint8_t *restrict p = work_p;

  memset(p + (thread->begin << 6), 0x00, (thread->end - thread->begin) << 6);
}

void kslice_clear_addr(void *p)
{
  work_p = p;
  run_threaded(clear_worker, work_cl, 0);
}

void kslice_clear(int s)
{
  work_p = kslice_get_address(s);
  run_threaded(clear_worker, work_cl, 0);
}

static void or_worker(struct ThreadData *thread)
{
  uint64_t *restrict p = work_p;
  uint64_t *restrict q = work_q;

  for (uint64_t idx = thread->begin << 3, end = thread->end << 3; idx < end;
      idx++)
    p[idx] |= q[idx];
}

void kslice_or(int s1, int s2)
{
  work_p = kslice_get_address(s1);
  work_q = kslice_get_address(s2);

  run_threaded(or_worker, work_cl, 0);
}

static void or_not_worker(struct ThreadData *thread)
{
  uint64_t *restrict p = work_p;
  uint64_t *restrict q = work_q;

  for (uint64_t idx = thread->begin << 3, end = thread->end << 3; idx < end;
      idx++)
    p[idx] |= ~q[idx];
}

void kslice_or_not(int s1, int s2)
{
  work_p = kslice_get_address(s1);
  work_q = kslice_get_address(s2);

  run_threaded(or_not_worker, work_cl, 0);
}

static void and_worker(struct ThreadData *thread)
{
  uint64_t *restrict p = work_p;
  uint64_t *restrict q = work_q;

  for (uint64_t idx = thread->begin << 3, end = thread->end << 3; idx < end;
      idx++)
    p[idx] &= q[idx];
}

void kslice_and(int s1, int s2)
{
  work_p = kslice_get_address(s1);
  work_q = kslice_get_address(s2);

  run_threaded(and_worker, work_cl, 0);
}

static void and_not_worker(struct ThreadData *thread)
{
  uint64_t *restrict p = work_p;
  uint64_t *restrict q = work_q;

  for (uint64_t idx = thread->begin << 3, end = thread->end << 3; idx < end;
      idx++)
    p[idx] &= ~q[idx];
}

void kslice_and_not(int s1, int s2)
{
  work_p = kslice_get_address(s1);
  work_q = kslice_get_address(s2);

  run_threaded(and_not_worker, work_cl, 0);
}

static void not_and_worker(struct ThreadData *thread)
{
  uint64_t *restrict p = work_p;
  uint64_t *restrict q = work_q;

  for (uint64_t idx = thread->begin << 3, end = thread->end << 3; idx < end;
      idx++)
    p[idx] = ~p[idx] & q[idx];
}

void kslice_not_and(int s1, int s2)
{
  work_p = kslice_get_address(s1);
  work_q = kslice_get_address(s2);

  run_threaded(not_and_worker, work_cl, 0);
}

void nor_worker(struct ThreadData *thread)
{
  uint64_t *restrict p = work_p;
  uint64_t *restrict q = work_q;

  for (uint64_t idx = thread->begin << 3, end = thread->end << 3; idx < end;
      idx++)
    p[idx] = ~(p[idx] | q[idx]);
}

void kslice_nor(int s1, int s2)
{
  work_p = kslice_get_address(s1);
  work_q = kslice_get_address(s2);

  run_threaded(nor_worker, work_cl, 0);
}

static void create_name(char *str, int s, int stm, const char *name, int n)
{
  int wk = KKSquare[s][0], bk = KKSquare[s][1];
  sprintf(str, "%s_%s.%c%c%c%c.%c.%d", g_tablename, name, 'a' + (wk & 7),
      '1' + (wk >> 3), 'a'+ (bk & 7), '1' + (bk >> 3), "wb"[stm], n);
}

void kslice_write_addr(void *p, int slice, int stm, const char *name, int n)
{
  char str[128];
  create_name(str, slice, stm, name, n);
  FILE *F = fopen(str, "wb");
  if (!F) {
    fprintf(stderr, "Could not open %s for writing.\n", str);
    exit(EXIT_FAILURE);
  }
  write_data(F, p, kslice_cache_lines << 6);
  fclose(F);
}

void kslice_write(int s, int slice, int stm, const char *name, int n)
{
  kslice_write_addr(kslice_get_address(s), slice, stm, name, n);
}

void kslice_read(int s, int slice, int stm, const char *name, int n)
{
  char str[128];
  create_name(str, slice, stm, name, n);
  FILE *F = fopen(str, "rb");
  if (!F) {
    fprintf(stderr, "Could not open %s for reading.\n", str);
    exit(EXIT_FAILURE);
  }
  read_data(F, kslice_get_address(s), kslice_cache_lines << 6);
  fclose(F);
}

void kslice_delete(int slice, int stm, const char *name, int n)
{
  char str[128];
  create_name(str, slice, stm, name, n);
  remove(str);
}

void kslice_sub_write_addr(void *p, int slice, int stm, const char *name)
{
  char str[128];
  create_name(str, slice, stm, name, 0);
  FILE *F = fopen(str, "wb");
  if (!F) {
    fprintf(stderr, "Could not open %s for writing.\n", str);
    exit(EXIT_FAILURE);
  }
  write_data(F, p, sub_size[stm]);
  fclose(F);
}

void kslice_sub_read(int s, int slice, int stm, const char *name)
{
  char str[128];
  create_name(str, slice, stm, name, 0);
  FILE *F = fopen(str, "rb");
  if (!F) {
    fprintf(stderr, "Could not open %s for reading.\n", str);
    exit(EXIT_FAILURE);
  }
  read_data(F, kslice_sub_get_base(s), sub_size[stm]);
  fclose(F);
}

void kslice_sub_and_not(int s1, int s2, int stm)
{
  work_p = kslice_sub_get_base(s1);
  work_q = kslice_sub_get_base(s2);

  run_threaded(and_not_worker, work_sub_cl[stm], 0);
}

static void count_worker(struct ThreadData *thread)
{
  uint64_t cnt = 0, *restrict p = work_p;

  for (uint64_t idx = thread->begin << 3, end = thread->end << 3; idx < end;
      idx++)
    cnt += popcnt(p[idx]);

  thread->cnt += cnt;
}

uint64_t kslice_count_addr(void *p)
{
  work_p = p;

  for (int t = 0; t < g_num_threads; t++)
    g_thread_data[t].cnt = 0;

  run_threaded(count_worker, work_clc, 0);

  // Count 1s in the last cache line up to kslice_size.
  uint64_t *restrict q = p;
  uint64_t cnt = 0, idx = (kslice_cache_lines - 1) << 3;
  for (; idx < (kslice_size >> 6); idx++)
    cnt += popcnt(q[idx]);
  uint64_t last = q[idx];
  last &= (1ULL << (kslice_size & 0x3f)) - 1;
  cnt += popcnt(last);

  for (int t = 0; t < g_num_threads; t++)
    cnt += g_thread_data[t].cnt;

  return cnt;
}

uint64_t kslice_count(int s)
{
  return kslice_count_addr(kslice_get_address(s));
}
