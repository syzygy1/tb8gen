/*
  Copyright (c) 2026 Ronald de Man

  This file is distributed under the terms of the GNU GPL, version 2.
*/

#include <assert.h>
#include <errno.h>
#include <inttypes.h>
#include <string.h>
#include <stdlib.h>
#include <sys/stat.h>

#include "defs.h"
#include "index.h"
#include "movegen.h"
#include "kslice.h"
#include "probe.h"
#include "tb8gen.h"
#include "threads.h"
#include "types.h"
#include "util.h"

uint8_t *kslice_buf[20];
uint8_t *kslice_sub_buf[19];
size_t kslice_size, kslice_sub_size[MAX_SETS];
size_t sub_offset[MAX_SETS];
size_t sub_size[2];
int8_t kslice_slot[463];
uint64_t kslice_cache_lines;
uint64_t kslice_read_count;

static uint64_t *work_cl, *work_clc;
static uint64_t *work_sub_cl[2];
static bool kslice_in_use[19];

static constexpr Bitboard LOWER  = 0x80c0e0f0f8fcfeffull;

static Bitboard around_bb(Bitboard b)
{
  static constexpr Bitboard FILE_A = 0x0101010101010101ull;
  static constexpr Bitboard FILE_H = 0x8080808080808080ull;

  Bitboard east  = (b & ~FILE_H) << 1;
  Bitboard west  = (b & ~FILE_A) >> 1;
  Bitboard north = b << 8;
  Bitboard south = b >> 8;

  Bitboard ne = (b & ~FILE_H) << 9;
  Bitboard nw = (b & ~FILE_A) << 7;
  Bitboard se = (b & ~FILE_H) >> 7;
  Bitboard sw = (b & ~FILE_A) >> 9;

  return b | east | west | north | south | ne | nw | se | sw;
}

static int kslice(int stm, int k1, int k2)
{
  return stm == WHITE ? KKMap[k2][InvTriangle[k1]]
                      : KKMap[InvTriangle[k1]][k2];
}

void kslice_iter_init(struct KSliceIterator *iter, int stm)
{
  *iter = (struct KSliceIterator){ .todo = 0, .stm = stm, .k1 = -1 };
}

bool kslice_iter_next(struct KSliceIterator *iter, int *s)
{
  if (!iter->todo) {
    if (++iter->k1 == 10)
      return false;
    iter->reserved = 0;
    iter->todo = ~KingMask[InvTriangle[iter->k1]];
    if (iter->k1 >= 6)
      iter->todo &= LOWER;
  }

  int k2 = pop_lsb(&iter->todo);
  Bitboard b =  KingMask[k2] & ~KingMask[InvTriangle[iter->k1]]
              & ~iter->reserved;
  iter->in_slices = iter->k1 < 6 ? b : b & LOWER;
  *s = kslice(iter->stm, iter->k1, k2);
  iter->releasing = false;
  return true;
}

bool kslice_iter_in(struct KSliceIterator *iter, int *in)
{
  if (!iter->in_slices) {
    return false;
  }

  int k2 = pop_lsb(&iter->in_slices);
  iter->reserved |= bit(k2);
  int s = kslice(iter->stm, iter->k1, k2);
  kslice_reserve(s);
  *in = s;
  return true;
}

bool kslice_iter_out(struct KSliceIterator *iter, int *out)
{
  if (!iter->releasing) {
    iter->out_slices = iter->reserved & ~around_bb(iter->todo);
    iter->releasing = true;
  } else {
    kslice_release(iter->release_slice);
  }

  if (!iter->out_slices)
    return false;

  int k2 = pop_lsb(&iter->out_slices);
  iter->reserved ^= bit(k2);
  *out = iter->release_slice = kslice(iter->stm, iter->k1, k2);
  return true;
}

// Convert number of bits to number of bytes rounded up to cache lines.
INLINE size_t bits_to_aligned(size_t size)
{
  size = (size + 7) >> 3;
  return (size + 0x3f) & ~0x3f;
}

uint8_t *alloc_kslice(void)
{
  size_t size = bits_to_aligned(kslice_size);
  uint8_t *p = alloc_huge(size);
  if (!p)
    out_of_mem();
  return p;
}

void kslice_setup(void)
{
  for (int i = 0; i < 20; i++)
    kslice_buf[i] = alloc_kslice();
  kslice_slot[0] = 19;
  for (int i = 0; i < 462; i++)
    kslice_slot[i + 1] = -1;
  kslice_cache_lines = bits_to_aligned(kslice_size) >> 6;
  work_cl = create_work(g_total_work, kslice_cache_lines, 0);
  work_clc = create_work(g_total_work, kslice_cache_lines - 1, 0);
  sub_size[0] = sub_size[1] = 0;
  for (int i = 0; i < ii.numsets; i++) {
    int stm = g_pos.pt[ii.first[i]] >> 3;
    sub_offset[i] = sub_size[stm];
    sub_size[stm] += bits_to_aligned(kslice_sub_size[i]);
  }
  uint64_t size = max(sub_size[WHITE], sub_size[BLACK]);
  for (int i = 0; i < 19; i++) {
    kslice_sub_buf[i] = alloc_huge(size);
    if (!kslice_sub_buf[i])
      out_of_mem();
  }
  work_sub_cl[WHITE] = create_work(g_total_work, sub_size[WHITE] >> 6, 0);
  work_sub_cl[BLACK] = create_work(g_total_work, sub_size[BLACK] >> 6, 0);
}

void kslice_free_buffers(void)
{
  for (int i = 0; i < 19; i++) {
    if (kslice_buf[i]) {
      free(kslice_buf[i]);
      kslice_buf[i] = nullptr;
    }
    if (kslice_sub_buf[i]) {
      free(kslice_sub_buf[i]);
      kslice_sub_buf[i] = nullptr;
    }
  }
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

void kslice_write_addr(void *p, int slice, int stm, const char *name, int n,
    uint64_t num)
{
  char str[128];
  create_name(str, slice, stm, name, n);
  FILE *F = fopen(str, "wb");
  if (!F) {
    fprintf(stderr, "Could not open %s for writing.\n", str);
    exit(EXIT_FAILURE);
  }
  if (num > 0) {
    file_write(&num, sizeof(num), F);
    write_data(F, p, kslice_cache_lines << 6);
  }
  fclose(F);
}

void kslice_write(int s, int slice, int stm, const char *name, int n,
    uint64_t num)
{
  kslice_write_addr(kslice_get_address(s), slice, stm, name, n, num);
}

bool kslice_test(int slice, int stm, const char *name, int n)
{
  char str[128];
  create_name(str, slice, stm, name, n);
  struct stat st;
  if (stat(str, &st) < 0) {
    if (errno == ENOENT)
      return false;
    fprintf(stderr, "Error trying to access %s.\n", str);
    exit(EXIT_FAILURE);
  }
  return st.st_size != 0;
}

bool kslice_read(int s, int slice, int stm, const char *name, int n)
{
  char str[128];
  create_name(str, slice, stm, name, n);

  FILE *F = fopen(str, "rb");
  if (!F && errno != ENOENT) {
    fprintf(stderr, "Could not open %s for reading.\n", str);
    exit(EXIT_FAILURE);
  }
  bool non_empty = F != 0;
  if (F) {
    struct stat st;
    fstat(fileno(F), &st);
    non_empty = st.st_size != 0;
  }
  if (non_empty) {
    file_read(&kslice_read_count, sizeof(uint64_t), F);
    read_data(F, kslice_get_address(s), kslice_cache_lines << 6);
  } else
    kslice_clear(s);
  if (F) fclose(F);
  return non_empty;
}

void kslice_delete(int slice, int stm, const char *name, int n)
{
  char str[128];
  create_name(str, slice, stm, name, n);
  remove(str);
}

void kslice_sub_clear(int s, int stm)
{
  if (sub_size[stm] == 0) return;

  work_p = kslice_sub_get_base(s);
  run_threaded(clear_worker, work_sub_cl[stm], 0);
}

void kslice_sub_write_addr(void *p, int slice, int stm, const char *name,
    uint64_t num)
{
  char str[128];
  create_name(str, slice, stm, name, -1);
  FILE *F = fopen(str, "wb");
  if (!F) {
    fprintf(stderr, "Could not open %s for writing.\n", str);
    exit(EXIT_FAILURE);
  }
  if (num > 0)
    write_data(F, p, sub_size[stm]);
  fclose(F);
}

void kslice_sub_read(int s, int slice, int stm, const char *name)
{
  char str[128];
  create_name(str, slice, stm, name, -1);
  FILE *F = fopen(str, "rb");
  if (!F) {
    fprintf(stderr, "Could not open %s for reading.\n", str);
    exit(EXIT_FAILURE);
  }
  struct stat st;
  fstat(fileno(F), &st);
  if (st.st_size > 0)
    read_data(F, kslice_sub_get_base(s), sub_size[stm]);
  else
    kslice_sub_clear(s, stm);
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

uint64_t kslice_sub_count_addr(void *p, int stm)
{
  if (sub_size[stm] == 0)
    return 0;

  work_p = p;

  uint64_t cnt = 0;
  for (int t = 0; t < g_num_threads; t++)
    g_thread_data[t].cnt = 0;
  run_threaded(count_worker, work_sub_cl[stm], 0);
  for (int t = 0; t < g_num_threads; t++)
    cnt += g_thread_data[t].cnt;

  return cnt;
}
