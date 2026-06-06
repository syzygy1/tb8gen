/*
  Copyright (c) 2026 Ronald de Man

  This file is distributed under the terms of the GNU GPL, version 2.
*/

#define _POSIX_C_SOURCE 200809L
#include <assert.h>
#include <errno.h>
#include <inttypes.h>
#include <stdio.h>
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
size_t kslice_sizes[2];
size_t kslice_sub_sizes[2][MAX_SETS];
size_t sub_offset[MAX_SETS];
size_t sub_size[2];
int8_t kslice_slot[463];
static uint64_t kslice_cache_lines[2];
uint64_t kslice_read_cost;
uint64_t kslice_read_count;

static struct Work *work_cl[2];
static struct Work *work_sub_cl[2]; // FIXME x2?
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
    iter->todo = ~king_mask(InvTriangle[iter->k1]);
    if (iter->k1 >= 6)
      iter->todo &= LOWER;
  }

  int k2 = pop_lsb(&iter->todo);
  Bitboard b =  king_mask(k2) & ~king_mask(InvTriangle[iter->k1])
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
  size_t size = bits_to_aligned(kslice_sizes[0]);
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
  kslice_cache_lines[0] = bits_to_aligned(kslice_sizes[0]) >> 6;
  kslice_cache_lines[1] = bits_to_aligned(kslice_sizes[1]) >> 6;
  work_cl[0] = create_work(calc_work_units(kslice_cache_lines[0], 1, 16),
      kslice_cache_lines[0], 0);
  work_cl[0]->schedule = WORK_STATIC;
  work_cl[1] = create_work(calc_work_units(kslice_cache_lines[1], 1, 16),
      kslice_cache_lines[1], 0);
  work_cl[1]->schedule = WORK_STATIC;
  sub_size[0] = sub_size[1] = 0;
  for (int i = 0; i < ri.numsets; i++) {
    int stm = g_pos.pt[ri.first[i]] >> 3;
    sub_offset[i] = sub_size[stm];
    sub_size[stm] += bits_to_aligned(kslice_sub_size[i]);
  }
  uint64_t size = max(sub_size[WHITE], sub_size[BLACK]);
  for (int i = 0; i < 19; i++) {
    kslice_sub_buf[i] = alloc_huge(size);
    if (!kslice_sub_buf[i])
      out_of_mem();
  }
  work_sub_cl[WHITE] = create_work(calc_work_units(sub_size[WHITE] >> 6, 1, 16),
      sub_size[WHITE] >> 6, 0);
  work_sub_cl[BLACK] = create_work(calc_work_units(sub_size[BLACK] >> 6, 1, 16),
      sub_size[BLACK] >> 6, 0);
  work_sub_cl[WHITE]->schedule = WORK_STATIC;
  work_sub_cl[BLACK]->schedule = WORK_STATIC;
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

static void run_cl(void (*func)(struct ThreadData *), int s)
{
  run_threaded(func, work_cl[s >= 441]);
}

static void run_sub_cl(void (*func)(struct ThreadData *), int stm)
{
  run_threaded(func, work_sub_cl[stm]);
}

static void set_worker(struct ThreadData *thread)
{
  uint8_t *restrict p = work_p;

  memset(p + (thread->begin << 6), 0xff, (thread->end - thread->begin) << 6);
}

void kslice_set_addr(void *p, int s)
{
  work_p = p;
  run_cl(set_worker, s);
}

void kslice_set(int s)
{
  work_p = kslice_get_address(s);
  run_cl(set_worker, s);
}

static void clear_worker(struct ThreadData *thread)
{
  uint8_t *restrict p = work_p;

  memset(p + (thread->begin << 6), 0x00, (thread->end - thread->begin) << 6);
}

void kslice_clear_addr(void *p, int s)
{
  work_p = p;
  run_cl(clear_worker, s);
}

void kslice_clear(int s)
{
  work_p = kslice_get_address(s);
  run_cl(clear_worker, s);
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

  run_cl(or_worker, max(s1, s2));
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

  run_cl(or_not_worker, max(s1, s2));
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

  run_cl(and_worker, max(s1, s2));
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

  run_cl(and_not_worker, max(s1, s2));
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

  run_cl(not_and_worker, max(s1, s2));
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

  run_cl(nor_worker, max(s1, s2));
}

uint64_t kslice_write_addr(void *p, int slice, int stm, const char *name, int n,
    uint64_t num)
{
  char str[64];
  create_name(str, slice, stm, name, n);
  FILE *F = file_open_write(str);
  if (num > 0) {
    kslice_clear_tail_addr(p, slice);
    file_write(&num, sizeof num, F);
    write_data_cache_aligned(F, p, kslice_cache_lines[slice >= 441] << 6);
  }
  uint64_t bytes = ftell(F);
  fclose(F);
  file_rename(str);
  return bytes;
}

uint64_t kslice_write(int s, int slice, int stm, const char *name, int n,
    uint64_t num)
{
  return kslice_write_addr(kslice_get_address(s), slice, stm, name, n, num);
}

bool kslice_test(int slice, int stm, const char *name, int n)
{
  char str[64];
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
  char str[64];
  create_name(str, slice, stm, name, n);

  FILE *F = fopen(str, "rb");
  if (!F && errno != ENOENT) {
    fprintf(stderr, "Could not open %s for reading.\n", str);
    exit(EXIT_FAILURE);
  }
  bool non_empty = false;
  kslice_read_cost = 0;
  if (F) {
    struct stat st;
    fstat(fileno(F), &st);
    kslice_read_cost = st.st_size;
    non_empty = st.st_size != 0;
  }
  if (non_empty) {
    file_read(&kslice_read_count, sizeof(uint64_t), F);
    read_data(F, kslice_get_address(s), kslice_cache_lines[slice >= 441] << 6);
    kslice_clear_tail(s);
  } else
    kslice_clear(s);
  if (F) fclose(F);
  return non_empty;
}

void kslice_read_or(int s, int slice, int stm, const char *name, int n)
{
  char str[64];
  create_name(str, slice, stm, name, n);

  FILE *F = fopen(str, "rb");
  if (!F && errno != ENOENT) {
    fprintf(stderr, "Could not open %s for reading.\n", str);
    exit(EXIT_FAILURE);
  }
  if (!F)
    return;
  struct stat st;
  fstat(fileno(F), &st);
  if (st.st_size > 0) {
    uint64_t num;
    file_read(&num, sizeof(uint64_t), F);
//    kslice_read_count += num;
    kslice_read_cost += st.st_size;
    read_data_or(F, kslice_get_address(s),
        kslice_cache_lines[slice >= 441] << 6);
    kslice_clear_tail(s);
  }
  fclose(F);
}

void kslice_delete(int slice, int stm, const char *name, int n)
{
  char str[64];
  create_name(str, slice, stm, name, n);
  remove(str);
}

void kslice_sub_clear_addr(void *p, int stm)
{
  if (sub_size[stm] == 0) return;

  work_p = p;
  run_sub_cl(clear_worker, stm);
}

void kslice_sub_clear(int s, int stm)
{
  kslice_sub_clear_addr(kslice_sub_get_base(s), stm);
}

void kslice_sub_write_addr(void *p, int slice, int stm, const char *name,
    uint64_t num)
{
  char str[64];
  create_name(str, slice, stm, name, -1);
  FILE *F = file_open_write(str);
  if (num > 0) {
    file_write(&num, sizeof num, F);
    write_data_cache_aligned(F, p, sub_size[stm]);
  }
  fclose(F);
  file_rename(str);
}

void kslice_sub_read(int s, int slice, int stm, const char *name)
{
  char str[64];
  create_name(str, slice, stm, name, -1);
  FILE *F = file_open_read(str);
  struct stat st;
  fstat(fileno(F), &st);
  if (st.st_size > 0) {
    file_read(&kslice_read_count, 8, F);
    read_data(F, kslice_sub_get_base(s), sub_size[stm]);
  } else
    kslice_sub_clear(s, stm);
  fclose(F);
}

bool kslice_test_count(int s, int stm, const char *name, int n, uint64_t *num)
{
  char str[64];
  create_name(str, s, stm, name, n);
  FILE *F = fopen(str, "rb");
  if (!F)
    return false;
  struct stat st;
  fstat(fileno(F), &st);
  if (st.st_size > 0)
    file_read(num, 8, F);
  else
    *num = 0;
  fclose(F);
  return true;
}

uint64_t kslice_size_count(int s, int stm, const char *name, int n,
    uint64_t *num)
{
  char str[64];
  create_name(str, s, stm, name, n);
  FILE *F = file_open_read(str);
  struct stat st;
  fstat(fileno(F), &st);
  if (st.st_size > 0)
    file_read(num, 8, F);
  else
    *num = 0;
  fclose(F);
  return st.st_size;
}

void kslice_sub_and_not(int s1, int s2, int stm)
{
  work_p = kslice_sub_get_base(s1);
  work_q = kslice_sub_get_base(s2);

  run_sub_cl(and_not_worker, stm);
}

INLINE void clear_tail(void *p, size_t num_bits, size_t num_words)
{
  uint64_t *restrict q = p;
  size_t idx = num_bits >> 6;
  int r = num_bits & 63;
  if (r)
    q[idx++] &= (1ULL << r) - 1;

  for (; idx < num_words; idx++)
    q[idx] = 0;
}

void kslice_clear_tail_addr(void *p, int s)
{
  uint8_t *restrict q = p;
  clear_tail(q, kslice_sizes[s >= 441], kslice_cache_lines[s >= 441] << 3);
}

void kslice_clear_tail(int s)
{
  kslice_clear_tail_addr(kslice_get_address(s), s);
}

static void count_worker(struct ThreadData *thread)
{
  uint64_t cnt = 0, *restrict p = work_p;

  for (uint64_t idx = thread->begin << 3, end = thread->end << 3; idx < end;
      idx++)
    cnt += popcnt(p[idx]);

  thread->cnt += cnt;
}

uint64_t kslice_count_addr(void *p, int s)
{
  kslice_clear_tail_addr(p, s);

  work_p = p;

  for (int t = 0; t < g_num_threads; t++)
    g_thread_data[t].cnt = 0;

  run_cl(count_worker, s);

  uint64_t cnt = 0;
  for (int t = 0; t < g_num_threads; t++)
    cnt += g_thread_data[t].cnt;

  return cnt;
}

uint64_t kslice_count(int s)
{
  return kslice_count_addr(kslice_get_address(s), s);
}

uint64_t kslice_sub_count_addr(void *p, int stm)
{
  if (sub_size[stm] == 0)
    return 0;

  work_p = p;

  uint64_t cnt = 0;
  for (int t = 0; t < g_num_threads; t++)
    g_thread_data[t].cnt = 0;
  run_sub_cl(count_worker, stm);
  for (int t = 0; t < g_num_threads; t++)
    cnt += g_thread_data[t].cnt;

  return cnt;
}
