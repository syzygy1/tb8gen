/*
  Copyright (c) 2026 Ronald de Man

  This file is distributed under the terms of the GNU GPL, version 2.
*/

#define _POSIX_C_SOURCE 200809L
#include <assert.h>
#include <errno.h>
#include <inttypes.h>
#include <stdbit.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#include "defs.h"
#include "index.h"
#include "movegen.h"
#include "kslicep.h"
#include "probe.h"
#include "tb8gen.h"
#include "threads.h"
#include "types.h"
#include "util.h"

#include "kslice_common.c"

uint8_t KK16Square[240][16][2];
uint8_t kk_to_slice[256];
uint8_t *k16slice_buf[12];
uint8_t *k16slice_sub_buf[11];
size_t kslice_size, kslice_alloc_size;
size_t k16slice_alloc_size;
size_t kslice_sub_size[MAX_SETS], kslice_sub_alloc_size[MAX_SETS];
size_t sub_offset[MAX_SETS], psub_offset[MAX_SETS];
size_t sub_size[2];
int8_t k16slice_slot[241];
uint64_t kslice_cache_lines, k16slice_cache_lines;
uint64_t k16slice_read_count[16];
uint64_t k16slice_read_cost;

static struct Work *work_cl;
static struct Work *work_cl16;
static struct Work *work_sub_cl[2];
static bool k16slice_in_use[11];

static uint16_t kmask[16];

// Each bit represents 2x2 squares.
static uint16_t around_bb(uint16_t b)
{
  static constexpr uint16_t FILE_AB = 0x1111;
  static constexpr uint16_t FILE_GH = 0x8888;

  uint16_t east  = (b & ~FILE_GH) << 1;
  uint16_t west  = (b & ~FILE_AB) >> 1;
  uint16_t north = b << 4;
  uint16_t south = b >> 4;

  uint16_t ne = (b & ~FILE_GH) << 5;
  uint16_t nw = (b & ~FILE_AB) << 3;
  uint16_t se = (b & ~FILE_GH) >> 3;
  uint16_t sw = (b & ~FILE_AB) >> 5;

  return b | east | west | north | south | ne | nw | se | sw;
}

INLINE int pop_lsb_u16(uint16_t *v)
{
  int s = stdc_trailing_zeros(*v);
  *v &= *v - 1;
  return s;
}

static int k16slice(int stm, int k1, int k2)
{
  return stm == WHITE ? kk_to_slice[(k1 << 4) | k2]
                      : kk_to_slice[(k2 << 4) | k1];
}

void k16slice_iter_init(struct K16SliceIterator *iter, int stm)
{
  *iter = (struct K16SliceIterator){ .todo = 0, .stm = stm, .k1 = -1 };
}

bool k16slice_iter_next(struct K16SliceIterator *iter, int *s)
{
  if (!iter->todo) {
    if (++iter->k1 == 16)
      return false;
    iter->reserved = 0;
    iter->todo = ~(1 << iter->k1);
  }

  int k2 = pop_lsb_u16(&iter->todo);
  iter->in_slices = kmask[k2] & ~(1 << iter->k1) & ~iter->reserved;
  *s = k16slice(iter->stm, iter->k1, k2);
  iter->releasing = false;
  return true;
}

bool k16slice_iter_in(struct K16SliceIterator *iter, int *in)
{
  if (!iter->in_slices) {
    return false;
  }

  int k2 = pop_lsb_u16(&iter->in_slices);
  iter->reserved |= 1 << k2;
  int s = k16slice(iter->stm, iter->k1, k2);
  k16slice_reserve(s);
  *in = s;
  return true;
}

bool k16slice_iter_out(struct K16SliceIterator *iter, int *out)
{
  if (!iter->releasing) {
    iter->out_slices = iter->reserved & ~around_bb(iter->todo);
    iter->releasing = true;
  } else {
    k16slice_release(iter->release_slice);
  }

  if (!iter->out_slices)
    return false;

  int k2 = pop_lsb_u16(&iter->out_slices);
  iter->reserved ^= 1 << k2;
  *out = iter->release_slice = k16slice(iter->stm, iter->k1, k2);
  return true;
}

// Convert number of bits to number of bytes rounded up to cache lines.
INLINE size_t bits_to_aligned(size_t size)
{
  size = (size + 7) >> 3;
  return (size + 0x3f) & ~0x3f;
}

uint8_t *alloc_k16slice(void)
{
  uint8_t *p = alloc_huge(k16slice_alloc_size);
  if (!p)
    out_of_mem();
  return p;
}

void kslice_alloc_buffers(void)
{
  for (int i = 0; i < 12; i++)
    if (!k16slice_buf[i])
      k16slice_buf[i] = alloc_k16slice();
  uint64_t size = max(sub_size[WHITE], sub_size[BLACK]);
  for (int i = 0; i < 11; i++)
    if (!k16slice_sub_buf[i]) {
      k16slice_sub_buf[i] = alloc_huge(size);
      if (!k16slice_sub_buf[i])
        out_of_mem();
    }
}

void kslice_free_buffers(void)
{
  // We keep k16slice_buf[11] (slice -1) around.
  for (int i = 0; i < 11; i++) {
    if (k16slice_buf[i]) {
      free(k16slice_buf[i]);
      k16slice_buf[i] = nullptr;
    }
    if (k16slice_sub_buf[i]) {
      free(k16slice_sub_buf[i]);
      k16slice_sub_buf[i] = nullptr;
    }
  }
}

void kslice_setup(void)
{
  int s = 0;
  for (int s1 = 0; s1 < 16; s1++)
    for (int s2 = 0; s2 < 16; s2++) {
      if (s2 == s1) continue;
      uint32_t v = _pdep_u32((s2 << 4) | s1, 0b110110110110);
      for (int r = 0; r < 16; r++) {
        uint32_t w = v | _pdep_u32(r, 0b001001001001);
        KK16Square[s][r][0] = w & 0x3f;
        KK16Square[s][r][1] = w >> 6;
        kk_to_slice[(s2 << 4) + s1] = s;
      }
      s++;
    }

  kslice_size = ri.sizes[0];
  kslice_alloc_size = bits_to_aligned(kslice_size);
  k16slice_alloc_size = 16 * bits_to_aligned(kslice_size);
  for (int i = 0; i < 16; i++)
    kmask[i] = around_bb(1 << i);
  kslice_alloc_buffers();
  k16slice_slot[0] = 11;
  for (int i = 0; i < 240; i++)
    k16slice_slot[i + 1] = -1;
  kslice_cache_lines = bits_to_aligned(kslice_size) >> 6;
  k16slice_cache_lines = 16 * kslice_cache_lines;
  work_cl = create_work(g_total_work, kslice_cache_lines, 0);
  work_cl16 = create_work(g_total_work, k16slice_cache_lines, 0);

  sub_size[0] = sub_size[1] = 0;
  for (int i = 0; i < ri.numsets; i++) {
    int stm = g_set_pt[i] >> 3;
    sub_offset[i] = sub_size[stm];
    kslice_sub_alloc_size[i] = bits_to_aligned(kslice_sub_size[i]);
    sub_size[stm] += 16 * kslice_sub_alloc_size[i];
  }
  // FIXME: it seems psub_offset[i] == sub_offset[i]...
  size_t psize = 0;
  for (int i = 0; i < ri.numsets; i++)
    if ((g_set_pt[i] >> 3) == WHITE) {
      psub_offset[i] = psize; // sub_size[BLACK];
      psize += 16 * kslice_sub_alloc_size[i];
//      sub_size[BLACK] += 16 * kslice_sub_alloc_size[i];
    }
  uint64_t size = max(sub_size[WHITE], sub_size[BLACK]);
  for (int i = 0; i < 11; i++) {
    k16slice_sub_buf[i] = alloc_huge(size);
    if (!k16slice_sub_buf[i])
      out_of_mem();
  }
  work_sub_cl[WHITE] = create_work(g_total_work, sub_size[WHITE] >> 6, 0);
  work_sub_cl[BLACK] = create_work(g_total_work, sub_size[BLACK] >> 6, 0);
}

void kslice_cleanup(void)
{
  for (int i = 0; i < 12; i++)
    if (k16slice_buf[i]) {
      free(k16slice_buf[i]);
      k16slice_buf[i] = nullptr;
    }
  for (int i = 0; i < 11; i++)
    if (k16slice_sub_buf[i]) {
      free(k16slice_sub_buf[i]);
      k16slice_sub_buf[i] = nullptr;
    }
}

void k16slice_reserve(int s)
{
  assert(k16slice_slot[s + 1] < 0);
  for (int i = 0; i < 11; i++) {
    if (k16slice_in_use[i]) continue;
    k16slice_in_use[i] = true;
    k16slice_slot[s + 1] = i;
    return;
  }
  assert(false);
}

void k16slice_release(int s)
{
  assert(k16slice_slot[s + 1] >= 0);
  k16slice_in_use[k16slice_slot[s + 1]] = false;
  k16slice_slot[s + 1] = -1;
}

static void *work_p, *work_q;

void kslice_set_addr(void *p)
{
  work_p = p;
  run_threaded(set_worker, work_cl);
}

void k16slice_set_addr(void *p)
{
  work_p = p;
  run_threaded(set_worker, work_cl16);
}

void k16slice_set(int s)
{
  work_p = k16slice_get_address(s);
  run_threaded(set_worker, work_cl16);
}

void k16slice_clear_addr(void *p)
{
  work_p = p;
  run_threaded(clear_worker, work_cl16);
}

void k16slice_clear(int s)
{
  work_p = k16slice_get_address(s);
  run_threaded(clear_worker, work_cl16);
}

void k16slice_not_addr(void *p)
{
  work_p = p;
  run_threaded(not_worker, work_cl16);
}

void k16slice_not(int s)
{
  work_p = k16slice_get_address(s);
  run_threaded(not_worker, work_cl16);
}

void k16slice_or_addr(void *p, void *q)
{
  work_p = p;
  work_q = q;

  run_threaded(or_worker, work_cl16);
}

void k16slice_or(int s1, int s2)
{
  work_p = k16slice_get_address(s1);
  work_q = k16slice_get_address(s2);

  run_threaded(or_worker, work_cl16);
}

void k16slice_ornot(int s1, int s2)
{
  work_p = k16slice_get_address(s1);
  work_q = k16slice_get_address(s2);

  run_threaded(ornot_worker, work_cl16);
}

void k16slice_and(int s1, int s2)
{
  work_p = k16slice_get_address(s1);
  work_q = k16slice_get_address(s2);

  run_threaded(and_worker, work_cl16);
}

void k16slice_andnot(int s1, int s2)
{
  work_p = k16slice_get_address(s1);
  work_q = k16slice_get_address(s2);

  run_threaded(andnot_worker, work_cl16);
}

void k16slice_notand(int s1, int s2)
{
  work_p = k16slice_get_address(s1);
  work_q = k16slice_get_address(s2);

  run_threaded(notand_worker, work_cl16);
}

void k16slice_nor(int s1, int s2)
{
  work_p = k16slice_get_address(s1);
  work_q = k16slice_get_address(s2);

  run_threaded(nor_worker, work_cl16);
}

uint64_t k16slice_write_addr(void *p, int slice, int stm, const char *name,
    int n, uint64_t num[16])
{
  char str[64];
  create_name(str, slice, stm, name, n);
  FILE *F = file_open_write(str);
  uint64_t cnt = 0;
  if (num) {
    for (int r = 0; r < 16; r++)
      cnt += num[r];
    if (cnt > 0) {
      file_write(num, 16 * 8, F);
      write_data_cache_aligned(F, p, k16slice_cache_lines << 6);
    }
  }
  else {
    uint64_t tmp[16] = { 0 };
    file_write(tmp, 16 * 8, F);
    write_data_cache_aligned(F, p, k16slice_cache_lines << 6);
  }
  uint64_t bytes = ftell(F);
  fclose(F);
  file_rename(str);
  return bytes;
}

uint64_t k16slice_write(int s, int slice, int stm, const char *name, int n,
    uint64_t num[16])
{
  return k16slice_write_addr(k16slice_get_address(s), slice, stm, name, n, num);
}

bool k16slice_test(int slice, int stm, const char *name, int n)
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

bool k16slice_test_count(int s, int stm, const char *name, int n,
    uint64_t num[16])
{
  char str[64];
  create_name(str, s, stm, name, n);
  FILE *F = fopen(str, "rb");
  if (!F)
    return false;
  struct stat st;
  fstat(fileno(F), &st);
  if (st.st_size > 0)
    file_read(num, 16 * 8, F);
  else
    memset(num, 0, 16 * 8);
  fclose(F);
  return true;
}

bool k16slice_read(int s, int slice, int stm, const char *name, int n)
{
  char str[64];
  create_name(str, slice, stm, name, n);

  FILE *F = fopen(str, "rb");
  if (!F && errno != ENOENT) {
    fprintf(stderr, "Could not open %s for reading.\n", str);
    exit(EXIT_FAILURE);
  }
  bool non_empty = false;
  if (F) {
    struct stat st;
    fstat(fileno(F), &st);
    k16slice_read_cost = st.st_size;
    non_empty = st.st_size != 0;
  }
  if (non_empty) {
    file_read(k16slice_read_count, sizeof k16slice_read_count, F);
    read_data(F, k16slice_get_address(s), k16slice_cache_lines << 6);
  } else {
    for (int r = 0; r < 16; r++)
      k16slice_read_count[r] = 0;
    k16slice_clear(s);
  }
  if (F) fclose(F);
  return non_empty;
}

void k16slice_read_or(int s, int slice, int stm, const char *name, int n)
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
    if (fseek(F, 16 * 8, SEEK_SET) != 0) {
      fprintf(stderr, "fseek() failed.\n");
      exit(EXIT_FAILURE);
    }
    k16slice_read_cost += st.st_size;
    read_data_or(F, k16slice_get_address(s), k16slice_cache_lines << 6);
  }
  fclose(F);
}

void k16slice_read_andnot(int s, int slice, int stm, const char *name, int n)
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
    if (fseek(F, 16 * 8, SEEK_SET) != 0) {
      fprintf(stderr, "fseek() failed.\n");
      exit(EXIT_FAILURE);
    }
    k16slice_read_cost += st.st_size;
    read_data_andnot(F, k16slice_get_address(s), k16slice_cache_lines << 6);
  }
  fclose(F);
}

void k16slice_delete(int slice, int stm, const char *name, int n)
{
  char str[64];
  create_name(str, slice, stm, name, n);
  remove(str);
}

void k16slice_sub_clear_addr(void *p, int stm)
{
  if (sub_size[stm] == 0) return;

  work_p = p;
  run_threaded(clear_worker, work_sub_cl[stm]);
}

void k16slice_sub_clear(int s, int stm)
{
  if (sub_size[stm] == 0) return;

  work_p = k16slice_sub_get_base(s);
  run_threaded(clear_worker, work_sub_cl[stm]);
}

void k16slice_sub_write_addr(void *p, int slice, int stm, const char *name,
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

void k16slice_sub_read(int s, int slice, int stm, const char *name)
{
  char str[64];
  create_name(str, slice, stm, name, -1);
  FILE *F = file_open_read(str);
  struct stat st;
  fstat(fileno(F), &st);
  if (st.st_size > 0) {
    uint64_t dummy;
    file_read(&dummy, 8, F);
    read_data(F, k16slice_sub_get_base(s), sub_size[stm]);
  } else
    k16slice_sub_clear(s, stm);
  fclose(F);
}

bool k16slice_sub_test_count(int s, int stm, const char *name, int n,
    uint64_t *num)
{
  char str[64];
  create_name(str, s, stm, name, n);
  FILE *F = fopen(str, "rb");
  if (!F)
    return false;
  struct stat st;
  fstat(fileno(F), &st);
  if (st.st_size > 0)
    file_read(num, sizeof(uint64_t), F);
  else
    *num = 0;
  fclose(F);
  return true;
}

void k16slice_sub_or_addr(void *p, void *q, int stm)
{
  if (sub_size[stm] == 0) return;

  work_p = p;
  work_q = q;
  run_threaded(or_worker, work_sub_cl[stm]);
}

void k16slice_sub_andnot(int s1, int s2, int stm)
{
  if (sub_size[stm] == 0) return;

  work_p = k16slice_sub_get_base(s1);
  work_q = k16slice_sub_get_base(s2);
  run_threaded(andnot_worker, work_sub_cl[stm]);
}

void k16slice_clear_tail_addr(void *p)
{
  uint8_t *restrict q = p;
  for (int r = 0; r < 16; r++, q += kslice_alloc_size)
    clear_tail(q, kslice_size, kslice_cache_lines << 3);
}

void k16slice_clear_tail(int s)
{
  k16slice_clear_tail_addr(k16slice_get_address(s));
}

static uint64_t kslice_count_addr(void *p)
{
  work_p = p;

  for (int t = 0; t < g_num_threads; t++)
    g_thread_data[t].cnt = 0;

  run_threaded(count_worker, work_cl);

  uint64_t cnt = 0;
  for (int t = 0; t < g_num_threads; t++)
    cnt += g_thread_data[t].cnt;

  return cnt;
}

uint64_t k16slice_count_addr(void *p, uint64_t count[16])
{
  k16slice_clear_tail_addr(p);

  uint64_t cnt = 0;
  uint8_t *q = p;
  for (int r = 0; r < 16; r++) {
    cnt += count[r] = kslice_count_addr(q);
    q += kslice_alloc_size;
  }
  return cnt;
}

uint64_t k16slice_count(int s, uint64_t count[16])
{
  return k16slice_count_addr(k16slice_get_address(s), count);
}

uint64_t k16slice_sub_count_addr(void *p, int stm)
{
  if (sub_size[stm] == 0)
    return 0;

  work_p = p;

  uint64_t cnt = 0;
  for (int t = 0; t < g_num_threads; t++)
    g_thread_data[t].cnt = 0;
  run_threaded(count_worker, work_sub_cl[stm]);
  for (int t = 0; t < g_num_threads; t++)
    cnt += g_thread_data[t].cnt;

  return cnt;
}
