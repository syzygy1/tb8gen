/*
  Copyright (c) 2026 Ronald de Man

  This file is distributed under the terms of the GNU GPL, version 2.
*/

#ifndef KSLICE_H
#define KSLICE_H

#include <assert.h>
#include <inttypes.h>
#include <stdatomic.h>

#include "defs.h"
#include "types.h"

struct KSliceIterator {
  Bitboard reserved, todo, in_slices, out_slices;
  int stm;
  int k1;
  int release_slice;
  bool releasing;
};

extern uint8_t *kslice_buf[20];
extern uint8_t *kslice_sub_buf[19];
extern size_t sub_offset[MAX_SETS];
extern int8_t kslice_slot[463];

extern size_t kslice_sizes[2];
#define kslice_size (kslice_sizes[0])
extern size_t kslice_sub_sizes[2][MAX_SETS];
#define kslice_sub_size (kslice_sub_sizes[0])
//extern size_t total_kslice_size;
//extern size_t reflection_offset;

extern uint64_t kslice_read_cost;
extern uint64_t kslice_read_count;

void kslice_iter_init(struct KSliceIterator *iter, int stm);
bool kslice_iter_next(struct KSliceIterator *iter, int *s);
bool kslice_iter_in(struct KSliceIterator *iter, int *in);
bool kslice_iter_out(struct KSliceIterator *iter, int *out);

INLINE void kslice_bit_flip(uint8_t *restrict p, uint64_t idx)
{
  p[idx >> 3] ^= 1 << (idx & 7);
}

INLINE void kslice_bit_set(uint8_t *restrict p, uint64_t idx)
{
  p[idx >> 3] |= 1 << (idx & 7);
}

INLINE void kslice_bit_set_atomic(uint8_t *restrict p, uint64_t idx)
{
  _Atomic uint8_t *restrict q = (_Atomic uint8_t *)p;

  atomic_fetch_or_explicit(&q[idx >> 3], 1 << (idx & 7), memory_order_relaxed);
}

INLINE void kslice_bit_clear(uint8_t *restrict p, uint64_t idx)
{
  p[idx >> 3] &= ~(1 << (idx & 7));
}

INLINE bool kslice_bit_test(uint8_t *restrict p, uint64_t idx)
{
  return p[idx >> 3] & (1 << (idx & 7));
}

INLINE uint8_t *kslice_get_address(int s)
{
  assert(kslice_slot[s + 1] >= 0);
  return kslice_buf[kslice_slot[s + 1]];
}

INLINE uint8_t *kslice_sub_get_base(int s)
{
  assert(kslice_slot[s + 1] >= 0);
//  return kslice_sub_buf[kslice_slot[s + 1]];
  return s == -1 ? kslice_buf[19] : kslice_sub_buf[kslice_slot[s + 1]];
}

INLINE uint8_t *kslice_sub_get_address(int s, int set)
{
  return kslice_sub_get_base(s) + sub_offset[set];
}

void kslice_setup(void);
void kslice_cleanup(void);
void kslice_alloc_buffers(int n);
void kslice_free_buffers(void);
void kslice_reserve(int s);
void kslice_release(int s);
void kslice_set(int s);
void kslice_set_addr(void *p, int s);
void kslice_clear(int s);
void kslice_clear_addr(void *p, int s);
void kslice_not(int s);
void kslice_not_addr(void *p, int s);
void kslice_or(int s1, int s2);
void kslice_or_addr(void *p, void *q, int s);
void kslice_or_not(int s1, int s2);
void kslice_and(int s1, int s2);
void kslice_and_not(int s1, int s2);
void kslice_and_not_addr(void *p, void *q, int s);
void kslice_not_and(int s1, int s2);
void kslice_nor(int s1, int s2);
void kslice_split_addr(void *p, void *q, int s);
uint64_t kslice_write(int s, int slice, int stm, const char *name, int n,
    uint64_t num);
uint64_t kslice_write_addr(void *p, int slice, int stm, const char *name, int n,
    uint64_t num);
bool kslice_test(int slice, int stm, const char *name, int n);
bool kslice_read(int s, int slice, int stm, const char *name, int n);
bool kslice_read_addr(void *p, int slice, int stm, const char *name, int n);
void kslice_read_or(int s, int slice, int stm, const char *name, int n);
void kslice_delete(int slice, int stm, const char *name, int n);
void kslice_sub_write_addr(void *p, int slice, int stm, const char *name,
    uint64_t cnt);
void kslice_sub_read(int s, int slice, int stm, const char *name);
void kslice_sub_and_not(int s1, int s2, int stm);
void kslice_clear_tail(int s);
void kslice_clear_tail_addr(void *p, int s);
uint64_t kslice_count(int s);
uint64_t kslice_count_addr(void *p, int s);
uint64_t kslice_sub_count_addr(void *p, int stm);
void kslice_sub_clear(int s, int stm);
void kslice_sub_clear_addr(void *p, int stm);
bool kslice_test_count(int s, int stm, const char *name, int n, uint64_t *num);
uint64_t kslice_size_count(int s, int stm, const char *name, int n,
    uint64_t *num);

#endif
