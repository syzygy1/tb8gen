/*
  Copyright (c) 2026 Ronald de Man

  This file is distributed under the terms of the GNU GPL, version 2.
*/

#ifndef KSLICE_H
#define KSLICE_H

#include <assert.h>
#include <inttypes.h>
#include <stdatomic.h>

#include "types.h"

struct KSliceManager {
  int16_t kslice;
  int16_t in[11];
  int16_t out[11];
};

extern uint8_t *kslice_buf[20];
extern int8_t kslice_slot[463];

extern size_t kslice_size;
extern int16_t KKMap[64][64];

extern uint64_t kslice_cache_lines;

// FIXME: make atomic
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

// FIXME: make atomic
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

struct KSliceManager *kslice_get_manager(int stm, int i);
void kslice_setup(void);
void kslice_cleanup(void);
uint8_t *get_kslice_address(int s);
void kslice_reserve(int s);
void kslice_release(int s);
void kslice_set(int s);
void kslice_clear(int s);
void kslice_or(int s1, int s2);
void kslice_or_not(int s1, int s2);
void kslice_and(int s1, int s2);
void kslice_and_not(int s1, int s2);
void kslice_not_and(int s1, int s2);
void kslice_write(int s, int slice, int stm, const char *name, int n);
void kslice_write_addr(void *p, int slice, int stm, const char *name, int n);
void kslice_read(int s, int slice, int stm, const char *name, int n);
void kslice_delete(int slice, int stm, const char *name, int n);
uint64_t kslice_count(int s);
uint64_t kslice_count_addr(void *p);

#endif
