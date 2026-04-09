/*
  Copyright (c) 2026 Ronald de Man

  This file is distributed under the terms of the GNU GPL, version 2.
*/

#ifndef KSLICEP_H
#define KSLICEP_H

#include <assert.h>
#include <inttypes.h>
#include <stdatomic.h>
#include <x86intrin.h>

#include "defs.h"
#include "types.h"

struct K16SliceIterator {
  uint16_t reserved, todo, in_slices, out_slices;
  int stm;
  int k1;
  int release_slice;
  bool releasing;
};

extern uint8_t KK16Square[240][16][2];
extern uint8_t kk_to_slice[256];

extern uint8_t *k16slice_buf[12];
extern uint8_t *k16slice_sub_buf[11];
extern size_t sub_offset[MAX_SETS], psub_offset[MAX_SETS];
extern int8_t k16slice_slot[241];

extern size_t kslice_size, kslice_sub_size[MAX_SETS];
extern size_t kslice_alloc_size, k16slice_alloc_size;
extern size_t kslice_sub_alloc_size[MAX_SETS];
extern size_t k16slice_cache_lines;

extern uint64_t k16slice_read_count[16];

void k16slice_iter_init(struct K16SliceIterator *iter, int stm);
bool k16slice_iter_next(struct K16SliceIterator *iter, int *s);
bool k16slice_iter_in(struct K16SliceIterator *iter, int *in);
bool k16slice_iter_out(struct K16SliceIterator *iter, int *out);

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

INLINE uint8_t *k16slice_get_address(int s)
{
  assert(k16slice_slot[s + 1] >= 0);
  return k16slice_buf[k16slice_slot[s + 1]];
}

INLINE uint8_t *k16slice_sub_get_base(int s)
{
  assert(k16slice_slot[s + 1] >= 0);
//  return kslice_sub_buf[kslice_slot[s + 1]];
  return s == -1 ? k16slice_buf[11] : k16slice_sub_buf[k16slice_slot[s + 1]];
}

INLINE uint8_t *k16slice_sub_get_address(int s, int set)
{
  return k16slice_sub_get_base(s) + sub_offset[set];
}

INLINE uint8_t *k16slice_psub_get_address(int s, int set)
{
  return k16slice_sub_get_base(s) + psub_offset[set];
}

INLINE size_t k16offset(const uint8_t *sq)
{
  return _pext_u32(sq[0] + (sq[1] << 8), 0x909) * kslice_alloc_size;
}

INLINE uint8_t *kslice_get_address(const uint8_t *sq)
{
  int s = kk_to_slice[_pext_u32(sq[0] + (sq[1] << 8), 0b11011000110110)];
  return k16slice_get_address(s) + k16offset(sq);
}

INLINE uint8_t *kslice_sub_get_address(const uint8_t *sq, int set)
{
  int s = kk_to_slice[_pext_u32(sq[0] + (sq[1] << 8), 0b11011000110110)];
  return k16slice_sub_get_address(s, set) + k16offset(sq);
}

INLINE uint8_t *kslice_psub_get_address(const uint8_t *sq, int set)
{
  int s = kk_to_slice[_pext_u32(sq[0] + (sq[1] << 8), 0b11011000110110)];
  return k16slice_psub_get_address(s, set) + k16offset(sq);
}

INLINE uint8_t *kslice_get_address_scratch(const uint8_t *sq)
{
  return k16slice_get_address(-1) + k16offset(sq);
}

void kslice_setup(void);
void kslice_cleanup(void);
void kslice_alloc_buffers(void);
void kslice_free_buffers(void);
uint8_t *alloc_k16slice(void);
void k16slice_reserve(int s);
void k16slice_release(int s);
void k16slice_set(int s);
void k16slice_set_addr(void *p);
void kslice_set_addr(void *p);
void k16slice_clear(int s);
void k16slice_clear_addr(void *p);
void k16slice_or(int s1, int s2);
void k16slice_or_addr(void *p, void *q);
void k16slice_or_not(int s1, int s2);
void k16slice_and(int s1, int s2);
void k16slice_and_not(int s1, int s2);
void k16slice_not_and(int s1, int s2);
void k16slice_nor(int s1, int s2);
void k16slice_write(int s, int slice, int stm, const char *name, int n,
    uint64_t num[16]);
void k16slice_write_addr(void *p, int slice, int stm, const char *name, int n,
    uint64_t num[16]);
bool k16slice_test(int slice, int stm, const char *name, int n);
bool k16slice_read(int s, int slice, int stm, const char *name, int n);
void k16slice_delete(int slice, int stm, const char *name, int n);
void k16slice_sub_write_addr(void *p, int slice, int stm, const char *name,
    uint64_t cnt);
void k16slice_sub_read(int s, int slice, int stm, const char *name);
void k16slice_sub_or_addr(void *p, void *q, int stm);
void k16slice_sub_and_not(int s1, int s2, int stm);
uint64_t k16slice_count(int s, uint64_t num[16]);
uint64_t k16slice_count_addr(void *p, uint64_t num[16]);
uint64_t k16slice_sub_count_addr(void *p, int stm);
uint64_t k16slice_psub_count_addr(void *p);
void k16slice_sub_clear_addr(void *p, int stm);
void k16slice_sub_clear(int s, int stm);

#endif
