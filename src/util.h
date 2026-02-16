#ifndef UTIL_H
#define UTIL_H

#include <stdbit.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#endif

#include "types.h"

#ifdef _WIN32
typedef HANDLE map_t;
typedef HANDLE FD;
#define FD_ERR INVALID_HANDLE_VALUE
#define SEP_STR ";"

#else
typedef size_t map_t;
typedef int FD;
#define FD_ERR -1
#define SEP_STR ":"

#endif

INLINE bool is_little_endian(void)
{
  return __STDC_ENDIAN_NATIVE__ == __STDC_ENDIAN_LITTLE__;
}

INLINE uint64_t from_le_u64(uint64_t v)
{
  return is_little_endian() ? v : __builtin_bswap64(v);
}

INLINE uint32_t from_le_u32(uint32_t v)
{
  return is_little_endian() ? v : __builtin_bswap32(v);
}

INLINE uint16_t from_le_u16(uint16_t v)
{
  return is_little_endian() ? v : __builtin_bswap16(v);
}

INLINE uint64_t from_be_u64(uint64_t v)
{
  return is_little_endian() ? __builtin_bswap64(v) : v;
}

INLINE uint32_t from_be_u32(uint32_t v)
{
  return is_little_endian() ? __builtin_bswap32(v) : v;
}

INLINE uint32_t to_be_u32(uint32_t v)
{
  return is_little_endian() ? __builtin_bswap32(v) : v;
}

INLINE uint32_t to_le_u32(uint32_t v)
{
  return is_little_endian() ? v : __builtin_bswap32(v);
}

INLINE uint64_t to_le_u64(uint64_t v)
{
  return is_little_endian() ? v : __builtin_bswap64(v);
}

INLINE uint64_t read_be_u64(const void *p)
{
  uint64_t v;
  memcpy(&v, p, sizeof v);
  return from_be_u64(v);
}

INLINE uint64_t read_be_u32(const void *p)
{
  uint32_t v;
  memcpy(&v, p, sizeof v);
  return from_be_u32(v);
}

INLINE void write_be_u32(void *p, uint32_t v)
{
  uint32_t w = to_be_u32(v);
  memcpy(p, &w, sizeof w);
}

INLINE void write_le_u64(void *p, uint64_t v)
{
  uint64_t w = to_le_u64(v);
  memcpy(p, &w, sizeof w);
}

INLINE void write_le_u32(void *p, uint32_t v)
{
  uint32_t w = to_le_u32(v);
  memcpy(p, &w, sizeof w);
}

INLINE uint64_t read_le_u64(const void *p)
{
  uint64_t v;
  memcpy(&v, p, sizeof v);
  return from_le_u64(v);
}

INLINE uint32_t read_le_u32(const void *p)
{
  uint32_t v;
  memcpy(&v, p, sizeof v);
  return from_le_u32(v);
}

INLINE uint16_t read_le_u16(const void *p)
{
  uint16_t v;
  memcpy(&v, p, sizeof v);
  return from_le_u16(v);
}

INLINE void write_u32(FILE *F, uint32_t v)
{
  fputc(v & 0xff, F);
  fputc((v >> 8) & 0xff, F);
  fputc((v >> 16) & 0xff, F);
  fputc((v >> 24) & 0xff, F);
}

INLINE void write_u16(FILE *F, uint16_t v)
{
  fputc(v & 0xff, F);
  fputc((v >> 8) & 0xff, F);
}

INLINE void write_u8(FILE *F, uint8_t v)
{
  fputc(v, F);
}

[[noreturn]] void out_of_mem(void);

FD open_file(const char *name);
void close_file(FD fd);

size_t file_size(FD fd);

void *map_file(int fd, bool shared, map_t *map);
void unmap_file(const void *data, map_t map);

void *alloc_aligned(uint64_t size, uintptr_t alignment);
void *alloc_huge(uint64_t size);

void write_u32(FILE *F, uint32_t v);
void write_u16(FILE *F, uint16_t v);
void write_u8(FILE *F, uint8_t v);

void file_read(void *ptr, size_t size, FILE *F);
void file_write(void *ptr, size_t size, FILE *F);

//void write_bits(FILE *F, uint32_t bits, int n);

void copy_data(FILE *F, FILE *G, size_t num);
void write_data(FILE *F, uint8_t *src, size_t size);
void read_data(FILE *F, uint8_t *dst, size_t size);

#endif
