#ifndef UTIL_H
#define UTIL_H

#include <stdbit.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

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

INLINE void write_u64(FILE *F, uint64_t v)
{
  write_u32(F, v);
  write_u32(F, v >> 32);
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

void write_data_fd(int fd, const void *data, ssize_t count);
void copy_data_fd(int fd_in, int fd_out, ssize_t count);

void make_dir(const char *pathname);
void change_dir(const char *pathname);

uint64_t llrand(void);

void *alloc_aligned(uint64_t size, uintptr_t alignment);
void *alloc_huge(uint64_t size);

void file_read(void *ptr, size_t size, FILE *F);
void file_write(void *ptr, size_t size, FILE *F);
void report_io(void);

FILE *file_open_read(const char *name);
FILE *file_open_write(const char *name);
void file_rename(const char *name);
bool file_exists(const char *name);
bool dir_exists(int n, int stm, const char *name);
void create_empty(const char *name);

void copy_data(FILE *F, FILE *G, size_t num);
void write_data(FILE *F, void *src, size_t size);
void write_data_cache_aligned(FILE *F, void *src, size_t size);
void write_data_transform_u8(FILE *F, void *src, size_t size, uint8_t *v);
void write_data_transform_u16(FILE *F, void *src, size_t size,
    uint16_t *v);
void write_data_transform_to_u8_u8(FILE *F, void *src, size_t size,
    uint8_t *v);
void write_data_transform_to_u8_u16(FILE *F, void *src, size_t size,
    uint8_t *v);
void write_data_as_u8_u8(FILE *F, void *src, size_t size);
void write_data_as_u8_u16(FILE *F, void *src, size_t size);
void read_data(FILE *F, void *dst, size_t size);
void read_data_transform_u8(FILE *F, void *dst, size_t size, uint8_t *v);
void read_data_transform_u16(FILE *F, void *dst, size_t size, uint16_t *v);
void read_data_transform_to_u8_u16(FILE *F, void *dst, size_t size,
    uint8_t *v);
void read_data_or(FILE *F, void *dst, size_t size);
void read_data_and(FILE *F, void *dst, size_t size);
void read_data_andnot(FILE *F, void *dst, size_t size);
void read_data_transform_or_u8(FILE *F, void *dst, size_t size, uint8_t *v);
void read_data_transform_or_u8_to_u16(FILE *F, void *dst, size_t size,
    uint16_t *v);

void create_dir(int n, int stm, const char *name);
void delete_dir(int n, const char *name);
void create_name(char *str, int s, int stm, const char *name, int n);
void create_name_10(char *str, int k, int stm, const char *name);
#ifdef HAS_PAWNS
void create_name_r(char *str, int s, int r, int stm, const char *name, int n);
void create_name_p(char *s, int sq, int stm, const char *name);
void create_name_sq(char *s, int s1, int s2, int stm, const char *name, int n);
#endif

void show_progress(const char *str, int k, int total, bool final);

#endif
