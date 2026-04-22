/*
  Copyright (c) 2011-2013, 2018, 2024, 2025, 2026 Ronald de Man

  This file is distributed under the terms of the GNU GPL, version 2.
*/

#define _GNU_SOURCE
#include <stdlib.h>

#include <assert.h>
#include <stdatomic.h>
#include <stddef.h>
#include <stdio.h>
#include <stdint.h>
#include <unistd.h>
#ifndef _WIN32
#include <unistd.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>
#include <fcntl.h>
#endif

#include <zstd.h>

#include "defs.h"
#include "threads.h"
#include "util.h"

#ifdef HAS_PAWNS
#include "kslicep.h"
#endif

static constexpr size_t HUGEPAGE_SIZE = 2 * 1024 * 1024;

FD open_file(const char *name)
{
#ifndef _WIN32
  return open(name, O_RDONLY);
#else
  return CreateFile(name, GENERIC_READ, FILE_SHARE_READ, nullptr, OPEN_EXISTING,
      FILE_FLAG_RANDOM_ACCESS, nullptr);
#endif
}

void close_file(FD fd)
{
#ifndef _WIN32
  close(fd);
#else
  CloseHandle(fd);
#endif
}

size_t file_size(FD fd)
{
#ifndef _WIN32
  struct stat statbuf;
  fstat(fd, &statbuf);
  return statbuf.st_size;
#else
  DWORD sizeLow, sizeHigh;
  sizeLow = GetFileSize(fd, &sizeHigh);
  return ((size_t)sizeHigh << 32) | sizeLow;
#endif
}

void *map_file(FD fd, bool shared, map_t *map)
{
#ifndef _WIN32
  *map = file_size(fd);
#ifdef __linux__
  void *data = mmap(nullptr, *map, PROT_READ,
      shared ? MAP_SHARED : MAP_PRIVATE | MAP_POPULATE, fd, 0);
#else
  void *data = mmap(nullptr, statbuf.st_size, PROT_READ, MAP_SHARED, fd, 0);
#endif
#ifdef MADV_RANDOM
  madvise(data, *map, MADV_RANDOM);
#endif
  if (data == MAP_FAILED) {
    fprintf(stderr, "mmap() failed.\n");
    exit(EXIT_FAILURE);
  }
  return data;

#else
  DWORD sizeLow, sizeHigh;
  sizeLow = GetFileSize(fd, &sizeHigh);
  *map = CreateFileMapping(fd, nullptr, PAGE_READONLY, sizeHigh, sizeLow,
      nullptr);
  if (!*map) {
    fprintf(stderr, "CreateFileMapping() failed.\n");
    exit(EXIT_FAILURE);
  }
  return MapViewOfFile(*map, FILE_MAP_READ, 0, 0, 0);

#endif
}

void unmap_file(const void *data, map_t map)
{
  if (!data) return;

#ifndef _WIN32
  munmap((void *)data, map);

#else
  UnmapViewOfFile(data);
  CloseHandle(map);

#endif
}

void write_data_fd(int fd, const void *data, ssize_t count)
{
  const uint8_t *ptr = data;
  while (count > 0) {
    ssize_t n = write(fd, ptr, count);
    if (n <= 0) {
      fprintf(stderr, "Error writing data.\n");
      exit(EXIT_FAILURE);
    }
    ptr += n;
    count -= n;
  }
}

void copy_data_fd(int fd_in, int fd_out, ssize_t count)
{
  while (count > 0) {
    ssize_t n = copy_file_range(fd_in, nullptr, fd_out, nullptr, count, 0);
    if (n <= 0) {
      fprintf(stderr, "Error copying data.\n");
      exit(EXIT_FAILURE);
    }
    count -= n;
  }
}

[[noreturn]] void out_of_mem(void)
{
  fprintf(stderr, "Could not allocate sufficient memory.\n");
  exit(EXIT_FAILURE);
}

void *alloc_aligned(uint64_t size, uintptr_t alignment)
{
  void *ptr;

#ifndef _WIN32
  posix_memalign(&ptr, alignment, size);
  if (!ptr) out_of_mem();

#else
  ptr = malloc(size + alignment - 1);
  if (!ptr) out_of_mem();
  ptr = (void *)((uintptr_t)(ptr + alignment - 1) & ~(alignment - 1));

#endif

  return ptr;
}

void *alloc_huge(uint64_t size)
{
  void *ptr;

#ifndef _WIN32
  posix_memalign(&ptr, HUGEPAGE_SIZE, size);

#else
  ptr = malloc(size);

#endif

  if (!ptr) out_of_mem();

#ifdef MADV_HUGEPAGE
  madvise(ptr, size, MADV_HUGEPAGE);
#endif

  return ptr;
}

void make_dir(const char *pathname)
{
  if (mkdir(pathname, 0755) == 0)
    return;

  if (errno == EEXIST) {
    struct stat st;
    if (stat(pathname, &st) == 0 && S_ISDIR(st.st_mode))
      return;
  }

  fprintf(stderr, "Could not create directory %s.\n", pathname);
  exit(EXIT_FAILURE);
}

void change_dir(const char *pathname)
{
  if (chdir(pathname) < 0) {
    fprintf(stderr, "Could not enter directory %s.\n", pathname);
    exit(EXIT_FAILURE);
  }
}

// adapted from glibc, random.c:
// TODO: replace with something better

static int myrand(void)
{
  static int tbl[31] = {
    -1726662223, 379960547, 1735697613, 1040273694, 1313901226,
    1627687941, -179304937, -2073333483, 1780058412, -1989503057,
    -615974602, 344556628, 939512070, -1249116260, 1507946756,
    -812545463, 154635395, 1388815473, -1926676823, 525320961,
    -1009028674, 968117788, -123449607, 1284210865, 435012392,
    -2017506339, -911064859, -370259173, 1132637927, 1398500161,
    -205601318
  };

  static int f = 3;
  static int r = 0;

  int result;

  tbl[f] = (unsigned)tbl[f] + (unsigned)tbl[r];
  result = (tbl[f] >> 1) & 0x7fffffff;

  f++;
  if (f >= 31) f = 0;
  r++;
  if (r >= 31) r= 0;

  return result;
}

uint64_t llrand(void)
{
  uint64_t rand1 = myrand(), rand2 = myrand();
  return (rand1 << 24) + rand2;
}

//static constexpr size_t COPYSIZE = 20*1024*1024ULL;
static constexpr size_t COPYSIZE = 50*1024ULL;

static size_t compress_bound;

static LOCK_T cmprs_mutex;

static FILE *cmprs_F;
static void *cmprs_ptr;
static size_t cmprs_size;
static size_t cmprs_idx;
static void *cmprs_v;
static int cmprs_type;

enum { COPY, U8U8, U16U16, U16U8, COPY_U16U8 };

struct CompressFrame {
  size_t cmprs_chunk;
  size_t chunk;
  size_t idx;
  uint8_t data[];
};

static constexpr size_t HEADER_SIZE = offsetof(struct CompressFrame, data);

struct CompressState {
  void *buffer;
  struct CompressFrame *frame;
  ZSTD_CCtx *c_ctx;
  ZSTD_DCtx *d_ctx;
};

static struct CompressState cmprs_state[COMPRESSION_THREADS];

static void init(void)
{
  static int initialised = 0;

  if (!initialised) {
    initialised = 1;
    LOCK_INIT(cmprs_mutex);
    compress_bound = ZSTD_compressBound(COPYSIZE);
    for (int i = 0; i < COMPRESSION_THREADS; i++) {
      cmprs_state[i].buffer = malloc(COPYSIZE);
      cmprs_state[i].frame = malloc(HEADER_SIZE + compress_bound);
      cmprs_state[i].c_ctx = ZSTD_createCCtx();
      cmprs_state[i].d_ctx = ZSTD_createDCtx();
    }
    create_compression_threads();
  }
}

static uint64_t total_read = 0, total_written = 0;

void file_read(void *ptr, size_t size, FILE *F)
{
  if (fread(ptr, 1, size, F) != size) {
    fprintf(stderr, "Error reading data.\n");
    exit(EXIT_FAILURE);
  }
  total_read += size;
}

void file_write(void *ptr, size_t size, FILE *F)
{
  if (fwrite(ptr, 1, size, F) != size) {
    fprintf(stderr, "Error writing data.\n");
    exit(EXIT_FAILURE);
  }
  total_written += size;
}

void report_io(void)
{
  printf("total bytes written = %lu\n", total_written);
  printf("total bytes read = %lu\n", total_read);
}

static size_t compress(struct CompressState *state, void *dst, void *src,
    size_t chunk)
{
  return ZSTD_compressCCtx(state->c_ctx, dst, compress_bound, src, chunk,
      ZSTD_LEVEL);
}

static void decompress(struct CompressState *state, void *dst, size_t chunk,
    void *src, size_t compressed)
{
  ZSTD_decompressDCtx(state->d_ctx, dst, chunk, src, compressed);
}

void copy_data(FILE *F, FILE *G, size_t size)
{
  init();

  uint8_t *buffer = cmprs_state[0].buffer;

  while (size) {
    size_t chunk = min(COPYSIZE, size);
    file_read(buffer, chunk, G);
    file_write(buffer, chunk, F);
    size -= chunk;
  }
}

static void write_data_worker(int t)
{
  struct CompressState *state = &cmprs_state[t];

  FILE *F = cmprs_F;
  uint8_t *restrict src = cmprs_ptr;
  while (true) {
    LOCK(cmprs_mutex);
    size_t idx = cmprs_idx;
    size_t chunk = min(COPYSIZE, cmprs_size - idx);
    cmprs_idx += chunk;
    UNLOCK(cmprs_mutex);
    if (chunk == 0)
      break;
    uint8_t *buf;
    switch (cmprs_type) {
    case COPY:
      buf = src + idx;
      break;
    case U8U8:
      uint8_t *restrict v8 = cmprs_v;
      buf = state->buffer;
      for (size_t i = 0; i < chunk; i++)
        buf[i] = v8[src[idx + i]];
      break;
    case U16U16:
      uint16_t *restrict v16 = cmprs_v;
      uint16_t *restrict buf16 = state->buffer;
      for (size_t i = 0; i < chunk / 2; i++)
        buf16[i] = v16[((uint16_t *)(src + idx))[i]];
      buf = (uint8_t *)buf16;
      break;
    case U16U8:
      v8 = cmprs_v;
      buf = state->buffer;
      for (size_t i = 0; i < chunk; i++)
        buf[i] = v8[((uint16_t *)src + idx)[i]];
      break;
    case COPY_U16U8:
      buf = state->buffer;
      for (size_t i = 0; i < chunk; i++)
        buf[i] = ((uint16_t *)src + idx)[i];
      break;
    }
    size_t cmprs_chunk = compress(state, state->frame->data, buf, chunk);
    state->frame->cmprs_chunk = cmprs_chunk;
    state->frame->chunk = chunk;
    state->frame->idx = idx;
    file_write(state->frame, cmprs_chunk + HEADER_SIZE, F);
  }
}

void write_data_transform_u8(FILE *F, void *src, size_t size, uint8_t *v)
{
  init();

  cmprs_F = F;
  cmprs_ptr = src;
  cmprs_size = size;
  cmprs_idx = 0;
  cmprs_v = v;
  cmprs_type = U8U8;
  run_compression(write_data_worker);
}

void write_data_transform_u16(FILE *F, void *src, size_t size, uint16_t *v)
{
  init();

  cmprs_F = F;
  cmprs_ptr = src;
  cmprs_size = size;
  cmprs_idx = 0;
  cmprs_v = v;
  cmprs_type = U16U16;
  run_compression(write_data_worker);
}

void write_data_transform_to_u8_u8(FILE *F, void *src, size_t size,
    uint8_t *v)
{
  write_data_transform_u8(F, src, size, v);
}

void write_data_transform_to_u8_u16(FILE *F, void *src, size_t size,
    uint8_t *v)
{
  init();

  cmprs_F = F;
  cmprs_ptr = src;
  cmprs_size = size;
  cmprs_idx = 0;
  cmprs_v = v;
  cmprs_type = U16U8;
  run_compression(write_data_worker);
}

void write_data(FILE *F, void *src, size_t size)
{
  init();

  cmprs_F = F;
  cmprs_ptr = src;
  cmprs_size = size;
  cmprs_idx = 0;
  cmprs_v = nullptr;
  cmprs_type = COPY;
  run_compression(write_data_worker);
}

void write_data_as_u8_u8(FILE *F, void *src, size_t size)
{
  write_data(F, src, size);
}

void write_data_as_u8_u16(FILE *F, void *src, size_t size)
{
  init();

  cmprs_F = F;
  cmprs_ptr = src;
  cmprs_size = size;
  cmprs_idx = 0;
  cmprs_type = COPY_U16U8;
  run_compression(write_data_worker);
}

static void read_data_worker(int t)
{
  struct CompressState *state = &cmprs_state[t];

  FILE *F = cmprs_F;
  uint8_t *restrict dst = cmprs_ptr;
  while (true) {
    size_t cmprs_chunk;
    LOCK(cmprs_mutex);
    if (cmprs_size == 0) {
      UNLOCK(cmprs_mutex);
      break;
    }
    file_read(&cmprs_chunk, sizeof(size_t), F);
    file_read(&state->frame->chunk, cmprs_chunk + HEADER_SIZE - sizeof(size_t),
        F);
    size_t chunk = state->frame->chunk;
    if (chunk > cmprs_size) {
      fprintf(stderr, "Error in read_data_worker.\n");
      exit(EXIT_FAILURE);
    }
    cmprs_size -= chunk;
    UNLOCK(cmprs_mutex);
    size_t idx = state->frame->idx;
    if (cmprs_type == COPY) {
      decompress(state, dst + idx, chunk, state->frame->data, cmprs_chunk);
      continue;
    }
    decompress(state, state->buffer, chunk, state->frame->data, cmprs_chunk);
    switch (cmprs_type) {
    case U8U8:
      uint8_t *restrict v8 = cmprs_v;
      uint8_t *restrict buf = state->buffer;
      for (size_t i = 0; i < chunk; i++)
        dst[idx + i] = v8[buf[i]];
      break;
    case U16U16:
      uint8_t *restrict v16 = cmprs_v;
      uint16_t *restrict buf16 = state->buffer;
      for (size_t i = 0; i < chunk / 2; i++)
        ((uint16_t *)(dst + idx))[i] = v16[buf16[i]];
      break;
    case U16U8:
      v8 = cmprs_v;
      buf16 = state->buffer;
      for (size_t i = 0; i < chunk / 2; i++)
        dst[idx / 2 + i] = v8[buf16[i]];
      break;
    }
  }
}

void read_data(FILE *F, void *dst, uint64_t size)
{
  init();

  cmprs_F = F;
  cmprs_ptr = dst;
  cmprs_size = size;
  cmprs_type = COPY;
  run_compression(read_data_worker);
}

void read_data_transform_u8(FILE *F, void *dst, size_t size, uint8_t *v)
{
  init();

  cmprs_F = F;
  cmprs_ptr = dst;
  cmprs_size = size;
  cmprs_v = v;
  cmprs_type = U8U8;
  run_compression(read_data_worker);
}

void read_data_transform_u16(FILE *F, void *dst, size_t size, uint16_t *v)
{
  init();

  cmprs_F = F;
  cmprs_ptr = dst;
  cmprs_size = size;
  cmprs_v = v;
  cmprs_type = U16U16;
  run_compression(read_data_worker);
}

void read_data_transform_to_u8_u16(FILE *F, void *dst, size_t size, uint8_t *v)
{
  init();

  cmprs_F = F;
  cmprs_ptr = dst;
  cmprs_size = size;
  cmprs_v = v;
  cmprs_type = U16U8;
  run_compression(read_data_worker);
}

void create_dir(int n, int stm, const char *name)
{
  char pathname[128];

  if (n >= 0)
    sprintf(pathname, "%d/%s/%c/", n, name, "wb"[stm]);
  else
    sprintf(pathname, "%s/%c/", name, "wb"[stm]);
  for (char *p = pathname + 1; *p; p++)
    if (*p == '/') {
      *p = 0;
      make_dir(pathname);
      *p = '/';
    }
}

#ifndef HAS_PAWNS
void create_name(char *str, int s, int stm, const char *name, int n)
{
  int wk = KKSquare[s][0], bk = KKSquare[s][1];
  if (n >= 0)
    sprintf(str, "%d/%s/%c/%c%c%c%c", n, name, "wb"[stm], 'a' + (wk & 7),
        '1' + (wk >> 3), 'a'+ (bk & 7), '1' + (bk >> 3));
  else
    sprintf(str, "%s/%c/%c%c%c%c", name, "wb"[stm], 'a' + (wk & 7),
        '1' + (wk >> 3), 'a'+ (bk & 7), '1' + (bk >> 3));
}
#else
void create_name_sq(char *str, int s1, int s2, int stm, const char *name, int n)
{
  if (n >= 0)
    sprintf(str, "%d/%s/%c/%c%c%c%c", n, name, "wb"[stm], 'a' + (s1 & 7),
        '1' + (s1 >> 3), 'a'+ (s2 & 7), '1' + (s2 >> 3));
  else
    sprintf(str, "%s/%c/%c%c%c%c", name, "wb"[stm], 'a' + (s1 & 7),
        '1' + (s1 >> 3), 'a'+ (s2 & 7), '1' + (s2 >> 3));
}

void create_name_r(char *str, int s, int r, int stm, const char *name, int n)
{
  int wk = KK16Square[s][r][0], bk = KK16Square[s][r][1];
  create_name_sq(str, wk, bk, stm, name, n);
}

void create_name(char *str, int s, int stm, const char *name, int n)
{
  create_name_r(str, s, 0, stm, name, n);
}

void create_name_p(char *str, int sq, int stm, const char *name)
{
  sprintf(str, "%s/%c/%c%c", name, "wb"[stm], 'a' + (sq & 7), '1' + (sq >> 3));
}
#endif

void create_name_10(char *str, int k, int stm, const char *name)
{
  int sq = InvTriangle[k];
  sprintf(str, "%s/%c/%c%c", name, "wb"[stm], 'a' + (sq & 7), '1' + (sq >> 3));
}
