/*
  Copyright (c) 2011-2013, 2018, 2024, 2025 Ronald de Man

  This file is distributed under the terms of the GNU GPL, version 2.
*/

#define _POSIX_C_SOURCE 200112L
#include <stdlib.h>

#include <assert.h>
#include <stdatomic.h>
#include <stddef.h>
#include <stdio.h>
#include <stdint.h>
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

//static constexpr size_t COPYSIZE = 20*1024*1024ULL;
static constexpr size_t COPYSIZE = 50*1024ULL;

static size_t compress_bound;

static LOCK_T cmprs_mutex;

static FILE *cmprs_F;
static void *cmprs_ptr;
static size_t cmprs_size;
static size_t cmprs_idx;

struct CompressFrame {
  size_t cmprs_chunk;
  size_t chunk;
  size_t idx;
  uint8_t data[];
};

static constexpr size_t HEADER_SIZE = offsetof(struct CompressFrame, data);

struct CompressState {
  uint8_t *buffer;
  struct CompressFrame *frame;
#ifdef USE_ZSTD
  ZSTD_CCtx *c_ctx;
  ZSTD_DCtx *d_ctx;
#endif
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

void file_read(void *ptr, size_t size, FILE *F)
{
  if (fread(ptr, 1, size, F) != size) {
    fprintf(stderr, "Error reading data from disk.\n");
    exit(EXIT_FAILURE);
  }
}

void file_write(void *ptr, size_t size, FILE *F)
{
  if (fwrite(ptr, 1, size, F) != size) {
    fprintf(stderr, "Error writing data to disk.\n");
    exit(EXIT_FAILURE);
  }
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
    buf = src + idx;
    size_t cmprs_chunk = compress(state, state->frame->data, buf, chunk);
    state->frame->cmprs_chunk = cmprs_chunk;
    state->frame->chunk = chunk;
    state->frame->idx = idx;
    file_write(state->frame, cmprs_chunk + HEADER_SIZE, F);
  }
}

void write_data(FILE *F, uint8_t *src, size_t size)
{
  init();

  cmprs_F = F;
  cmprs_ptr = src;
  cmprs_size = size;
  cmprs_idx = 0;
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
    decompress(state, dst + idx, chunk, state->frame->data, cmprs_chunk);
  }
}

void read_data(FILE *F, uint8_t *dst, uint64_t size)
{
  init();

  cmprs_F = F;
  cmprs_ptr = dst;
  cmprs_size = size;
  run_compression(read_data_worker);
}
