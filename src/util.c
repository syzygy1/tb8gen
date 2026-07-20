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
#include <time.h>
#include <unistd.h>
#ifndef _WIN32
#include <unistd.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>
#include <fcntl.h>
#endif

#if defined(__AVX512F__) || defined(__AVX2__)
#include <immintrin.h>
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

static constexpr size_t MIN_COPYSIZE = 50*1024ULL;
static constexpr size_t MAX_COPYSIZE = 512*1024ULL;
static constexpr int TARGET_COMPRESSION_CHUNKS_PER_THREAD = 4;

static size_t compress_bound;

static LOCK_T cmprs_mutex;

static FILE *cmprs_F;
static void *cmprs_ptr;
static size_t cmprs_size;
static size_t cmprs_chunk_size;
static size_t cmprs_chunk_alignment;
static size_t cmprs_chunks_issued;
static size_t cmprs_idx;
static void *cmprs_v;
static int cmprs_type;

enum {
  COPY, U8U8, U16U16, U16U8, COPY_U16U8,
  U8_OR, U8_ANDNOT, U8_AND, U8U8_OR, U8U16_OR
};

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
    compress_bound = ZSTD_compressBound(MAX_COPYSIZE);
    for (int i = 0; i < COMPRESSION_THREADS; i++) {
      cmprs_state[i].buffer = alloc_aligned(MAX_COPYSIZE, 64);
      cmprs_state[i].frame = malloc(HEADER_SIZE + compress_bound);
      cmprs_state[i].c_ctx = ZSTD_createCCtx();
      cmprs_state[i].d_ctx = ZSTD_createDCtx();
    }
    create_compression_threads();
  }
}

static size_t align_compression_chunk(size_t chunk, size_t remaining)
{
  size_t alignment = cmprs_chunk_alignment;

  if (alignment > 1) {
    if (chunk < remaining) {
      chunk = (chunk + alignment - 1) & ~(alignment - 1);
      if (chunk > remaining)
        chunk = remaining & ~(alignment - 1);
    }
    if (chunk == 0)
      chunk = remaining;
  }

  return chunk;
}

static size_t next_compression_chunk(void)
{
  size_t remaining = cmprs_size - cmprs_idx;
  if (!remaining)
    return 0;

  size_t chunk;
  if (cmprs_chunks_issued < COMPRESSION_THREADS) {
    chunk = MIN_COPYSIZE;
  } else if (remaining > cmprs_chunk_size * COMPRESSION_THREADS * 2) {
    chunk = cmprs_chunk_size;
  } else {
    size_t chunks = COMPRESSION_THREADS * TARGET_COMPRESSION_CHUNKS_PER_THREAD;
    chunk = (remaining + chunks - 1) / chunks;
    if (chunk < MIN_COPYSIZE)
      chunk = MIN_COPYSIZE;
    if (chunk > cmprs_chunk_size)
      chunk = cmprs_chunk_size;
  }

  if (chunk > remaining)
    chunk = remaining;

  cmprs_chunks_issued++;
  return align_compression_chunk(chunk, remaining);
}

static void prepare_write_data(FILE *F, void *src, size_t size,
    size_t alignment)
{
  cmprs_F = F;
  cmprs_ptr = src;
  cmprs_size = size;
  cmprs_chunk_size = MAX_COPYSIZE;
  cmprs_chunk_alignment = alignment;
  cmprs_chunks_issued = 0;
  cmprs_idx = 0;
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

FILE *file_open_read(const char *name)
{
  FILE *F = fopen(name, "rb");
  if (!F) {
    fprintf(stderr, "Could not open %s for reading.\n", name);
    exit(EXIT_FAILURE);
  }
  return F;
}

FILE *file_open_write(const char *name)
{
  char *str = malloc(strlen(name) + 5);
  if (!str)
    out_of_mem();
  strcat(strcpy(str, name), ".tmp");
  FILE *F = fopen(str, "wb");
  if (!F) {
    fprintf(stderr, "Could not open %s for writing.\n", str);
    exit(EXIT_FAILURE);
  }
  free(str);
  return F;
}

void file_rename(const char *name)
{
  char *str = malloc(strlen(name) + 5);
  if (!str)
    out_of_mem();
  strcat(strcpy(str, name), ".tmp");
  if (rename(str, name) < 0) {
    fprintf(stderr, "Could not rename %s into %s.\n", str, name);
    exit(EXIT_FAILURE);
  }
  free(str);
}

bool file_exists(const char *name)
{
  struct stat st;
  return stat(name, &st) == 0;
}

void create_empty(const char *name)
{
  int fd = creat(name, 0666);
  if (fd >= 0)
    close(fd);
}

INLINE size_t compress(struct CompressState *state, void *dst, void *src,
    size_t chunk)
{
  return ZSTD_compressCCtx(state->c_ctx, dst, compress_bound, src, chunk,
      ZSTD_LEVEL);
}

INLINE void decompress(struct CompressState *state, void *dst, size_t chunk,
    void *src, size_t compressed)
{
  ZSTD_decompressDCtx(state->d_ctx, dst, chunk, src, compressed);
}

static void decompress_logical(struct CompressState *state, uint8_t *dst,
    size_t compressed, int op)
{
  ZSTD_inBuffer in = { state->frame->data, compressed, 0 };
  uint8_t *buf = state->buffer;

  for (;;) {
    ZSTD_outBuffer out = { buf, 64 * 1024, 0 };

    size_t res = ZSTD_decompressStream(state->d_ctx, &out, &in);
    if (ZSTD_isError(res)) {
      fprintf(stderr, "zstd error: %s\n", ZSTD_getErrorName(res));
      exit(EXIT_FAILURE);
    }

    if (op == U8_OR || op == U8_ANDNOT || op == U8_AND)
      assert((out.pos & 0x3f) == 0);

    if (op == U8_OR) {

#ifdef __AVX512F__

      for (size_t off = 0; off < out.pos; off += 64) {
        __m512i d = _mm512_load_si512((__m512i *)(buf + off));

        if (_mm512_test_epi64_mask(d, d) == 0)
          continue;

        uint8_t *p = dst + off;

        __m512i old = _mm512_load_si512((__m512i *)p);
        _mm512_store_si512((__m512i *)p, _mm512_or_si512(old, d));
      }

#elifdef __AVX2__

      for (size_t off = 0; off < out.pos; off += 64) {
        __m256i d0 = _mm256_load_si256((__m256i *)(buf + off));
        __m256i d1 = _mm256_load_si256((__m256i *)(buf + off + 32));

        if (_mm256_testz_si256(d0, d0) && _mm256_testz_si256(d1, d1))
          continue;

        uint8_t *p = dst + off;

        __m256i old0 = _mm256_load_si256((__m256i *)(p));
        __m256i old1 = _mm256_load_si256((__m256i *)(p + 32));

        _mm256_store_si256((__m256i *)(p), _mm256_or_si256(old0, d0));
        _mm256_store_si256((__m256i *)(p + 32), _mm256_or_si256(old1, d1));
      }

#else

      for (size_t off = 0; off < out.pos; off += 8) {
        uint64_t d = *(uint64_t *)(buf + off);

        if (!d) continue;

        uint64_t *p = (uint64_t *)(dst + off);
        *p |= d;
      }

#endif

    } else if (op == U8_ANDNOT) {

#ifdef __AVX512F__

      for (size_t off = 0; off < out.pos; off += 64) {
        __m512i d = _mm512_load_si512((__m512i *)(buf + off));

        if (_mm512_test_epi64_mask(d, d) == 0)
          continue;

        uint8_t *p = dst + off;

        __m512i old = _mm512_load_si512((__m512i *)p);
        _mm512_store_si512((__m512i *)p, _mm512_andnot_si512(d, old));
      }

#elifdef __AVX2__

      for (size_t off = 0; off < out.pos; off += 64) {
        __m256i d0 = _mm256_load_si256((__m256i *)(buf + off));
        __m256i d1 = _mm256_load_si256((__m256i *)(buf + off + 32));

        if (_mm256_testz_si256(d0, d0) && _mm256_testz_si256(d1, d1))
          continue;

        uint8_t *p = dst + off;

        __m256i old0 = _mm256_load_si256((__m256i *)(p));
        __m256i old1 = _mm256_load_si256((__m256i *)(p + 32));

        _mm256_store_si256((__m256i *)(p), _mm256_andnot_si256(d0, old0));
        _mm256_store_si256((__m256i *)(p + 32), _mm256_andnot_si256(d1, old1));
      }

#else

      for (size_t off = 0; off < out.pos; off += 8) {
        uint64_t d = *(uint64_t *)(buf + off);

        if (!d) continue;

        uint64_t *p = (uint64_t *)(dst + off);
        *p &= ~d;
      }

#endif

    } else if (op == U8_AND) {

#ifdef __AVX512F__

      for (size_t off = 0; off < out.pos; off += 64) {
        __m512i d = _mm512_load_si512((__m512i *)(buf + off));

        uint8_t *p = dst + off;

        __m512i old = _mm512_load_si512((__m512i *)p);
        _mm512_store_si512((__m512i *)p, _mm512_and_si512(old, d));
      }

#elifdef __AVX2__

      for (size_t off = 0; off < out.pos; off += 64) {
        __m256i d0 = _mm256_load_si256((__m256i *)(buf + off));
        __m256i d1 = _mm256_load_si256((__m256i *)(buf + off + 32));

        uint8_t *p = dst + off;

        __m256i old0 = _mm256_load_si256((__m256i *)(p));
        __m256i old1 = _mm256_load_si256((__m256i *)(p + 32));

        _mm256_store_si256((__m256i *)(p), _mm256_and_si256(old0, d0));
        _mm256_store_si256((__m256i *)(p + 32), _mm256_and_si256(old1, d1));
      }

#else

      for (size_t off = 0; off < out.pos; off += 8) {
        uint64_t d = *(uint64_t *)(buf + off);

        uint64_t *p = (uint64_t *)(dst + off);
        *p &= d;
      }

#endif

    } else if (op == U8U8_OR) {

      uint8_t *restrict v8 = cmprs_v;

      for (size_t off = 0; off < out.pos; off++)
        dst[off] |= v8[buf[off]];

    } else { // U8U16_OR

      uint16_t *restrict v16 = cmprs_v;

      for (size_t off = 0; off < out.pos; off++)
        ((uint16_t *)dst)[off] |= v16[buf[off]];

      dst += out.pos; // hack
    }

    dst += out.pos;

    if (res == 0)
      break;
  }

  if (in.pos != in.size) {
    fprintf(stderr, "Trailing compressed data in decompress_or().\n");
    exit(EXIT_FAILURE);
  }
}

void copy_data(FILE *F, FILE *G, size_t size)
{
  init();

  uint8_t *buffer = cmprs_state[0].buffer;

  while (size) {
    size_t chunk = min(MAX_COPYSIZE, size);
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
    size_t chunk = next_compression_chunk();
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

  prepare_write_data(F, src, size, 1);
  cmprs_v = v;
  cmprs_type = U8U8;
  run_compression(write_data_worker);
}

void write_data_transform_u16(FILE *F, void *src, size_t size, uint16_t *v)
{
  init();

  prepare_write_data(F, src, size, 2);
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

  prepare_write_data(F, src, size, 1);
  cmprs_v = v;
  cmprs_type = U16U8;
  run_compression(write_data_worker);
}

void write_data(FILE *F, void *src, size_t size)
{
  init();

  prepare_write_data(F, src, size, 1);
  cmprs_v = nullptr;
  cmprs_type = COPY;
  run_compression(write_data_worker);
}

void write_data_cache_aligned(FILE *F, void *src, size_t size)
{
  init();

  prepare_write_data(F, src, size, 64);
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

  prepare_write_data(F, src, size, 1);
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
    if (cmprs_type >= U8_OR) {
      // Frame offsets describe the byte-oriented source stream.  U8U16_OR
      // expands every source byte to one u16, so its destination offset is
      // twice the frame offset.  Using the source offset here misaligns all
      // frames after the first one and corrupts the reconstructed u16 table.
      size_t dst_idx = cmprs_type == U8U16_OR ? 2 * idx : idx;
      decompress_logical(state, dst + dst_idx, cmprs_chunk, cmprs_type);
      continue;
    }
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
      uint16_t *restrict v16 = cmprs_v;
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

void read_data_or(FILE *F, void *dst, uint64_t size)
{
  init();

  cmprs_F = F;
  cmprs_ptr = dst;
  cmprs_size = size;
  cmprs_type = U8_OR;
  run_compression(read_data_worker);
}

void read_data_and(FILE *F, void *dst, uint64_t size)
{
  init();

  cmprs_F = F;
  cmprs_ptr = dst;
  cmprs_size = size;
  cmprs_type = U8_AND;
  run_compression(read_data_worker);
}

void read_data_andnot(FILE *F, void *dst, uint64_t size)
{
  init();

  cmprs_F = F;
  cmprs_ptr = dst;
  cmprs_size = size;
  cmprs_type = U8_ANDNOT;
  run_compression(read_data_worker);
}

void read_data_transform_or_u8(FILE *F, void *dst, uint64_t size, uint8_t *v)
{
  init();

  cmprs_F = F;
  cmprs_ptr = dst;
  cmprs_size = size;
  cmprs_v = v;
  cmprs_type = U8U8_OR;
  run_compression(read_data_worker);
}

void read_data_transform_or_u8_to_u16(FILE *F, void *dst, uint64_t size,
    uint16_t *v)
{
  init();

  cmprs_F = F;
  cmprs_ptr = dst;
  cmprs_size = size;
  cmprs_v = v;
  cmprs_type = U8U16_OR;
  run_compression(read_data_worker);
}

void create_dir(int n, int stm, const char *name)
{
  char pathname[64];

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

bool dir_exists(int n, int stm, const char *name)
{
  char pathname[64];

  if (n >= 0)
    sprintf(pathname, "%d/%s/%c/", n, name, "wb"[stm]);
  else
    sprintf(pathname, "%s/%c/", name, "wb"[stm]);
  return file_exists(pathname);
}

void delete_dir(int n, const char *name)
{
  char dir[64];
  if (n >= 0)
    sprintf(dir, "%d/%s", n, name);
  else
    strcpy(dir, name);

  if (!file_exists(dir))
    return;

  char str[80];
  for (int stm = 0; stm < 2; stm++) {
    sprintf(str, "%s/%c/done", dir, "wb"[stm]);
    unlink(str);
    sprintf(str, "%s/%c", dir, "wb"[stm]);
    rmdir(str);
  }
  rmdir(dir);
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

static double now_sec(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

static void fmt_duration(char buf[32], double t)
{
  uint32_t secs = t;
  uint32_t days  = secs / 86400;
  secs %= 86400;
  uint32_t hours = secs / 3600;
  secs %= 3600;
  uint32_t mins  = secs / 60;
  secs %= 60;

  if (days)
    snprintf(buf, 32, "%ud %02uh %02um %02us", days, hours, mins, secs);
  else if (hours)
    snprintf(buf, 32, "%uh %02um %02us", hours, mins, secs);
  else if (mins)
    snprintf(buf, 32, "%um %02us", mins, secs);
  else
    snprintf(buf, 32, "%us", secs);
}

void show_progress(const char *phase, int k, int total, bool final)
{
  static double t_init, t0, last = 0.0;
  static int first_k;
  static bool init = true, first = true, disable = false;

  if (disable)
    return;

  double t = now_sec();

  if (first) {
    t0 = t;
    first_k = k;
    first = false;
    if (init) {
      if (!isatty(STDERR_FILENO))
        disable = true;
      t_init = t0;
      init = false;
    }
  }
  else if (!final && t - last < 0.25)
    return;

  last = t;

  char ebuf[32];
  fmt_duration(ebuf, t - t_init);

  if (!final) {
    char etabuf[32];
    int done = k - first_k;
    int remaining = total - k;
    if (done > 0 && remaining > 0)
      fmt_duration(etabuf, (t - t0) * remaining / done);
    fprintf(stderr, "\r\033[K%s  %s  %3u/%u  [%s left]", ebuf, phase, k, total,
        done > 0 && remaining > 0 ? etabuf : "--");
  }
  else {
    fprintf(stderr, "\r\033[k%s  %s\033[K\n", ebuf, phase);
    first = true;
  }

  fflush(stderr);
}
