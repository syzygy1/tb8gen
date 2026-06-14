/*
  Copyright (c) 2011-2013, 2025 Ronald de Man

  This file is distributed under the terms of the GNU GPL, version 2.
*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include "defs.h"
#include "checksum.h"
#include "threads.h"
#include "util.h"
#include "hash/citycrc.h"
#include "hash/xxhash.h"

#define CHUNK (1ULL << 24)

static uint64_t checksum1[2];
static uint64_t checksum2[2];
static uint64_t *results = nullptr;
static XXH128_hash_t *xxh_results = nullptr;
static char *data;
static size_t fsize;
static bool checksum_found;
static bool checksum_match;
static struct Work *work = nullptr;

static_assert(sizeof checksum1 == sizeof(XXH128_hash_t));

static void calc_cityhash_worker(struct ThreadData *thread)
{
  uint64_t idx;
  uint64_t end = thread->end;

  for (idx = thread->begin; idx < end; idx++) {
    uint64_t start = idx * CHUNK;
    uint64_t chunk = CHUNK;
    if (start + chunk > fsize)
      chunk = fsize - start;
    CityHashCrc256(data + start, chunk, &results[4 * idx]);
  }
}

static void calc_xxhash_worker(struct ThreadData *thread)
{
  uint64_t idx;
  uint64_t end = thread->end;

  for (idx = thread->begin; idx < end; idx++) {
    uint64_t start = idx * CHUNK;
    uint64_t chunk = CHUNK;
    if (start + chunk > fsize)
      chunk = fsize - start;
    xxh_results[idx] = XXH3_128bits(data + start, chunk);
  }
}

static void calc_cityhash(const char *name)
{
  FD fd = open_file(name);
  data = map_file(fd, false, &fsize);
  size_t orig_size = fsize;
  if ((fsize & 0x3f) == 0x10) {
    fsize &= ~0x3fULL;
    memcpy(checksum1, data + fsize, 16);
    checksum_found = true;
  } else {
    if (fsize & 0x3f) {
      printf("Size of %s is not a multiple of 64.\n", name);
      exit(EXIT_FAILURE);
    }
    checksum_found = false;
  }

  int chunks = (fsize + CHUNK - 1) / CHUNK;
  results = malloc(32 * chunks);
  if (!work)
    work = create_work(g_total_work, chunks, 0);
  else
    work_refill_units(work, g_total_work, chunks, 0);
  run_threaded(calc_cityhash_worker, work);
  CityHashCrc128((char *)results, 32 * chunks, checksum2);
  unmap_file(data, orig_size);
  free(results);

  if (checksum_found)
    checksum_match = (memcmp(checksum1, checksum2, sizeof checksum1) == 0);
}

static void calc_xxhash(const char *name)
{
  FD fd = open_file(name);
  data = map_file(fd, false, &fsize);
  size_t orig_size = fsize;
  if ((fsize & 0x3f) == 0x10) {
    fsize &= ~0x3fULL;
    memcpy(checksum1, data + fsize, 16);
    checksum_found = true;
  } else {
    if (fsize & 0x3f) {
      printf("Size of %s is not a multiple of 64.\n", name);
      exit(EXIT_FAILURE);
    }
    checksum_found = false;
  }

  int chunks = (fsize + CHUNK - 1) / CHUNK;
  xxh_results = malloc(chunks * sizeof(XXH128_hash_t));
  if (!work)
    work = create_work(g_total_work, chunks, 0);
  else
    work_refill_units(work, g_total_work, chunks, 0);
  run_threaded(calc_xxhash_worker, work);
  unmap_file(data, orig_size);
  XXH128_hash_t tmp = XXH3_128bits(xxh_results, chunks * sizeof(XXH128_hash_t));
  memcpy(checksum2, &tmp, sizeof tmp);
  free(results);

  if (checksum_found)
    checksum_match = (memcmp(checksum1, checksum2, sizeof checksum1) == 0);
}

void print_checksum(const char *name, char *sum)
{
  FD fd = open_file(name);
  data = map_file(fd, true, &fsize);
  if ((fsize & 0x3f) == 0x10) {
    memcpy(checksum1, data + (fsize & ~0x3fULL), 16);
  } else {
    printf("No checksum found.\n");
    exit(EXIT_FAILURE);
  }
  unmap_file(data, fsize);

  int i;
  static char nibble[16] = "0123456789abcdef";
  uint8_t *c = (uint8_t *)checksum1;

  for (i = 0; i < 16; i++) {
    sum[2 * i] = nibble[c[i] >> 4];
    sum[2 * i + 1] = nibble[c[i] & 0x0f];
  }
  sum[32] = 0;
}

void calc_checksum(const char *name, bool new)
{
  if (new)
    calc_xxhash(name);
  else
    calc_cityhash(name);
}

void add_checksum(const char *name, bool new)
{
  calc_checksum(name, new);
  if (checksum_found) {
    printf("%s checksum already present.\n",
        checksum_match ? "Matching" : "Non-matching");
    exit(EXIT_FAILURE);
  }
  FILE *F = fopen(name, "ab");
  if (!F) {
    printf("Could not open %s for appending.\n", name);
    exit(EXIT_FAILURE);
  }
  fwrite(checksum2, 16, 1, F);
  fclose(F);
}

void add_xxhash(const char *name)
{
  add_checksum(name, true);
}

void add_cityhash(const char *name)
{
  add_checksum(name, false);
}

void verify_checksum(const char *name, bool new)
{
  printf("%s: ", name);
  calc_checksum(name, new);
  if (!checksum_found) {
    printf("No checksum present.\n");
    exit(EXIT_FAILURE);
  }
  if (!checksum_match)
    printf("FAIL!\n");
  else
    printf("OK!\n");
}
