/*
  Copyright (c) 2011-2013, 2018, 2025, 2026 Ronald de Man

  This file is distributed under the terms of the GNU GPL, version 2.
*/

#ifndef THREADS_H
#define THREADS_H

#include <stdalign.h>
#include <stdatomic.h>

#include "defs.h"
#include "probe.h"
#include "types.h"

#if defined(__STDC_NO_THREADS__) || !__has_include(<threads.h>)
#include "c11threads.h"
#else
#include <threads.h>
#endif

#define LOCK_T mtx_t
#define LOCK_INIT(x) mtx_init(&(x), mtx_plain)
#define LOCK(x) mtx_lock(&(x))
#define UNLOCK(x) mtx_unlock(&(x))

struct ThreadData {
  alignas(64) size_t begin, end;
  uint64_t cnt;
  int thread_id;
  int affinity;
};

enum WorkSchedule {
  WORK_DYNAMIC,
  WORK_STATIC
};

struct Work {
  uint64_t *range;
  int units;
  int capacity;
  uint64_t size;
  uint64_t mask;
  uint64_t min_chunk;
  int factor;
  enum WorkSchedule schedule;
};

extern int g_num_threads;
extern struct ThreadData *g_thread_data;
extern bool g_thread_affinity;
extern int g_total_work;
extern struct timeval g_start_time, g_cur_time;

void init_threads(void);
void run_threaded(void (*func)(struct ThreadData *), struct Work *work,
    bool report_time);
void run_single(void (*func)(struct ThreadData *), struct Work *work,
    bool report_time);
void fill_work(int n, uint64_t size, uint64_t mask, uint64_t *w);
void fill_work_offset(int n, uint64_t size, uint64_t mask, uint64_t *w,
    uint64_t offset);
int calc_work_units(uint64_t size, int factor, uint64_t min_chunk);
uint64_t *alloc_work(int n);
void work_init(struct Work *work, uint64_t size, uint64_t mask,
    enum WorkSchedule schedule, int factor, uint64_t min_chunk);
void work_init_units(struct Work *work, int units, uint64_t size, uint64_t mask,
    enum WorkSchedule schedule);
void work_refill(struct Work *work, uint64_t size);
void work_refill_units(struct Work *work, int units, uint64_t size,
    uint64_t mask);
void work_free(struct Work *work);
struct Work *create_work(int n, uint64_t size, uint64_t mask);

void create_compression_threads(void);
void run_compression(void (*func)(int t));

#endif
