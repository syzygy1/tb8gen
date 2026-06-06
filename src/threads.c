/*
  Copyright (c) 2011-2013, 2018, 2025, 2026 Ronald de Man

  This file is distributed under the terms of the GNU GPL, version 2.
*/

#ifdef __linux__
#define _GNU_SOURCE
#include <sched.h>
#endif

#include <stdatomic.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#include "defs.h"
#include "types.h"
#include "threads.h"
#include "util.h"

#if defined(__STDC_NO_THREADS__) || !__has_include(<threads.h>)
#include "c11threads_win32.c"
#endif

struct ThreadData *g_thread_data;

int g_num_threads;
int g_total_work;
bool g_thread_affinity;

static void setaffinity(int i)
{
#ifdef __linux__
  cpu_set_t cpuset;
  CPU_ZERO(&cpuset);
  CPU_SET(i, &cpuset);
  if (sched_setaffinity(0, sizeof(cpuset), &cpuset))
    perror("sched_setaffinity");
#else
  (void)i;
#endif
}

static thrd_t *threads;

struct Queue {
  void (*func)(struct ThreadData *);
  struct Work *work;
  alignas(64) _Atomic int counter;
} Queue;

alignas(64) static struct Queue queue;

struct timeval g_start_time, g_cur_time;

void fill_work(int n, uint64_t size, uint64_t mask, uint64_t *w)
{
  w[0] = 0;
  w[n] = size;

  for (size_t i = 1; i < n; i++)
    w[i] = ((i * size) / (uint64_t)n) & ~mask;
}

void fill_work_offset(int n, uint64_t size, uint64_t mask, uint64_t *w,
    uint64_t offset)
{
  fill_work(n, size, mask, w);
  for (int i = 0; i <= n; i++)
    w[i] += offset;
}

int calc_work_units(uint64_t size, int factor, uint64_t min_chunk)
{
  int units = g_num_threads * factor;
  uint64_t max_units = min_chunk
      ? (size + min_chunk - 1) / min_chunk
      : (uint64_t)units;

  if (max_units < (uint64_t)units)
    units = max_units;
  if (units < 1)
    units = 1;

  return units;
}

uint64_t *alloc_work(int n)
{
  return malloc((n + 1) * sizeof(uint64_t));
}

static void ensure_work_capacity(struct Work *work, int units)
{
  if (work->capacity < units) {
    free(work->range);
    work->range = alloc_work(units);
    if (!work->range)
      out_of_mem();
    work->capacity = units;
  }
}

void work_init(struct Work *work, uint64_t size, uint64_t mask,
    enum WorkSchedule schedule, int factor, uint64_t min_chunk)
{
  *work = (struct Work){
    .mask = mask,
    .min_chunk = min_chunk,
    .factor = factor,
    .schedule = schedule
  };
  work_refill(work, size);
}

void work_init_units(struct Work *work, int units, uint64_t size, uint64_t mask,
    enum WorkSchedule schedule)
{
  *work = (struct Work){ .schedule = schedule };
  work_refill_units(work, units, size, mask);
}

void work_refill(struct Work *work, uint64_t size)
{
  int units = calc_work_units(size, work->factor, work->min_chunk);
  work_refill_units(work, units, size, work->mask);
}

void work_refill_units(struct Work *work, int units, uint64_t size,
    uint64_t mask)
{
  if (units < 1)
    units = 1;
  ensure_work_capacity(work, units);
  fill_work(units, size, mask, work->range);
  work->units = units;
  work->size = size;
  work->mask = mask;
}

void work_free(struct Work *work)
{
  free(work->range);
  *work = (struct Work){0};
}

struct Work *create_work(int n, uint64_t size, uint64_t mask)
{
  struct Work *work = calloc(1, sizeof(*work));
  if (!work)
    out_of_mem();
  work_init_units(work, n, size, mask, WORK_DYNAMIC);

  return work;
}

uint64_t *create_work_offset(int n, uint64_t size, uint64_t mask,
    uint64_t offset)
{
  uint64_t *w = alloc_work(n);
  fill_work_offset(n, size, mask, w, offset);

  return w;
}

alignas(64) struct Pool {
  mtx_t mtx;
  cnd_t cv_work, cv_done;
  size_t generation;
  size_t helpers_done;
  size_t helpers_target;
} pool, cpool;

static int worker(void *arg);

#if 0
static void update_time(bool report_time)
{
  int secs, usecs;
  struct timeval stop_time;

  gettimeofday(&stop_time, nullptr);
  secs = stop_time.tv_sec - g_cur_time.tv_sec;
  usecs = stop_time.tv_usec - g_cur_time.tv_usec;
  if (usecs < 0) {
    usecs += 1000000;
    secs--;
  }
  if (report_time)
    printf("time taken = %3d:%02d.%03d\n", secs / 60, secs % 60, usecs/1000);
  g_cur_time = stop_time;
}
#endif

void init_threads(void)
{
  assume(g_num_threads >= 1); // to get rid of some warnings

  mtx_init(&pool.mtx, mtx_plain);
  pool.generation = 0;

  mtx_init(&cpool.mtx, mtx_plain);
  cpool.generation = 0;

  g_thread_data = alloc_aligned(g_num_threads * sizeof(*g_thread_data), 64);

  for (int i = 0; i < g_num_threads; i++) {
    g_thread_data[i].thread_id = i;
    g_thread_data[i].affinity = -1;
  }

  threads = malloc(g_num_threads * sizeof(*threads));

  for (int i = 0; i < g_num_threads - 1; i++) {
    int rc = thrd_create(&threads[i], worker, (void *)&(g_thread_data[i]));
    if (rc != thrd_success) {
      fprintf(stderr, "ERROR: thrd_create() returned %d\n", rc);
      exit(EXIT_FAILURE);
    }
  }
  threads[g_num_threads - 1] = thrd_current();

  if (g_thread_affinity) {
    for (int i = 0; i < g_num_threads; i++)
      g_thread_data[i].affinity = i;
    setaffinity(g_thread_data[g_num_threads - 1].affinity);
  }
}

static int worker(void *arg)
{
  struct ThreadData *thread = arg;
  size_t seen_generation = 0;

  while (true) {
    mtx_lock(&pool.mtx);
    while (pool.generation == seen_generation)
      cnd_wait(&pool.cv_work, &pool.mtx);

    seen_generation = pool.generation;
    mtx_unlock(&pool.mtx);

    struct Work *work = queue.work;
    int total = work->units;

    if (work->schedule == WORK_STATIC) {
      int w = atomic_fetch_add_explicit(&queue.counter, 1,
          memory_order_relaxed);
      if (w < total) {
        thread->begin = work->range[w];
        thread->end = work->range[w + 1];
        queue.func(thread);
      }
    } else while (true) {
      int w = atomic_fetch_add_explicit(&queue.counter, 1,
          memory_order_relaxed);
      if (w >= total) break;
      thread->begin = work->range[w];
      thread->end = work->range[w + 1];
      queue.func(thread);
    }

    mtx_lock(&pool.mtx);
    pool.helpers_done++;
    if (pool.helpers_done == pool.helpers_target)
      cnd_signal(&pool.cv_done);
    mtx_unlock(&pool.mtx);
  }

  return 0;
}

#if 0
#ifdef HAS_PAWNS
int group_worker(void *arg)
{
  struct ThreadData *thread = (struct ThreadData *)arg;
  struct GroupData *g = thread->group;
  int t = thread->thread_id;

  if (t != g_num_threads - 1) {
    if (thread->affinity >= 0)
      setaffinity(thread->affinity);
  }

  do {
    barrier_wait(&(g->barrier));

    int total = g->queue.total;

    while (true) {
      int w = atomic_fetch_add_explicit(&(g->queue.counter), 1,
            memory_order_relaxed);
      if (w >= total) break;
      thread->begin = g->queue.work[w];
      thread->end = g->queue.work[w + 1];
      g->queue.func(thread);
    }

    barrier_wait(&(g->barrier));
  } while (t != g_num_threads - 1);

  return 0;
}

void run_group(struct GroupData *g, void (*func)(struct ThreadData *),
    uint64_t *work, int report_time)
{
  int secs, usecs;
  struct timeval stop_time;

  g->queue.func = func;
  g->queue.work = work;
  g->queue.total = g_total_work;
  g->queue.counter = 0;

  worker((void *)&(g_thread_data[g_num_threads - 1]));

  gettimeofday(&stop_time, nullptr);
  secs = stop_time.tv_sec - g_cur_time.tv_sec;
  usecs = stop_time.tv_usec - g_cur_time.tv_usec;
  if (usecs < 0) {
    usecs += 1000000;
    secs--;
  }
  if (report_time)
    printf("time taken = %3d:%02d.%03d\n", secs / 60, secs % 60, usecs/1000);
  g_cur_time = stop_time;
}

#endif
#endif

void run_threaded(void (*func)(struct ThreadData *), struct Work *work)
{
  int total_work = work->units;
  if (total_work <= 1) {
    if (total_work == 1) {
      struct ThreadData *thread = &g_thread_data[g_num_threads - 1];
      thread->begin = work->range[0];
      thread->end = work->range[1];
      func(thread);
    }
    return;
  }

  queue.func = func;
  queue.work = work;
  queue.counter = 0;
  int total = work->units;
  struct ThreadData *thread = &g_thread_data[g_num_threads - 1];
  int helpers = g_num_threads - 1;

  mtx_lock(&pool.mtx);
  pool.helpers_done = 0;
  pool.helpers_target = helpers;
  atomic_store_explicit(&queue.counter, 0, memory_order_relaxed);
  pool.generation++;
  cnd_broadcast(&pool.cv_work);
  mtx_unlock(&pool.mtx);

  if (work->schedule == WORK_STATIC) {
    int w = atomic_fetch_add_explicit(&queue.counter, 1,
        memory_order_relaxed);
    if (w < total) {
      thread->begin = work->range[w];
      thread->end = work->range[w + 1];
      queue.func(thread);
    }
  } else while (true) {
    int w = atomic_fetch_add_explicit(&queue.counter, 1,
        memory_order_relaxed);
    if (w >= total) break;
    thread->begin = work->range[w];
    thread->end = work->range[w + 1];
    queue.func(thread);
  }

  mtx_lock(&pool.mtx);
  while (pool.helpers_done != pool.helpers_target)
    cnd_wait(&pool.cv_done, &pool.mtx);
  mtx_unlock(&pool.mtx);
}

static struct ThreadData cmprs_data[COMPRESSION_THREADS];
static void (*cmprs_func)(int t);

static thrd_t cmprs_threads[COMPRESSION_THREADS];

static int cmprs_worker(void *arg)
{
  struct ThreadData *thread = arg;
  int t = thread->thread_id;
  size_t seen_generation = 0;

  while (true) {
    mtx_lock(&cpool.mtx);
    while (cpool.generation == seen_generation)
      cnd_wait(&cpool.cv_work, &cpool.mtx);

    seen_generation = cpool.generation;
    mtx_unlock(&cpool.mtx);

    cmprs_func(t);

    mtx_lock(&cpool.mtx);
    cpool.helpers_done++;
    if (cpool.helpers_done == COMPRESSION_THREADS - 1)
      cnd_signal(&cpool.cv_done);
    mtx_unlock(&cpool.mtx);
  }

  return 0;
}

void create_compression_threads(void)
{
  for (int i = 0; i < COMPRESSION_THREADS; i++)
    cmprs_data[i].thread_id = i;

  for (int i = 0; i < COMPRESSION_THREADS - 1; i++) {
    int rc = thrd_create(&cmprs_threads[i], cmprs_worker, &cmprs_data[i]);
    if (rc != thrd_success) {
      fprintf(stderr, "ERROR: thrd_create() returned %d\n", rc);
      exit(EXIT_FAILURE);
    }
  }
}

void run_compression(void (*func)(int t))
{
  cmprs_func = func;

  mtx_lock(&cpool.mtx);
  cpool.helpers_done = 0;
  cpool.generation++;
  cnd_broadcast(&cpool.cv_work);
  mtx_unlock(&cpool.mtx);

  func(COMPRESSION_THREADS - 1);

  mtx_lock(&cpool.mtx);
  while (cpool.helpers_done != COMPRESSION_THREADS - 1)
    cnd_wait(&cpool.cv_done, &cpool.mtx);
  mtx_unlock(&cpool.mtx);
}
