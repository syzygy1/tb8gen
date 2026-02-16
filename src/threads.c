/*
  Copyright (c) 2011-2013, 2018, 2025 Ronald de Man

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

typedef struct {
  int needed;
  int called;
  mtx_t mutex;
  cnd_t cond;
} barrier_t;

static int barrier_init(barrier_t *barrier, int needed)
{
  barrier->needed = needed;
  barrier->called = 0;
  mtx_init(&barrier->mutex, mtx_plain);
  cnd_init(&barrier->cond);
  return 0;
}

static int barrier_destroy(barrier_t *barrier)
{
  mtx_destroy(&barrier->mutex);
  cnd_destroy(&barrier->cond);
  return 0;
}

static int barrier_wait(barrier_t *barrier)
{
  mtx_lock(&barrier->mutex);
  barrier->called++;
  if (barrier->called == barrier->needed) {
    barrier->called = 0;
    cnd_broadcast(&barrier->cond);
  } else {
    cnd_wait(&barrier->cond,&barrier->mutex);
  }
  mtx_unlock(&barrier->mutex);
  return 0;
}

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
static barrier_t barrier;

struct Queue {
  void (*func)(struct ThreadData *);
  uint64_t *work;
  _Atomic int counter;
  int total;
} Queue;

static struct Queue queue;

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

uint64_t *alloc_work(int n)
{
  return malloc((n + 1) * sizeof(uint64_t));
}

uint64_t *create_work(int n, uint64_t size, uint64_t mask)
{
  uint64_t *w = alloc_work(n);
  fill_work(n, size, mask, w);

  return w;
}

uint64_t *create_work_offset(int n, uint64_t size, uint64_t mask,
    uint64_t offset)
{
  uint64_t *w = alloc_work(n);
  fill_work_offset(n, size, mask, w, offset);

  return w;
}

int worker(void *arg);

void init_threads(void)
{
  assume(g_num_threads >= 1); // to get rid of some warnings

  g_thread_data = alloc_aligned(g_num_threads * sizeof(*g_thread_data), 64);

  for (int i = 0; i < g_num_threads; i++) {
    g_thread_data[i].thread_id = i;
    g_thread_data[i].affinity = -1;
//    g_thread_data[i].p = g_pos;
  }

  threads = malloc(g_num_threads * sizeof(*threads));
  barrier_init(&barrier, g_num_threads);

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

int worker(void *arg)
{
  struct ThreadData *thread = (struct ThreadData *)arg;
  int t = thread->thread_id;

  if (t != g_num_threads - 1) {
    if (thread->affinity >= 0)
      setaffinity(thread->affinity);
  }

  do {
    barrier_wait(&barrier);

    int total = queue.total;

    while (1) {
      int w = atomic_fetch_add_explicit(&queue.counter, 1,
            memory_order_relaxed);
      if (w >= total) break;
      thread->begin = queue.work[w];
      thread->end = queue.work[w + 1];
      queue.func(thread);
    }

    barrier_wait(&barrier);
  } while (t != g_num_threads - 1);

  return 0;
}

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

    while (1) {
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

void run_threaded(void (*func)(struct ThreadData *), uint64_t *work,
    bool report_time)
{
  int secs, usecs;
  struct timeval stop_time;

  queue.func = func;
  queue.work = work;
  queue.total = g_total_work;
  queue.counter = 0;

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

void run_single(void (*func)(struct ThreadData *), uint64_t *work,
    bool report_time)
{
  int secs, usecs;
  struct timeval stop_time;
  struct ThreadData *thread = &(g_thread_data[g_num_threads - 1]);

  thread->begin = work[0];
  thread->end = work[g_total_work];
  func(thread);

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

static struct ThreadData cmprs_data[COMPRESSION_THREADS];
static void (*cmprs_func)(int t);

static thrd_t cmprs_threads[COMPRESSION_THREADS];
static barrier_t cmprs_barrier;

static int cmprs_worker(void *arg)
{
  struct ThreadData *thread = arg;
  int t = thread->thread_id;

  do {
    barrier_wait(&cmprs_barrier);

    cmprs_func(t);

    barrier_wait(&cmprs_barrier);
  } while (t != COMPRESSION_THREADS - 1);

  return 0;
}

void create_compression_threads(void)
{
  for (int i = 0; i < COMPRESSION_THREADS; i++)
    cmprs_data[i].thread_id = i;

  barrier_init(&cmprs_barrier, COMPRESSION_THREADS);

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
  cmprs_worker(&cmprs_data[COMPRESSION_THREADS - 1]);
}
