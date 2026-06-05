/*
  Copyright (c) 2026 Ronald de Man

  This file is distributed under the terms of the GNU GPL, version 2.
*/

#include <getopt.h>
#include <inttypes.h>
#include <stdarg.h>
#include <stdatomic.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <unistd.h>

#include "defs.h"
#include "generate.h"
#include "index.h"
#include "join.h"
#include "kslice.h"
#include "merge.h"
#include "movegen.h"
#include "probe.h"
#include "stats.h"
#include "tb8gen.h"
#include "threads.h"
#include "types.h"
#include "util.h"

#define TBPATH "RTBPATH"
#define STATSDIR "RTBSTATSDIR"
#define WORKDIR "TB8DIR"

Position g_pos;
bool g_only_generate, g_use_rans, symmetric, used_rans = false;
bool g_cleanup;
bool one_sided, wins_only;
int one_sided_stm;
char *g_tablename;
char *g_output_dir;
struct Work *work_g, *work_capt[MAX_SETS];
struct MaxFen mf;

const char *typename[3] = { "wdl", "dtm", "dtz" };

static struct option options[] = {
  { "threads", 1, nullptr, 't' },
  { "stats", 0, nullptr, 's' },
  { "path", 1, nullptr, 'p' },
  { "rans", 0, nullptr, 'r' },
  { "layout", 1, nullptr, 'l' },
  { "workdir", 1, nullptr, 'w' },
  { 0 }
};

int main(int argc, char **argv)
{
  int val, lindex;
  uint8_t pcs[16];
  uint8_t pt[8];
  int layout = -1;

  const char *path = getenv(TBPATH);
  const char *workdir = getenv(WORKDIR);
  g_num_threads = 1;

  while ((val = getopt_long(argc, argv, "at:gp:rl:cw:", options, &lindex)) != -1)
    switch (val) {
    case 'a':
      g_thread_affinity = true;
      break;
    case 't':
      g_num_threads = atoi(optarg);
      break;
    case 'g':
      g_only_generate = true;
      break;
    case 'p':
      path = optarg;
      break;
    case 'r':
      g_use_rans = true;
      break;
    case 'l':
      layout = atoi(optarg);
      break;
    case 'c':
      g_cleanup = true;
      break;
    case 'w':
      workdir = optarg;
      break;
    }

  if (optind >= argc) {
    printf("No tablebase specified.\n");
    exit(EXIT_FAILURE);
  }
  g_tablename = argv[optind];
  g_output_dir = getcwd(nullptr, 0);
  if (!g_output_dir) {
    fprintf(stderr, "Could not determine current working directory.\n");
    exit(EXIT_FAILURE);
  }

  init_movegen();
  init_ranking();
  init_tablebases(path);
  init_threads();

  for (int i = 0; i < 16; i++)
    pcs[i] = 0;

  // TODO make more robust: check that each side respects order K,Q,R,B,N,P
  int color = 0, k = 0;
  for (char *s = g_tablename; *s; s++)
    if (*s == 'v' && !color)
      color = 8;
    else {
      for (int i = PAWN; i <= KING; i++)
        if (*s == PieceChar[i]) {
          pcs[i | color]++;
          pt[k++] = i | color;
          break;
        }
    }
  int numpcs = k;

  if (layout < 0 || layout > 2)
    layout = numpcs <= 6 ? 0 : numpcs == 7 ? 1 : 2;

  if (!color) exit(EXIT_FAILURE);

  if (numpcs < 3) {
    fprintf(stderr, "Need at least 3 pieces.\n");
    exit(EXIT_FAILURE);
  }

  if (pcs[WKING] != 1 || pcs[BKING] != 1) {
    fprintf(stderr, "Need one white king and one black king.\n");
    exit(EXIT_FAILURE);
  }

  if (pcs[WPAWN] || pcs[BPAWN]) {
    fprintf(stderr, "Can't handle pawns.\n");
    exit(EXIT_FAILURE);
  }

  symmetric = true;
  for (int i = PAWN; i <= KING; i++)
    symmetric = symmetric && pcs[i] == pcs[i + 8];

  g_num_threads = min(max(g_num_threads, 1), MAX_THREADS);
  printf("number of threads = %d\n", g_num_threads);

  // TODO: increase work units per thread as number of pieces increases.
  g_total_work = g_num_threads > 1 ? 1 * (g_num_threads + 0) : 1;

  g_pos.num = numpcs;

  static const int piece_order[16] = {
    0, 0, 3, 5, 7, 9, 1, 0,
    0, 0, 4, 6, 8, 10, 2, 0
  };

  for (int i = 0; i < numpcs; i++)
    for (int j = i + 1; j < numpcs; j++)
      if (piece_order[pt[i]] > piece_order[pt[j]])
        Swap(pt[i], pt[j]);

  for (int i = 0; i < numpcs; i++)
    g_pos.pt[i] = pt[i];

  k = 0;
  for (int i = 2; i < numpcs; i++)
    if (!(pt[i] & 0x08))
      g_pos.pcs[WHITE][k++] = i;
  g_pos.pcs[WHITE][k] = -1;

  k = 0;
  for (int i = 2; i < numpcs; i++)
    if (pt[i] & 0x08)
      g_pos.pcs[BLACK][k++] = i;
  g_pos.pcs[BLACK][k] = -1;

  // Initialize main IdxInfo struct.
  uint8_t mult[MAX_SETS] = { 0 };
  k = 0;
  for (int i = 2; i < numpcs;) {
    int j = i;
    for (; i < numpcs && pt[i] == pt[j]; i++)
      pc_to_set[i] = k;
    mult[k] = i - j;
    k++;
  }
  ri = rank_info[rank_mult(mult)];
  kslice_sizes[0] = ri.sizes[0];
  kslice_sizes[1] = ri.sizes[1];

#if 0
  {
    // Test the perfect ranking functions.
    uint8_t sq[MAX_PIECES];
    sq[0] = 0;
    sq[1] = 63;
    Bitboard occ = bit(sq[0]) | bit(sq[1]);
    uint64_t size = ri.sizes[1];
    printf("Testing perfecting indexing on %lu positions.\n", size);
    struct timespec start, end;
    timespec_get(&start, TIME_UTC);
    for (uint64_t idx = 0; idx < size; idx++) {
      unrank_reflection(idx, sq, occ, &ri);
      uint64_t rk = rank_reflection(sq, occ, ri.numsets, &ri);
      if (rk != idx) {
        printf("%lu != %lu\n", idx, rk);
        unrank_reflection(idx, sq, occ, &ri);
        uint64_t rk = rank_reflection(sq, occ, ri.numsets, &ri);
        printf("%lu != %lu\n", idx, rk);
        exit(EXIT_FAILURE);
      }
    }
    timespec_get(&end, TIME_UTC);
    printf("Test passed!\n");
    double elapsed = (double)(end.tv_sec - start.tv_sec)
        + 1e-9 * (double)(end.tv_nsec - start.tv_nsec);
    printf("elapsed: %.6f s\n", elapsed);
    printf("rate:    %.3f Miter/s\n", size / elapsed / 1e6);
    exit(0);
  }
#endif

  // Initialize RankInfo structs for running through positions with
  // a captured piece.
#if 0
  for (k = 0; k < ii.numsets; k++) {
    capt_ri[k] = ri;
    capt_ri[k].mult[k]--;
    if (capt_ri[k].mult[k] == 0) {
      for (int i = k + 1; i < ii.numsets; i++) {
        capt_ri[k].first[i - 1] = capt_ri[k].first[i];
        capt_ri[k].mult[i - 1] = capt_ri[k].mult[i];
        capt_ri[k].last[i - 1] = capt_ri[k].last[i];
      }
      capt_ii[k].numsets--;
    }
    calc_factors(&capt_ii[k]);
    kslice_sub_size[k] = capt_ii[k].size;
  }
#endif
  for (k = 0; k < ri.numsets; k++) {
    uint8_t mult2[MAX_SETS];
    memcpy(mult2, mult, sizeof mult2);
    mult2[k]--;
    if (mult2[k] == 0)
      for (int i = k + 1; i < ri.numsets; i++)
        mult[i - 1] = mult[i];
    *capt_ri = rank_info[rank_mult(mult2)];
    kslice_sub_size[k] = capt_ri[k].sizes[0];
  }

  kslice_setup();

  // Align work units on cache lines of 64 x 8 = 512 positions.
  work_g = create_work(g_total_work, kslice_sizes[0], 0x1ff);
  for (int i = 0; i < ri.numsets; i++)
    work_capt[i] = create_work(g_total_work, capt_ri[i].sizes[0], 0x1ff);
  init_generation_work();
  init_merge_work();

  if (workdir && *workdir)
    change_dir(workdir);

  make_dir(g_tablename);
  change_dir(g_tablename);

  generate();

  delete_intermediate_slices();

  kslice_free_buffers(); // Free memory but keep slice "-1".

#if 0
  printf("\n########## %s ##########\n", g_tablename);
  print_stats(WHITE);
  print_stats(BLACK);
#endif
  printf("\n");

  // Estimate sizes of different DTZ formats.
  double ewh = entropy_one_sided(WHITE);
  double ebl = entropy_one_sided(BLACK);
  double elo_w = entropy_loss_only(WHITE);
  double elo_b = entropy_loss_only(BLACK);
  double elo = elo_w + (symmetric ? 0.0 : elo_b);
  double ewi_w = entropy_win_only(WHITE);
  double ewi_b = entropy_win_only(WHITE);
  double ewi = ewi_w + (symmetric ? 0.0 : ewi_b);

  printf("entropy wtm  = %lf\n", ewh);
  printf("entropy btm  = %lf\n", ebl);
  if (!symmetric) {
    printf("entropy loss = %lf (%lf + %lf)\n", elo, elo_w, elo_b);
    printf("entropy win  = %lf (%lf + %lf)\n\n", ewi, ewi_w, ewi_b);
  } else {
    printf("entropy loss = %lf\n", elo);
    printf("entropy win  = %lf\n\n", ewi);
  }

  // Add 0.1% of overhead for having a table for each side.
  // Perhaps this should be higher?
  if (!symmetric) {
    elo *= 1.001;
    ewi *= 1.001;
  }

  // Determine the DTZ format to use.
  one_sided = !symmetric && min(ewh, ebl) < min(elo, ewi);
  wins_only = ewi <= elo;
  one_sided_stm = ewh <= ebl ? WHITE : BLACK;

  printf("DTZ format: %s only.\n\n",
        one_sided ? one_sided_stm == WHITE ? "white" : "black"
      : wins_only ? "wins" : "losses");

  if (file_exists("maxfens")) {
    FILE *F = file_open_read("maxfens");
    file_read(&mf, sizeof mf, F);
    fclose(F);
  } else {
    for (int stm = 0; stm < 2; stm++) {
      int n, m;
      for (n = MAX_STATS / 2 - 3; n > DRAW_RULE; n--)
        if (g_stats[stm][stats_n(n)])
          break;
      for (m = MAX_STATS / 2 - 3; m > n; m--)
        if (g_stats[stm ^ 1][MAX_STATS - 1 - m])
          break;
      mf.dtz[stm][1] = m > n ? 2 * m : n > DRAW_RULE ? 2 * n + 1 : -1;
      for (n = DRAW_RULE; n >= 1; n--)
        if (g_stats[stm][stats_n(n)])
          break;
      for (m = DRAW_RULE; m > n; m--)
        if (g_stats[stm ^ 1][MAX_STATS - 1 - m])
          break;
      mf.dtz[stm][0] = m > n ? 2 * m : n >= 1 ? 2 * n + 1 : -1;
    }
  }

  merge(WHITE);
  merge(BLACK);

  cleanup_generation();

  // Read out the files in "stats".
  collect_stats(WHITE);
  collect_stats(BLACK);

  size_t stats_file_len = strlen(g_output_dir) + strlen(g_tablename) + 6;
  char *stats_file = malloc(stats_file_len);
  if (!stats_file)
    out_of_mem();
  snprintf(stats_file, stats_file_len, "%s/%s.txt", g_output_dir, g_tablename);
  FILE *F = file_open_write(stats_file);
  fprintf(F, "########## %s ##########\n", g_tablename);
  print_stats(F, WHITE);
  print_stats(F, BLACK);
  fprintf(F, "\n");
  print_max_fens(F, &mf);
  fclose(F);
  file_rename(stats_file);
  free(stats_file);

  printf("\n########## %s ##########\n", g_tablename);
  print_stats(stdout, WHITE);
  print_stats(stdout, BLACK);
  printf("\n");
  print_max_fens(stdout, &mf);

  kslice_cleanup();

  if (!g_only_generate) {

    switch (layout) {
      case 0:
        join_slices(pcs, pt);
        break;
      case 1:
        join_slices_10();
        break;
      case 2:
        join_slices_462();
        break;
    }

  }

  if (g_cleanup) {
    change_dir("..");
    rmdir(g_tablename);
  }

  report_io();

  return 0;
}
