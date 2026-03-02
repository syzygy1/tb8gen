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

#include "defs.h"
#include "generate.h"
#include "index.h"
#include "kslice.h"
#include "merge.h"
#include "movegen.h"
#include "probe.h"
#include "stats.h"
#include "tb8gen.h"
#include "threads.h"
#include "types.h"

#if 0
#include "compress.h"
#include "permute.h"
#endif

#define TBPATH "RTBPATH"
#define STATSDIR "RTBSTATSDIR"

Position g_pos;
bool g_only_generate, g_use_rans, symmetric;
bool one_sided, wins_only;
int one_sided_stm;
char *g_tablename;
uint64_t *work_g, *work_capt[MAX_SETS];

static struct option options[] = {
  { "threads", 1, nullptr, 't' },
  { "stats", 0, nullptr, 's' },
  { "path", 1, nullptr, 'p' },
//  { "rans", 0, nullptr, 'r' },
  { 0 }
};

int main(int argc, char **argv)
{
  int val, lindex;
  uint8_t pcs[16];
  uint8_t pt[8];

  const char *path = getenv(TBPATH);
  g_num_threads = 1;

  while ((val = getopt_long(argc, argv, "at:gp:r", options, &lindex)) != -1)
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
    }

  if (optind >= argc) {
    printf("No tablebase specified.\n");
    exit(EXIT_FAILURE);
  }
  g_tablename = argv[optind];

  init_tablebases(path);
  init_movegen();
  init_tables();
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

  g_num_threads = max(g_num_threads, 1);
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
  for (int i = 0; i < numpcs; i++)
    if (!(pt[i] & 0x08))
      g_pos.pcs[WHITE][k++] = i;
  g_pos.pcs[WHITE][k] = -1;

  k = 0;
  for (int i = 0; i < numpcs; i++)
    if (pt[i] & 0x08)
      g_pos.pcs[BLACK][k++] = i;
  g_pos.pcs[BLACK][k] = -1;

  // Initialize main IdxInfo struct.
  k = 0;
  for (int i = 2; i < numpcs;) {
    int j = i;
    for (; i < numpcs && pt[i] == pt[j]; i++)
      pc_to_set[i] = k;
    ii.first[k] = j;
    ii.mult[k] = i - j;
    ii.last[k] = i - 1;
    k++;
  }
  ii.numsets = k;
  calc_factors(&ii);
  kslice_size = ii.size;

  // Initialize IdxInfo structs for running through positions with
  // a captured piece.
  for (k = 0; k < ii.numsets; k++) {
    capt_ii[k] = ii;
    capt_ii[k].mult[k]--;
    calc_factors(&capt_ii[k]);
    kslice_sub_size[k] = capt_ii[k].size;
  }

  kslice_setup();

  // Align work units on cache lines of 64 x 8 = 512 positions.
  work_g = create_work(g_total_work, kslice_size, 0x1ff);
  for (int i = 0; i < ii.numsets; i++)
    work_capt[i] = create_work(g_total_work, capt_ii[i].size, 0x1ff);

  make_dir(g_tablename);
  change_dir(g_tablename);

  generate();

  kslice_free_buffers(); // Free memory but keep slice "-1".

  for (int stm = 0; stm < 2; stm++) {
    // Remove some double counting.
    g_stats[stm][2] -= g_stats[stm][1];
    g_stats[stm][3 + DRAW_RULE] -= g_stats[stm][2 + DRAW_RULE];
    uint64_t tot = 0;
    for (int i = 0; i < MAX_STATS; i++)
      tot += g_stats[stm][i];
    g_stats[stm][MAX_STATS / 2 + 1] = 462 * kslice_size - tot;
  }

  printf("\n########## %s ##########\n", g_tablename);
  print_stats(WHITE);
  print_stats(BLACK);
  printf("\n");

  // Estimate sizes of different DTZ formats.
  double ewh, ebl, elo, ewi;
  ewh = entropy_one_sided(WHITE);
  ebl = entropy_one_sided(BLACK);
  elo = entropy_loss_only(WHITE) + (symmetric? 0.0 : entropy_loss_only(BLACK));
  ewi = entropy_win_only(WHITE) + (symmetric ? 0.0 : entropy_win_only(BLACK));

  printf("entropy wtm  = %lf\n", ewh);
  printf("entropy btm  = %lf\n", ebl);
  printf("entropy loss = %lf\n", elo);
  printf("entropy win  = %lf\n\n", ewi);

  // Add 1% of overhead for having a table for each side.
  // Perhaps this should be higher.
  if (!symmetric) {
    elo *= 1.01;
    ewi *= 1.01;
  }

  // Determine the DTZ format to use.
  one_sided = !symmetric && min(ewh, ebl) < min(elo, ewi);
  wins_only = ewi < elo;
  one_sided_stm = ewh < ebl ? WHITE : BLACK;

  printf("DTZ format: %s only.\n\n",
        one_sided ? one_sided_stm == WHITE ? "white" : "black"
      : wins_only ? "wins" : "losses");

  merge(WHITE);
  merge(BLACK);

  // Read out the files in "stats".
  collect_stats(WHITE);
  collect_stats(BLACK);

  printf("########## %s ##########\n", g_tablename);
  print_stats(WHITE);
  print_stats(BLACK);
  printf("\n");

  printf("entropy wtm  = %lf\n", entropy_one_sided(WHITE));
  printf("entropy btm  = %lf\n", entropy_one_sided(BLACK));
  printf("entropy loss = %lf\n",
      entropy_loss_only(WHITE) + (symmetric? 0.0 : entropy_loss_only(BLACK)));
  printf("entropy win  = %lf\n\n",
      entropy_win_only(WHITE) + (symmetric ? 0.0 : entropy_win_only(BLACK)));

  kslice_cleanup();

  report_io();
}
