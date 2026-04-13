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
#include "generatep.h"
#include "indexp.h"
#include "joinp.h"
#include "kslicep.h"
#include "merge.h"
#include "movegen.h"
#include "probe.h"
#include "stats.h"
#include "tb8genp.h"
#include "threads.h"
#include "types.h"

#define TBPATH "RTBPATH"
#define STATSDIR "RTBSTATSDIR"

Position g_pos;
bool flipped = false;
bool g_only_generate, g_use_rans, symmetric, used_rans = false;
bool one_sided, wins_only;
int one_sided_stm;
char *g_tablename;
uint64_t *work_g, *work_g16, *work_capt[MAX_SETS];

const char *name[3] = { "wdl", "dtm", "dtz" };

static struct option options[] = {
  { "threads", 1, nullptr, 't' },
  { "stats", 0, nullptr, 's' },
  { "path", 1, nullptr, 'p' },
  { "rans", 0, nullptr, 'r' },
  { "layout", 1, nullptr, 'l' },
  { 0 }
};

int main(int argc, char **argv)
{
  int val, lindex;
  uint8_t pcs[16];
  uint8_t pt[8];
  int layout = -1;

  const char *path = getenv(TBPATH);
  g_num_threads = 1;

  while ((val = getopt_long(argc, argv, "at:gp:rl:", options, &lindex)) != -1)
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

#if 1
  if (layout < 0 || layout > 2)
    layout = numpcs <= 6 ? 1 : numpcs == 7 ? 1 : 2;
#endif

  if (!color) exit(EXIT_FAILURE);

  if (numpcs < 3) {
    fprintf(stderr, "Need at least 3 pieces.\n");
    exit(EXIT_FAILURE);
  }

  if (pcs[WKING] != 1 || pcs[BKING] != 1) {
    fprintf(stderr, "Need one white king and one black king.\n");
    exit(EXIT_FAILURE);
  }

  if (pcs[WPAWN] + pcs[BPAWN] != 1) {
    fprintf(stderr, "Expecting exactly one pawn.\n");
    exit(EXIT_FAILURE);
  }

  g_num_threads = min(max(g_num_threads, 1), MAX_THREADS);
  printf("number of threads = %d\n", g_num_threads);

  // TODO: increase work units per thread as number of pieces increases.
  g_total_work = g_num_threads > 1 ? 1 * (g_num_threads + 0) : 1;

  g_pos.num = numpcs;

  // Flip all piece colors if the pawn is white.
  if (pcs[WPAWN] > 0) {
    for (int i = 0; i < numpcs; i++)
      pt[i] ^= 0x08;
    flipped = true;
  }

  // Move the wK and bK to the front, followed by the black pawn.
  static const int piece_order[16] = {
    0, 0, 4, 6, 8, 10, 1, 0,
    0, 3, 5, 7, 9, 11, 2, 0
  };

  for (int i = 0; i < numpcs; i++)
    for (int j = i + 1; j < numpcs; j++)
      if (piece_order[pt[i]] > piece_order[pt[j]])
        Swap(pt[i], pt[j]);

  for (int i = 0; i < numpcs; i++)
    g_pos.pt[i] = pt[i];

  k = 0;
  for (int i = 3; i < numpcs; i++)
    if (!(pt[i] & 0x08))
      g_pos.pcs[WHITE][k++] = i;
  g_pos.pcs[WHITE][k] = -1;

  k = 0;
  for (int i = 3; i < numpcs; i++)
    if (pt[i] & 0x08)
      g_pos.pcs[BLACK][k++] = i;
  g_pos.pcs[BLACK][k] = -1;

  // Initialize main IdxInfo struct.
  k = 0;
  for (int i = 3; i < numpcs;) {
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
  work_g16 = create_work(g_total_work, k16slice_alloc_size << 3, 0x1ff);
  for (int i = 0; i < ii.numsets; i++)
    work_capt[i] = create_work(g_total_work, capt_ii[i].size, 0x1ff);

  make_dir(g_tablename);
  change_dir(g_tablename);
layout=0;
  for (int p = 0; p < 24; p++) {

    g_pos.sq[2] = InvPawnFlip[0][p];
    char str[128];
    sprintf(str, "%c%c", 'a' + (g_pos.sq[2] & 7), '1' + (g_pos.sq[2] >> 3));
    make_dir(str);
    change_dir(str);

    memset(g_stats, 0, sizeof g_stats);

    generate();

    uint64_t num_kslices = 0;
    for (int s = 0; s < 240; s++)
      for (int r = 0; r < 16; r++) {
        g_pos.sq[0] = KK16Square[s][r][0];
        g_pos.sq[1] = KK16Square[s][r][1];
        num_kslices += !is_broken(&g_pos);
      }

    for (int stm = 0; stm < 2; stm++) {
      // Remove some double counting.
      g_stats[stm][3] -= g_stats[stm][1] + g_stats[stm][2];
      g_stats[stm][DRAW_RULE + 5] -=
        g_stats[stm][DRAW_RULE + 3] - g_stats[stm][DRAW_RULE + 4];
      uint64_t total = 0;
      for (int i = 0; i < MAX_STATS; i++)
        total += g_stats[stm][i];
      g_stats[stm][MAX_STATS / 2 + 2] = num_kslices * kslice_size - total;
    }

    printf("\n########## %s - %d ##########\n", g_tablename, p);
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

    // Add 0.1% of overhead for having a table for each side.
    // Perhaps this should be higher?
    if (!symmetric) {
      elo *= 1.001;
      ewi *= 1.001;
    }

    // FIXME: store DTZ format for each pawn slice on disk for later use
    // by join_slices().

    // Determine the DTZ format to use.
    one_sided = !symmetric && min(ewh, ebl) <= min(elo, ewi);
    wins_only = ewi <= elo;
    one_sided_stm = ewh <= ebl ? WHITE : BLACK;

    printf("DTZ format: %s only.\n\n",
        one_sided ? one_sided_stm == WHITE ? "white" : "black"
        : wins_only ? "wins" : "losses");

    merge(WHITE);
    merge(BLACK);

    // We could now delete all files except those in stats/ and merged/,
    // or let merge() delete those.

    // Read out the files in "stats".
    collect_stats(WHITE);
    collect_stats(BLACK);

    printf("\n########## %s ##########\n", g_tablename);
    print_stats(WHITE);
    print_stats(BLACK);
    printf("\n");

    kslice_cleanup();

    // only compress wdl/dtz wtm/btm slices.
    switch (layout) {
    case 1:
      compress_slice_p();
      break;
    case 2:
      compress_slice_pk();
      break;
    }

    change_dir("..");
  }

  switch (layout) {
  case 1:
    join_slices_p();
    break;
  case 2:
    join_slices_pk();
    break;
  }

  report_io();

  return 0;
}
