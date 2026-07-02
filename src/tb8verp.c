/*
  Copyright (c) 2026 Ronald de Man

  This file is distributed under the terms of the GNU GPL, version 2.
*/

#include <getopt.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "defs.h"
#include "index.h"
#include "kslicep.h"
#include "movegen.h"
#include "probe.h"
#include "stats.h"
#include "tb8gen.h"
#include "threads.h"
#include "types.h"
#include "util.h"
#include "verify.h"

#define TBPATH "RTBPATH"
#define WORKDIR "TB8DIR"

Position g_pos;
char pawnstr[24][3];
int8_t g_sets[2][8];
int8_t g_piece_set[2][8];
uint8_t g_set_pt[8];
bool flipped = false;
bool g_cleanup, symmetric;
bool one_sided, wins_only;
int one_sided_stm;
char *g_tablename;
char *g_output_dir;
struct Work *work_g;
struct MaxFen mf;

const char *typename[3] = { "wdl", "dtm", "dtz" };

static struct option options[] = {
  { "threads", 1, nullptr, 't' },
  { "path", 1, nullptr, 'p' },
  { "workdir", 1, nullptr, 'w' },
  { 0 }
};

int main(int argc, char **argv)
{
  int val, lindex;
  uint8_t pcs[16] = { 0 };
  uint8_t pt[8];

  const char *path = getenv(TBPATH);
  const char *workdir = getenv(WORKDIR);
  g_num_threads = 1;

  while ((val = getopt_long(argc, argv, "at:p:cw:", options, &lindex)) != -1)
    switch (val) {
    case 'a':
      g_thread_affinity = true;
      break;
    case 't':
      g_num_threads = atoi(optarg);
      break;
    case 'p':
      path = optarg;
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

  if (pcs[WPAWN] + pcs[BPAWN] != 1) {
    fprintf(stderr, "Expecting exactly one pawn.\n");
    exit(EXIT_FAILURE);
  }

  g_num_threads = min(max(g_num_threads, 1), MAX_THREADS);
  printf("number of threads = %d\n", g_num_threads);

  g_pos.num = numpcs;

  if (pcs[WPAWN] > 0) {
    for (int i = 0; i < numpcs; i++)
      pt[i] ^= 0x08;
    flipped = true;
  }

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

  // Initialize main RankInfo struct.
  uint8_t mult[MAX_SETS] = { 0 };
  k = 0;
  for (int i = 3; i < numpcs;) {
    int j = i;
    for (; i < numpcs && pt[i] == pt[j]; i++)
      pc_to_set[i] = k;
    mult[k] = i - j;
    k++;
  }
  ri = rank_info_61[rank_mult(mult)];
  kslice_size = ri.sizes[0];

  // Initialize RankInfo structs for running through positions with
  // a captured piece.
  for (k = 0; k < ri.numsets; k++) {
    capt_ri[k] = ri;
    capt_ri[k].mult[k]--;
    calc_factors(&capt_ri[k], 61);
    kslice_sub_size[k] = capt_ri[k].sizes[0];
  }

  for (int i = 0; i < ri.numsets; i++)
    g_set_pt[i] = g_pos.pt[ri.first[i]];

  memset(g_piece_set, -1, sizeof g_piece_set);
  for (int i = 0; i < ri.numsets; i++)
    g_piece_set[g_set_pt[i] >> 3][g_set_pt[i] & 7] = i;

  for (int stm = 0; stm < 2; stm++) {
    int k = 0;
    for (int i = 0; i < ri.numsets; i++)
      if ((g_set_pt[i] >> 3) == stm)
        g_sets[stm][k++] = i;
    g_sets[stm][k] = -1;
  }

  kslice_setup();
  kslice_alloc_buffers();

  // Align work units on cache lines of 64 x 8 = 512 positions.
  work_g = create_work(g_total_work, kslice_size, 0x1ff);
  init_verification_work();

  for (int i = 0; i < 24; i++) {
    pawnstr[i][0] = 'a' + (i / 6);
    pawnstr[i][1] = '1' + (((i % 6) + 1) ^ (flipped ? 7 : 0));
  }

  if (workdir && *workdir)
    change_dir(workdir);

  char verify_dir[64];
  snprintf(verify_dir, sizeof verify_dir, "v%s", g_tablename);
  make_dir(verify_dir);
  change_dir(verify_dir);

  verify();

  if (g_cleanup) {
    change_dir("..");
    rmdir(verify_dir);
  }

  report_io();

  return 0;
}
