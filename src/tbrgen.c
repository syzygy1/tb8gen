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
#include "index.h"
//#include "join.h"
#include "movegen.h"
#include "probe.h"
#include "rgenerate.h"
#include "rstats.h"
#include "rwdl.h"
//#include "tb8gen.h"
#include "threads.h"
#include "types.h"
#include "util.h"

#define TBPATH "RTBPATH"
#define STATSDIR "RTBSTATSDIR"
#define WORKDIR "TB8DIR"

bool g_only_generate, g_use_rans, symmetric, used_rans = false;
bool g_cleanup;
bool one_sided, wins_only;
int one_sided_stm;
char *g_tablename;
char *g_output_dir;
struct MaxFen mf;

size_t table_size;
size_t table_diagonal;
size_t table_sub_size[MAX_SETS];
size_t kslice_sizes[2];

uint8_t *g_table[2];

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
    layout = numpcs <= 5 ? 0 : numpcs == 6 ? 1 : 2;

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

  // Initialize main RankInfo struct.
  uint8_t mult[MAX_SETS] = { 0 };
  k = 0;
  for (int i = 2; i < numpcs; k++) {
    int j = i;
    while (i < numpcs && pt[i] == pt[j])
      i++;
    mult[k] = i - j;
  }
  ri = rank_info_62[rank_mult(mult)];
  kslice_sizes[0] = ri.sizes[0];
  kslice_sizes[1] = ri.sizes[1];
  table_size = 462 * ri.sizes[0];
  table_diagonal = 441 * ri.sizes[0];

  g_table[0] = alloc_huge((table_size + 63) & ~0x3f);
  g_table[1] = alloc_huge((table_size + 63) & ~0x3f);
  if (!g_table[0] || !g_table[1])
    out_of_mem();

  // Initialize RankInfo structs for running through positions with
  // a captured piece.
  for (k = 0; k < ri.numsets; k++) {
    capt_ri[k] = ri;
    capt_ri[k].mult[k]--;
    calc_factors(&capt_ri[k], 62);
    table_sub_size[k] = 462 * capt_ri[k].sizes[0];
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

//  if (workdir && *workdir)
//    change_dir(workdir);

//  make_dir(g_tablename);
//  change_dir(g_tablename);

  generate();

  collect_stats(WHITE);
  collect_stats(BLACK);

  // Estimate sizes of different DTZ formats.
  double ewh = entropy_one_sided(WHITE);
  double ebl = entropy_one_sided(BLACK);
  double elo_w = entropy_loss_only(WHITE);
  double elo_b = entropy_loss_only(BLACK);
  double elo = elo_w + (symmetric ? 0.0 : elo_b);
  double ewi_w = entropy_win_only(WHITE);
  double ewi_b = entropy_win_only(BLACK);
  double ewi = ewi_w + (symmetric ? 0.0 : ewi_b);

  printf("\nentropy wtm  = %.1lf\n", ewh);
  printf("entropy btm  = %.1lf\n", ebl);
  if (!symmetric) {
    printf("entropy loss = %.1lf (%.1lf + %.1lf)\n", elo, elo_w, elo_b);
    printf("entropy win  = %.1lf (%.1lf + %.1lf)\n\n", ewi, ewi_w, ewi_b);
  } else {
    printf("entropy loss = %.1lf\n", elo);
    printf("entropy win  = %.1lf\n\n", ewi);
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
  printf("\n########## %s ##########\n", g_tablename);
  print_stats(stdout, WHITE);
  print_stats(stdout, BLACK);
  printf("\n");
  print_max_fens(stdout, &mf);

  calc_stats_checksums();

  if (!g_only_generate) {
    reset_bloss_captures_for_wdl();
    rjoin_wdl(pcs, pt, layout);
    fix_bloss_after_wdl(WHITE);
    fix_bloss_after_wdl(BLACK);
  }

#if 0
  if (file_exists("maxfens")) {
    FILE *F = file_open_read("maxfens");
    file_read(&mf, sizeof mf, F);
    fclose(F);
  } else {
    for (int stm = 0; stm < 2; stm++) {
      int n, m;
      for (n = MAX_STATS / 2 - 4; n > DRAW_RULE; n--)
        if (g_stats[stm][stats_n(n)])
          break;
      for (m = MAX_STATS / 2 - 4; m > n; m--)
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
#endif

#if 0
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
#endif

#if 0
  printf("\n########## %s ##########\n", g_tablename);
  print_stats(stdout, WHITE);
  print_stats(stdout, BLACK);
  printf("\n");
#endif

#if 0
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
#endif

  return 0;
}
