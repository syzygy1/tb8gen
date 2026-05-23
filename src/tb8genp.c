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
#include "generatep.h"
#include "index.h"
#include "joinp.h"
#include "kslicep.h"
#include "merge.h"
#include "movegen.h"
#include "probe.h"
#include "stats.h"
#include "tb8genp.h"
#include "threads.h"
#include "types.h"
#include "util.h"

#define TBPATH "RTBPATH"
#define STATSDIR "RTBSTATSDIR"
#define WORKDIR "TB8DIR"

Position g_pos;
bool flipped = false;
bool g_only_generate, g_use_rans, symmetric, used_rans = false;
bool g_cleanup;
bool one_sided, wins_only;
int one_sided_stm;
char *g_tablename;
char *g_output_dir;
struct Work *work_g, *work_g16, *work_capt[MAX_SETS];
struct DtzFormat dtz_format[24];
struct MaxFen mf, mmf;

const char *typename[3] = { "wdl", "dtm", "dtz" };

char pawnstr[24][3];

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
  int layout = 0;

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
  init_unrank();
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

  if (layout <= 0 || layout > 2)
    layout = numpcs <= 7 ? 1 : numpcs == 8 ? 1 : 2;

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
    if (capt_ii[k].mult[k] == 0) {
      for (int i = k + 1; i < ii.numsets; i++) {
        capt_ii[k].first[i - 1] = capt_ii[k].first[i];
        capt_ii[k].mult[i - 1] = capt_ii[k].mult[i];
        capt_ii[k].last[i - 1] = capt_ii[k].last[i];
      }
      capt_ii[k].numsets--;
    }
    calc_factors(&capt_ii[k]);
    kslice_sub_size[k] = capt_ii[k].size;
  }

  kslice_setup();

  // Align work units on cache lines of 64 x 8 = 512 positions.
  work_g = create_work(g_total_work, kslice_size, 0x1ff);
  work_g16 = create_work(g_total_work, k16slice_alloc_size << 3, 0x1ff);
  for (int i = 0; i < ii.numsets; i++)
    work_capt[i] = create_work(g_total_work, capt_ii[i].size, 0x1ff);
  init_generation_work();

  for (int i = 0; i < 24; i++) {
    pawnstr[i][0] = 'a' + (i / 6);
    pawnstr[i][1] = '1' + (((i % 6) + 1) ^ (flipped ? 7 : 0));
  }

  if (workdir && *workdir)
    change_dir(workdir);

  make_dir(g_tablename);
  change_dir(g_tablename);
  make_dir("stats");

  for (int q = 0; q < 24; q++) {
    g_pos.sq[2] = InvPawnFlip[0][q];
    make_dir(pawnstr[q]);
    change_dir(pawnstr[q]);

    memset(g_stats, 0, sizeof g_stats);

    generate();

    delete_intermediate_slices();

    // Estimate sizes of different DTZ formats.
    double ewh, ebl, elo, ewi;
    ewh = entropy_one_sided(WHITE);
    ebl = entropy_one_sided(BLACK);
    elo = entropy_loss_only(WHITE) + (symmetric? 0.0 : entropy_loss_only(BLACK));
    ewi = entropy_win_only(WHITE) + (symmetric ? 0.0 : entropy_win_only(BLACK));

    printf("\nentropy wtm  = %lf\n", ewh);
    printf("entropy btm  = %lf\n", ebl);
    printf("entropy loss = %lf\n", elo);
    printf("entropy win  = %lf\n\n", ewi);

    // Add 0.1% of overhead for having a table for each side.
    // Perhaps this should be higher?
    if (!symmetric) {
      elo *= 1.001;
      ewi *= 1.001;
    }

    // Determine the DTZ format to use.
    one_sided = !symmetric && min(ewh, ebl) <= min(elo, ewi);
//    one_sided = true;
    wins_only = ewi <= elo;
    one_sided_stm = ewh <= ebl ? WHITE : BLACK;

    printf("DTZ format: %s only.\n\n",
        one_sided ? one_sided_stm == WHITE ? "white" : "black"
        : wins_only ? "wins" : "losses");

    dtz_format[q] = (struct DtzFormat){
      .one_sided = one_sided, .wins_only = wins_only,
      .one_sided_stm = one_sided_stm
    };

    memset(&mf, 0, sizeof mf);

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

    merge(WHITE);
    merge(BLACK);

    cleanup_generation();

    // Read out the files in "stats".
    collect_stats(WHITE);
    collect_stats(BLACK);

    fprintf(stdout, "\n########## %s - %s ##########\n", g_tablename,
        pawnstr[q]);
    print_stats(stdout, WHITE);
    print_stats(stdout, BLACK);
    fprintf(stdout, "\n");
    print_max_fens(stdout, &mf);

    // Keep track of max DTZ and corresponding FENs across pawn slices.
    for (int stm = 0; stm < 2; stm++)
      for (int i = 0; i < 2; i++)
        if (mf.dtz[stm][i] > mmf.dtz[stm][i]) {
          mmf.dtz[stm][i] = mf.dtz[stm][i];
          strcpy(mmf.fen[stm][i], mf.fen[stm][i]);
        }

    char str[64];
    sprintf(str, "../stats/%s", pawnstr[q]);
    FILE *F = file_open_write(str);
    write_data(F, g_stats, sizeof g_stats);
    fclose(F);
    file_rename(str);

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

    if (!g_cleanup) continue;

    if (g_pos.sq[2] >= 16 && g_pos.sq[2] < 40)
      rmdir(pawnstr[q - 1]);
    else {
      rmdir(pawnstr[q - 2]);
      rmdir(pawnstr[q - 1]);
      rmdir(pawnstr[q]);
    }
  }

  switch (layout) {
  case 1:
    join_slices_p();
    break;
  case 2:
    join_slices_pk();
    break;
  }

  memset(g_stats, 0, sizeof g_stats);
  uint64_t tmp[2][MAX_STATS];
  for (int q = 0; q < 24; q++) {
    char str[64];
    sprintf(str, "stats/%s", pawnstr[q]);
    FILE *F = file_open_read(str);
    read_data(F, tmp, sizeof tmp);
    fclose(F);
    for (int stm = 0; stm < 2; stm++)
      for (int i = 0; i < MAX_STATS; i++)
        g_stats[stm][i] += tmp[stm][i];
  }

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
  print_max_fens(F, &mmf);
  fclose(F);
  file_rename(stats_file);
  free(stats_file);

  printf("\n########## %s ##########\n", g_tablename);
  print_stats(stdout, WHITE);
  print_stats(stdout, BLACK);
  printf("\n");
  print_max_fens(stdout, &mmf);

  if (g_cleanup) {
    char str[64];
    for (int q = 0; q < 24; q++) {
      int psq = InvPawnFlip[0][q];
      for (int k1 = 0; k1 < 64; k1++) {
        if (k1 == psq) continue;
        for (int k2 = 0; k2 < 64; k2++) {
          if (k2 == psq || (king_mask(k1) & bit(k2))) continue;
          for (int stm = 0; stm < 2; stm++) {
            sprintf(str, "%s/", pawnstr[q]);
            create_name_sq(str + strlen(str), k1, k2, stm, "stats", -1);
            unlink(str);
          }
        }
      }
      sprintf(str, "%s/stats", pawnstr[q]);
      delete_dir(-1, str);
      sprintf(str, "%s/merge_info.w", pawnstr[q]);
      unlink(str);
      sprintf(str, "%s/merge_info.b", pawnstr[q]);
      unlink(str);
      sprintf(str, "%s/maxfens", pawnstr[q]);
      unlink(str);
      sprintf(str, "%s/generate_info", pawnstr[q]);
      unlink(str);
      sprintf(str, "stats/%s", pawnstr[q]);
      unlink(str);
      rmdir(pawnstr[q]);
    }
    rmdir("stats");
    change_dir("..");
    rmdir(g_tablename);
  }

  report_io();

  return 0;
}
