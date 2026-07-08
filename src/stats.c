/*
  Copyright (c) 2026 Ronald de Man

  This file is distributed under the terms of the GNU GPL, version 2.
*/

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "compress.h"
#include "defs.h"
#include "stats.h"
#include "types.h"
#include "util.h"
#include "hash/xxhash.h"

#ifdef HAS_PAWNS
#include "generatep.h"
#include "kslicep.h"

XXH128_hash_t dtz_partial_checksum[24];
#endif

uint64_t g_stats[2][MAX_STATS];
XXH128_hash_t wdl_checksum, dtz_checksum[2];

#include "stats_common.c"

void collect_stats(int stm)
{
  char str[64];
  uint64_t tmp[MAX_STATS];

  uint64_t *stats = g_stats[stm];

  memset(stats, 0, sizeof g_stats[stm]);

#ifndef HAS_PAWNS

  for (int s = 0; s < 462; s++) {
    create_name(str, s, stm, "stats", -1);
    FILE *F = file_open_read(str);
    read_data(F, tmp, sizeof tmp);
    fclose(F);
    for (int i = 0; i < MAX_STATS; i++)
      stats[i] += tmp[i];
  }

#else

  for (int s = 0; s < 240; s++)
    for (int r = 0; r < 16; r++) {
      g_slice.sq[0] = KK16Square[s][r][0];
      g_slice.sq[1] = KK16Square[s][r][1];

      if (is_broken(&g_slice))
        continue;

      create_name_sq(str, g_slice.sq[0], g_slice.sq[1], stm, "stats", -1);
      FILE *F = file_open_read(str);
      read_data(F, tmp, sizeof tmp);
      fclose(F);
      for (int i = 0; i < MAX_STATS; i++)
        stats[i] += tmp[i];
    }

#endif
}
