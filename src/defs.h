/*
  Copyright (c) 2011-2013, 2026 Ronald de Man

  This file is distributed under the terms of the GNU GPL, version 2.
*/

#ifndef DEFS_H
#define DEFS_H

static constexpr int MAX_PIECES = 8;
static constexpr int MAX_SETS = MAX_PIECES - 2;
static constexpr int MAX_STATS = 2560;
static constexpr int DRAW_RULE = 2 * 50;
static constexpr int MAX_VAL = (MAX_STATS / 2 - DRAW_RULE) / 2;
static constexpr int MAX_THREADS = 32;

static constexpr size_t MIN_STREAMING_STORES_SIZE = (size_t)1024*1024;

#endif
