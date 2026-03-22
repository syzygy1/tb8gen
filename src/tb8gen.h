/*
  Copyright (c) 2026 Ronald de Man

  This file is distributed under the terms of the GNU GPL, version 2.
*/

#ifndef TB8GEN_H
#define TB8GEN_H

enum { WL_BOTH = 0, WL_WTM, WL_BTM, W_ONLY, L_ONLY };
enum { LT_PIECE = 0, LT_PIECE_K, LT_PIECE_KK, LT_PAWN_FILE, LT_PAWN_RANK };

extern bool one_sided, wins_only;
extern int one_sided_stm;

extern size_t kslice_size, kslice_sub_size[MAX_SETS];
extern char *g_tablename;
extern bool symmetric;
extern Position g_pos;
extern uint64_t *work_g, *work_capt[MAX_SETS];
extern uint64_t g_stats[2][MAX_STATS];
extern bool g_use_rans;
extern bool used_rans;
extern const char *name[3];

static constexpr int g_num_pawns = 0;

#endif
