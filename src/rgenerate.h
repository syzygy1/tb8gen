/*
  Copyright (c) 2026 Ronald de Man

  This file is distributed under the terms of the GNU GPL, version 2.
*/

#ifndef RGENERATE_H
#define RGENERATE_H

#include <inttypes.h>

#include "defs.h"
#include "threads.h"
#include "types.h"

extern uint8_t *g_table[2];
extern size_t table_size, table_sub_size[MAX_SETS];
extern struct Work work_g_dynamic, work_g_static;

struct RIdxState {
  alignas(64) Bitboard occ[8];
  Bitboard bb[8];
  uint32_t sub[MAX_SETS + 1];
  uint8_t sq[2];
};

INLINE Bitboard ridx_state_inc(struct RIdxState *is, const struct RankInfo *ri)
{
  uint32_t *const restrict sub = is->sub;
  int i = ri->numsets;

  for (;;) {
    i--;
    if (++sub[i + 1] < ri->factor[i]) break;
    sub[i + 1] = 0;
    if (i == 0) {
      sub[0]++;
      is->sq[0] = KKSquare[sub[0]][0];
      is->sq[1] = KKSquare[sub[0]][1];
      is->bb[0] = is->occ[0] = bit(is->sq[0]) | bit(is->sq[1]);
      break;
    }
  }

  Bitboard occ = is->occ[i];
  for (; i < ri->numsets; i++) {
    occ |= is->bb[i + 1] = unrank_binomial(is->sub[i + 1], ri->mult[i], occ);
    is->occ[i + 1] = occ;
  }
  return occ;
}

INLINE Bitboard ridx_state_add(struct RIdxState *is, uint64_t v,
    const struct RankInfo *restrict ri)
{
  uint32_t *const restrict sub = is->sub;
  int i = ri->numsets;

  for (;;) {
    uint64_t s = (uint64_t)sub[(--i) + 1] + v;
    uint32_t f = ri->factor[i];

    if (s < f) {
      sub[i + 1] = s;
      break;
    }
    v = divmod_recip(s, f, ri->recip[i], &sub[i + 1]);

    if (i == 0) {
      sub[0] += v;
      is->sq[0] = KKSquare[sub[0]][0];
      is->sq[1] = KKSquare[sub[0]][1];
      is->bb[0] = is->occ[0] = bit(is->sq[0]) | bit(is->sq[1]);
      break;
    }
  }

  Bitboard occ = is->occ[i];
  for (; i < ri->numsets; i++) {
    occ |= is->bb[i + 1] = unrank_binomial(is->sub[i + 1], ri->mult[i], occ);
    is->occ[i + 1] = occ;
  }
  return occ;
}

INLINE bool ridx_state_legal(const struct RIdxState *is, int stm, Bitboard occ)
{
  int ksq = is->sq[stm ^ 1];
  for (int i = 0; g_sets[stm][i] >= 0; i++) {
    int k = g_sets[stm][i];
    Bitboard b = non_king_piece_attacks(g_set_pt[k], ksq, occ);
    if (b & is->bb[k + 1])
      return false;
  }
  return true;
}

INLINE void ridx_state_to_sq(const struct RIdxState *is, uint8_t *restrict sq,
    const struct RankInfo *ri)
{
  sq[0] = is->sq[0];
  sq[1] = is->sq[1];
  for (int k = 0; k < ri->numsets; k++) {
    Bitboard b = is->bb[k + 1];
    int i = ri->first[k];
    while (b)
      sq[i++] = pop_lsb(&b);
  }
}

INLINE int stats_n(int n)
{
  return 2 + n + 2 * (n > DRAW_RULE);
}

void generate(void);
void reset_bloss_captures_for_wdl(void);
void fix_bloss_after_wdl(int stm);
Bitboard ridx_state_init(struct RIdxState *is, uint64_t idx,
    const struct RankInfo *ri);
//void init_generation_work(void);
//void delete_intermediate_slices(void);
//void cleanup_generation(void);

#endif
