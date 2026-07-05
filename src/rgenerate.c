/*
  Copyright (c) 2026 Ronald de Man

  This file is distributed under the terms of the GNU GPL, version 2.
*/

#include <inttypes.h>
#include <stdint.h>
#include <stdlib.h>
#include <immintrin.h>

#include "defs.h"
#include "index.h"
#include "movegen.h"
#include "probe.h"
#include "rgenerate.h"
#include "threads.h"
#include "types.h"
#include "util.h"

#if 0
static constexpr uint8_t LOSS_IN_0  = 0x01;
static constexpr uint8_t CAPT_LOSS  = 0x02;
static constexpr uint8_t CAPT_BLOSS = 0;
static constexpr uint8_t CAPT_DRAW  = 0;
static constexpr uint8_t CAPT_CWIN  = 0;
static constexpr uint8_t CAPT_WIN   = 0xfe;
static constexpr uint8_t ILLEGAL    = 0xff;
#endif

struct Work work_g_dynamic, work_g_static;
static struct Work work_capt_dynamic[MAX_SETS];

static constexpr uint64_t GENERATE_MIN_CHUNK = 1ULL << 9;
static constexpr int GENERATE_DYNAMIC_FACTOR = 4;

void init_generation_work(void)
{
  work_init(&work_g_dynamic, table_size, 0x3f, WORK_DYNAMIC,
      GENERATE_DYNAMIC_FACTOR, GENERATE_MIN_CHUNK);
  work_init(&work_g_static, table_size, 0x3f, WORK_STATIC, 1,
      GENERATE_MIN_CHUNK);

  for (int k = 0; k < ri.numsets; k++)
    work_init(&work_capt_dynamic[k], table_sub_size[k], 0x3f, WORK_DYNAMIC,
        GENERATE_DYNAMIC_FACTOR, GENERATE_MIN_CHUNK);
}

struct RIdxState {
  alignas(64) Bitboard occ[8];
  Bitboard bb[8];
  uint32_t sub[MAX_SETS + 1];
  uint8_t sq[2];
};

Bitboard ridx_state_init(struct RIdxState *is, uint64_t idx,
    const struct RankInfo *ri)
{
  for (int k = ri->numsets - 1; k >= 0; k--)
    idx = divmod_recip(idx, ri->factor[k], ri->recip[k], &is->sub[k + 1]);
  is->sub[0] = idx;
  is->sq[0] = KKSquare[idx][0];
  is->sq[1] = KKSquare[idx][1];
  is->bb[0] = is->occ[0] = bit(is->sq[0]) | bit(is->sq[1]);
  for (int i = 0; i < ri->numsets; i++) {
    is->bb[i + 1] = unrank_binomial(is->sub[i + 1], ri->mult[i], is->occ[i]);
    is->occ[i + 1] = is->occ[i] | is->bb[i + 1];
  }
  return is->occ[ri->numsets];
}

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

INLINE void set_max_atomic(uint8_t *restrict p, uint8_t v)
{
  _Atomic uint8_t *restrict q = (_Atomic uint8_t *)p;
  uint8_t old = atomic_load_explicit(q, memory_order_relaxed);
  while (    old < v
         && !atomic_compare_exchange_weak_explicit(q, &old, v,
           memory_order_relaxed, memory_order_relaxed)) ;
}

INLINE void mark_king_uncaptures(int stm, int k, uint8_t *restrict p,
    Bitboard occ, struct RIdxState *is, uint8_t v)
{
  alignas(64) Bitboard bb[8];
  uint8_t ksq[2] = { is->sq[0], is->sq[1] };

  // The stm king uncaptures, so we need to place a piece where the king was.
  is->bb[k + 1] |= bit(ksq[stm]);

  Bitboard b = king_attacks(ksq[stm]) & ~king_attacks(ksq[stm ^ 1]) & ~occ;
  while (b) {
    ksq[stm] = pop_lsb(&b);
    is->bb[0] = bit(ksq[0]) | bit(ksq[1]);
    int s = KKMap[ksq[0]][ksq[1]];
    int t = KK_transform[ksq[0]][ksq[1]];
    Bitboard *set_bb = is->bb;
    if (t) {
      transform_set_bb(t, is->bb, bb);
      set_bb = bb;
    }
    uint64_t idx = rank_bb_from(set_bb, s, 0, set_bb[0], &ri);
    set_max_atomic(p + idx, v);
    if (is->sub[0] < 441 && s >= 441) {
      __m512i x = _mm512_load_si512((__m512i *)set_bb);
      x = flip_main_8xbb(x);
      _mm512_store_si512((__m512i *)bb, x);
      idx = rank_bb_from(bb, s, 0, bb[0], &ri);
      set_max_atomic(p + idx, v);
    }
  }

  // Restore the set occupancy bitboards.
  is->bb[0] = is->occ[0];
  is->bb[k + 1] ^= bit(is->sq[stm]);
}

// Uncapture a piece in set k by a piece in set j.
INLINE void mark_uncaptures(int stm, const int pc, int k,
   uint8_t *restrict p, Bitboard occ, struct RIdxState *is, uint8_t v)
{
  int j = g_piece_set[stm][pc];
  if (j < 0) return;

  int l = min(j, k);
  uint64_t idx0 = is->sub[0];
  for (int i = 0; i < l; i++)
    idx0 = idx0 * ri.factor[i] + is->sub[i + 1];

  Bitboard bb = is->bb[j + 1];
  while (bb) {
    Bitboard to_bb = bb & -bb;
    is->bb[j + 1] ^= to_bb;
    is->bb[k + 1] ^= to_bb;
    bb ^= to_bb;
    Bitboard b = piece_moves(pc, lsb(to_bb), occ);
    while (b) {
      Bitboard from_bb = b & -b;
      is->bb[j + 1] ^= from_bb;
      uint64_t idx = rank_bb_from(is->bb, idx0, l, is->occ[l], &ri);
      set_max_atomic(p + idx, v);
      is->bb[j + 1] ^= from_bb;
      b ^= from_bb;
    }
    is->bb[j + 1] ^= to_bb;
    is->bb[k + 1] ^= to_bb;
  }
}

INLINE void mark_uncapture_king(int stm, const int pc, int ksq,
    uint8_t *restrict p, Bitboard occ, struct RIdxState *is)
{
  int k = g_piece_set[stm][pc];
  if (k < 0) return;

  uint64_t idx0 = is->sub[0];
  for (int i = 0; i < k; i++)
    idx0 = idx0 * ri.factor[i] + is->sub[i + 1];

  Bitboard b = piece_moves(pc, ksq, occ);
  while (b) {
    Bitboard from_bb = b & -b;
    is->bb[k + 1] ^= from_bb;
    uint64_t idx = rank_bb_from(is->bb, idx0, k, is->occ[k], &ri);
    p[idx] = 0xff;
    is->bb[k + 1] ^= from_bb;
    b ^= from_bb;
  }
}

INLINE void mark_king_unmoves(int stm, uint8_t *restrict p, Bitboard occ,
    struct RIdxState *is, uint8_t v)
{
  alignas(64) Bitboard bb[8];
  uint8_t ksq[2] = { is->sq[0], is->sq[1] };

  Bitboard b = king_attacks(ksq[stm]) & ~king_attacks(ksq[stm ^ 1]) & ~occ;
  while (b) {
    ksq[stm] = pop_lsb(&b);
    is->bb[0] = bit(ksq[0]) | bit(ksq[1]);
    int s = KKMap[ksq[0]][ksq[1]];
    int t = KK_transform[ksq[0]][ksq[1]];
    Bitboard *set_bb = is->bb;
    if (t) {
      transform_set_bb(t, is->bb, bb);
      set_bb = bb;
    }
    uint64_t idx = rank_bb_from(set_bb, s, 0, set_bb[0], &ri);
    set_max_atomic(p + idx, v);
    if (is->sub[0] < 441 && s >= 441) {
      __m512i x = _mm512_load_si512((__m512i *)set_bb);
      x = flip_main_8xbb(x);
      _mm512_store_si512((__m512i *)bb, x);
      idx = rank_bb_from(bb, s, 0, bb[0], &ri);
      set_max_atomic(p + idx, v);
    }
  }

  is->bb[0] = is->occ[0];
}

INLINE void mark_unmoves(int stm, const int pc, uint8_t *restrict const p,
    Bitboard occ, struct RIdxState *is, uint8_t v)
{
  int k = g_piece_set[stm][pc];
  if (k < 0) return;

  uint64_t idx0 = is->sub[0];
  for (int i = 0; i < k; i++)
    idx0 = idx0 * ri.factor[i] + is->sub[i + 1];

  Bitboard bb = is->bb[k + 1];
  while (bb) {
    Bitboard from_bb = bb & -bb;
    is->bb[k + 1] ^= from_bb;
    bb ^= from_bb;
    Bitboard b = piece_moves(pc, lsb(from_bb), occ);
    while (b) {
      Bitboard to_bb = b & -b;
      is->bb[k + 1] ^= to_bb;
      uint64_t idx = rank_bb_from(is->bb, idx0, k, is->occ[k], &ri);
      set_max_atomic(p + idx, v);
      is->bb[k + 1] ^= to_bb;
      b ^= to_bb;
    }
    is->bb[k + 1] ^= from_bb;
  }
}

static uint8_t loss_byte;

INLINE bool check_king_moves(int stm, uint8_t *restrict p, Bitboard occ,
    struct RIdxState *is)
{
  alignas(64) Bitboard bb[8];
  uint8_t ksq[2] = { is->sq[0], is->sq[1] };

  Bitboard b = king_attacks(ksq[stm]) & ~king_attacks(ksq[stm ^ 1]) & ~occ;
  while (b) {
    ksq[stm] = pop_lsb(&b);
    is->bb[0] = bit(ksq[0]) | bit(ksq[1]);
    int s = KKMap[ksq[0]][ksq[1]];
    int t = KK_transform[ksq[0]][ksq[1]];
    Bitboard *set_bb = is->bb;
    if (t) {
      transform_set_bb(t, is->bb, bb);
      set_bb = bb;
    }
    uint64_t idx = rank_bb_from(set_bb, s, 0, set_bb[0], &ri);
    if (p[idx] < loss_byte) {
      is->bb[0] = is->occ[0];
      return false;
    }
  }

  is->bb[0] = is->occ[0];
  return true;
}

INLINE bool check_moves(int stm, const int pc, uint8_t *restrict p,
    Bitboard occ, struct RIdxState *is)
{
  int k = g_piece_set[stm][pc];
  if (k < 0) return true;

  uint64_t idx0 = is->sub[0];
  for (int i = 0; i < k; i++)
    idx0 = idx0 * ri.factor[i] + is->sub[i + 1];

  Bitboard bb = is->bb[k + 1];
  while (bb) {
    Bitboard from_bb = bb & -bb;
    is->bb[k + 1] ^= from_bb;
    bb ^= from_bb;

    Bitboard b = piece_moves(pc, lsb(from_bb), occ);
    while (b) {
      Bitboard to_bb = b & -b;
      is->bb[k + 1] ^= to_bb;
      uint64_t idx = rank_bb_from(is->bb, idx0, k, is->occ[k], &ri);
      is->bb[k + 1] ^= to_bb;
      if (p[idx] < loss_byte) {
        is->bb[k + 1] ^= from_bb;
        return false;
      }
      b ^= to_bb;
    }
    is->bb[k + 1] ^= from_bb;
  }

  return true;
}

#if 0
INLINE uint8_t check_loss(int stm, uint8_t *restrict p, Bitboard occ,
    struct RIdxState *is)
{
  int k;
  uint8_t v = 0xff;
  if ((k = g_piece_set[stm][QUEEN]) >= 0) {
    v = check_moves(stm, QUEEN, k, p, occ, is);
    if (v == 0) return 0;
  }
  if ((k = g_piece_set[stm][ROOK]) >= 0) {
    v = best(v, check_moves(stm, ROOK, k, p, occ, is));
    if (v == 0) return 0;
  }
  if ((k = g_piece_set[stm][BISHOP]) >= 0) {
    v = best(v, check_moves(stm, BISHOP, k, p, occ, is));
    if (v == 0) return 0;
  }
  if ((k = g_piece_set[stm][KNIGHT]) >= 0) {
    v = best(v, check_moves(stm, KNIGHT, k, p, occ, is));
    if (v == 0) return 0;
  }
  // king unmoves
}
#endif

static constexpr uint8_t wdl_to_capt_val[5] = {
//  CAPT_WIN, CAPT_CWIN, CAPT_DRAW, CAPT_BLOSS, CAPT_LOSS
  254, 253 - DRAW_RULE - 1, 127, 2 + DRAW_RULE, 2
};

static void clear_table_worker(struct ThreadData *thread)
{
  uint8_t *restrict p = g_table[g_pos.stm];

  size_t begin = thread->begin, end = thread->end;

#ifdef __AVX512F__

  __m512i z = _mm512_setzero_si512();
  if (end - begin >= MIN_STREAMING_STORES_SIZE) {
    for (size_t idx = begin; idx < end; idx += 64)
      _mm512_stream_si512((__m512i *)(p + idx), z);
  } else {
    for (size_t idx = begin; idx < end; idx += 64)
      _mm512_store_si512((__m512i *)(p + idx), z);
  }

#elifdef __AVX2__

  __m256i z = _mm256_setzero_si256();
  if (end - begin >= MIN_STREAMING_STORES_SIZE) {
    for (size_t idx = begin; idx < end; idx += 64) {
      _mm256_stream_si256((__m256i *)(p + idx), z);
      _mm256_stream_si256((__m256i *)(p + idx + 32), z);
    }
  } else {
    for (size_t idx = begin; idx < end; idx += 64) {
      _mm256_store_si256((__m256i *)(p + idx), z);
      _mm256_store_si256((__m256i *)(p + idx + 32), z);
    }
  }

#else

  memset(p + begin, 0, end - begin);

#endif
}

static void clear_table(void)
{
  for (int stm = 0; stm < 2; stm++) {
    g_pos.stm = stm;
    run_threaded(clear_table_worker, &work_g_static);
  }
}

static int work_set;

// TODO: keep track of cursedness
// or just do a quick iteration?
static void calc_capt_worker(struct ThreadData *thread)
{
  struct RIdxState is;
  int k = work_set;

  Position pos = g_pos;
  int stm = pos.stm ^ 1;
  int n = --pos.num;
  int m = ri.last[k];
  pos.pt[m] = pos.pt[n];

  uint8_t *restrict p = g_table[stm];

  Bitboard occ = ridx_state_init(&is, thread->begin, &capt_ri[k]);
  for (uint64_t idx = thread->begin, end = thread->end; idx < end;
      idx++, occ = ridx_state_inc(&is, &capt_ri[k]))
  {
    if (!ridx_state_legal(&is, pos.stm, occ))
      continue;
    pos.occ = occ;
    ridx_state_to_sq(&is, pos.sq, &capt_ri[k]);
    pos.sq[m] = pos.sq[n];
    uint8_t v = wdl_to_capt_val[probe_wdl(&pos, -2, 2) + 2];
    mark_king_uncaptures(stm, k, p, occ, &is, v);
    mark_uncaptures(stm, KNIGHT, k, p, occ, &is, v);
    mark_uncaptures(stm, BISHOP, k, p, occ, &is, v);
    mark_uncaptures(stm, ROOK  , k, p, occ, &is, v);
    mark_uncaptures(stm, QUEEN , k, p, occ, &is, v);
  }
}

static void calc_captures(void)
{
  for (int k = 0; k < ri.numsets; k++) {
    g_pos.stm = g_set_pt[k] >> 3;
    work_set = k;
    run_threaded(calc_capt_worker, &work_capt_dynamic[k]);
  }
}

static void calc_illegal_worker_tmpl(struct ThreadData *thread, const int pc)
{
  struct RIdxState is;
  int stm = g_pos.stm;
  int k = g_piece_set[stm][pc];
  uint8_t *restrict p = g_table[stm];

  Bitboard occ = ridx_state_init(&is, thread->begin, &capt_ri[k]);
  for (uint64_t idx = thread->begin, end = thread->end; idx < end;
      idx++, occ = ridx_state_inc(&is, &capt_ri[k]))
  {
    int ksq = is.sq[stm ^ 1];
    mark_uncapture_king(stm, pc, ksq, p, occ, &is);
  }
}

static void calc_illegal_knight_worker(struct ThreadData *thread)
{
  calc_illegal_worker_tmpl(thread, KNIGHT);
}

static void calc_illegal_bishop_worker(struct ThreadData *thread)
{
  calc_illegal_worker_tmpl(thread, BISHOP);
}

static void calc_illegal_rook_worker(struct ThreadData *thread)
{
  calc_illegal_worker_tmpl(thread, ROOK);
}

static void calc_illegal_queen_worker(struct ThreadData *thread)
{
  calc_illegal_worker_tmpl(thread, QUEEN);
}

static void calc_illegal(void)
{
  for (int stm = 0; stm < 2; stm++) {
    int k;
    g_pos.stm = stm;
    if ((k = g_piece_set[stm][KNIGHT]) >= 0)
      run_threaded(calc_illegal_knight_worker, &work_capt_dynamic[k]);
    if ((k = g_piece_set[stm][BISHOP]) >= 0)
      run_threaded(calc_illegal_bishop_worker, &work_capt_dynamic[k]);
    if ((k = g_piece_set[stm][ROOK]) >= 0)
      run_threaded(calc_illegal_rook_worker, &work_capt_dynamic[k]);
    if ((k = g_piece_set[stm][QUEEN]) >= 0)
      run_threaded(calc_illegal_queen_worker, &work_capt_dynamic[k]);
  }
}

// Mark potential mates.
static void calc_mates_worker(struct ThreadData *thread)
{
  uint8_t *restrict table_w = g_table[WHITE];
  uint8_t *restrict table_b = g_table[BLACK];

#ifdef __AVX512BW__

  const __m512i v_ff   = _mm512_set1_epi8((char)0xff);
  const __m512i v_zero = _mm512_setzero_si512();
  const __m512i v_one  = _mm512_set1_epi8(1);

  for (uint64_t idx = thread->begin, end = thread->end; idx < end; idx += 64) {
    __m512i w = _mm512_load_si512((__m512 *)(table_w + idx));
    __m512i b = _mm512_load_si512((__m512 *)(table_b + idx));

    __mmask64 w_is_ff = _mm512_cmpeq_epi8_mask(w, v_ff);
    __mmask64 b_is_ff = _mm512_cmpeq_epi8_mask(b, v_ff);
    __mmask64 w_is_0  = _mm512_cmpeq_epi8_mask(w, v_zero);
    __mmask64 b_is_0  = _mm512_cmpeq_epi8_mask(b, v_zero);

    __mmask64 set_w = b_is_ff & w_is_0;
    __mmask64 set_b = w_is_ff & b_is_0;

    w = _mm512_mask_mov_epi8(w, set_w, v_one);
    b = _mm512_mask_mov_epi8(b, set_b, v_one);

    _mm512_store_si512((__m512 *)(table_w + idx), w);
    _mm512_store_si512((__m512 *)(table_b + idx), b);
  }

#else

  for (uint64_t idx = thread->begin, end = thread->end; idx < end; idx++)
    if (table_w[idx] == 0xff) {
      if (table_b[idx] == 0) {
        table_b[idx] = 1;
      }
    } else if (table_b[idx] == 0xff) {
      if (table_w[idx] == 0) {
        table_w[idx] = 1;
      }
    }

#endif
}

static void calc_mates(void)
{
  run_threaded(calc_mates_worker, &work_g_static);
}

alignas(64) static uint8_t action[256];
alignas(64) static uint8_t mark_val[256];

static _Atomic bool not_finished;

static void iter(struct ThreadData *thread)
{
  struct RIdxState is;
  int stm = g_pos.stm;
  uint8_t *restrict table = g_table[stm];
  uint8_t *restrict table_opp = g_table[stm ^ 1];

  bool finished = true;

  uint64_t last = thread->begin;
  ridx_state_init(&is, last, &ri);
  for (uint64_t idx = last, end = thread->end; idx < end; idx++) {
    int t = table[idx];
    int w = action[t];
    if (!w) continue;
    Bitboard occ = ridx_state_add(&is, idx - last, &ri);
    last = idx;
    switch (w) {
    case 1:
      if (   !check_moves(stm, QUEEN , table_opp, occ, &is)
          || !check_moves(stm, ROOK  , table_opp, occ, &is)
          || !check_moves(stm, BISHOP, table_opp, occ, &is)
          || !check_moves(stm, KNIGHT, table_opp, occ, &is)
          || !check_king_moves(stm, table_opp, occ, &is))
      {
        table[idx] = 0;  // Reset to unresolved.
        break;
      }
      [[fallthrough]];
    case 2:
      uint8_t v = mark_val[t];
      int prev = stm ^ 1;
      mark_king_unmoves(prev, table_opp, occ, &is, v);
      mark_unmoves(prev, KNIGHT, table_opp, occ, &is, v);
      mark_unmoves(prev, BISHOP, table_opp, occ, &is, v);
      mark_unmoves(prev, ROOK  , table_opp, occ, &is, v);
      mark_unmoves(prev, QUEEN , table_opp, occ, &is, v);
      finished = false;
      break;
#if 0
    case 3:
      capt_bloss();
      break;
#endif
    default:
      unreachable();
    }
  }

  if (!finished)
    atomic_store_explicit(&not_finished, true, memory_order_relaxed);
}

// UNRESOLVED = 0
// MATE/L0 = 1
// L1 = 2
// L2 = 3
// ...
// L100 = 101
// CAPT_BLOSS/L101 = 102
// ...
// L124 = 125
// CAPT_DRAW = 126
// W124 = 127
// ...
// W101 = 150
// PAWN_CWIN = 151
// CAPT_CWIN = 152
// W100 = 153
// ...
// W2 = 251
// W1 = 252
// PAWN_WIN = 253
// CAPT_WIN = 254
// ILLEGAL = 255

uint16_t map[MAX_STATS];
// 0 -> ILLEGAL
// 1 -> CAPT_WIN
// 2 -> PAWN_WIN
// 3 -> W1
// ...
// 102 -> W100
// 103 -> CAPT_CWIN
// 104 -> PAWN_CWIN
// 105 -> W101
// ...
// MAX_STATS / 2 + 1 -> CAPT_DRAW
// MAX_STATS / 2 + 2 -> DRAW
// ...
// MAX_STATS - 3 -> L2
// MAX_STATS - 2 -> L1
// MAX_STATS - 1 -> L0

// ply = 0
// L0_w -> W1_b
// L1_w -> W2_b
//
// ply = 1
// L0_b -> W1_w
// L1_b -> W2_w
// W1_b -> L2_w
// W2_b -> L3_w
//
// ply = 2
// W1_w -> L2_b
// W2_w -> L3_b
// L2_w -> W3_b
// L3_w -> W4_b
//
// ...
//
// ply = 98
// W97_w -> L98_b
// W98_w -> L99_b
// L98_w -> W99_b
// L99_w -> W100_b
//
// ply = 99
// L98_b -> W99_w
// L99_b -> W100_w
// W99_b -> L100_w
// W100_b -> L101_w
//
// ply = 100
// W99_w -> L100_b
// W100_w -> L101_b
// L100_w -> W101_b
// L101_w -> W102_b
//
// WIN_CAPT to be added to W1_*
// CWIN_CAPT to be added to W101_*
// BLOSS_CAPT == L101_*
// LOSS_CAPT == L1_*

static int epoch = 0;

INLINE uint8_t win_to_byte(int n)
{
  if (epoch == 0)
    return 253 - n;
  else
    return 253 - n;
}

INLINE uint8_t loss_to_byte(int n)
{
  if (epoch == 0)
    return 1 + n;
  else
    return 1 + n;
}

static void iterate(void)
{
  memset(&mark_val, 0, sizeof mark_val);

  mark_val[255] = 1;
  mark_val[254] = mark_val[253] = loss_to_byte(2);
  mark_val[1] = win_to_byte(1);
  for (int n = 1; n <= DRAW_RULE; n++) {
    mark_val[win_to_byte(n)] = loss_to_byte(n + 1);
    mark_val[loss_to_byte(n)] = win_to_byte(n + 1);
  }

  int stm = WHITE;
  int n = 0;
  while (true) {
    // White to move, n is even
    // W(n-1)_w -> L(n  )_b
    // W(n  )_w -> L(n+1)_b
    // L(n  )_w -> W(n+1)_b
    // L(n+1)_w -> W(n+2)_b

    loss_byte =  n == 0 ? 255
               : win_to_byte(n);

    memset(&action, 0, sizeof action);
    if (n > 0) {
      action[win_to_byte(n - 1)] = 2;
      action[win_to_byte(n)] = 2;
      if (n == 2) {
        action[254] = action[253] = 2;
      }
    }
    action[loss_to_byte(n)] = 1;
    action[loss_to_byte(n + 1)] = 1;

    g_pos.stm = stm;
    atomic_store_explicit(&not_finished, false, memory_order_relaxed);
    printf("iteration %d\n", n);
    run_threaded(iter, &work_g_dynamic);
//    if (!atomic_load_explicit(&not_finished, memory_order_relaxed)) break;
    stm ^= 1;
    n++;

    // Black to move, n is odd
    // L(n-1)_w -> W(n  )_b
    // L(n  )_w -> W(n+1)_b
    // W(n  )_w -> L(n+1)_b
    // W(n+1)_w -> L(n+2)_b

    loss_byte =  n == 1 ? 255
               : win_to_byte(n - 1);

    memset(&action, 0, sizeof action);
    action[loss_to_byte(n - 1)] = 1;
    action[loss_to_byte(n)]     = 1;
    action[win_to_byte(n)]      = 2;
    action[win_to_byte(n + 1)]  = 2;
    if (n == 1)
      action[254] = action[253] = 2;

    g_pos.stm = stm;
//    atomic_store_explicit(&not_finished, false, memory_order_relaxed);
    printf("iteration %d\n", n);
    run_threaded(iter, &work_g_dynamic);
    if (!atomic_load_explicit(&not_finished, memory_order_relaxed)) break;
    stm ^= 1;
    n++;
  }
}

void generate(void)
{
  init_generation_work();

  printf("Initializing table.\n");
  clear_table();
  printf("Probing subtables.\n");
  calc_captures();
  printf("Calculating illegal positions.\n");
  calc_illegal();
  printf("Calculating cnadidate mate positions.\n");
  calc_mates();
  iterate();
}
