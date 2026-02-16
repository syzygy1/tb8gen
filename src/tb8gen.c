/*
  Copyright (c) 2011-2017, 2025, 2026 Ronald de Man

  This file is distributed under the terms of the GNU GPL, version 2.
*/

#include <getopt.h>
#include <inttypes.h>
#include <stdatomic.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <x86intrin.h>

#include "defs.h"
#include "kslice.h"
#include "movegen.h"
#include "probe.h"
#include "threads.h"
#include "types.h"

#if 0
#include "compress.h"
#include "reduce.h"
#include "stats.h"
#include "permute.h"
#endif

#define TBPATH "RTBPATH"
#define STATSDIR "RTBSTATSDIR"

static constexpr int DRAW_RULE  = 2 * 50;
static constexpr int MAX_SETS = MAX_PIECES - 2;

int16_t KKMap[64][64];
uint8_t MirrorMask[64];
bool FlipTest[64][64];

size_t kslice_size;
Position g_pos;

bool g_only_generate, g_use_rans, symmetric;
char *g_tablename;

struct IdxInfo {
  int numpcs;    // total number of pieces, including kings
  int numsets;   // number of sets of like pieces, excluding kings.
  int norm[MAX_PIECES];
  uint64_t factor[MAX_PIECES];
  int first[MAX_SETS];          // index of first piece of each set
  int mult[MAX_SETS];           // number of like pieces in each set
  uint32_t subfactor[MAX_SETS]; // total number of placements for a set
};

static struct IdxInfo ii, capt_ii[MAX_SETS];

static uint64_t *work_g, *work_capt[MAX_SETS];

void init_tables(void)
{
  static constexpr Bitboard A1D1D4 = 0x080c0e0f;
  static constexpr Bitboard A1D4   = 0x08040201;
  static constexpr Bitboard LOWER  = 0x80c0e0f0f8fcfeff;

  for (int s = 0; s < 64; s++)
    MirrorMask[s] = ((s & 0x04) ? 0x07 : 0x00) | ((s & 0x20) ? 0x38 : 0x00);

  for (int i = 0; i < 64; i++)
    for (int j = 0; j < 64; j++) {
      int s1 = i ^ MirrorMask[i];
      int s2 = j ^ MirrorMask[i];
      if (!(bit(s1) & A1D1D4) || ((bit(s1) & A1D4) && !(bit(s2) & LOWER))) {
        FlipTest[i][j] = true;
        s1 = FlipDiag[s1];
        s2 = FlipDiag[s2];
      }
      KKMap[i][j] = KKIdx[Triangle[s1]][s2];
    }
}

// To be optimized at some point.
INLINE void sort_squares(int n, uint8_t *restrict p)
{
  for (int i = 0; i < n; i++)
    for (int j = i + 1; j < n; j++)
      if (p[i] > p[j])
        Swap(p[i], p[j]);
}

INLINE int rank_among_free(uint8_t sq, Bitboard occ)
{
  return sq - popcnt(occ & ((1ULL << sq) - 1));
}

INLINE Bitboard unrank_binomial(uint64_t idx, int n, uint8_t *restrict p,
    Bitboard occ)
{
  if (n == 0)
    return occ;

  Bitboard b = ~occ;
  for (int i = n - 1; i > 0; i--) {
    int r = i;
    while (idx >= Binomial[i + 1][r + 1])
      r++;
    idx -= Binomial[i + 1][r];
    Bitboard b1 = _pdep_u64(1ULL << r, b);
    p[i] = lsb(b1);
    occ |= b1;
  }
  Bitboard b1 = _pdep_u64(1ULL << idx, b);
  p[0] = lsb(b1);
  occ |= b1;

  return occ;
}

// We expect a normalized position.
INLINE uint64_t sq_to_idx(uint8_t *restrict sq)
{
  uint64_t idx = 0;
  Bitboard occ = bit(sq[0]) | bit(sq[1]);

  int n = ii.numpcs;
  for (int k = 2; k < n;) {
    int t = k + ii.norm[k];
    sort_squares(ii.norm[k], &sq[k]);
    size_t s = 0;
    Bitboard occ2 = occ;
    for (int i = k; i < t; i++) {
      int rank = rank_among_free(sq[i], occ);
      occ2 |= bit(sq[i]);
      s += Binomial[i - k + 1][rank];
    }
    idx += s * ii.factor[k];
    occ = occ2;
    k = t;
  }

  return idx;
}

INLINE uint64_t sq_to_idx_mirror(uint8_t *restrict sq)
{
  uint8_t p[MAX_PIECES];

  for (int i = 2; i < ii.numpcs; i++)
    p[i] = FlipDiag[p[i]];

  return sq_to_idx(p);
}

INLINE Bitboard idx_to_sq_unpack(uint32_t *sub, uint8_t *sq, struct IdxInfo *ii)
{
  Bitboard occ = bit(sq[0]) | bit(sq[1]);
  for (int i = 0; i < ii->numsets; i++)
    occ = unrank_binomial(sub[i], ii->mult[i], sq + ii->first[i], occ);
  return occ;
}

static void idx_to_sq_init(uint64_t idx, uint32_t *sub, struct IdxInfo *ii)
{
  for (int i = 2, j = 0; i < ii->numpcs; i += ii->norm[i], j++) {
    sub[j] = idx / ii->factor[i];
    idx -= sub[j] * ii->factor[i];
  }
}

static Bitboard idx_to_sq(uint32_t *sub, uint8_t *restrict sq)
{
  return idx_to_sq_unpack(sub, sq, &ii);
}

INLINE void idx_to_sq_inc(uint32_t *sub, struct IdxInfo *ii)
{
  for (int i = ii->numsets - 1; ++sub[i] >= ii->subfactor[i]; i--)
    sub[i] = 0;
}

static Bitboard capt_idx_to_sq(uint32_t *sub, uint8_t *restrict sq,
    const int k)
{
  return idx_to_sq_unpack(sub, sq, &capt_ii[k]);
}


static void calc_factors(struct IdxInfo *ii)
{
  for (int i = 0; i < 8; i++)
    ii->factor[i] = 0;

  for (int i = 0, n = 62; i < ii->numsets; i++) {
    ii->subfactor[i] = subfactor(ii->mult[i], n);
    n -= ii->mult[i];
  }

  uint64_t f = 1;
  for (int i = ii->numsets - 1; i >= 0; i--) {
    ii->factor[ii->first[i]] = f;
    f *= ii->subfactor[i];
  }
  ii->factor[0] = f;
  ii->subfactor[0]++;
}

INLINE void normalize(uint8_t *restrict sq, uint8_t *restrict sq2)
{
  for (int i = 0; i < MAX_PIECES; i++)
    sq2[i] = sq[i] ^ MirrorMask[sq[0]];

  if (FlipTest[sq[0]][sq[1]])
    for (int i = 0; i < MAX_PIECES; i++)
      sq2[i] = FlipDiag[sq2[i]];
}

INLINE void mirror_diagonal(uint8_t *restrict sq)
{
  for (int i = 0; i < MAX_PIECES; i++)
    sq[i] = FlipDiag[sq[i]];
}

INLINE void mark_king_unmoves(int stm, Bitboard occ, uint8_t *restrict sq,
    int s)
{
  uint8_t sq2[MAX_PIECES];
  Bitboard b = king_attacks(sq[stm]) & ~(occ | king_attacks(sq[stm ^ 1]));
  while (b) {
    sq[stm] = pop_lsb(&b);
    normalize(sq, sq2);
    int s2 = KKMap[sq2[0]][sq2[1]];
    uint8_t *p = kslice_get_address(s2);
    kslice_bit_set_atomic(p, sq_to_idx(sq2));
    if (s < 441 && s2 >= 441) {
      mirror_diagonal(sq2);
      kslice_bit_set_atomic(p, sq_to_idx(sq2));
    }
  }
}

INLINE void mark_unmoves(int k, uint8_t *restrict p, Bitboard occ,
    uint8_t *restrict sq)
{
  uint8_t sq2[MAX_PIECES];
  Bitboard b = non_king_piece_moves(g_pos.pt[k], sq[k], occ);
  while (b) {
    sq[k] = pop_lsb(&b);
    for (int i = 0; i < MAX_PIECES; i++)
      sq2[i] = sq[i];
    uint64_t idx = sq_to_idx(sq2);
    kslice_bit_set_atomic(p, idx);
  }
}

static int work_slice, work_set;
static int work_lower, work_upper;
static const bool *work_v;

static void uncapture_pieces_worker(struct ThreadData *thread)
{
  Position pos = g_pos;
  uint32_t sub[MAX_SETS];
  int stm = pos.stm;
  pos.stm ^= 1;
  int n = --pos.num;
  int k = work_set;
  int s = work_slice;
  int lower = work_lower;
  int upper = work_upper;
  const bool *v = work_v;

  int m = capt_ii[k].first[k] + capt_ii[k].mult[k];
  pos.pt[m] = pos.pt[n];

  uint8_t *restrict p = kslice_get_address(s);

  idx_to_sq_init(thread->begin, sub, &capt_ii[k]);
  for (uint64_t idx = thread->begin, end = thread->end; idx < end;
      idx++, idx_to_sq_inc(sub, &capt_ii[k]))
  {
    pos.occ = capt_idx_to_sq(sub, pos.sq, k);
    pos.sq[m] = pos.sq[n];
    if (opp_king_attacked(&pos) || !v[probe_wdl(&pos, lower, upper) + 2])
      continue;
    pos.sq[m] = pos.sq[stm];
    mark_king_unmoves(stm, pos.occ, pos.sq, s);
    pos.sq[stm] = pos.sq[m];
    // Uncapture by non-king pieces
    for (int i = 1; pos.pcs[stm][i] >= 0; i++) {
      int j = pos.pcs[stm][i];
      pos.sq[m] = pos.sq[j];
      mark_unmoves(j, p, pos.occ, pos.sq);
      pos.sq[j] = pos.sq[m];
    }
  }
}

// Mark stm positions with a capture that results in wdl == v.
static void uncapture_pieces(int stm, int s, int lower, int upper,
    const bool *v)
{
  work_slice = s;
  work_lower = lower;
  work_upper = upper;
  work_v = v;

  g_pos.stm = stm;
  g_pos.sq[0] = KKSquare[s][0];
  g_pos.sq[1] = KKSquare[s][1];

  for (int k = 0; k < ii.numsets; k++) {
    int m = capt_ii[k].first[k] + capt_ii[k].mult[k];
    if ((g_pos.pt[m] >> 3) != (stm ^ 1))
      continue;
    work_set = k;
    run_threaded(uncapture_pieces_worker, work_capt[k], 0);
  }
}

static void predecessors_worker(struct ThreadData *thread)
{
  uint32_t sub[MAX_SETS];
  Position pos = g_pos;
  int stm = pos.stm;
  int s = work_slice;

  uint8_t *restrict p = kslice_get_address(-1);
  uint8_t *restrict q = kslice_get_address(s);

  idx_to_sq_init(thread->begin, sub, &ii);
  for (int64_t idx = thread->begin, end = thread->end; idx < end;
      idx++, idx_to_sq_inc(sub, &ii))
  {
    if (!kslice_bit_test(p, idx))
      continue;
    Bitboard occ = idx_to_sq(sub, pos.sq);
    uint8_t tmp = pos.sq[stm];
    mark_king_unmoves(stm, occ, pos.sq, s);
    pos.sq[stm] = tmp;
    for (int i = 1; pos.pcs[stm][i] >= 0; i++) {
      int j = pos.pcs[stm][i];
      uint8_t tmp = pos.sq[j];
      mark_unmoves(j, q, occ, pos.sq);
      pos.sq[j] = tmp;
    }
  }
}

// Calculate stm predecessors of stm^1 positions in scratch.
static void predecessors(int stm, int s)
{
  work_slice = s;
  g_pos.stm = stm;
  g_pos.sq[0] = KKSquare[s][0];
  g_pos.sq[1] = KKSquare[s][1];

  run_threaded(predecessors_worker, work_g, 0);
}

INLINE bool check_king_moves(int stm, Bitboard occ, uint8_t *restrict sq)
{
  uint8_t sq2[MAX_PIECES];

  Bitboard b = king_attacks(sq[stm]) & ~(occ | king_attacks(sq[stm ^ 1]));
  while (b) {
    sq[stm] = pop_lsb(&b);
    normalize(sq, sq2);
    int s2 = KKMap[sq2[0]][sq2[1]];
    uint8_t *p = kslice_get_address(s2);
    if (kslice_bit_test(p, sq_to_idx(sq2)))
      return false;
  } 

  return true;
}

INLINE bool check_moves(int k, uint8_t *restrict p, Bitboard occ,
    uint8_t *restrict sq)
{
  uint8_t sq2[MAX_PIECES];
  Bitboard b = non_king_piece_moves(g_pos.pt[k], sq[k], occ);
  while (b) {
    sq[k] = pop_lsb(&b);
    for (int i = 0; i < MAX_PIECES; i++)
      sq2[i] = sq[i];
    uint64_t idx = sq_to_idx(sq2);
    if (kslice_bit_test(p, idx))
      return false;
  }

  return true;
}

static void check_successors_worker(struct ThreadData *thread)
{
  uint32_t sub[MAX_SETS];
  Position pos = g_pos;
  int stm = pos.stm;
  int s = work_slice;

  uint8_t *restrict p = kslice_get_address(-1);
  uint8_t *restrict q = kslice_get_address(s);

  idx_to_sq_init(thread->begin, sub, &ii);
  for (int64_t idx = thread->begin, end = thread->end; idx < end;
      idx++, idx_to_sq_inc(sub, &ii))
  {
    if (!kslice_bit_test(p, idx))
      continue;
    Bitboard occ = idx_to_sq(sub, pos.sq);
    for (int i = 1; pos.pcs[stm][i] >= 0; i++) {
      int j = pos.pcs[stm][i];
      uint8_t tmp = pos.sq[j];
      bool v = check_moves(j, q, occ, pos.sq);
      pos.sq[j] = tmp;
      if (!v) goto clear_bit;
    }
    uint8_t tmp = pos.sq[stm];
    bool v = check_king_moves(stm, occ, pos.sq);
    pos.sq[stm] = tmp;
    if (v) continue;
clear_bit:
    kslice_bit_flip(p, idx);
  }
}

// Verify stm positions in bmap as loss against stm^1 positions.
static void check_successors(int stm, int s)
{
  work_slice = s;
  g_pos.stm = stm;
  g_pos.sq[0] = KKSquare[s][0];
  g_pos.sq[1] = KKSquare[s][1];

  run_threaded(check_successors_worker, work_g, 0);
}

static void calc_mask_L0_worker(struct ThreadData *thread)
{
  uint32_t sub[MAX_SETS];
  Position pos = g_pos;
  pos.stm = BLACK;

  uint8_t *restrict mask_w = kslice_buf[0];
  uint8_t *restrict mask_b = kslice_buf[1];
  uint8_t *restrict L0_w = kslice_buf[2];
  uint8_t *restrict L0_b = kslice_buf[3];

  idx_to_sq_init(thread->begin, sub, &ii);

  for (uint64_t idx = thread->begin, end = thread->end; idx < end;
      idx++, idx_to_sq_inc(sub, &ii))
  {
    pos.occ = idx_to_sq(sub, pos.sq);
    bool chk_b = opp_king_attacked(&pos); // Black gives check?
    pos.stm = WHITE;
    bool chk_w = opp_king_attacked(&pos); // White gives check?
    if (chk_b) {
      // If black gives check, the btm position is illegal.
      kslice_bit_flip(mask_b, idx);
      // Test whether the wtm position is mate, i.e. loss in 0.
      if (!chk_w && !has_legal_moves(&pos) && !has_legal_caps(&pos))
        kslice_bit_set(L0_w, idx);
    }
    pos.stm = BLACK;
    if (chk_w) {
      kslice_bit_flip(mask_w, idx);
      if (!chk_b && !has_legal_moves(&pos) && !has_legal_caps(&pos))
        kslice_bit_set(L0_b, idx);
    }
  }
}

// The initials WINS table contains the ILLEGAL positions.
static void calc_mask_L0(void)
{
  uint64_t broken_w = 0, broken_b = 0, loss0_w = 0, loss0_b = 0;

  uint8_t *restrict mask_w = kslice_buf[0];
  uint8_t *restrict mask_b = kslice_buf[1];
  uint8_t *restrict L0_w = kslice_buf[2];
  uint8_t *restrict L0_b = kslice_buf[3];

  for (int s = 0; s < 462; s++) {
    memset(mask_w, 0xff, kslice_cache_lines << 6);
    memset(mask_b, 0xff, kslice_cache_lines << 6);
    memset(L0_w, 0, kslice_cache_lines << 6);
    memset(L0_b, 0, kslice_cache_lines << 6);

    g_pos.sq[0] = KKSquare[s][0];
    g_pos.sq[1] = KKSquare[s][1];

    run_threaded(calc_mask_L0_worker, work_g, 0);

    broken_w += kslice_size - kslice_count_addr(mask_w);
    broken_b += kslice_size - kslice_count_addr(mask_b);
    loss0_w += kslice_count_addr(L0_w);
    loss0_b += kslice_count_addr(L0_b);
    kslice_write_addr(mask_w, s, WHITE, "mask", 0);
    kslice_write_addr(mask_b, s, BLACK, "mask", 0);
    kslice_write_addr(L0_w, s, WHITE, "L", 0);
    kslice_write_addr(L0_b, s, BLACK, "L", 0);
  }

  printf("broken_w = %lu\n", broken_w);
  printf("broken_b = %lu\n", broken_b);
  printf("l0_w = %lu\n", loss0_w);
  printf("l0_b = %lu\n", loss0_b);
}

static void check_loss_in_1_worker(struct ThreadData *thread)
{
  uint32_t sub[MAX_SETS];
  Position pos = g_pos;

  uint8_t *p = kslice_get_address(work_slice);

  idx_to_sq_init(thread->begin, sub, &ii);
  for (uint64_t idx = thread->begin, end = thread->end; idx < end;
      idx++, idx_to_sq_inc(sub, &ii))
  {
    if (kslice_bit_test(p, idx)) {
      pos.occ = idx_to_sq(sub, pos.sq);
      if (has_legal_moves(&pos) || !has_legal_caps(&pos))
        kslice_bit_flip(p, idx);
    }
  } 
}

// From the positions with no non-losing captures, remove the positions with
// legal moves and the positions with no legal captures. What remains is L1.
static void check_loss_in_1(int stm, int s)
{
  work_slice = s;
  g_pos.sq[0] = KKSquare[s][0];
  g_pos.sq[1] = KKSquare[s][1];
  g_pos.stm = stm;

  run_threaded(check_loss_in_1_worker, work_g, 0);
}

static void calc_W1_L1(int stm, bool *more_w, bool *more_l)
{
  static const bool v1[5] = { true, false, false, false, false };
  static const bool v2[5] = { true, true, true, true, false };
  uint64_t w1 = 0, l1 = 0;
  uint64_t cw = 0;

  // Create the CAPT_WIN bitmap.
  for (int i = 0; i < 462; i++) {
    struct KSliceManager *mgr = kslice_get_manager(stm, i);
    for (int j = 0; mgr->in[j] >= 0; j++) {
      kslice_reserve(mgr->in[j]);
      kslice_clear(mgr->in[j]);
    }

    uncapture_pieces(stm, mgr->kslice, -2, -1, v1);

    for (int j = 0; mgr->out[j] >= 0; j++) {
      int s = mgr->out[j];
      cw += kslice_count(s);
      kslice_read(-1, s, stm, "mask", 0); // everything but ILLEGAL
      kslice_or_not(s, -1);
      kslice_write(s, s, stm, "capt_win", 0); // ILLEGAL + CAPT_WIN
      kslice_release(s);
    }
  }

  // Create the W1 bitmap.
  for (int i = 0; i < 462; i++) {
    struct KSliceManager *mgr = kslice_get_manager(stm, i);
    for (int j = 0; mgr->in[j] >= 0; j++) {
      kslice_reserve(mgr->in[j]);
      kslice_clear(mgr->in[j]);
    }

    kslice_read(-1, mgr->kslice, stm ^ 1, "L", 0);
    predecessors(stm, mgr->kslice);

    for (int j = 0; mgr->out[j] >= 0; j++) {
      int s = mgr->out[j];
      kslice_read(-1, s, stm, "capt_win", 0);
      kslice_or(s, -1);
      kslice_read(-1, s, stm, "mask", 0); // ILLEGAL positions
      kslice_and(s, -1);
      w1 += kslice_count(s);
      kslice_write(s, s, stm, "W", 1); // exactly W1
      kslice_and_not(-1, s);  // remove W1 from mask
      kslice_write(-1, s, stm, "mask", 0); // all but W1 and ILLEGAL
      kslice_release(s);
    }
  }

  // Create the L1 bitmap.
  for (int i = 0; i < 462; i++) {
    struct KSliceManager *mgr = kslice_get_manager(stm, i);
    for (int j = 0; mgr->in[j] >= 0; j++) {
      kslice_reserve(mgr->in[j]);
      kslice_clear(mgr->in[j]);
    }

    // Mark all positions with a non-losing capture.
    uncapture_pieces(stm, mgr->kslice, 1, 2, v2);

    for (int j = 0; mgr->out[j] >= 0; j++) {
      int s = mgr->out[j];
      kslice_write(s, s, stm, "capt_bloss", 0);
      kslice_read(-1, s, stm, "mask", 0);
      kslice_not_and(s, -1);
      check_loss_in_1(stm, s);
      l1 += kslice_count(s);
      kslice_write(s, s, stm, "L", 1);
      kslice_release(s);
    }
  }

  *more_w = l1 != 0;
  *more_l = w1 != 0;

  printf("l1_%c = %lu\n", "wb"[stm], l1);
  printf("w1_%c = %lu (%lu)\n", "wb"[stm], w1, cw);
}

static bool calc_L_n(int stm, int n)
{
  uint64_t cnt = 0;

  // Calculate potential losses in n = predecessors(W(n-1))
  for (int i = 0; i < 462; i++) {
    struct KSliceManager *mgr = kslice_get_manager(stm, i);
    for (int j = 0; mgr->in[j] >= 0; j++) {
      kslice_reserve(mgr->in[j]);
      kslice_clear(mgr->in[j]);
    }

    kslice_read(-1, mgr->kslice, stm ^ 1, "W", n - 1);
    predecessors(stm, mgr->kslice);

    for (int j = 0; mgr->out[j] >= 0; j++) {
      int s = mgr->out[j];
      // Currently, we remove illegal positions, known wins, and positions
      // with a drawing capture.
      // Known wins other than CAPT_WIN will anyway not verify as a
      // potential loss because they have a move into a lost position.
      // Illegal positions can be moved by checking legality.
      // CAPT_WIN and CAPT_DRAW(/CWIN/BLOSS) can be checked by probing.
      // So we could either not mask at all or mask with capt_bloss.
      //
      // If we mask out capt_bloss, then why not include illegal?
      // Reason: illegal might explode capt_bloss
      //
      // Idea: create bitmap of all positions with only losing captures?
      // -> if position has capture and is not in this bitmap, then not
      // losing.
      // Depending on material, this is more or less efficient.
      // If stm is badly losing, then probably no drawing caps, so
      // capt_bloss would be small.
      // if stm is not badly losing, the captures are probably draw, so
      // few positions with only losing captures.
      //
      // Another idea: create smaller (n-1)-piece bitmaps for positions after
      // capture. Then probe those bitmaps.
      // These bitmaps either per K-slice (1/462) or per wK/bK-slice (1/8).
      kslice_read(-1, s, stm, "mask", 0);
      kslice_and(s, -1);
      kslice_read(-1, s, stm, n <= DRAW_RULE ? "capt_bloss" : "capt_draw", 0);
      kslice_and_not(s, -1);
      kslice_write(s, s, stm, "PL", 0);
      kslice_release(s);
    }
  }

  // Verify potential losses.
  for (int i = 0; i < 462; i++) {
    struct KSliceManager *mgr = kslice_get_manager(stm, i);
    for (int j = 0; mgr->in[j] >= 0; j++) {
      int s = mgr->in[j];
      kslice_reserve(s);
      // FIXME: only read k-slices if PL(s) is not empty.
      // mask MUST include W_in_<=N, but we could work with deltas
      kslice_read(s, s, stm ^ 1, "mask", 0);
    }

    int s = mgr->kslice;
    kslice_read(-1, s, stm, "PL", 0);
    check_successors(stm, s);
    cnt += kslice_count(-1);
    kslice_write(-1, s, stm, "L", n);
    kslice_delete(s, stm, "PL", 0);

    for (int j = 0; mgr->out[j] >= 0; j++)
      kslice_release(mgr->out[j]);
  }

  printf("l%d_%c = %lu\n", n, "wb"[stm], cnt);
  return cnt != 0;
}

static bool calc_W_n(int stm, int n)
{
  uint64_t cnt = 0;

  // Calculate wins in n = predecessors(L(n-1))
  for (int i = 0; i < 462; i++) {
    struct KSliceManager *mgr = kslice_get_manager(stm, i);
    for (int j = 0; mgr->in[j] >= 0; j++) {
      kslice_reserve(mgr->in[j]);
      kslice_clear(mgr->in[j]);
    }

    kslice_read(-1, mgr->kslice, stm ^ 1, "L", n - 1);
    predecessors(stm, mgr->kslice);

    for (int j = 0; mgr->out[j] >= 0; j++) {
      int s = mgr->out[j];
      // We are currently removing illegal positions and known faster wins.
      // We can check easily for illegal positions.
      // But we MUST remove W_in_<=(N-1)
      // Again, we could work with deltas.
      kslice_read(-1, s, stm, "mask", 0);
      kslice_and(s, -1);
      cnt += kslice_count(s);
      kslice_write(s, s, stm, "W", n);
      kslice_and_not(-1, s);
      kslice_write(-1, s, stm, "mask", 0);
      kslice_release(s);
    }
  }

  printf("w%d_%c = %lu\n", n, "wb"[stm], cnt);
  return cnt != 0;
}

static bool calc_L101(int stm, bool more)
{
  static const bool v1[5] = { true, true, true, false, false };
  uint64_t cnt = 0;

  // Calculate positions with a capture to a draw or better.
  for (int i = 0; i < 462; i++) {
    struct KSliceManager *mgr = kslice_get_manager(stm, i);
    for (int j = 0; mgr->in[j] >= 0; j++) {
      kslice_reserve(mgr->in[j]);
      kslice_clear(mgr->in[j]);
    }

    uncapture_pieces(stm, mgr->kslice, 0, 1, v1);

    for (int j = 0; mgr->out[j] >= 0; j++) {
      int s = mgr->out[j];
      kslice_write(s, s, stm, "capt_draw", 0);
      if (!more) {
        kslice_read(-1, s, stm, "capt_bloss", 0);
        kslice_and_not(-1, s);
        kslice_read(s, s, stm, "mask", 0);
        kslice_and(s, -1);
        kslice_write(s, s, stm, "PL", 0);
      }
      kslice_release(s);
    }
  }

  // Positions with a capture into a blessed loss are in L101 if
  // they have no better capture and all moves are to positions that
  // are illegal or lose in at most 100 ply.

  if (more) {
    // Calculate potential losses in 101 = predecessors(W100)
    for (int i = 0; i < 462; i++) {
      struct KSliceManager *mgr = kslice_get_manager(stm, i);
      for (int j = 0; mgr->in[j] >= 0; j++) {
        kslice_reserve(mgr->in[j]);
        kslice_clear(mgr->in[j]);
      }

      kslice_read(-1, mgr->kslice, stm ^ 1, "W", 100);
      predecessors(stm, mgr->kslice);

      for (int j = 0; mgr->out[j] >= 0; j++) {
        int s = mgr->out[j];
        kslice_read(-1, s, stm, "capt_bloss", 0);
        kslice_or(s, -1);
        kslice_read(-1, s, stm, "mask", 0);
        kslice_and(s, -1);
        kslice_read(-1, s, stm, "capt_draw", 0);
        kslice_and_not(s, -1);
        kslice_write(s, s, stm, "PL", 0);
        kslice_release(s);
      }
    }
  }

  // Verify potential losses.
  for (int i = 0; i < 462; i++) {
    struct KSliceManager *mgr = kslice_get_manager(stm, i);
    for (int j = 0; mgr->in[j] >= 0; j++) {
      int s = mgr->in[j];
      kslice_reserve(s);
      kslice_read(s, s, stm ^ 1, "mask", 0);
    }

    int s = mgr->kslice;
    kslice_read(-1, s, stm, "PL", 0);
    check_successors(stm, s);
    cnt += kslice_count(-1);
    kslice_write(-1, s, stm, "L", 101);
    kslice_delete(s, stm, "PL", 0);

    for (int j = 0; mgr->out[j] >= 0; j++)
      kslice_release(mgr->out[j]);
  }

  printf("l%d_%c = %lu\n", 101, "wb"[stm], cnt);
  return cnt != 0;
}

static bool calc_W101(int stm, bool more)
{
  static const bool v1[5] = { false, true, false, false, false };
  uint64_t cnt = 0;

  // Calculate positions with a capture to a cursed win.
  for (int i = 0; i < 462; i++) {
    struct KSliceManager *mgr = kslice_get_manager(stm, i);
    for (int j = 0; mgr->in[j] >= 0; j++) {
      kslice_reserve(mgr->in[j]);
      kslice_clear(mgr->in[j]);
    }

    uncapture_pieces(stm, mgr->kslice, -2, 0, v1);

    for (int j = 0; mgr->out[j] >= 0; j++) {
      int s = mgr->out[j];
      kslice_write(s, s, stm, "capt_cwin", 0);
      if (!more) {
        kslice_read(-1, s, stm, "mask", 0);
        kslice_and(s, -1);
        cnt += kslice_count(s);
        kslice_write(s, s, stm, "W", 101);
        kslice_and_not(-1, s);
        kslice_write(-1, s, stm, "mask", 0);
      }
      kslice_release(s);
    }
  }

  if (more) {
    // Calculate predecessors of L100.
    for (int i = 0; i < 462; i++) {
      struct KSliceManager *mgr = kslice_get_manager(stm, i);
      for (int j = 0; mgr->in[j] >= 0; j++) {
        kslice_reserve(mgr->in[j]);
        kslice_clear(mgr->in[j]);
      }

      kslice_read(-1, mgr->kslice, stm ^ 1, "L", 100);
      predecessors(stm, mgr->kslice);

      for (int j = 0; mgr->out[j] >= 0; j++) {
        int s = mgr->out[j];
        kslice_read(-1, s, stm, "capt_cwin", 0);
        kslice_or(s, -1);                   // add CAPT_CWIN positions
        kslice_read(-1, s, stm, "mask", 0);
        kslice_and(s, -1);
        cnt += kslice_count(s);
        kslice_write(s, s, stm, "W", 101);
        kslice_and_not(-1, s);
        kslice_write(-1, s, stm, "mask", 0);
        kslice_release(s);
      }
    }
  }

  printf("w%d_%c = %lu\n", 101, "wb"[stm], cnt);
  return cnt != 0;
}

void merge(int stm, int s)
{
}

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
  ii.numpcs = numpcs;
  k = 0;
  for (int i = 2; i < numpcs; i += ii.norm[i]) {
    for (int j = i; j < numpcs && pt[j] == pt[i]; j++)
      ii.norm[i]++;
    ii.first[k] = i;
    ii.mult[k] = ii.norm[i];
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
  }

  kslice_size = ii.factor[0];
  kslice_setup();

  // Align work units on cache lines of 64 x 8 = 512 positions.
  work_g = create_work(g_total_work, kslice_size, 0x1ff);
  for (int i = 0; i < ii.numsets; i++)
    work_capt[i] = create_work(g_total_work, capt_ii[i].factor[0], 0x1ff);

  calc_mask_L0();

  bool more_ww = true, more_wb = true, more_lw = true, more_lb = true;
  bool more_wb_next, more_ww_next;

  calc_W1_L1(WHITE, &more_wb, &more_lb);
  calc_W1_L1(BLACK, &more_ww, &more_lw);

  for (int n = 2; n <= DRAW_RULE; n++) {
    more_wb_next = more_lw && calc_L_n(WHITE, n);
    more_ww_next = more_lb && calc_L_n(BLACK, n);

    more_lb = more_ww && calc_W_n(WHITE, n);
    more_lw = more_wb && calc_W_n(BLACK, n);

    more_wb = more_wb_next;
    more_ww = more_ww_next;
  }

  more_wb_next = calc_L101(WHITE, more_lw);
  more_ww_next = calc_L101(BLACK, more_lb);

  more_lb = calc_W101(WHITE, more_ww);
  more_lw = calc_W101(BLACK, more_wb);

  more_wb = more_wb_next;
  more_ww = more_ww_next;

  for (int n = DRAW_RULE + 2; more_ww || more_wb || more_lw || more_lb; n++) {
    more_wb_next = more_lw && calc_L_n(WHITE, n);
    more_ww_next = more_lb && calc_L_n(BLACK, n);

    more_lb = more_ww && calc_W_n(WHITE, n);
    more_lw = more_wb && calc_W_n(BLACK, n);

    more_wb = more_wb_next;
    more_ww = more_ww_next;
  }

  kslice_cleanup();

  // Not yet implemented:
  for (int s = 0; s < 462; s++) {
    merge(WHITE, s);
    merge(BLACK, s);
  }
}
