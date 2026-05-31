/*
  Copyright (c) 2026 Ronald de Man

  This file is distributed under the terms of the GNU GPL, version 2.
*/

#include <stdbit.h>
#include <stdint.h>

#include "defs.h"
#include "index.h"
#include "probe.h"
#include "types.h"

struct IdxInfo ii, capt_ii[MAX_SETS];
int pc_to_set[MAX_PIECES];
Bitboard Unrank2[62 * 61 / 2], Unrank3[62 * 61 * 60 / 6];

void init_unrank(void)
{
  int idx = 0;
  for (int s1 = 1; s1 < 62; s1++)
    for (int s0 = 0; s0 < s1; s0++)
      Unrank2[idx++] = bit(s0) | bit(s1);

  idx = 0;
  for (int s2 = 2; s2 < 62; s2++)
    for (int s1 = 1; s1 < s2; s1++)
      for (int s0 = 0; s0 < s1; s0++)
        Unrank3[idx++] = bit(s0) | bit(s1) | bit(s2);
}

uint64_t sq_to_idx(uint8_t *restrict sq)
{
  Bitboard occ = bit(sq[0]) | bit(sq[1]);
  if (has_pawns)
    occ |= bit(sq[2]);

  return sq_to_idx_helper(sq, 0, occ, &ii);
}

uint64_t capt_sq_to_idx(uint8_t *restrict sq, int k)
{
  Bitboard occ = bit(sq[0]) | bit(sq[1]);
  if (has_pawns)
    occ |= bit(sq[2]);

  return sq_to_idx_helper(sq, 0, occ, &capt_ii[k]);
}

static int find_partition(int len, uint8_t mult[]);

void calc_factors(struct IdxInfo *ii)
{
  for (int i = 0, n = has_pawns ? 61 : 62; i < ii->numsets; i++) {
    ii->factor[i] = Binomial[ii->mult[i]][n];
    ii->recip[i] = recip(ii->factor[i]);
    n -= ii->mult[i];
  }

  uint64_t f = 1;
  for (int i = ii->numsets - 1; i >= 0; i--)
    f *= ii->factor[i];
  ii->size = f;

  uint8_t mult[6];
  for (int i = 0; i < ii->numsets; i++)
    mult[i] = ii->mult[i];
  ii->part_id = find_partition(ii->numsets, mult);
}

void idx_state_init(struct IdxState *is, uint64_t idx, uint8_t *restrict sq,
    const struct IdxInfo *ii)
{
  for (int k = ii->numsets - 1; k >= 0; k--)
    idx = divmod_recip(idx, ii->factor[k], ii->recip[k], &is->sub[k]);
  is->n = 0;
  is->occ[0] = has_pawns ? bit(sq[0]) | bit(sq[1]) | bit(sq[2])
                         : bit(sq[0]) | bit(sq[1]);
}

// Perfect indexing code for positions with Ks on the diagonal starts here.

static constexpr Bitboard LOWER_DIAG_MASK = UINT64_C(0x0080c0e0f0f8fcfe);
static constexpr Bitboard MAIN_DIAG_MASK = UINT64_C(0x8040201008040201);

static uint8_t partition[30][7] = {
  { 0 },
  { 1 },
  { 2 }, { 1, 1},
  { 3 }, { 2, 1}, { 1, 1, 1 },
  { 4 }, { 3, 1}, { 2, 2 }, { 2, 1, 1 }, { 1, 1, 1, 1, 1 },
  { 5 }, { 4, 1}, { 3, 2 }, { 3, 1, 1 }, { 2, 2, 1 }, { 2, 1, 1, 1},
         { 1, 1, 1, 1, 1},
  { 6 }, { 5, 1}, { 4, 2 }, { 4, 1, 1 }, { 3, 3 }, { 3, 2, 1 }, { 3, 1, 1, 1},
         { 2, 2, 2 }, { 2, 2, 1, 1 }, { 2, 1, 1, 1, 1 }, { 1, 1, 1, 1, 1, 1}
};

static int8_t next_partition[30][8];
static uint8_t transition_id[30][8];

uint64_t reflection_size[30];

// Per transition, one case per number of 2-orbits filled.
struct TransitionCase {
  uint8_t d;
  uint8_t rem;
  uint64_t diag_tail;
  uint64_t diag_block;
  uint64_t broken_tail;
  uint64_t per_full_block;
  uint64_t block;
  uint64_t prefix;
};

// There are 201 valid cases.
static struct TransitionCase transition_cases[45][12][4];

INLINE Bitboard flip_main_diag(Bitboard b)
{
  Bitboard t;

  t  = (b ^ (b << 28)) & UINT64_C(0x0f0f0f0f00000000);
  b ^= t ^ (t >> 28);

  t  = (b ^ (b << 14)) & UINT64_C(0x3333000033330000);
  b ^= t ^ (t >> 14);

  t  = (b ^ (b << 7)) & UINT64_C(0x5500550055005500);
  b ^= t ^ (t >> 7);

  return b;
}

INLINE uint64_t binom(int n, int k)
{
  return (k < 0 || k > n) ? 0 : Binomial[k][n];
}

INLINE int popcnt32(uint32_t m)
{
  return stdc_count_ones(m);
}

// Fold (p,s) into the range 0...11.
// p = number of empty 2-orbits, s = number of empty 1-orbits.
INLINE int fold_ps(int p, int s)
{
  int used_p = 28 - p, used_s = 6 - s;
  return used_p * 7 - used_p * used_p + used_s;
}

static uint8_t unfold_ps[12][2];

static int find_partition(int len, uint8_t mult[])
{
  uint8_t m[6] = { 0 };

  for (int i = 0; i < len; i++)
    m[i] = mult[i];

  for (int i = 0; i < len; i++)
    for (int j = i + 1; j < len; j++)
      if (m[i] < m[j])
        Swap(m[i], m[j]);

  for (int i = 0; i < 30; i++)
    if (memcmp(m, partition[i], 6) == 0)
      return i;

  assume(0);
}

static uint64_t count_broken_residual_cases(int m, int p, int s, int d)
{
  int rem = m - 2 * d;
  int one_min = max(1, rem - s);
  int one_max = min(p - d, rem);
  uint64_t total = 0;
  for (int one = one_min; one <= one_max; one++) {
    int f = rem - one;
    total += (binom(p - d, one) * binom(s, f)) << (one - 1);
  }
  return total;
}

static uint64_t count_trivial_part(int part_id, int free)
{
  uint8_t *mult = partition[part_id];
  uint64_t r = 1;
  for (int i = 0; mult[i]; i++) {
    int m = mult[i];
    r *= binom(free, m);
    free -= m;
  }
  return r;
}

void init_perfect_ranker(void)
{
  int id = 0;
  for (int used_p = 0; 2 * used_p <= 5; used_p++)
    for (int used_s = 0; 2 * used_p + used_s <= 5; used_s++) {
      unfold_ps[id][0] = 28 - used_p;
      unfold_ps[id][1] = 6 - used_s;
      id++;
    }

  memset(next_partition, -1, sizeof next_partition);

  uint8_t mult[7];
  for (int p = 0; p < 30; p++) {
    int m = 0;
    for (int i = 0; partition[p][i]; i++) {
      if (partition[p][i] == m)
        continue;
      m = partition[p][i];
      int k = 0;
      for (int j = 0; partition[p][j]; j++)
        if (j != i)
          mult[k++] = partition[p][j];
      next_partition[p][m] = find_partition(k, mult);
    }
  }

  uint64_t count[30][4][7] = { 0 };

  for (int p = 25; p <= 28; p++)
    for (int s = 0; s <= 6; s++)
      count[0][p - 25][s] = 1;

  for (int id = 1; id < 30; id++) {
    uint8_t *part = partition[id];
    int m = part[0];
    int next_part = next_partition[id][m];
    for (int p = 25; p <= 28; p++)
      for (int s = 0; s <= 6; s++) {
        uint64_t total = 0;
        for (int d = 0; d <= min(p, m / 2); d++) {
          int rem = m - 2 * d;
          uint64_t per_full = 0;
          if (rem <= s)
            per_full += binom(s, rem) * count[next_part][p - d - 25][s - rem];
          uint64_t broken_cases = count_broken_residual_cases(m, p, s, d);
          if (broken_cases != 0)
            per_full += broken_cases * count_trivial_part(next_part,
                2 * p + s - m);
          total += binom(p, d) * per_full;
        }
        count[id][p - 25][s] = total;
      }
    reflection_size[id] = count[id][28 - 25][6];
  }

  id = 0;
  for (int part_id = 1; part_id < 30; part_id++) {
    for (int m = 1; m <= MAX_SETS; m++) {
      int next_part = next_partition[part_id][m];
      if (next_part < 0)
        continue;
      for (int ps = 0; ps < 12; ps++) {
        int p = unfold_ps[ps][0];
        int s = unfold_ps[ps][1];
        uint64_t prefix = 0;
        for (int d = 0; d <= min(p, m / 2); d++) {
          int rem = m - 2 * d;
          uint64_t diag_tail = 0;
          uint64_t diag_block = 0;
          if (rem <= s) {
            diag_tail = count[next_part][p - d - 25][s - rem];
            diag_block = binom(s, rem) * diag_tail;
          }
          uint64_t broken_cases = count_broken_residual_cases(m, p, s, d);
          uint64_t broken_tail =
            broken_cases == 0 ? 0 : count_trivial_part(next_part, 2 * p + s - m);
          uint64_t per_full = diag_block + broken_cases * broken_tail;
          uint64_t block = binom(p, d) * per_full;
          transition_cases[id][ps][d] = (struct TransitionCase) {
            .d = d,
            .rem = rem,
            .diag_tail = diag_tail,
            .diag_block = diag_block,
            .broken_tail = broken_tail,
            .per_full_block = per_full,
            .block = block,
            .prefix = prefix
          };
          prefix += block;
        }
      }
      transition_id[part_id][m] = id++;
    }
  }
}

static uint64_t rank_combination(Bitboard subset, Bitboard universe)
{
  uint64_t dense = _pext_u64(subset, universe);

  uint64_t r = 0;
  for (int j = 1; dense; j++)
    r += Binomial[j][pop_lsb(&dense)];
  return r;
}

static uint64_t rank_trivial_from(uint8_t *restrict sq, int k, Bitboard occ,
    const struct IdxInfo *ii)
{
  uint64_t idx = 0;
  for (; k < ii->numsets; k++) {
    size_t s;
    int i = ii->first[k];
    int m = ii->mult[k];
    if (m == 1) {
      s = rank_among_free(sq[i], occ);
      occ |= bit(sq[i]);
    } else if (m == 2) {
      Bitboard b = ~occ;
      Bitboard b1 = bit(sq[i]) | bit(sq[i + 1]);
      occ |= b1;
      b1 = _pext_u64(b1, b);
      s = pop_lsb(&b1);
      s += Binomial[2][lsb(b1)];
    } else {
      Bitboard b = ~occ, b1 = 0;
      for (int j = 0; j < m; j++)
        b1 |= bit(sq[i + j]);
      occ |= b1;
      b1 = _pext_u64(b1, b);
      s = 0;
      for (int j = 1; b1; j++)
        s += Binomial[j][pop_lsb(&b1)];
    }
    idx = idx * ii->factor[k] + s;
  }
  return idx;
}

INLINE uint64_t count_broken_residual_before(int rem, int p, int s, int one)
{
  uint64_t total = 0;
  int one_min = max(1, rem - s);
  for (int oo = one_min; oo < one; oo++) {
    int f = rem - oo;
    total += binom(p, oo) * binom(s, f) * (1ull << (oo - 1));
  }
  return total;
}

uint64_t rank_reflection(uint8_t *restrict sq, Bitboard occ,
    const struct IdxInfo *ii)
{
  int part_id = ii->part_id;
  Bitboard pair_mask = LOWER_DIAG_MASK;
  Bitboard diag_mask = MAIN_DIAG_MASK & ~occ;
  int p = 28, s = 6;

  uint64_t rank = 0;
  for (int k = 0; k < ii->numsets; k++) {
    int m = ii->mult[k];
    Bitboard bb = 0;
    for (int i = 0; i < m; i++)
      bb |= bit(sq[ii->first[k] + i]);
    occ |= bb;
    int tid = transition_id[part_id][m];

    Bitboard mirror = flip_main_diag(bb);
    Bitboard full_mask = bb & mirror & pair_mask;
    Bitboard one_mask = (bb ^ mirror) & pair_mask;
    int d = popcnt(full_mask);    // Number of 2-orbits fully filled.

    const struct TransitionCase *c = &transition_cases[tid][fold_ps(p, s)][d];

    rank += c->prefix;
    rank += rank_combination(full_mask, pair_mask) * c->per_full_block;
    pair_mask &= ~full_mask;
    p -= d;

    if (!one_mask) {
      rank += rank_combination(bb, diag_mask) * c->diag_tail;
      diag_mask &= ~bb;
      s = popcnt(diag_mask);
      part_id = next_partition[part_id][m];
      continue;
    }

    int one = popcnt(one_mask);     // Number of 2-orbits half filled.
    int f = popcnt(bb & diag_mask); // Number of 1-orbits filled.
    rank += c->diag_block;
    uint64_t r = count_broken_residual_before(c->rem, p, s, one);

    uint64_t rone = rank_combination(one_mask, pair_mask);
    r += (rone * binom(s, f) + rank_combination(bb, diag_mask)) << (one - 1);
    rank += r * c->broken_tail;

    // Canonical orientation: orient_mask <= bitwise complement within oo bits.
    uint32_t orient_mask = _pext_u64(bb, one_mask);
    uint32_t mask = (1u << one) - 1u;
    uint32_t comp = (~orient_mask) & mask;
    uint32_t canon = orient_mask < comp ? orient_mask : comp;
    // Among 2^oo orientations paired with complements, canonical ones are
    // exactly 0 .. 2^(oo-1)-1 after min(x,~x) for this ordering.
    assert(canon < (1u << (one - 1)));
    rank += canon * c->broken_tail;

    if (comp < orient_mask) {
      for (int i = 2; i < MAX_PIECES; i++)
        sq[i] = FlipDiag[sq[i]];
      bb = flip_main_diag(bb);
    }
    return rank + rank_trivial_from(sq, k + 1, occ | bb, ii);
  }
  return rank;
}

static Bitboard unrank_combination(uint64_t idx, Bitboard universe, int m)
{
  if (m == 0)
    return 0;

  Bitboard dense;
  if (m == 1)
    dense = bit(idx);
  else if (m == 2)
    dense = Unrank2[idx];
  else if (m == 3)
    dense = Unrank3[idx];
  else {
    dense = 0;
    int r = popcnt(universe) - 1;
    for (int i = m - 1; i > 0; i--) {
      while (idx < Binomial[i + 1][r])
        r--;
      idx -= Binomial[i + 1][r];
      dense |= bit(r);
      r--;
    }
    dense |= bit(idx);
  }
  return _pdep_u64(dense, universe);
}

static void unrank_trivial(uint8_t *restrict sq, uint64_t idx, int g,
    Bitboard occ, const struct IdxInfo *ii)
{
  uint32_t sub[MAX_SETS];
  for (int k = ii->numsets - 1; k >= g; k--)
    idx = divmod_recip(idx, ii->factor[k], ii->recip[k], &sub[k]);
  for (; g < ii->numsets; g++)
    occ = unrank_binomial(sub[g], ii->mult[g], sq + ii->first[g], occ);
}

void unrank_reflection(uint64_t idx, uint8_t *restrict sq, Bitboard occ,
    const struct IdxInfo *ii)
{
  int part_id = ii->part_id;
  Bitboard pair_mask = LOWER_DIAG_MASK;
  Bitboard diag_mask = MAIN_DIAG_MASK & ~occ;
  int p = 28, s = 6;

  for (int k = 0; k < ii->numsets; k++) {
    int m = ii->mult[k];
    int tid = transition_id[part_id][m];

    Bitboard bb = 0;
    const struct TransitionCase *c = &transition_cases[tid][fold_ps(p, s)][0];
    if (idx >= c->block) {
      idx -= c->block;
      for (;;) {
        c++;
        p--;
        if (idx < c->block) break;
        idx -= c->block;
      }

      uint64_t rfull = idx / c->per_full_block;
      idx %= c->per_full_block;
      Bitboard full_lower = unrank_combination(rfull, pair_mask, c->d);
      bb = full_lower | flip_main_diag(full_lower);
      pair_mask &= ~full_lower;
    }

    // Unbroken symmetry case.
    if (idx < c->diag_block) {
      uint64_t rsing = idx / c->diag_tail;
      idx %= c->diag_tail;
      bb |= unrank_combination(rsing, diag_mask, c->rem);
      s -= c->rem;
      diag_mask &= ~bb;
      occ |= bb;
      for (int i = ii->first[k]; bb; i++)
        sq[i] = pop_lsb(&bb);
      part_id = next_partition[part_id][m];
      continue;
    }

    idx -= c->diag_block;
    uint64_t q = idx / c->broken_tail;
    idx %= c->broken_tail;

    int one = max(1, c->rem - s);
    for (;; one++) {
      int f = c->rem - one;
      uint64_t one_cases = (binom(p, one) * binom(s, f)) << (one - 1);
      if (q < one_cases)
        break;
      q -= one_cases;
    }
    int f = c->rem - one;
    uint64_t orient = q & ((1ull << (one - 1)) - 1);
    q >>= one - 1;
    uint64_t rsing = q % binom(s, f);
    q /= binom(s, f);

    bb |= unrank_combination(rsing, diag_mask, f);

    Bitboard one_mask = unrank_combination(q, pair_mask, one);
    Bitboard orient_sparse = _pdep_u64(orient, one_mask);
    Bitboard lower_one = one_mask & orient_sparse;
    Bitboard upper_one = flip_main_diag(one_mask & ~orient_sparse);
    bb |= lower_one | upper_one;
    occ |= bb;
    for (int i = ii->first[k]; bb; i++)
      sq[i] = pop_lsb(&bb);
    unrank_trivial(sq, idx, k + 1, occ, ii);
    return;
  }
}
