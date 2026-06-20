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

struct RankInfo ri, capt_ri[MAX_SETS];
int pc_to_set[MAX_PIECES];
Bitboard Unrank2[62 * 61 / 2], Unrank3[62 * 61 * 60 / 6];
uint32_t Binomial[8][64];
uint64_t MirrorMask[64];
bool FlipTest[64][64];
uint8_t KK_transform[64][64];

const int16_t KKIdx[10][64] = {
  {  -1,  -1,  -1,   0,   1,   2,   3,   4,
     -1,  -1,  -1,   5,   6,   7,   8,   9,
     10,  11,  12,  13,  14,  15,  16,  17,
     18,  19,  20,  21,  22,  23,  24,  25,
     26,  27,  28,  29,  30,  31,  32,  33,
     34,  35,  36,  37,  38,  39,  40,  41,
     42,  43,  44,  45,  46,  47,  48,  49,
     50,  51,  52,  53,  54,  55,  56,  57  },
  {  58,  -1,  -1,  -1,  59,  60,  61,  62,
     63,  -1,  -1,  -1,  64,  65,  66,  67,
     68,  69,  70,  71,  72,  73,  74,  75,
     76,  77,  78,  79,  80,  81,  82,  83,
     84,  85,  86,  87,  88,  89,  90,  91,
     92,  93,  94,  95,  96,  97,  98,  99,
    100, 101, 102, 103, 104, 105, 106, 107,
    108, 109, 110, 111, 112, 113, 114, 115  },
  { 116, 117,  -1,  -1,  -1, 118, 119, 120,
    121, 122,  -1,  -1,  -1, 123, 124, 125,
    126, 127, 128, 129, 130, 131, 132, 133,
    134, 135, 136, 137, 138, 139, 140, 141,
    142, 143, 144, 145, 146, 147, 148, 149,
    150, 151, 152, 153, 154, 155, 156, 157,
    158, 159, 160, 161, 162, 163, 164, 165,
    166, 167, 168, 169, 170, 171, 172, 173  },
  { 174,  -1,  -1,  -1, 175, 176, 177, 178,
    179,  -1,  -1,  -1, 180, 181, 182, 183,
    184,  -1,  -1,  -1, 185, 186, 187, 188,
    189, 190, 191, 192, 193, 194, 195, 196,
    197, 198, 199, 200, 201, 202, 203, 204,
    205, 206, 207, 208, 209, 210, 211, 212,
    213, 214, 215, 216, 217, 218, 219, 220,
    221, 222, 223, 224, 225, 226, 227, 228  },
  { 229, 230,  -1,  -1,  -1, 231, 232, 233,
    234, 235,  -1,  -1,  -1, 236, 237, 238,
    239, 240,  -1,  -1,  -1, 241, 242, 243,
    244, 245, 246, 247, 248, 249, 250, 251,
    252, 253, 254, 255, 256, 257, 258, 259,
    260, 261, 262, 263, 264, 265, 266, 267,
    268, 269, 270, 271, 272, 273, 274, 275,
    276, 277, 278, 279, 280, 281, 282, 283  },
  { 284, 285, 286, 287, 288, 289, 290, 291,
    292, 293,  -1,  -1,  -1, 294, 295, 296,
    297, 298,  -1,  -1,  -1, 299, 300, 301,
    302, 303,  -1,  -1,  -1, 304, 305, 306,
    307, 308, 309, 310, 311, 312, 313, 314,
    315, 316, 317, 318, 319, 320, 321, 322,
    323, 324, 325, 326, 327, 328, 329, 330,
    331, 332, 333, 334, 335, 336, 337, 338  },
  {  -1,  -1, 339, 340, 341, 342, 343, 344,
     -1,  -1, 345, 346, 347, 348, 349, 350,
     -1,  -1, 441, 351, 352, 353, 354, 355,
     -1,  -1,  -1, 442, 356, 357, 358, 359,
     -1,  -1,  -1,  -1, 443, 360, 361, 362,
     -1,  -1,  -1,  -1,  -1, 444, 363, 364,
     -1,  -1,  -1,  -1,  -1,  -1, 445, 365,
     -1,  -1,  -1,  -1,  -1,  -1,  -1, 446  },
  {  -1,  -1,  -1, 366, 367, 368, 369, 370,
     -1,  -1,  -1, 371, 372, 373, 374, 375,
     -1,  -1,  -1, 376, 377, 378, 379, 380,
     -1,  -1,  -1, 447 ,381, 382, 383, 384,
     -1,  -1,  -1,  -1 ,448, 385, 386, 387,
     -1,  -1,  -1,  -1 , -1, 449, 388, 389,
     -1,  -1,  -1,  -1 , -1,  -1, 450, 390,
     -1,  -1,  -1,  -1 , -1,  -1,  -1, 451  },
  { 452, 391, 392, 393 ,394, 395, 396, 397,
     -1,  -1,  -1,  -1 ,398, 399, 400, 401,
     -1,  -1,  -1,  -1 ,402, 403, 404, 405,
     -1,  -1,  -1,  -1 ,406, 407, 408, 409,
     -1,  -1,  -1,  -1 ,453, 410, 411, 412,
     -1,  -1,  -1,  -1 , -1, 454, 413, 414,
     -1,  -1,  -1,  -1 , -1,  -1, 455, 415,
     -1,  -1,  -1,  -1 , -1,  -1,  -1, 456  },
  { 457, 416, 417, 418 ,419, 420, 421, 422,
     -1, 458, 423, 424 ,425, 426, 427, 428,
     -1,  -1,  -1,  -1 , -1, 429, 430, 431,
     -1,  -1,  -1,  -1 , -1, 432, 433, 434,
     -1,  -1,  -1,  -1 , -1, 435, 436, 437,
     -1,  -1,  -1,  -1 , -1, 459, 438, 439,
     -1,  -1,  -1,  -1 , -1,  -1, 460, 440,
     -1,  -1,  -1,  -1 , -1,  -1,  -1, 461  }
};

uint8_t KKSquare[462][2];
int16_t KKMap[64][64];

// Regular RankInfo structs used for generation.

struct RankInfo rank_info_61[32];
struct RankInfo rank_info_62[64];
struct RankInfo rank_info_63[64];

static void unrank_mult(int idx, uint8_t mult[MAX_SETS])
{
  int k = 0;
  uint32_t b = idx;
  while (b) {
    int s = stdc_first_trailing_one(b);
    b >>= s;
    mult[k++] = s;
  }
  while (k < MAX_SETS)
    mult[k++] = 0;
}

int rank_mult(uint8_t mult[MAX_SETS])
{
  int idx = 0, s = 0;
  for (int k = 0; k < MAX_SETS && mult[k]; k++) {
    s += mult[k];
    idx |= bit(s - 1);
  }
  return idx;
}

static_assert(MAX_SETS == 6);

static uint8_t partition[30][7] = {
  { 0 },
  { 1 },
  { 2 }, { 1, 1},
  { 3 }, { 2, 1}, { 1, 1, 1 },
  { 4 }, { 3, 1}, { 2, 2 }, { 2, 1, 1 }, { 1, 1, 1, 1 },
  { 5 }, { 4, 1}, { 3, 2 }, { 3, 1, 1 }, { 2, 2, 1 }, { 2, 1, 1, 1},
         { 1, 1, 1, 1, 1},
  { 6 }, { 5, 1}, { 4, 2 }, { 4, 1, 1 }, { 3, 3 }, { 3, 2, 1 }, { 3, 1, 1, 1},
         { 2, 2, 2 }, { 2, 2, 1, 1 }, { 2, 1, 1, 1, 1 }, { 1, 1, 1, 1, 1, 1}
};

static int find_partition(int len, uint8_t mult[])
{
  uint8_t m[MAX_SETS] = { 0 };

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

static uint8_t unfold_ps[12][2];

enum {
  ID, FLIP_V, FLIP_H, ROT180,
  FLIP_MAIN, ROT270, ROT90, FLIP_ANTI
};

// a1d4 -> ID        FLIP_MAIN
// h1e4 -> FLIP_H    ROT90
// h8e5 -> ROT180    FLIP_ANTI
// a8d5 -> FLIP_V    ROT270

static int transform_sq(int s, int t)
{
  int f = s & 7, r = s >> 3;

  switch (t) {
  case ID:        return  r      * 8 + f;
  case FLIP_V:    return (7 - r) * 8 + f;
  case FLIP_H:    return  r      * 8 + (7 - f);
  case ROT180:    return (7 - r) * 8 + (7 - f);
  case FLIP_MAIN: return  f      * 8 + r;
  case ROT90:     return  f      * 8 + (7 - r);
  case ROT270:    return (7 - f) * 8 + r;
  case FLIP_ANTI: return (7 - f) * 8 + (7 - r);
  }
  __builtin_unreachable();
}

void init_ranking(void)
{
  for (int i = 0; i < 10; i++)
    for (int j = 0; j < 64; j++)
      if (KKIdx[i][j] >= 0) {
        KKSquare[KKIdx[i][j]][0] = InvTriangle[i];
        KKSquare[KKIdx[i][j]][1] = j;
      }

  static constexpr Bitboard A1D1D4 = 0x080c0e0full;
  static constexpr Bitboard A1D4   = 0x08040201ull;
  static constexpr Bitboard LOWER  = 0x80c0e0f0f8fcfeffull;

  for (int s = 0; s < 64; s++)
    MirrorMask[s] =  ((s & 0x04) ? 0x0707070707070707ull : 0)
                   | ((s & 0x20) ? 0x3838383838383838ull : 0);

  for (int wk = 0; wk < 64; wk++)
    for (int bk = 0; bk < 64; bk++) {
      int s1 = wk ^ (MirrorMask[wk] & 0xff);
      int s2 = bk ^ (MirrorMask[wk] & 0xff);
      if (!(bit(s1) & A1D1D4) || ((bit(s1) & A1D4) && !(bit(s2) & LOWER))) {
        FlipTest[wk][bk] = true;
        s1 = ((s1 & 7) << 3) | (s1 >> 3);
        s2 = ((s2 & 7) << 3) | (s2 >> 3);
      }
      KKMap[wk][bk] = KKIdx[Triangle[s1]][s2];
    }

  for (int wk = 0; wk < 64; wk++) {
    int s, t;
    for (t = 0; t < 8; t++) {
      s = transform_sq(wk, t);
      if (bit(s) & A1D1D4)
        break;
    }
    assert(t < 8);
    if (bit(s) & A1D4) {
      assert(t < 4);
      for (int bk = 0; bk < 64; bk++) {
        int s1 = transform_sq(bk, t);
        if (bit(s1) & LOWER)
          KK_transform[wk][bk] = t;
        else
          KK_transform[wk][bk] = t + 4;
      }
    } else {
      for (int bk = 0; bk < 64; bk++)
        KK_transform[wk][bk] = t;
    }
  }

  // Binomial[k][n] = Bin(n, k)
  for (int j = 0; j < 64; j++)
    Binomial[0][j] = 1;
  for (int i = 1; i < 8; i++)
    for (int j = 1; j < 64; j++)
      Binomial[i][j] = Binomial[i - 1][j - 1] + Binomial[i][j - 1];

  int idx = 0;
  for (int s1 = 1; s1 < 62; s1++)
    for (int s0 = 0; s0 < s1; s0++)
      Unrank2[idx++] = bit(s0) | bit(s1);

  idx = 0;
  for (int s2 = 2; s2 < 62; s2++)
    for (int s1 = 1; s1 < s2; s1++)
      for (int s0 = 0; s0 < s1; s0++)
        Unrank3[idx++] = bit(s0) | bit(s1) | bit(s2);

  // Initialize tables for perfect ranking when kings are on the diagonal.

  int id = 0;
  for (int used_p = 0; 2 * used_p <= 5; used_p++)
    for (int used_s = 0; 2 * used_p + used_s <= 5; used_s++) {
      unfold_ps[id][0] = 28 - used_p;
      unfold_ps[id][1] = 6 - used_s;
      id++;
    }

  int8_t next_partition[30][7];
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
        for (int d = 0; d <= min(p - 25, m / 2); d++) {
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
  }

  uint8_t transition_id[30][8] = { 0 };

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
        for (int d = 0; d <= min(p - 25, m / 2); d++) {
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

  // Initialize RankInfo[] structs.

  for (int id = 0; id < 64; id++) {
    struct RankInfo *ri = &rank_info_62[id];
    uint8_t *m = ri->mult;
    unrank_mult(id, m);
    int n = 0;
    while (n < MAX_PIECES && m[n])
      n++;
    ri->numsets = n;
    int part_id = find_partition(n, m);
    for (int k = 0, p = part_id; k < n; k++) {
      ri->transition_id[k] = transition_id[p][m[k]];
      p = next_partition[p][m[k]];
    }
    calc_factors(ri, 62);
    ri->sizes[1] = count[part_id][28 - 25][6];
    int i = 2;
    for (int k = 0; k < ri->numsets; k++) {
      ri->first[k] = i;
      ri->last[k] = i + m[k] - 1;
      i += m[k];
    }
  }

  for (int id = 0; id < 32; id++) {
    struct RankInfo *ri = &rank_info_61[id];
    uint8_t *m = ri->mult;
    unrank_mult(id, m);
    int n = 0;
    while (n < MAX_PIECES && m[n])
      n++;
    ri->numsets = n;
    calc_factors(ri, 61);
    int i = 3;
    for (int k = 0; k < ri->numsets; k++) {
      ri->first[k] = i;
      ri->last[k] = i + m[k] - 1;
      i += m[k];
    }
  }

  for (int id = 0; id < 64; id++) {
    struct RankInfo *ri = &rank_info_63[id];
    uint8_t *m = ri->mult;
    unrank_mult(id, m);
    int n = 0;
    while (n < MAX_PIECES && m[n])
      n++;
    ri->numsets = n;
    calc_factors(ri, 63);
  }
}

void transform_set_bb(int t, Bitboard *restrict set_bb,
    Bitboard *restrict set_bb2)
{
  __m512i x = _mm512_load_si512((__m512i *)set_bb);
  switch (t) {
  case FLIP_V:
    x = flip_vertical_8xbb(x);
    break;
  case FLIP_H:
    x = flip_horizontal_8xbb(x);
    break;
  case ROT270:
    x = rotate270_8xbb(x);
    break;
  case FLIP_MAIN:
    x = flip_main_8xbb(x);
    break;
  case FLIP_ANTI:
    x = flip_anti_8xbb(x);
    break;
  case ROT180:
    x = rotate180_8xbb(x);
    break;
  case ROT90:
    x = rotate90_8xbb(x);
    break;
  default:
    unreachable();
  }
  _mm512_store_si512((__m512i *)set_bb2, x);
}

INLINE uint64_t rank_bb_from(Bitboard *restrict set_bb, int k, Bitboard occ,
    const struct RankInfo *ri)
{
  uint64_t idx = 0;
  for (; k < ri->numsets; k++) {
    Bitboard b = _pext_u64(set_bb[k + 1], ~occ);
    occ |= set_bb[k + 1];
    size_t s = 0;
    for (int j = 1; b; j++)
      s += Binomial[j][pop_lsb(&b)];
    idx = idx * ri->factor[k] + s;
  }
  return idx;
}

// Directly operating on struct IdxState2 might be even better:
// Bitboard b = _pext_u64(is->set_bb[k + 1], ~is->set_bb[k]);
// No need to keep track of occ.
uint64_t rank_bb(Bitboard *restrict set_bb, const struct RankInfo *ri)
{
  return rank_bb_from(set_bb, 0, set_bb[0], ri);
}

uint64_t sq_to_idx(uint8_t *restrict sq)
{
  Bitboard occ = bit(sq[0]) | bit(sq[1]);
  if (has_pawns)
    occ |= bit(sq[2]);

  return sq_to_idx_helper(sq, 0, occ, &ri);
}

uint64_t capt_sq_to_idx(uint8_t *restrict sq, int k)
{
  Bitboard occ = bit(sq[0]) | bit(sq[1]);
  if (has_pawns)
    occ |= bit(sq[2]);

  return sq_to_idx_helper(sq, 0, occ, &capt_ri[k]);
}

void calc_factors(struct RankInfo *ri, int n)
{
  uint64_t f = 1;
  for (int i = 0; i < ri->numsets; i++) {
    ri->factor[i] = Binomial[ri->mult[i]][n];
    if (ri->factor[i] != 1)
      ri->recip[i] = recip(ri->factor[i]);
    f *= ri->factor[i];
    n -= ri->mult[i];
  }
  ri->sizes[0] = f;
}

void idx_state_init2(struct IdxState2 *is, uint64_t idx, uint8_t *restrict sq,
    const struct RankInfo *ri)
{
  for (int k = ri->numsets - 1; k >= 0; k--)
    idx = divmod_recip(idx, ri->factor[k], ri->recip[k], &is->sub[k]);
  is->n = 0;
  is->sq[0] = sq[0];
  is->sq[1] = sq[1];
  if (has_pawns)
    is->sq[2] = sq[2];
  is->occ[0] = has_pawns ? bit(sq[0]) | bit(sq[1]) | bit(sq[2])
                         : bit(sq[0]) | bit(sq[1]);
}

void idx_state_init(struct IdxState *is, uint64_t idx, uint8_t *restrict sq,
    const struct RankInfo *ri)
{
  for (int k = ri->numsets - 1; k >= 0; k--)
    idx = divmod_recip(idx, ri->factor[k], ri->recip[k], &is->sub[k]);
  is->n = 0;
  is->occ[0] = has_pawns ? bit(sq[0]) | bit(sq[1]) | bit(sq[2])
                         : bit(sq[0]) | bit(sq[1]);
}

// Perfect indexing code for positions with Ks on the diagonal starts here.

static constexpr Bitboard LOWER_DIAG_MASK = UINT64_C(0x0080c0e0f0f8fcfe);
static constexpr Bitboard MAIN_DIAG_MASK = UINT64_C(0x8040201008040201);

static uint64_t rank_combination(Bitboard subset, Bitboard universe)
{
  uint64_t dense = _pext_u64(subset, universe);

  uint64_t r = 0;
  for (int j = 1; dense; j++)
    r += Binomial[j][pop_lsb(&dense)];
  return r;
}

uint64_t rank_trivial_from(uint8_t *restrict sq, int k, Bitboard occ,
    const uint8_t *restrict first, const struct RankInfo *ri)
{
  uint64_t idx = 0;
  for (; k < ri->numsets; k++) {
    size_t s;
    int i = first[k];
    int m = ri->mult[k];
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
    idx = idx * ri->factor[k] + s;
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
    const uint8_t *restrict first, const struct RankInfo *ri)
{
  Bitboard pair_mask = LOWER_DIAG_MASK;
  Bitboard diag_mask = MAIN_DIAG_MASK & ~occ;
  int p = 28, s = 6;

  uint64_t rank = 0;
  for (int k = 0; k < ri->numsets; k++) {
    int m = ri->mult[k];
    Bitboard bb = 0;
    for (int i = 0; i < m; i++)
      bb |= bit(sq[first[k] + i]);
    occ |= bb;
    int tid = ri->transition_id[k];

    Bitboard mirror = flip_main_diag(bb);
    Bitboard full_mask = bb & mirror  & pair_mask;
    Bitboard one_mask = (bb ^ mirror) & pair_mask;
    int d = popcnt(full_mask);    // Number of 2-orbits fully filled.

    const struct TransitionCase *c = &transition_cases[tid][fold_ps(p, s)][d];

    rank += c->prefix;
    rank += rank_combination(full_mask, pair_mask) * c->per_full_block;
    pair_mask &= ~full_mask;
    p -= d;

    if (one_mask) {
      // The current set breaks symmetry.
      int one = popcnt(one_mask);     // Number of 2-orbits half filled.
      int f = popcnt(bb & diag_mask); // Number of 1-orbits filled.
      rank += c->diag_block;
      uint64_t r = count_broken_residual_before(c->rem, p, s, one);

      uint64_t rone = rank_combination(one_mask, pair_mask);
      r += (rone * binom(s, f) + rank_combination(bb, diag_mask)) << (one - 1);
      rank += r * c->broken_tail;

      // Canonical orientation: orient_mask <= bitwise complement within oo bits.
      uint32_t orient_mask = _pext_u64(bb, one_mask);
      uint32_t comp = (~orient_mask) & ((1u << one) - 1u);
      uint32_t canon = orient_mask < comp ? orient_mask : comp;
      // Among 2^oo orientations paired with complements, canonical ones are
      // exactly 0 .. 2^(oo-1)-1 after min(x,~x) for this ordering.
      assert(canon < (1u << (one - 1)));
      rank += canon * c->broken_tail;

      if (comp < orient_mask) {
        mirror_diagonal(sq);
        occ = flip_main_diag(occ);
      }
      return rank + rank_trivial_from(sq, k + 1, occ, first, ri);
    }

    rank += rank_combination(bb, diag_mask) * c->diag_tail;
    diag_mask &= ~bb;
    s = popcnt(diag_mask);
  }
  return rank;
}

// Again directly operating on IdxState2 seems possible and even
// advantageous. pair_mask and diag_mask will mask out previously
// set bits.
uint64_t rank_bb_ref(Bitboard *set_bb, const struct RankInfo *ri)
{
  Bitboard occ = set_bb[0];
  Bitboard pair_mask = LOWER_DIAG_MASK;
  Bitboard diag_mask = MAIN_DIAG_MASK & ~occ;
  int p = 28, s = 6;

  uint64_t rank = 0;
  for (int k = 0; k < ri->numsets; k++) {
    Bitboard bb = set_bb[k + 1];
    occ |= bb;
    int tid = ri->transition_id[k];

    // Perhaps execute flip_main_8xbb() at the start, so no need to
    // repeat?
    Bitboard mirror = flip_main_diag(bb);
    Bitboard full_mask = bb & mirror  & pair_mask;
    Bitboard one_mask = (bb ^ mirror) & pair_mask;
    int d = popcnt(full_mask);    // Number of 2-orbits fully filled.

    const struct TransitionCase *c = &transition_cases[tid][fold_ps(p, s)][d];

    rank += c->prefix;
    rank += rank_combination(full_mask, pair_mask) * c->per_full_block;
    pair_mask &= ~full_mask;
    p -= d;

    if (one_mask) {
      // The current set breaks symmetry.
      int one = popcnt(one_mask);     // Number of 2-orbits half filled.
      int f = popcnt(bb & diag_mask); // Number of 1-orbits filled.
      rank += c->diag_block;
      uint64_t r = count_broken_residual_before(c->rem, p, s, one);

      uint64_t rone = rank_combination(one_mask, pair_mask);
      r += (rone * binom(s, f) + rank_combination(bb, diag_mask)) << (one - 1);
      rank += r * c->broken_tail;

      // Canonical orientation: orient_mask <= bitwise complement within oo bits.
      uint32_t orient_mask = _pext_u64(bb, one_mask);
      uint32_t comp = (~orient_mask) & ((1u << one) - 1u);
      uint32_t canon = orient_mask < comp ? orient_mask : comp;
      // Among 2^oo orientations paired with complements, canonical ones are
      // exactly 0 .. 2^(oo-1)-1 after min(x,~x) for this ordering.
      assert(canon < (1u << (one - 1)));
      rank += canon * c->broken_tail;

      if (comp < orient_mask) {
        alignas(64) Bitboard mirror_bb[8];
        set_bb[7] = occ;
        __m512i x = _mm512_load_si512((__m512i *)set_bb);
        x = flip_main_8xbb(x);
        _mm512_load_si512((__m512i *)mirror_bb);
        occ = mirror_bb[7];
        return rank + rank_bb_from(mirror_bb, k + 1, occ, ri);
      }
      return rank + rank_bb_from(set_bb, k + 1, occ, ri);
    }

    rank += rank_combination(bb, diag_mask) * c->diag_tail;
    diag_mask &= ~bb;
    s = popcnt(diag_mask);
  }
  return rank;
}

uint64_t sq_to_idx_ref(uint8_t *restrict sq)
{
  Bitboard occ = bit(sq[0]) | bit(sq[1]);

  return rank_reflection(sq, occ, ri.first, &ri);
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

static Bitboard unrank_trivial(uint8_t *restrict sq, uint64_t idx, int g,
    Bitboard occ, const struct RankInfo *ri)
{
  uint32_t sub[MAX_SETS];
  for (int k = ri->numsets - 1; k >= g; k--)
    idx = divmod_recip(idx, ri->factor[k], ri->recip[k], &sub[k]);
  for (; g < ri->numsets; g++)
    occ = unrank_binomial(sub[g], ri->mult[g], sq + ri->first[g], occ);
  return occ;
}

Bitboard unrank_reflection(uint64_t idx, uint8_t *restrict sq, Bitboard occ,
    const struct RankInfo *ri)
{
  Bitboard pair_mask = LOWER_DIAG_MASK;
  Bitboard diag_mask = MAIN_DIAG_MASK & ~occ;
  int p = 28, s = 6;

  for (int k = 0; k < ri->numsets; k++) {
    int tid = ri->transition_id[k];

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
      for (int i = ri->first[k]; bb; i++)
        sq[i] = pop_lsb(&bb);
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
    for (int i = ri->first[k]; bb; i++)
      sq[i] = pop_lsb(&bb);
    return unrank_trivial(sq, idx, k + 1, occ, ri);
  }
  return occ;
}
