/*
  Copyright (c) 2025 Ronald de Man

  This file is distributed under the terms of the GNU GPL, version 2.
*/

#include <assert.h>
#include <inttypes.h>
#include <stdbit.h>
#include <stdio.h>
#include <stdlib.h>

#include "defs.h"
#include "rans.h"
#include "types.h"
#include "util.h"

void free_code_rans(struct RansCode *c)
{
  free(c);
}

// Based on from https://github.com/rygorous/ryg_rans.
static void normalize_freqs(struct RansCode *c, int64_t *freq, int num_syms)
{
  uint64_t total = 0;
  for (int i = 0; i < num_syms; i++)
    total += freq[i];

  uint64_t sum = c->cum_freq[0] = 0;
  for (int i = 0; i < num_syms; i++) {
    sum += freq[i];
    c->cum_freq[i + 1] = sum * (1 << 16) / total;
  }
  for (int i = num_syms + 1; i < 4096; i++)
    c->cum_freq[i] = c->cum_freq[num_syms];

  assert(c->cum_freq[4095] == (uint16_t)(1 << 16));

  for (int i = 0; i < num_syms; i++)
    if (freq[i] && c->cum_freq[i + 1] == c->cum_freq[i]) {
      uint32_t best_freq = UINT32_MAX;
      int best_steal = -1;
      for (int j = 0; j < 4095; j++) {
        uint32_t fr = (uint16_t)(c->cum_freq[j + 1] - c->cum_freq[j]);
        if (fr > 1 && fr < best_freq) {
          best_freq = fr;
          best_steal = j;
        }
      }
      assert(best_steal >= 0);
      if (best_steal < i) {
        for (int j = best_steal + 1; j <= i; j++)
          c->cum_freq[j]--;
      } else {
        for (int j = i + 1; j <= best_steal; j++)
          c->cum_freq[j]++;
      }
    }

  for (int i = 0; i < 4095; i++)
    c->d.freq[i] = c->cum_freq[i + 1] - c->cum_freq[i];
  c->d.freq[4095] = 0;

  assert(c->cum_freq[0] == 0 && c->cum_freq[4095] == (uint16_t)(1 << 16));
  for (int i = 0; i < num_syms; i++) {
    if (freq[i] == 0)
      assert(c->d.freq[i] == 0);
    else
      assert(c->d.freq[i] > 0);
  }
}

static void sort_freqs(struct RansCode *c)
{
  for (int i = 0; i < 4096; i++)
    c->map[i] = i;
  for (int i = 0; i < 4096; i++)
    for (int j = i + 1; j < 4096; j++)
      if (c->d.freq[i] < c->d.freq[j]) {
        Swap(c->d.freq[i], c->d.freq[j]);
        Swap(c->map[i], c->map[j]);
      }
  for (int i = 0; i < 4096; i++)
    c->inv[c->map[i]] = i;
  uint32_t sum = 0;
  for (int i = 0; i < 4096; i++) {
    c->cum_freq[i] = sum;
    sum += c->d.freq[i];
  }
}

// Again based on https://github.com/rygorous/ryg_rans but with a smaller
// alias data structure.
void make_alias_table(struct RansDecode *d, struct RansCode *c)
{
  uint32_t sum = 1u << 16;
  uint16_t tgt_sum = sum / 4096;

  uint16_t remaining[4096];
  uint16_t alias[4096];
  uint8_t divider[4096];
  for (int i = 0; i < 4096; i++) {
    remaining[i] = d->freq[i];
    divider[i] = tgt_sum;
    alias[i] = i;
  }

  int cur_large = 0, cur_small = 0;
  for (;;) {
    while (cur_large < 4096 && remaining[cur_large] <= tgt_sum)
      cur_large++;
    if (cur_large == 4096) break;
    while (remaining[cur_small] >= tgt_sum)
      cur_small = (cur_small + 1) & 4095;
    alias[cur_small] = cur_large;
    divider[cur_small] = remaining[cur_small];
    remaining[cur_large] -= tgt_sum - remaining[cur_small];
    remaining[cur_small] = tgt_sum;
  }

  uint16_t assigned[4096];
  for (int b = 0; b < 4096; b++) {
    if (c)
      for (int i = 0; i < divider[b]; i++)
        c->alias_remap[c->cum_freq[b] + i] = b * tgt_sum + i;
    assigned[b] = divider[b];
    d->entry[b].alias_div = (alias[b] << 4) | (divider[b] & 15);
  }
  for (int b = 0; b < 4096; b++) {
    int s = alias[b];
    d->entry[b].start = s != b ? assigned[s] - divider[b] : 0;
    if (c)
      for (int i = divider[b]; i < tgt_sum; i++)
        c->alias_remap[c->cum_freq[s] + assigned[s] + i - divider[b]] =
          b * tgt_sum + i;
    assigned[s] += tgt_sum - divider[b];
  }

  // sanity check
  for (int i = 0; i < 4096; i++)
    assert(assigned[i] == d->freq[i]);
}

struct RansCode *create_code_rans(int64_t *freq, int num_syms)
{
  struct RansCode *c = malloc(sizeof *c);
  normalize_freqs(c, freq, num_syms);
  sort_freqs(c);
  make_alias_table(&(c->d), c);
  return c;
}

// Encode value v satisfying v_min <= v <= v_max.
INLINE void rans_enc_val_uni(RansState *r, uint32_t v, uint32_t v_min,
    uint32_t v_max, uint8_t **p)
{
  if (v_min == v_max) return;
  v -= v_min;
  uint32_t range = v_max - v_min + 1;

  // approximately uniform distribution
  uint32_t cum_freq = (v << 16) / range;
  uint32_t freq = ((v + 1) << 16) / range - cum_freq;

  uint64_t x_max = (uint64_t)freq << 47;
  uint64_t x = *r;

  if (x >= x_max) {
    *p -= 4;
    write_le_u32(*p, (uint32_t)x);
    x >>= 32;
  }

  *r = ((x / freq) << 16) + (x % freq) + cum_freq;
}

// Encode value v satisfying v_min <= v <= v_max. Use probability f >> 16
// for each value v < v_max and leave the rest for v == v_max.
INLINE void rans_enc_val_top(RansState *r, uint32_t v, uint32_t v_min,
   uint32_t v_max, uint32_t f, uint8_t **p)
{
  assert((v_max - v_min + 1) * f < 65536);

  if (v_min == v_max) return;

  uint32_t cum_freq = (v - v_min) * f;
  uint32_t freq = v == v_max ? 65536 - cum_freq : f;

  uint64_t x_max = (uint64_t)freq << 47;
  uint64_t x = *r;

  if (x >= x_max) {
    *p -= 4;
    write_le_u32(*p, (uint32_t)x);
    x >>= 32;
  }

  *r = ((x / freq) << 16) + (x % freq) + cum_freq;
}

// Decode value v in the range 0 <= v < range
INLINE uint32_t rans_dec_val_uni(RansState *r, uint32_t v_min, uint32_t v_max)
{
  if (v_min == v_max)
    return v_min;
  uint32_t range = v_max - v_min + 1;

  uint64_t x = *r;
  uint32_t w = x & 0xffff;

  uint32_t v = (w * range + range - 1) >> 16;
  uint32_t cum_freq = (v << 16) / range;
  uint32_t freq = ((v + 1) << 16) / range - cum_freq;
  assert(cum_freq <= w && w < cum_freq + freq);
  *r = freq * (x >> 16) + (w - cum_freq);

  return v + v_min;
}

INLINE uint16_t rans_dec_val_top(RansState *r, uint32_t v_min, uint32_t v_max,
    uint32_t f)
{
  if (v_min == v_max)
    return v_min;

  uint64_t x = *r;

  uint32_t w = x & 0xffff;
  uint32_t v = v_min + min(w / f, v_max - v_min);
  uint32_t cum_freq = (v - v_min) * f;
  uint32_t freq = v == v_max ? 65536 - cum_freq : f;
  assert(cum_freq <= w && w < cum_freq + freq);
  *r = freq * (x >> 16) + (w - cum_freq);

  return v;
}

size_t write_freq_table(FILE *F, struct RansCode *c, int num)
{
  uint8_t buf[10000], *p = buf + 10000;

  RansState rans;
  rans_enc_init(&rans, false);

  int bound[4096];
  bound[num] = 0;
  int R = 0;
  for (int s = num - 1; s >= 0; s--) {
    int P = s > 0 ? c->d.freq[s - 1] : INT32_MAX;
    R += c->d.freq[s];
    int v_max = min(R, P);
    bound[s] = max(1, max(bound[s + 1], v_max - c->d.freq[s]));
    if (bound[s] > bound[s + 1])
      continue;
    int v_min = max((R + num - s - 1) / (num - s), v_max - bound[s]);
    if (v_max * (v_max - v_min + 1) >= (num - s))
      rans_enc_val_uni(&rans, c->d.freq[s], v_min, v_max, &p);
    else
      rans_enc_val_top(&rans, c->d.freq[s], v_min, v_max,
          (v_max << 16) / (num - s), &p);
  }

  for (int s = num - 1; s > 0; s--) {
    int P = bound[s - 1];
    rans_enc_val_top(&rans, bound[s], 1, P, 1, &p);
  }
  rans_enc_val_uni(&rans, bound[0], 0, (1<<16) - ((1<<16) + num - 1) / num, &p);
  rans_enc_val_uni(&rans, num, 2, 4095, &p);
  rans_enc_flush(&rans, &p, p - sizeof(RansState));

  size_t numbytes = buf + 10000 - p;
  file_write(p, numbytes, F);

  return numbytes;
}

const uint8_t *read_freq_table(struct RansDecode *d, int *num_syms,
    const uint8_t *p)
{
  RansState rans;

  rans_dec_init(&rans, &p, p + sizeof(RansState));
  int num = rans_dec_val_uni(&rans, 2, 4095);
  *num_syms = num;
  rans_dec_renorm(&rans, &p);
  d->freq[0] = rans_dec_val_uni(&rans, 0,
      (1 << 16) - ((1 << 16) + num - 1) / num);
  for (int s = 1; s < num; s++) {
    rans_dec_renorm(&rans, &p);
    d->freq[s] = rans_dec_val_top(&rans, 1, d->freq[s - 1], 1);
  }
  int R = 1 << 16;
  for (int s = 0; s < num; s++) {
    int P = s > 0 ? d->freq[s - 1] : INT32_MAX;
    if (d->freq[s] > d->freq[s + 1] && s < num - 1) {
      d->freq[s] = min(R, P) - d->freq[s];
      R -= d->freq[s];
      continue;
    }
    rans_dec_renorm(&rans, &p);
    int v_max = min(R, P);
    int v_min = max((R + num - s - 1) / (num - s), v_max - d->freq[s]);
    if (v_max * (v_max - v_min + 1) >= (num - s))
      d->freq[s] = rans_dec_val_uni(&rans, v_min, v_max);
    else {
      d->freq[s] = rans_dec_val_top(&rans, v_min, v_max,
          (v_max << 16) / (num - s));
    }
    R -= d->freq[s];
  }

  return p;
}
