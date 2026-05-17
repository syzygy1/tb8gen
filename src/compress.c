/*
  Copyright (c) 2011-2017, 2025, 2026 Ronald de Man

  This file is distributed under the terms of the GNU GPL, version 2.
*/

#include <inttypes.h>
#include <math.h>
#include <stdbit.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "checksum.h"
#include "compress.h"
#include "defs.h"
#include "huffman.h"
#include "index.h"
#include "joinp.h"
#include "permute.h"
#include "probe.h"
#include "rans.h"
#include "tb8gen.h"
#include "threads.h"
#include "types.h"
#include "util.h"

static_assert(__STDC_ENDIAN_NATIVE__ == __STDC_ENDIAN_LITTLE__);

int g_compress_type;

static constexpr int MAX_NEW = 250;

static struct DtzMap *current_map;

// If DTZ Huffman coding loses both more than 2% against entropy and
// more than 1/64 bit per symbol, use rANS to avoid integer-code redundancy.
static constexpr double HUFFMAN_MAX_RELATIVE_REDUNDANCY = 0.02;
static constexpr double HUFFMAN_MAX_ABSOLUTE_REDUNDANCY = 1.0 / 64.0;

static double calc_entropy_bits(const int64_t *freq, int num_syms,
    uint64_t *total_freq)
{
  uint64_t total = 0;
  for (int i = 0; i < num_syms; i++)
    total += freq[i];

  *total_freq = total;
  if (total == 0) return 0.0;

  double entropy = 0.0;
  for (int i = 0; i < num_syms; i++)
    if (freq[i])
      entropy += freq[i] * log2((double)total / freq[i]);

  return entropy;
}

static uint64_t calc_huffman_bits(const struct HuffCode *c)
{
  uint64_t bits = 0;
  for (int i = 0; i < c->num_syms; i++)
    bits += c->length[i] * c->freq[i];

  return bits;
}

static bool huffman_misses_entropy(const struct HuffCode *c, int num_syms)
{
  uint64_t total;
  double entropy = calc_entropy_bits(c->freq, num_syms, &total);
  if (total == 0) return false;

  double huff_bits = calc_huffman_bits(c);
  double redundancy_per_symbol = (huff_bits - entropy) / total;
  double relative_redundancy = entropy > 0.0
      ? (huff_bits - entropy) / entropy
      : INFINITY;

  return redundancy_per_symbol > HUFFMAN_MAX_ABSOLUTE_REDUNDANCY
      && relative_redundancy > HUFFMAN_MAX_RELATIVE_REDUNDANCY;
}

static int64_t *release_huffman_freq(struct HuffCode *c)
{
  int64_t *freq = c->freq;
  free(c);
  return freq;
}

struct Symbol {
  uint16_t pattern[2];
  union {
    struct {
      uint8_t first, second;
    };
    uint16_t sym;
  };
  int len;
};

static struct Symbol symtable[MAXSYMB];
static int64_t pairfreq[MAXSYMB + 1][MAXSYMB + 1];
static uint16_t symcode[256][256];
static int replace[MAXSYMB];

static int num_syms;
static uint64_t num_blocks, real_num_blocks;
static uint32_t num_indices;
static uint16_t *sizetable;
static uint8_t *indextable;
static int blocksize;
static int idxbits;
static int num_vals;

static int num_ctrl, cur_ctrl_idx;

struct {
  int64_t freq; // maybe remove later; currently still used for dtz
  int s1, s2;
  int sym;
} newpairs[MAX_NEW + 25];

struct {
  int64_t freq;
  int s1, s2;
} paircands[MAX_NEW];

static uint8_t newtest[MAXSYMB][MAXSYMB];
static uint32_t (*countfirst)[MAX_NEW + 1][MAXSYMB + 1];
static uint32_t (*countsecond)[MAX_NEW + 1][MAXSYMB + 1];

static uint64_t *restrict work = nullptr, *restrict work_adj = nullptr;

static struct {
  void *data;
  size_t size;
} compress_state;

static uint8_t pairfirst[4][MAXSYMB], pairsecond[4][MAXSYMB];

static int t1[4][4], t2[4][4];
static int dcfreq[4][4];
static int wdl_vals[5];
//static uint8_t wdl_flags;
static uint8_t dc_to_val[4];

//static int freq_t1, freq_t2;
//static int64_t dcfreq_dtz;

static const size_t BUFFER_SIZE = 256 * 1024;
static void *out_buffer;

static int64_t (*countfreq_wdl)[9][16];
static void *countfreq_dtz;
typedef int64_t cf_u8[255][255];
typedef int64_t cf_u16[MAX_VAL][MAX_VAL];

struct BWState {
  FILE *F;
  uint8_t *buf;
  uint64_t bitbuf;
  size_t idx, cap;
  int bitcnt;
};

INLINE void init_bits(struct BWState *bw, FILE *F, uint8_t *buf, size_t cap)
{
  *bw = (struct BWState){.F = F, .buf = buf, .cap = cap};
}

INLINE void put_bits(struct BWState *bw, uint32_t bits, int n)
{
  bw->bitbuf = bw->bitbuf << n | bits;
  bw->bitcnt += n;
  if (bw->bitcnt >= 32) {
    write_be_u32(bw->buf + bw->idx, bw->bitbuf >> (bw->bitcnt - 32));
    bw->bitcnt -= 32;
    bw->idx += 4;
    if (bw->idx == bw->cap) {
      file_write(bw->buf, bw->cap, bw->F);
      bw->idx = 0;
    }
  }
}

INLINE void pad_bits(struct BWState *bw, int n)
{
  for (; n >= 32; n -= 32)
    put_bits(bw, 0, 32);
  put_bits(bw, 0, n);
}

INLINE void flush_bits(struct BWState *bw)
{
  file_write(bw->buf, bw->idx, bw->F);
}

void compress_alloc_wdl()
{
  countfirst = malloc(g_num_threads * sizeof(*countfirst));
  countsecond = malloc(g_num_threads * sizeof(*countsecond));
  countfreq_wdl = malloc(g_num_threads * sizeof(*countfreq_wdl));
  if (!countfirst || !countsecond || !countfreq_wdl)
    out_of_mem();
}

void compress_free_wdl()
{
  free(countfirst);
  free(countsecond);
  free(countfreq_wdl);
}

void compress_alloc_dtz(bool wide)
{
  countfirst = malloc(g_num_threads * sizeof(*countfirst));
  countsecond = malloc(g_num_threads * sizeof(*countsecond));
  countfreq_dtz = malloc(
      g_num_threads * (wide ? sizeof(cf_u16) : sizeof(cf_u8)));
  if (!countfirst || !countsecond || !countfreq_dtz)
    out_of_mem();
}

void compress_free_dtz(void)
{
  free(countfirst);
  free(countsecond);
  free(countfreq_dtz);
}

static void calc_symbol_tweaks(struct HuffCode *c)
{
  // Find longer symbols with shorter code as potential replacement at the
  // end of a block.

  for (int i = 0; i < num_syms; i++)
    replace[i] = i;

  for (int i = 0; i < num_syms; i++) {
    if (c->freq[i] == 0) continue;
    int l = c->length[i];
    int s = i;
    while (symtable[s].len > 1) {
      s = symtable[s].pattern[0];
      if (c->length[replace[s]] > l)
	replace[s] = i;
    }
  }
}

#define read_symbol(x) _Generic((x), \
  u8: symcode[x][(&(x))[1]], \
  u16: x \
);

#define write_symbol(x, y) _Generic((x), \
  u8: ((x) = (y)->first, (&(x))[1] = (y)->second ), \
  u16: (x) = (y)->sym \
);

#define T u8
#include "compress_tmpl.c"
#undef T

#define T u16
#include "compress_tmpl.c"
#undef T

void compress_init_wdl(bool vals[])
{
  int i;

  for (i = 0; i < 5; i++)
    wdl_vals[i] = vals[i];

  num_vals = 0;
  for (i = 0; i < 5; i++)
    if (wdl_vals[i])
      num_vals++;

  if (num_vals == 0) {
    wdl_vals[2] = true;
    num_vals = 1;
  }

  dc_to_val[0] = wdl_vals[0] ? 0 : 1;
  dc_to_val[1] = wdl_vals[2] ? 2 : dc_to_val[0];
  dc_to_val[2] = wdl_vals[2] ? 2 : (wdl_vals[0] ? 0 : (wdl_vals[3] ? 3 : 1));
  dc_to_val[3] = wdl_vals[4] ? 4 : dc_to_val[2];

  g_compress_type = 0;
  if (num_vals == 1)
    g_compress_type = 1;

  if (!out_buffer && !(out_buffer = malloc(BUFFER_SIZE)))
    out_of_mem();
}

static void count_pairs_wdl(struct ThreadData *thread)
{
  uint64_t idx = thread->begin;
  uint64_t end = thread->end;
  uint8_t *restrict data = compress_state.data;
  int t = thread->thread_id;

  if (idx == 0) idx = 1;
  int s1 = data[idx - 1];
  for (; idx < end; idx++) {
    int s2 = data[idx];
    countfreq_wdl[t][s1][s2]++;
    s1 = s2;
  }
}

static void remove_wdl_worker(struct ThreadData *thread)
{
  uint8_t *restrict data = compress_state.data;
  uint64_t size = compress_state.size;
  int s, t;

  static bool dc_map[5][16] = {
    { 1, 0, 0, 0, 0, 1, 1, 1, 1 },
    { 0, 1, 0, 0, 0, 1, 1, 1, 1 },
    { 0, 0, 1, 0, 0, 0, 1, 1, 1 },
    { 0, 0, 0, 1, 0, 0, 0, 1, 1 },
    { 0, 0, 0, 0, 1, 0, 0, 0, 1 }
  };

  for (uint64_t idx = thread->begin, end = thread->end; idx < end; ) {
    for (; idx < end; idx++)
      if (data[idx] >= 5) break;
    if (idx == end) break;
    // now idx points to a first dontcare of possibly more
    s = data[idx];
    uint64_t idx2;
    for (idx2 = idx + 1; idx2 < end; idx2++) {
      if (data[idx2] < 5) break;
      if (data[idx2] < s)
        s = data[idx2];
    }
    int f = 0;
    if (idx > 0 && dc_map[data[idx - 1]][s]) {
      t = data[idx - 1];
      f = 1;
      if (idx > 1 && data[idx - 2] == t)
        f = 2;
    }
    if (idx2 < size && dc_map[data[idx2]][s] && f != 2) {
      if (!f || (idx2 < size - 1 && dc_map[data[idx2]][data[idx2 + 1]]))
        t = data[idx2];
      f = 1;
    }
    if (f)
      for (; idx < idx2; idx++)
        data[idx] = t;
    else if (idx2 - idx > 5) {
      s = dc_to_val[s - 5];
      for (; idx < idx2; idx++)
        data[idx] = s;
    }
    idx = idx2;
  }
}

static void adjust_work_dontcares_wdl(uint64_t *restrict work1,
    uint64_t *restrict work2)
{
  uint64_t end = work1[g_total_work];
  uint8_t *restrict data = compress_state.data;

  work2[0] = work1[0];
  for (int i = 1; i < g_total_work; i++) {
    uint64_t idx = work1[i];
    if (idx < work2[i - 1]) {
      work2[i] = work2[i - 1];
      continue;
    }
    while (idx < end && (data[idx -1] >= 5 || data[idx] >= 5))
      idx++;
    work2[i] = idx;
  }
  work2[g_total_work] = work1[g_total_work];
}

static int d[5] = { 5, 5, 6, 7, 8 };

int64_t *construct_pairs_wdl(void *data, uint64_t size, int minfreq,
    int maxsymbols, int *nsyms)
{
  int k;

  if (!work) {
    work = alloc_work(g_total_work);
    work_adj = alloc_work(g_total_work);
  }

  compress_state.data = data;
  compress_state.size = size;
  fill_work(g_total_work, size, 0, work);

  for (int i = 0; i < 9; i++) {
    symtable[i].pattern[0] = i;
    symtable[i].len = 1;
    for (int j = 0; j < 256; j++)
      symcode[i][j] = i;
  }

  adjust_work_dontcares_wdl(work, work_adj);
  run_threaded(remove_wdl_worker, work_adj, 0);

  num_syms = 9;

  num_ctrl = num_syms - 1;
  cur_ctrl_idx = 0;

  if (maxsymbols > 0) {
    maxsymbols += num_syms;
    if (maxsymbols > MAXSYMB)
      maxsymbols = MAXSYMB;
  } else {
    maxsymbols = MAXSYMB + 1 - num_vals;
  }

  // FIXME: do this only once
  for (int i = 0; i < MAXSYMB; i++)
    for (int j = 0; j < MAXSYMB; j++)
      newtest[i][j] = 0;

  int num;

  for (int i = 0; i < num_syms; i++)
    for (int j = 0; j < num_syms; j++)
      pairfreq[i][j] = 0;

  for (int t = 0; t < g_num_threads; t++)
    for (int i = 0; i < num_syms; i++)
      for (int j = 0; j < num_syms; j++)
        countfreq_wdl[t][i][j] = 0;
  run_threaded(count_pairs_wdl, work, 0);
  for (int t = 0; t < g_num_threads; t++)
    for (int i = 0; i < num_syms; i++)
      for (int j = 0; j < num_syms; j++)
        pairfreq[i][j] += countfreq_wdl[t][i][j];

  while (num_syms < maxsymbols) {

    num = 0;
    for (int i = 0; i < 5; i++) {
      if (!wdl_vals[i]) continue;
      for (int j = 0; j < 5; j++) {
        if (!wdl_vals[j]) continue;
        int64_t pf = pairfreq[i][j];
        for (k = d[i]; k < 9; k++)
          pf += pairfreq[k][j];
        for (int l = d[j]; l < 9; l++)
          pf += pairfreq[i][l];
        for (k = d[i]; k < 9; k++)
          for (int l = d[j]; l < 9; l++)
            pf += pairfreq[k][l];
        if (   pf >= minfreq
            && (num < MAX_NEW || pf > paircands[MAX_NEW - 1].freq))
        {
          for (k = 0; k < num; k++)
            if (paircands[k].freq < pf) break;
          if (num < MAX_NEW) num++;
          for (int l = num - 1; l > k; l--)
            paircands[l] = paircands[l - 1];
          paircands[k].freq = pf;
          paircands[k].s1 = i;
          paircands[k].s2 = j;
        }
      }
    }

    for (int i = 0; i < 5; i++) {
      if (!wdl_vals[i]) continue;
      for (int j = 9; j < num_syms; j++) {
        if (1 + symtable[j].len > 256) continue;
        int64_t pf = pairfreq[i][j];
        for (k = d[i]; k < 9; k++)
          pf += pairfreq[k][j];
        if (    pf >= minfreq
            && (num < MAX_NEW || pf > paircands[MAX_NEW - 1].freq))
        {
          for (k = 0; k < num; k++)
            if (paircands[k].freq < pf) break;
          if (num < MAX_NEW - 1) num++;
          for (int l = num - 1; l > k; l--)
            paircands[l] = paircands[l - 1];
          paircands[k].freq = pf;
          paircands[k].s1 = i;
          paircands[k].s2 = j;
        }
      }
    }

    for (int i = 9; i < num_syms; i++) {
      if (symtable[i].len + 1 > 256) continue;
      for (int j = 0; j < 5; j++) {
        if (!wdl_vals[j]) continue;
        int64_t pf = pairfreq[i][j];
        for (int l = d[j]; l < 9; l++)
          pf += pairfreq[i][l];
        if (   pf >= minfreq
            && (num < MAX_NEW || pf > paircands[MAX_NEW - 1].freq))
        {
          for (k = 0; k < num; k++)
            if (paircands[k].freq < pf) break;
          if (num < MAX_NEW - 1) num++;
          for (int l = num - 1; l > k; l--)
            paircands[l] = paircands[l - 1];
          paircands[k].freq = pf;
          paircands[k].s1 = i;
          paircands[k].s2 = j;
        }
      }
    }

    for (int i = 9; i < num_syms; i++)
      for (int j = 9; j < num_syms; j++) {
        int64_t pf = pairfreq[i][j];
        if (   pf >= minfreq
            && (num < MAX_NEW || pf > paircands[MAX_NEW - 1].freq)
            && (symtable[i].len + symtable[j].len <= 256))
        {
          for (k = 0; k < num; k++)
            if (paircands[k].freq < pf) break;
          if (num < MAX_NEW - 1) num++;
          for (int l = num - 1; l > k; l--)
            paircands[l] = paircands[l - 1];
          paircands[k].freq = pf;
          paircands[k].s1 = i;
          paircands[k].s2 = j;
        }
      }

    // Estimate the number of skipped pairs to make sure they will be
    // considered in the next iteration (before running out of symbols).
    int skipped = 0;
    int64_t max = paircands[0].freq;
    int num2 = 0; // count new symbols
    int num3 = 0; // count new pairs
    for (int i = 0; i < num && num_syms + skipped < maxsymbols; i++) {
      while (max > paircands[i].freq * 2) {
        skipped += i;
        max /= 2;
      }
      int j, s1 = paircands[i].s1, s2 = paircands[i].s2;
      // Check for conflicts with more frequent candidates.
      for (j = 0; j < i; j++) {
        int t1 = paircands[j].s1;
        int t2 = paircands[j].s2;
        if (s2 == t1) break;
        if (s1 == t2) break;
        if (s2 < 5 && t1 < 5) {
          for (k = max(d[s2], d[t1]); k < 9; k++)
            if (pairfreq[s1][k] && pairfreq[k][t2]) goto lab;
          if (s1 < 5)
            for (int l = d[s1]; l < 9; l++)
              for (k = max(d[s2], d[t1]); k < 9; k++)
                if (pairfreq[l][k] && pairfreq[k][t2]) goto lab;
          if (t2 < 5)
            for (int l = d[t2]; l < 9; l++)
              for (k = max(d[s2], d[t1]); k < 9; k++)
                if (pairfreq[s1][k] && pairfreq[k][l]) goto lab;
          if (s1 < 5 && t2 < 5)
            for (int l = d[s1]; l < 9; l++)
              for (int m = d[t2]; m < 9; m++)
                for (k = max(d[s2], d[t1]); k < 9; k++)
                  if (pairfreq[l][k] && pairfreq[k][m]) goto lab;
        }
        if (t2 < 5 && s1 < 5) {
          for (k = max(d[t2], d[s1]); k < 9; k++)
            if (pairfreq[t1][k] && pairfreq[k][s2]) goto lab;
          if (t1 < 5)
            for (int l = d[t1]; l < 9; l++)
              for (k = max(d[t2], d[s1]); k < 9; k++)
                if (pairfreq[l][k] && pairfreq[k][s2]) goto lab;
          if (s2 < 5)
            for (int l = d[s2]; l < 9; l++)
              for (k = max(d[t2], d[s1]); k < 9; k++)
                if (pairfreq[t1][k] && pairfreq[k][l]) goto lab;
          if (t1 < 5 && s2 < 5)
            for (int l = d[t1]; l < 9; l++)
              for (int m = d[s2]; m < 9; m++)
                for (k = max(d[t2], d[s1]); k < 9; k++)
                  if (pairfreq[l][k] && pairfreq[k][m]) goto lab;
        }
        if (s1 < 5 && t1 < 5) {
          if (s2 == t2)
            for (k = max(d[s1], d[t1]); k < 9; k++)
              if (pairfreq[k][s2]) goto lab;
          if (s2 < 5 && t2 < 5)
            for (k = max(d[s1], d[t1]); k < 9; k++)
              for (int l = max(d[s2], d[t2]); l < 9; l++)
                if (pairfreq[k][l]) goto lab;
        }
        if (s2 < 5 && t2 < 5)
          if (s1 == t1)
            for (k = max(d[s2], d[t2]); k < 9; k++)
              if (pairfreq[s1][k]) goto lab;
      }
lab:
      if (j == i) { // no conflict
        int tmp = num3;
        newpairs[tmp].sym = num_syms;
        newpairs[tmp].s1 = s1;
        newpairs[tmp].s2 = s2;
        tmp++;
        if (s1 < 5)
          for (k = d[s1]; k < 9; k++)
            if (pairfreq[k][s2]) {
              newpairs[tmp].sym = num_syms;
              newpairs[tmp].s1 = k;
              newpairs[tmp].s2 = s2;
              tmp++;
            }
        if (s2 < 5)
          for (int l = d[s2]; l < 9; l++)
            if (pairfreq[s1][l]) {
              newpairs[tmp].sym = num_syms;
              newpairs[tmp].s1 = s1;
              newpairs[tmp].s2 = l;
              tmp++;
            }
        if (s1 < 5 && s2 < 5)
          for (k = d[s1]; k < 9; k++)
            for (int l = d[s2]; l < 9; l++)
              if (pairfreq[k][l]) {
                newpairs[tmp].sym = num_syms;
                newpairs[tmp].s1 = k;
                newpairs[tmp].s2 = l;
                tmp++;
              }
        if (tmp > MAX_NEW) break;
        num3 = tmp;
        num_syms++;
        num2++;
      } else
        skipped++;
    }

    num = num2;
    if (num == 0) break;

    for (int i = 0; i < num3;) {
      k = newpairs[i].sym;
      symtable[k].len =  symtable[newpairs[i].s1].len
                       + symtable[newpairs[i].s2].len;
      symtable[k].pattern[0] = newpairs[i].s1;
      symtable[k].pattern[1] = newpairs[i].s2;
      if (!cur_ctrl_idx)
        num_ctrl++;
      if (num_ctrl == 256) break;
      symtable[k].first = num_ctrl;
      symtable[k].second = cur_ctrl_idx;
      symcode[num_ctrl][cur_ctrl_idx] = k;
      for (i++; i < num3 && newpairs[i].sym == k; i++);
      cur_ctrl_idx++;
      if (cur_ctrl_idx == 256) cur_ctrl_idx = 0;
    }

    for (int i = 0; i < num_syms - num; i++)
      for (int j = num_syms - num; j < num_syms; j++)
        pairfreq[i][j] = 0;
    for (int i = num_syms - num; i < num_syms; i++)
      for (int j = 0; j < num_syms; j++)
        pairfreq[i][j] = 0;

    for (int i = 0; i < num3; i++)
      newtest[newpairs[i].s1][newpairs[i].s2] = i + 1;

    // Thread this later.
    for (int t = 0; t < g_num_threads; t++)
      for (int i = 0; i < num3; i++)
        for (int j = 0; j < num_syms; j++) {
          countfirst[t][i][j] = 0;
          countsecond[t][i][j] = 0;
        }

    adjust_work_replace_u8(work);
    run_threaded(replace_pairs_u8, work, 0);

    for (int t = 0; t < g_num_threads; t++)
      for (int i = 0; i < num3; i++)
        for (int j = 0; j < num_syms; j++) {
          pairfreq[j][newpairs[i].s1] -= countfirst[t][i][j];
          pairfreq[j][newpairs[i].sym] += countfirst[t][i][j];
          pairfreq[newpairs[i].s2][j] -= countsecond[t][i][j];
          pairfreq[newpairs[i].sym][j] += countsecond[t][i][j];
        }

    for (int i = 0; i < num3; i++) {
      pairfreq[newpairs[i].s1][newpairs[i].s2] = 0;
      newtest[newpairs[i].s1][newpairs[i].s2] = 0;
    }

  }

  int64_t *freq = setup_code_u8(data, size);

  // Map remaining don't cares to a value.
  for (k = 5; k < 9; k++)
    if (freq[k] > 0) {
      int i;
      for (i = 0; i < 5; i++)
        if (wdl_vals[i]) break;
      for (int j = i + 1; j < 5; j++)
        if (d[j] <= k && freq[j] > freq[i])
          i = j;
      symcode[k][0] = i;
      freq[i] += freq[k];
      freq[k] = 0;
    }

  // Remove unused symbols.
  uint16_t map[num_syms];
  int i;
  for (i = 0, k = 0; i < 5; i++)
    if (wdl_vals[i])
      map[i] = k++;
  for (int i = 9; i < num_syms; i++)
    map[i] = k++;
  for (i = 0, k = 0; i < 5; i++)
    if (wdl_vals[i]) {
      symtable[k] = symtable[i];
      for (int j = 0; j < 256; j++)
        symcode[i][j] = k;
      for (int l = 5; l < 9; l++)
        if (symcode[l][0] == i)
          for (int j = 0; j < 256; j++)
            symcode[l][j] = k;
      freq[k] = freq[i];
      k++;
    }
  for (i = 9; i < num_syms; i++) {
    symtable[k] = symtable[i];
    symcode[symtable[k].first][symtable[k].second] = k;
    symtable[k].pattern[0] = map[symtable[k].pattern[0]];
    symtable[k].pattern[1] = map[symtable[k].pattern[1]];
    freq[k] = freq[i];
    k++;
  }
  num_syms = k;

  *nsyms = num_syms;
  return freq;
}

void compress_init_dtz(struct DtzMap *map, bool wide)
{
  current_map = map;

  num_vals = map->max_num;

  for (int i = 0; i < num_vals; i++) {
    symtable[i].pattern[0] = i;
    symtable[i].len = 1;
    if (!wide)
      for (int j = 0; j < 256; j++)
        symcode[i][j] = i;
  }
  g_compress_type = 0;
  if (num_vals <= 1)
    g_compress_type = 1;

  if (!out_buffer && !(out_buffer = malloc(BUFFER_SIZE)))
    out_of_mem();
}

int64_t *construct_pairs_dtz(void *data, uint64_t size, int minfreq,
    int maxsymbols, int *nsyms, bool wide)
{
  int k;

  if (!work)
    work = alloc_work(g_total_work);
  if (!work_adj)
    work_adj = alloc_work(g_total_work);

  num_syms = num_vals;

  num_ctrl = num_vals - 1;
  cur_ctrl_idx = 0;
  compress_state.data = data;
  compress_state.size = size;

  fill_work(g_total_work, size, 0, work);

  if (maxsymbols > 0) {
    maxsymbols += num_syms;
    maxsymbols = min(maxsymbols, 4095);
  } else
    maxsymbols = 4095;

  if (!wide) {
    adjust_work_dontcares_dtz_u8(work, work_adj);
    run_threaded(remove_dtz_worker_u8, work_adj, 0);
  } else {
    adjust_work_dontcares_dtz_u16(work, work_adj);
    run_threaded(remove_dtz_worker_u16, work_adj, 0);
  }

  // first round of freq counting, to fill in dont cares
  for (int i = 0; i < num_syms; i++)
    for (int j = 0; j < num_syms; j++)
      pairfreq[i][j] = 0;

  if (!wide && num_syms > 255) {
    fprintf(stderr, "Ran out of symbols.\n");
    exit(EXIT_FAILURE);
  }

  if (!wide) {
    cf_u8 *countfreq = countfreq_dtz;
    memset(countfreq, 0, g_num_threads * sizeof *countfreq);
    run_threaded(count_pairs_dtz_u8, work, 0);
    for (int t = 0; t < g_num_threads; t++)
      for (int i = 0; i < num_syms; i++)
        for (int j = 0; j < num_syms; j++)
          pairfreq[i][j] += countfreq[t][i][j];
  } else {
    cf_u16 *countfreq = countfreq_dtz;
    memset(countfreq, 0, g_num_threads * sizeof *countfreq);
    run_threaded(count_pairs_dtz_u16, work, 0);
    for (int t = 0; t < g_num_threads; t++)
      for (int i = 0; i < num_syms; i++)
        for (int j = 0; j < num_syms; j++)
          pairfreq[i][j] += countfreq[t][i][j];
  }

  for (int j = 0; j < num_syms; j++)
    pairfirst[0][j] = pairsecond[0][j] = 0;
  for (int i = 1; i < num_syms; i++)
    for (int j = 0; j < num_syms; j++) {
      if (pairfreq[j][i] > pairfreq[j][pairfirst[0][j]])
        pairfirst[0][j] = i;
      if (pairfreq[i][j] > pairfreq[pairsecond[0][j]][j])
        pairsecond[0][j] = i;
    }

  dcfreq[0][0] = -1;
  for (int i = 0; i < num_syms; i++)
    for (int j = 0; j < num_syms; j++)
      if (pairfreq[i][j] > dcfreq[0][0]) {
        dcfreq[0][0] = pairfreq[i][j];
        t1[0][0] = i;
        t2[0][0] = j;
      }

  for (int i = 0; i < MAXSYMB; i++)
    for (int j = 0; j < MAXSYMB; j++)
      newtest[i][j] = 0;

  if (!wide) {
    adjust_work_dontcares_u8(work, work_adj);
    run_threaded(fill_dontcares_u8, work_adj, 0);
  } else {
    adjust_work_dontcares_u16(work, work_adj);
    run_threaded(fill_dontcares_u16, work_adj, 0);
  }

  for (int i = 0; i < num_syms; i++)
    for (int j = 0; j < num_syms; j++)
      pairfreq[i][j] = 0;

  if (!wide) {
    cf_u8 *countfreq = countfreq_dtz;
    memset(countfreq, 0, g_num_threads * sizeof *countfreq);
    run_threaded(count_pairs_dtz_u8, work, 0);
    for (int t = 0; t < g_num_threads; t++)
      for (int i = 0; i < num_syms; i++)
        for (int j = 0; j < num_syms; j++)
          pairfreq[i][j] += countfreq[t][i][j];
  } else {
    cf_u16 *countfreq = countfreq_dtz;
    memset(countfreq, 0, g_num_threads * sizeof *countfreq);
    run_threaded(count_pairs_dtz_u16, work, 0);
    for (int t = 0; t < g_num_threads; t++)
      for (int i = 0; i < num_syms; i++)
        for (int j = 0; j < num_syms; j++)
          pairfreq[i][j] += countfreq[t][i][j];
  }

  for (int i = 0; i < num_syms; i++)
    for (int j = 0; j < 256; j++)
      symcode[i][j] = i;

  while (num_syms < maxsymbols) {

    int num = 0;
    for (int i = 0; i < num_syms; i++)
      for (int j = 0; j < num_syms; j++)
        if (   pairfreq[i][j] >= minfreq
            && (num < MAX_NEW - 1 || pairfreq[i][j] > newpairs[num - 1].freq)
            && symtable[i].len + symtable[j].len <= 256)
        {
          for (k = 0; k < num; k++)
            if (newpairs[k].freq < pairfreq[i][j]) break;
          if (num < MAX_NEW - 1) num++;
          for (int l = num - 1; l > k; l--)
            newpairs[l] = newpairs[l - 1];
          newpairs[k].freq = pairfreq[i][j];
          newpairs[k].s1 = i;
          newpairs[k].s2 = j;
        }

    for (int i = 0; i < num_syms; i++)
      pairfirst[0][i] = pairsecond[0][i] = 0;

    // keep track of number of skipped pairs to make sure they'll be
    // considered in the next iteration (before running out of symbols)
    int skipped = 0; // just a rough estimate
    int64_t max = newpairs[0].freq;
    int i, j;
    for (i = 0, j = 0; i < num && num_syms + j + skipped <= maxsymbols; i++) {
      while (max > newpairs[i].freq * 2) {
        skipped += i;
        max /= 2;
      }
      if (!pairsecond[0][newpairs[i].s1] && !pairfirst[0][newpairs[i].s2])
        newpairs[j++] = newpairs[i];
      else
        skipped++;
      pairfirst[0][newpairs[i].s1] = 1;
      pairsecond[0][newpairs[i].s2] = 1;
    }
    num = j;

    if (num_syms + num > maxsymbols)
      num = maxsymbols - num_syms;

    for (i = 0; i < num; i++) {
      newpairs[i].sym = num_syms;
      symtable[num_syms].len =
        symtable[newpairs[i].s1].len + symtable[newpairs[i].s2].len;
      symtable[num_syms].pattern[0] = newpairs[i].s1;
      symtable[num_syms].pattern[1] = newpairs[i].s2;
      if (!wide) {
        if (!cur_ctrl_idx)
          num_ctrl++;
        if (num_ctrl == 256) break;
        symtable[num_syms].first = num_ctrl;
        symtable[num_syms].second = cur_ctrl_idx;
        symcode[num_ctrl][cur_ctrl_idx] = num_syms;
        cur_ctrl_idx++;
        if (cur_ctrl_idx == 256) cur_ctrl_idx = 0;
      } else {
        symtable[num_syms].sym = num_syms;
      }
      num_syms++;
    }

    if (i != num) {
      fprintf(stderr, "Ran short of symbols.\n");
      exit(EXIT_FAILURE);
    }

    if (num == 0) break;

    for (i = 0; i < num_syms - num; i++)
      for (int j = num_syms - num; j < num_syms; j++)
        pairfreq[i][j] = 0;
    for (; i < num_syms; i++)
      for (j = 0; j < num_syms; j++)
        pairfreq[i][j] = 0;

    for (int i = 0; i < num; i++)
      newtest[newpairs[i].s1][newpairs[i].s2] = i + 1;

    // thread this later
    for (int t = 0; t < g_num_threads; t++)
      for (int i = 0; i < num; i++)
        for (int j = 0; j < num_syms; j++) {
          countfirst[t][i][j] = 0;
          countsecond[t][i][j] = 0;
        }

    if (!wide) {
      adjust_work_replace_u8(work);
      run_threaded(replace_pairs_u8, work, 0);
    } else {
      adjust_work_replace_u16(work);
      run_threaded(replace_pairs_u16, work, 0);
    }

    for (int t = 0; t < g_num_threads; t++)
      for (int i = 0; i < num; i++)
        for (int j = 0; j < num_syms; j++) {
          pairfreq[j][newpairs[i].s1] -= countfirst[t][i][j];
          pairfreq[j][newpairs[i].sym] += countfirst[t][i][j];
          pairfreq[newpairs[i].s2][j] -= countsecond[t][i][j];
          pairfreq[newpairs[i].sym][j] += countsecond[t][i][j];
        }

    for (int i = 0; i < num; i++)
      pairfreq[newpairs[i].s1][newpairs[i].s2] = 0;

    for (int i = 0; i < num; i++)
      newtest[newpairs[i].s1][newpairs[i].s2] = 0;
  }

  int64_t *freq = !wide ? setup_code_u8(data, size)
                        : setup_code_u16(data, size);
  *nsyms = num_syms;

  return freq;
}

int64_t *construct_pairs(void *data, uint64_t size, int minfreq,
    int maxsymbols, int *nsyms, bool wide, bool wdl)
{
  return wdl ? construct_pairs_wdl(data, size, minfreq, maxsymbols, nsyms)
    : construct_pairs_dtz(data, size, minfreq, maxsymbols, nsyms, wide);
}

// Write out the final DTM tablebase file in the "old format" layout.
void write_final(struct tb_handle *F, FILE *G)
{
  if (    F->type == WDL
      || (F->type == DTZ && one_sided && !used_rans))
  {
    // The first 4 bytes form a magic number.
    write_u32(G, magic[F->type]);

    // The next byte holds the number of pieces and some flags.
    // split: separate tables for wtm and btm (always unless symmetric)
    // num_tables > 2: indicates a pawnful table, but this flag is not used.
    write_u8(G, (g_pos.num << 4) | (F->split ? 0x01: 0x00)
        | (F->num_tables > 2 ? 0x02 : 0x00));

    // A pawnless tablebase file contains 2 tables, or 1 if symmetric.
    // A pawnful tablebase file contains 8 tables, or 4 if symmetric (one per
    // file a/b/c/d).
  } else {
    // Files incompatible with old probing code get a new magic number.
    write_u32(G, magic2[F->type]);

    // Version 0 means the same table layout as the old format.
    write_u8(G, 0);

    // For non-WDL tables, we encode only wtm, btm, wins or losses.
    // For pawnful DTZ, this can vary per pawn slice/file/rank
    // -> to be looked into later.
    if (F->type != WDL) {
      uint8_t dist_format = !symmetric && !one_sided ? TWO_SIDED : 0;
      if (one_sided)
        dist_format |= WTM_OR_BTM | (one_sided_stm == WHITE ? WTM_ONLY : 0);
      else
        dist_format |= WIN_OR_LOSS | (wins_only ? WIN_ONLY : 0);
      write_u8(G, dist_format);
    }
  }

  // One "order" nibble for pawnless tables, two nibbles for pawnful tables.
  // For a non-split pawnful tablebase we pack the two nibbles per table in
  // one byte. Otherwise, we pack corresponding nibbles for wtm and btm
  // in one byte. The "numorder" value below is the number of encoded "order"
  // bytes per table.
  int numorder = 1;
  if (g_num_pawns > 0 && (F->perm[0][1] >> 4) != 0x0f)
    numorder = 2;

  // For each set of wtm/btm tables (1 or 6 sets of 2), write out the "order"
  // bytes and the permutations of piece types (low nibbles for one table,
  // high nibbles for the other table).
  if (F->split) {
    for (int i = 0; i < F->num_tables; i += 2) {
      for (int j = 0; j < numorder; j++)
        write_u8(G, (F->perm[i][j] >> 4) | (F->perm[i + 1][j] & 0xf0));
      for (int j = 0; j < g_pos.num; j++)
        write_u8(G, (F->perm[i][j] & 0x0f) | ((F->perm[i + 1][j] & 0x0f) << 4));
    }
  } else if (F->type == WDL) {
    for (int i = 0; i < F->num_tables; i++) {
      for (int j = 0; j < numorder; j++)
        write_u8(G, (F->perm[i][j] >> 4) | (F->perm[i][j] & 0xf0));
      for (int j = 0; j < g_pos.num; j++)
        write_u8(G, (F->perm[i][j] & 0x0f) | ((F->perm[i][j] & 0x0f) << 4));
    }
  } else {
    for (int i = 0; i < F->num_tables; i++) {
      for (int j = 0; j < numorder; j++)
        write_u8(G, (F->perm[i][j] >> 4));
      for (int j = 0; j < g_pos.num; j++)
        write_u8(G, (F->perm[i][j] & 0x0f));
    }
  }

  // Align on an even position in the file.
  if (ftell(G) & 0x01)
    write_u8(G, 0);

  // Write out information describing each table. If all positions in a table
  // have the same value, e.g. the wtm WDL table for KQvK (-> all wins),
  // we only write that single value.
  for (int i = 0; i < F->num_tables; i++) {
    if (!(F->flags[i] & 0x80)) {

      // Some flags. The most-significant bit indicates a single-value table.
      write_u8(G, F->flags[i]);

      // A compressed block of table i has 1 << blocksize bytes.
      write_u8(G, F->blocksize[i]);

      // The index of a position is divided by 1 << idxbits to find
      // a "main entry" in the table's index.
      write_u8(G, F->idxbits[i]);

      // A small number of non-existing blocks at the end of the index
      // to prevent an out-of-bounds error when accessing the index.
      write_u8(G, F->num_blocks[i] - F->real_num_blocks[i]);

      // The number of existing compressed blocks in the table.
      write_u32(G, F->real_num_blocks[i]);

      if (!(F->flags[i] & 0x40)) {
        // A representation of a canonical Huffman code.
        // For background, see: Alistair Moffat and Andrew Turpin,
        // "On the Implementation of Minimum Redundancy Prefix Codes",
        // IEEE Transactions on Communications, Vol. 45, No. 10, October 1997.
        // https://tinyurl.com/w5r92fw8
        struct HuffCode *c = F->code.huff[i];
        write_u8(G, c->max_len);
        write_u8(G, c->min_len);
        for (int j = c->min_len; j <= c->max_len; j++)
          write_u16(G, c->offset[j]);

        // Up to 4095 Re-Pair symbols. Each symbol is either the concatenation
        // of two other symbols (encoded as 2 x 12 bits = 3 bytes) or a single
        // WDL/DTZ value (encoded as a combination of 0xfff and 12 bits).
        struct Symbol *stable = F->symtable[i];
        write_u16(G, F->num_syms[i]);
        for (int j = 0; j < F->num_syms[i]; j++) {
          int k = c->map[j];
          if (stable[k].len == 1) {
            int s1 = stable[k].pattern[0];
            write_u8(G, s1 & 0xff);
            write_u8(G, (s1 >> 8) | 0xf0);
            write_u8(G, 0xff);
          } else {
            int s1 = c->inv[stable[k].pattern[0]];
            int s2 = c->inv[stable[k].pattern[1]];
            write_u8(G, s1 & 0xff);
            write_u8(G, (s1 >> 8) | ((s2 << 4) & 0xff));
            write_u8(G, (s2 >> 4));
          }
        }
      } else {
        // RANS. Write out the frequency table.
        struct RansCode *c = F->code.rans[i];
        size_t sz = write_freq_table(G, c, F->num_syms[i]);
        printf("num_syms = %d, size freq_table = %lu\n", F->num_syms[i], sz);

        // Up to 4095 Re-Pair symbols. Each symbol is either the concatenation
        // of two other symbols (encoded as 2 x 12 bits = 3 bytes) or a single
        // WDL/DTZ value (encoded as a combination of 0xfff and 12 bits).
        struct Symbol *stable = F->symtable[i];
        for (int j = 0; j < F->num_syms[i]; j++) {
          int k = c->map[j];
          if (stable[k].len == 1) {
            int s1 = stable[k].pattern[0];
            write_u8(G, s1 & 0xff);
            write_u8(G, (s1 >> 8) | 0xf0);
            write_u8(G, 0xff);
          } else {
            int s1 = c->inv[stable[k].pattern[0]];
            int s2 = c->inv[stable[k].pattern[1]];
            write_u8(G, s1 & 0xff);
            write_u8(G, (s1 >> 8) | ((s2 << 4) & 0xff));
            write_u8(G, (s2 >> 4));
          }
        }
      }
    } else {
      // The single value case.
      write_u8(G, F->flags[i]);
      write_u8(G, F->single_val[i]);
    }
    // Align on an even position.
    if (ftell(G) & 0x01)
      write_u8(G, 0);
  }

  // For DTZ tables, win, cursed win, blessed loss and loss distances are
  // mapped to the same range, each sorted by frequency. We have to store
  // these maps.
  if (F->type == DTZ) {
    for (int i = 0; i < F->num_tables; i++)
      if (F->flags[i] & 2) {
        if (!(F->flags[i] & 16)) {
          for (int j = 0; j < 4; j++) {
            write_u8(G, F->map[i]->num[j]);
            for (int k = 0; k < F->map[i]->num[j]; k++)
              write_u8(G, F->map[i]->map[j][k]);
          }
        } else {
          if (ftell(G) & 0x01)
            write_u8(G, 0);
          for (int j = 0; j < 4; j++) {
            write_u16(G, F->map[i]->num[j]);
            for (int k = 0; k < F->map[i]->num[j]; k++)
              write_u16(G, F->map[i]->map[j][k]);
          }
        }
      }
    if (ftell(G) & 0x01)
      write_u8(G, 0);
  }

  // Write the "main" entries of the index. Each entry consists of a 32-bit
  // value and a 16 bit value packed into 6 bytes. The 32-bit value of the
  // kth entry is the number of the block which contains the value of the
  // position with index I_k = k * (1 << indexbits) + (1 << (indexbits - 1)).
  // The 16-bit value is the offset within that compressed block at which
  // the value with this index is found. Note that each compressed block
  // contains at most 65536 values.
  for (int i = 0; i < F->num_tables; i++) {
    if (F->flags[i] & 0x80) continue;
    copy_data(G, F->H[i], 6 * F->num_indices[i]);
  }

  // For each table, write out a list of 16-bit values representing the
  // number of values stored in the compressed blocks. This is the 2nd level
  // index. To find the value of a position:
  // - Encode the position into an index I.
  // - Calculate k for which abs(I - I_k) is minimal.
  // - Use the main index to find the block and offset of the value of the
  //   position with index I_k.
  // - Use the 16-bit values in the 2nd level index, starting from the
  //   16-bit value corresponding to block which contains I_k, to find
  //   the block that contains I and the offset of I within that block.
  //
  //  The "indexbits" value is chosen such that about 8 16-bit values have
  //  to be added or subtracted on average.
  for (int i = 0; i < F->num_tables; i++) {
    if (F->flags[i] & 0x80) continue;
    copy_data(G, F->H[i], 2 * F->num_blocks[i]);
  }

  // Align to 64-byte boundary, so compressed blocks are aligned on cache
  // lines (unless blocks are just 32 bytes, which happens for some more
  // trivial tables).
  uint64_t idx = ftell(G);
  while (idx & 0x3f) {
    write_u8(G, 0);
    idx++;
  }

  // Copy the actual compressed data (which has previously been written
  // to temporary files).
  for (int i = 0; i < F->num_tables; i++) {
    if (F->flags[i] & 0x80) continue;
    uint64_t datasize = (uint64_t)F->real_num_blocks[i] << F->blocksize[i];
    copy_data(G, F->H[i], datasize);
    while (datasize & 0x3f) {
      write_u8(G, 0);
      datasize++;
    }
  }
}

static void compress_data(struct tb_handle *F, int num, FILE *G,
    void *data, uint64_t size, int minfreq, bool wide)
{
  int num_syms;
  int64_t *freq = construct_pairs(data, size, minfreq, 0, &num_syms, wide,
      F->type == WDL);
  bool rans = g_use_rans && F->type == DTZ;
  if (!rans) {
    struct HuffCode *c = create_code(freq, num_syms);

    if (F->type == DTZ && huffman_misses_entropy(c, num_syms)) {
      rans = true;
      freq = release_huffman_freq(c);
    } else {
      if (!wide)
        calc_block_sizes_u8(data, size, c, F->default_blocksize);
      else
        calc_block_sizes_u16(data, size, c, F->default_blocksize);

      F->code.huff[num] = c;
    }
  }

  if (rans) {
    used_rans = true;
    struct RansCode *c = create_code_rans(freq, num_syms);

    if (!wide)
      calc_block_sizes_rans_u8(data, size, c, F->default_blocksize);
    else
      calc_block_sizes_rans_u16(data, size, c, F->default_blocksize);

    F->code.rans[num] = c;
  }
  F->blocksize[num] = blocksize;
  F->idxbits[num] = idxbits;
  F->real_num_blocks[num] = real_num_blocks;
  F->num_blocks[num] = num_blocks;
  F->num_indices[num] = num_indices;
  F->num_syms[num] = num_syms;
  F->symtable[num] = malloc(sizeof(symtable));
  memcpy(F->symtable[num], symtable, sizeof(symtable));

  struct Symbol *stable = F->symtable[num];
  if (F->type == WDL) {
    F->flags[num] = 0x0c; // To get the same output as rtbgen produces.
  } else {
    F->map[num] = current_map;
    F->flags[num] =  current_map->stm | (1 << 2) | (1 << 3)
                   | ((int)current_map->wide << 4);
    int i, j, k;
    for (i = 0; i < 4; i++)
      if (current_map->num[i] == num_vals)
        break;
    for (j = 0; j < 4; j++)
      if (j != i) {
        for (k = 0; k < current_map->num[j]; k++)
          if (current_map->map[j][k] != current_map->map[i][k])
            break;
        if (k < current_map->num[j])
          break;
      }
    if (j == 4) {
      for (k = 0; k < num_vals; k++)
        stable[k].pattern[0] = current_map->map[i][k];
    } else {
      F->flags[num] |= (1 << 1);
    }
  }

  file_write(indextable, 6 * num_indices, G);
  file_write(sizetable, 2 * num_blocks, G);

  free(indextable);
  free(sizetable);

  if (!rans) {
    struct HuffCode *c = F->code.huff[num];
    if (!wide)
      write_ctb_data_u8(G, data, c, size, F->blocksize[num]);
    else
      write_ctb_data_u16(G, data, c, size, F->blocksize[num]);
  } else {
    struct RansCode *c = F->code.rans[num];
    F->flags[num] |= 0x40;
    if (!wide)
      write_ctb_data_rans_u8(G, data, c, size, F->blocksize[num]);
    else
      write_ctb_data_rans_u16(G, data, c, size, F->blocksize[num]);
  }
}

struct tb_handle *create_tb(const char *tablename, int type, int blocksize)
{
  struct tb_handle *F;

  F = calloc(1, sizeof(struct tb_handle));

  F->num_tables = 0;
  F->type = type;
  F->split = !symmetric;
  F->default_blocksize = blocksize;
  strcpy(F->name, tablename);

  return F;
}

void compress_data_single_valued(struct tb_handle *F, int num)
{
  int i;

  F->flags[num] = 0x8c;
  if (F->type == WDL) {
    for (i = 0; i < 5; i++)
      if (wdl_vals[i]) break;
    F->single_val[num] = i;
  } else {
    F->map[num] = current_map;
    F->single_val[num] = current_map->map[0][0];
  }
}

void compress_tb(struct tb_handle *F, int num, void *data, uint64_t tb_size,
    uint8_t *perm, int minfreq, bool wide)
{
  char name[64];
  char ext[16];

  if (num < 0)
    num = F->num_tables++;
  else if (F->num_tables < num + 1)
    F->num_tables = num + 1;

  for (int i = 0; i < g_pos.num; i++)
    F->perm[num][i] = perm[i];

  if (g_compress_type == 0) {
    sprintf(ext, ".%d", 1 + num);
    strcpy(name, F->name);
    strcat(name, suffix[F->type]);
    strcat(name, ext);
    FILE *G = file_open_write(name);
    compress_data(F, num, G, data, tb_size, minfreq, wide);
    fclose(G);
    file_rename(name);
  } else {
    compress_data_single_valued(F, num);
  }
}

void merge_tb(struct tb_handle *F)
{
  char str[64];
  char ext[32];

  size_t name_len = strlen(g_output_dir) + strlen(F->name)
      + strlen(suffix[F->type]) + 2;
  char *name = malloc(name_len);
  if (!name)
    out_of_mem();
  snprintf(name, name_len, "%s/%s%s", g_output_dir, F->name, suffix[F->type]);
  FILE *G = file_open_write(name);

  for (int i = 0; i < F->num_tables; i++) {
    if (F->flags[i] & 0x80) continue;
    sprintf(ext, ".%d", 1 + i);
    strcpy(str, F->name);
    strcat(str, suffix[F->type]);
    strcat(str, ext);
    F->H[i] = file_open_read(str);
  }

  write_final(F, G);

  for (int i = 0; i < F->num_tables; i++) {
    if (F->flags[i] & 0x80) continue;
    fclose(F->H[i]);
  }

  fclose(G);

  char *tmp = malloc(strlen(name) + 5);
  if (!tmp)
    out_of_mem();
  strcat(strcpy(tmp, name), ".tmp");
  add_checksum(tmp);
  file_rename(name);
  free(tmp);
  free(name);

  for (int i = 0; i < F->num_tables; i++) {
    if (F->flags[i] & 0x80) continue;
    sprintf(ext, ".%d", 1 + i);
    strcpy(str, F->name);
    strcat(str, suffix[F->type]);
    strcat(str, ext);
    unlink(str);
  }

  for (int i = 0; i < F->num_tables; i++) {
    if (F->flags[i] & 0x80) continue;
    if (F->flags[i] & 0x40)
      free_code_rans(F->code.rans[i]);
    else
      free_code(F->code.huff[i]);
    free(F->symtable[i]);
  }

  free(F);
}

static constexpr int default_blocksize[3] = { 6, 6, 10 };

void compress_data_slice(const char *name, int stm, int type, void *data,
    uint64_t tb_size, uint8_t *perm, int minfreq, bool wide, bool big)
{
  FILE *F = file_open_write(name);

  if (g_compress_type == 1) {
    write_u8(F, 0xff);
    if (type == WDL) {
      for (int i = 0; i < 5; i++)
        if (wdl_vals[i]) {
          write_u8(F, i);
          break;
        }
    } else {
      write_u8(F, current_map->map[0][0]);
#ifdef HAS_PAWNS
      if (type != WDL)
        write_u8(F, g_dist_format);
#endif
    }
    goto finished;
  }

  void *code;
  int num_syms;
  int64_t *freq = construct_pairs(data, tb_size, minfreq, 0, &num_syms, wide,
      type == WDL);
  bool rans = g_use_rans && type == DTZ;
  if (!rans) {
    struct HuffCode *c = create_code(freq, num_syms);
    if (type == DTZ && huffman_misses_entropy(c, num_syms)) {
      rans = true;
      freq = release_huffman_freq(c);
    } else {
      code = c;
      if (!wide)
        calc_block_sizes_u8(data, tb_size, code, default_blocksize[type]);
      else
        calc_block_sizes_u16(data, tb_size, code, default_blocksize[type]);
    }
  }

  if (rans) {
    used_rans = true;
    code = create_code_rans(freq, num_syms);
    if (!wide)
      calc_block_sizes_rans_u8(data, tb_size, code, default_blocksize[type]);
    else
      calc_block_sizes_rans_u16(data, tb_size, code, default_blocksize[type]);
  }

  bool mapped = false;
  if (type == DTZ) {
    int i, j, k;
    for (i = 0; i < 4; i++)
      if (current_map->num[i] == num_vals)
        break;
    for (j = 0; j < 4; j++)
      if (j != i) {
        for (k = 0; k < current_map->num[j]; k++)
          if (current_map->map[j][k] != current_map->map[i][k])
            break;
        if (k < current_map->num[j])
          break;
      }
    if (j == 4) {
      for (k = 0; k < num_vals; k++)
        symtable[k].pattern[0] = current_map->map[i][k];
    } else {
      mapped = true;
    }
  }

  // Mapped or not.
  if (type != WDL) {
    write_u8(F, !mapped ? 0 : !wide ? 1 : 2);
#ifdef HAS_PAWNS
    write_u8(F, g_dist_format);
#endif
  }

  // Piece permutation.
#ifndef HAS_PAWNS
  if (!big)
    for (int i = 0; i < ii.numsets; i++)
      write_u8(F, perm[i]);
  else
    for (int i = 0; i < ii.numsets + 1; i++)
      write_u8(F, perm[i]);
#else
  if (!big)
    for (int i = 0; i < ii.numsets + 1; i++)
      write_u8(F, perm[i]);
  else
    for (int i = 0; i < ii.numsets + 2; i++)
      write_u8(F, perm[i]);
#endif

  if (ftell(F) & 0x01)
    write_u8(F, 0);

  // Entropy coding method.
  write_u8(F, !rans ? 1 : 2);

  // Each compressed block is 1 << blocksize bytes.
  write_u8(F, blocksize);

  // The index of a position is divided by 1 << idxbits to find a
  // "main entry" in the table's index.
  write_u8(F, idxbits);

  // A small number of non-existing blocks at the end of the index to prevent
  // an out-of-bounds error when accessing the index.
  write_u8(F, num_blocks - real_num_blocks);

  // The number of existing compressed blocks in the table.
  write_u32(F, real_num_blocks);

  if (!rans) {
    // A representation of a canonical Huffman code.
    // For background, see: Alistair Moffat and Andrew Turpin,
    // "On the Implementation of Minimum Redundancy Prefix Codes",
    // IEEE Transactions on Communications, Vol. 45, No. 10, October 1997.
    // https://tinyurl.com/w5r92fw8
    struct HuffCode *c = code;
    write_u8(F, c->max_len);
    write_u8(F, c->min_len);
    for (int i = c->min_len; i <= c->max_len; i++)
      write_u16(F, c->offset[i]);

    // Up to 4095 Re-Pair symbols. Each symbol is either the concatenation
    // of two other symbols (encoded as 2 x 12 bits = 3 bytes) or a single
    // WDL/DTZ value (encoded as a combination of 0xfff and 12 bits).
    write_u16(F, num_syms);
    for (int i = 0; i < num_syms; i++) {
      int k = c->map[i];
      if (symtable[k].len == 1) {
        int s1 = symtable[k].pattern[0];
        write_u8(F, s1 & 0xff);
        write_u8(F, (s1 >> 8) | 0xf0);
        write_u8(F, 0xff);
      } else {
        int s1 = c->inv[symtable[k].pattern[0]];
        int s2 = c->inv[symtable[k].pattern[1]];
        write_u8(F, s1 & 0xff);
        write_u8(F, (s1 >> 8) | ((s2 << 4) & 0xff));
        write_u8(F, s2 >> 4);
      }
    }
  } else {
    // rANS. Write out the frequency table.
    struct RansCode *c = code;
    size_t ft_size = write_freq_table(F, c, num_syms);
    printf("num_syms = %d, freq_table size = %lu\n", num_syms, ft_size);

    // Up to 4095 Re-Pair symbols. Each symbol is either the concatenation
    // of two other symbols (encoded as 2 x 12 bits = 3 bytes) or a single
    // WDL/DTZ value (encoded as a combination of 0xfff and 12 bits).
    for (int i = 0; i < num_syms; i++) {
      int k = c->map[i];
      if (symtable[k].len == 1) {
        int s1 = symtable[k].pattern[0];
        write_u8(F, s1 & 0xff);
        write_u8(F, (s1 >> 8) | 0xf0);
        write_u8(F, 0xff);
      } else {
        int s1 = c->inv[symtable[k].pattern[0]];
        int s2 = c->inv[symtable[k].pattern[1]];
        write_u8(F, s1 & 0xff);
        write_u8(F, (s1 >> 8) | ((s2 << 4) & 0xff));
        write_u8(F, s2 >> 4);
      }
    }
  }

  if (ftell(F) & 0x01)
    write_u8(F, 0);

  if (mapped) {
    struct DtzMap *map = current_map;
    if (!wide) {
      for (int i = 0; i < 4; i++) {
        write_u8(F, map->num[i]);
        for (int j = 0; j < map->num[i]; j++)
          write_u8(F, map->map[i][j]);
      }
      if (ftell(F) & 0x01)
        write_u8(F, 0);
    } else {
      for (int i = 0; i < 4; i++) {
        write_u16(F, map->num[i]);
        for (int j = 0; j < map->num[i]; j++)
          write_u16(F, map->map[i][j]);
      }
    }
  }

  // Write the "main" entries of the index. Each entry consists of a 32-bit
  // value and a 16 bit value packed into 6 bytes. The 32-bit value of the
  // kth entry is the number of the block which contains the value of the
  // position with index I_k = k * (1 << indexbits) + (1 << (indexbits - 1)).
  // The 16-bit value is the offset within that compressed block at which
  // the value with this index is found. Note that each compressed block
  // contains at most 65536 values.
  file_write(indextable, 6 * num_indices, F);
  free(indextable);

  // For each table, write out a list of 16-bit values representing the
  // number of values stored in the compressed blocks. This is the 2nd level
  // index. To find the value of a position:
  // - Encode the position into an index I.
  // - Calculate k for which abs(I - I_k) is minimal.
  // - Use the main index to find the block and offset of the value of the
  //   position with index I_k.
  // - Use the 16-bit values in the 2nd level index, starting from the
  //   16-bit value corresponding to block which contains I_k, to find
  //   the block that contains I and the offset of I within that block.
  //
  //  The "indexbits" value is chosen such that about 8 16-bit values have
  //  to be added or subtracted on average.
  file_write(sizetable, 2 * num_blocks, F);
  free(sizetable);

  while (ftell(F) & 0x3f)
    write_u8(F, 0);

  if (!rans) {
    struct HuffCode *c = code;
    if (!wide)
      write_ctb_data_u8(F, data, c, tb_size, blocksize);
    else
      write_ctb_data_u16(F, data, c, tb_size, blocksize);
    free_code(c);
  } else {
    struct RansCode *c = code;
    if (!wide)
      write_ctb_data_rans_u8(F, data, c, tb_size, blocksize);
    else
      write_ctb_data_rans_u16(F, data, c, tb_size, blocksize);
    free_code_rans(c);
  }
  while (ftell(F) & 0x3f)
    write_u8(F, 0);

finished:
  fclose(F);
  file_rename(name);
}
