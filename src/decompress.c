/*
  Copyright (c) 2011-2013, 2018, 2025 Ronald de Man

  This file is distributed under the terms of the GNU GPL, version 2.
*/

#include <fcntl.h>
#include <math.h>
#include <stdbit.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "defs.h"
#include "decompress.h"
#include "probe.h"
#include "threads.h"
#include "util.h"

static uint64_t expand_symbol(uint8_t *restrict dst, int sym, uint64_t idx,
    uint64_t end, const uint8_t *sympat, const uint8_t *symlen)
{
  if (idx == end) return idx;
  uint32_t w = read_le_u32(sympat + 3 * sym);
  if (symlen[sym] == 0) {
    dst[idx++] = min(0xff, w & 0x0fff);
    return idx;
  }
  idx = expand_symbol(dst, w & 0x0fff, idx, end, sympat, symlen);
  idx = expand_symbol(dst, (w >> 12) & 0x0fff, idx, end, sympat, symlen);
  return idx;
}

static uint8_t *table_data;
static struct PairsData *precomp;
static int table_const_value;
static uint64_t table_size;
static struct Work work_decompression;

static void decompress_worker(struct ThreadData *thread)
{
  uint64_t idx = thread->begin;
  uint64_t end = thread->end;
  struct PairsData *d = precomp;
  uint8_t *restrict dst = table_data;

  if (!d) {
    int s = table_const_value;
    memset(dst + idx, s, end - idx);
    return;
  }

  int l;
  int m = d->min_len;
  const uint16_t *offset = d->offset;
  const uint64_t *base = d->base - m;
  const uint8_t *symlen = d->symlen;
  const uint8_t *sympat = d->sympat;
  int sym, bitcnt;

  uint32_t mainidx = idx >> d->idx_bits;
  int litidx = (int)(idx & ((1u << d->idx_bits) - 1))
              - (1u << (d->idx_bits - 1));
  uint32_t block = read_le_u32(d->index_table + 6 * mainidx);
  litidx += read_le_u16(d->index_table + 6 * mainidx + 4);

  if (litidx < 0) {
    while (litidx < 0)
      litidx += d->size_table[--block] + 1;
  } else {
    mainidx++;
    while (litidx > d->size_table[block])
      litidx -= d->size_table[block++] + 1;
  }

  uint64_t idx2 =  (1ULL << (d->idx_bits - 1))
                 + ((uint64_t)mainidx << d->idx_bits);
  if (litidx > 0) {
    idx += d->size_table[block++] + 1 - litidx;
    while (idx >= idx2) {
      idx2 += 1ULL << d->idx_bits;
      mainidx++;
    }
  }

  const uint8_t *data = d->data + ((uint64_t)block << d->block_size);
  while (idx < end) {
    int size = d->size_table[block] + 1;
    while (idx + size > idx2) {
      if (   read_le_u32(d->index_table + 6 * mainidx) != block
          || read_le_u16(d->index_table + 6 * mainidx + 4) != idx2 - idx)
      {
        fprintf(stderr, "ERROR in main index!!\n");
        exit(EXIT_FAILURE);
      }
      idx2 += 1ULL << d->idx_bits;
      mainidx++;
    }
    block++;

    uint64_t blockend = idx + size;
    if (blockend > table_size) blockend = table_size;

    if (d->compr_type == 0) { // Huffman
      uint32_t *ptr = (uint32_t *)data;
      uint64_t code = from_be_u64(*(uint64_t *)ptr);
      ptr += 2;
      bitcnt = 0;
      while (idx < blockend) {
        l = m;
        while (code < base[l]) l++;
        sym = from_le_u16(offset[l]) + ((code - base[l]) >> (64 - l));
        idx = expand_symbol(dst, sym, idx, blockend, sympat, symlen);
        code <<= l;
        bitcnt += l;
        if (bitcnt >= 32) {
          bitcnt -= 32;
          code |= (uint64_t)from_be_u32(*ptr++) << bitcnt;
        }
      }
      data += 1 << d->block_size;
    } else { // rANS
      const uint8_t *p = data;
      data += 1 << d->block_size;
      RansState rans;
      rans_dec_init(&rans, &p, data);
      uint64_t blockstart = idx;
      idx = blockend;
      while (p < data) {
        int sym = rans_dec_get(&rans, d->rans);
        rans_dec_renorm(&rans, &p);
        idx -= d->symlen[sym] + 1;
        expand_symbol(dst, sym, idx, blockend, sympat, symlen);
      }
      while (idx > blockstart) {
        int sym = rans_dec_get(&rans, d->rans);
        idx -= d->symlen[sym] + 1;
        expand_symbol(dst, sym, idx, blockend, sympat, symlen);
      }
      idx = blockend;
    }
  }
}

void decompress_table(struct TbTable2 *table, uint8_t *decomp_table,
   size_t size)
{
  table_data = decomp_table;
  table_size = size;
  precomp = table->precomp;
  if (!precomp) {
    struct TbTableConst *tbl = (struct TbTableConst *)table;
    table_const_value = tbl->const_value;
  }
  work_init(&work_decompression, size, 0x1ff, WORK_DYNAMIC, 4, 1ULL << 9);
  run_threaded(decompress_worker, &work_decompression);
}

#if 0
void print_decompression_info(struct EncInfo *ei)
{
  struct PairsData *d = ei->precomp;
  if (d->compr_type == 2)
    return;
  if (d->compr_type == 0)
    printf("min_len = %d; max_len = %d.\n", d->min_len, d->max_len);
  printf("block_size = %d bytes.\n", 1 << d->block_size);
  uint64_t size = d->tb_size;
  uint64_t *freq = calloc(4096, sizeof(uint64_t));
  uint64_t used_bits = 0, block = 0, symcnt = 0;
  if (d->compr_type == 0) { /* huffman */
    uint64_t *base = d->base - d->min_len;
    for (uint64_t idx = 0; idx < size; block++) {
      const uint8_t *ptr = d->data + (block << d->block_size);
      uint64_t bitbuf = read_be_u64(ptr), pending = bitbuf;
      ptr += 8;
      int sym, bitcnt = 0;
      uint64_t blockend = idx + d->size_table[block] + 1;
      // We use the "zero refill latency" trick from
      // https://dougallj.wordpress.com/2022/08/26/reading-bits-with-zero-refill-latency/
      while (idx < blockend) {
        int l = d->start[bitbuf >> (64 - STARTBITS)];
        while (bitbuf < base[l]) l++;
        sym = from_le_u16(d->offset[l]) + ((bitbuf - base[l]) >> (64 - l));
        idx += d->symlen[sym] + 1;
        freq[sym]++;
        symcnt++;
        bitbuf = pending << l;
        bitcnt += l;
        used_bits += l;
        pending = bitbuf | read_be_u64(ptr) >> (64 - bitcnt);
        ptr += (bitcnt >> 5) * sizeof(uint32_t);
        bitcnt &= 31;
      }
      idx = blockend;
    }
  } else { /* rANS */
    for (uint64_t idx = 0; idx < size; block++) {
      const uint8_t *p = d->data + (block << d->block_size);
      const uint8_t *end = p + (1 << d->block_size);
      uint64_t blockend = idx + d->size_table[block] + 1;
      RansState rans;
      rans_dec_init(&rans, &p, end);
      used_bits += 8 * (end - p);
      if (rans) used_bits += 64 - stdc_leading_zeros(rans);
      else used_bits++;
      while (p < end) {
        int sym = rans_dec_get(&rans, d->rans);
        rans_dec_renorm(&rans, &p);
        freq[sym]++;
        symcnt++;
        idx += d->symlen[sym] + 1;
      }
      while (idx < blockend) {
        int sym = rans_dec_get(&rans, d->rans);
        freq[sym]++;
        symcnt++;
        idx += d->symlen[sym] + 1;
      }
    }
  }
  printf("total blocks = %lu\n", block);
  uint64_t bits_per_block = block << (3 + d->block_size);
  printf("fill rate = %.2lf%%, unused bits/block = %.2lf\n",
      (double)used_bits/bits_per_block * 100.0,
      (double)(8 << d->block_size) - (double)used_bits / block);
  printf("num_syms = %u, sym_cnt = %lu, positions/symbol = %.3lf\n",
      d->num_syms, symcnt, (double)size/symcnt);
  int total = 0;
  for (int i = 0; i < d->num_syms; i++)
    total += d->symlen[i] + 1;
  printf("total length dictionary = %d bytes\n", total);
  double entropy = 0.0;
  for (int i = 0; i < d->num_syms; i++) {
    if (freq[i] == 0) continue;
    double p = (double)freq[i] / symcnt;
    entropy += -p * log(p);
  }
  entropy /= log(2.0);
  double actual = (double)bits_per_block / symcnt;
  printf("entropy = %.3lf bits/symbol\n", entropy);
  if (d->compr_type == 0) {
    double huff_entropy = (double)used_bits / symcnt;
    printf("huffman code = %.3lf bits/symbol (%.3lf%%)\n", huff_entropy,
        entropy / huff_entropy * 100.0);
  }
  printf("compressed blocks = %.3lf bits/symbol (%.3lf%%)\n", actual,
      entropy / actual * 100.0);
  double per_block = entropy * symcnt / (8 * block);
  printf("entropy per block = %.3lf bytes, %.3lf bits below maximum\n",
      per_block, (double)(8 << d->block_size) - 8.0 * per_block);
}

void print_freq_table(struct EncInfo *ei)
{
  struct PairsData *d = ei->precomp;
  if (d->compr_type != 1) return;

  uint32_t R = 1 << 16;
  uint32_t P = 1 << 16;
  int num = d->num_syms;
  for (int s = 0; s < num; s++) {
    uint32_t v_max = min(R, P);
    uint32_t v_min = (R + num - s - 1) / (num - s);
    uint32_t f = d->rans->freq[s];
    printf("%d: v_min = %u, v_max = %u, f = %u (%u), R = %u/%u\n",
        s, v_min, v_max, f, v_max - f, R, num - s - 1);
    R -= f;
    P = f;
  }
}
#endif
