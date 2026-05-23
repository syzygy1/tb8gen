/*
  Copyright (c) 2011-2017, 2025 Ronald de Man

  This file is distributed under the terms of the GNU GPL, version 2.
*/

#define NAME(f) EVALUATOR(f,T)

static int64_t *NAME(setup_code)(T *restrict data, uint64_t size)
{
  int s;

  int64_t *freq = calloc(MAXSYMB, sizeof(uint64_t));

  for (uint64_t idx = 0; idx < size; idx += symtable[s].len) {
    s = read_symbol(data[idx]);
    freq[s]++;
  }

  return freq;
}

// same as for dtz
void NAME(adjust_work_dontcares)(uint64_t *restrict work1,
    uint64_t *restrict work2)
{
  size_t idx;
  size_t end = work1[compress_work_units];
  T *restrict data = compress_state.data;
  int i;
  int num = num_vals;

  work2[0] = work1[0];
  for (i = 1; i < compress_work_units; i++) {
    idx = work1[i];
    if (idx < work2[i - 1]) {
      work2[i] = work2[i - 1];
      continue;
    }
    while (idx < end && data[idx] == num)
      idx++;
    work2[i] = idx;
  }
  work2[compress_work_units] = work1[compress_work_units];
}

// same as for dtz
void NAME(fill_dontcares)(struct ThreadData *thread)
{
  size_t idx, idx2;
  size_t end = thread->end;
  T *restrict data = compress_state.data;
  size_t size = compress_state.size;
  int s1;
  int k;
  int num = num_vals;

  idx = thread->begin;
  while (idx < end) {
    // find start of don't care sequence
    for (; idx < end; idx++)
      if (data[idx] == num)
	break;
    if (idx == end) break;
    // find end of sequence and determine max pair frequency
    for (idx2 = idx + 1; idx2 < end; idx2++)
      if (data[idx2] != num)
	break;
    int64_t max_val;
    if (idx2 - idx > 1)
      max_val = dcfreq[0][0];
    else {
      max_val = -1;
      if (idx > 0)
	max_val = pairfreq[data[idx - 1]][pairfirst[0][data[idx - 1]]];
      if (idx2 < size)
        max_val = max(max_val, pairfreq[pairsecond[0][data[idx2]]][data[idx2]]);
    }
    // now replace whenever we reach max pair frequency for the sequence
    int dc = 1;
    if (idx > 0) {
      s1 = data[idx - 1];
      k = pairfirst[0][s1];
      if (pairfreq[s1][k] == max_val) {
	data[idx] = k;
	s1 = k;
	dc = 0;
      }
    }
    idx2 = idx + 1;
    while (idx2 < end && data[idx2] == num) {
      if (dc) {
	data[idx2 - 1] = t1[0][0];
	s1 = t2[0][0];
	data[idx2] = s1;
	dc = 0;
      } else {
	k = pairfirst[0][s1];
	dc = 1;
	if (pairfreq[s1][k] == max_val) {
	  data[idx2] = k;
	  s1 = k;
	  dc = 0;
	}
      }
      idx2++;
    }
    if (dc && idx2 < size) {
      s1 = data[idx2];
      k = pairsecond[0][s1];
      if (pairfreq[k][s1] == max_val)
	data[idx2 - 1] = k;
    }
  }
}

static void NAME(count_pairs_dtz)(struct ThreadData *thread)
{
  int s1, s2;
  size_t idx = thread->begin;
  size_t end = thread->end;
  T *restrict data = compress_state.data;
  int t = thread->thread_id;
  NAME(cf) *countfreq = countfreq_dtz;

  if (idx == 0) idx = 1;
  if (idx >= end) return;
  s1 = data[idx - 1];
  for (; idx < end; idx++) {
    s2 = data[idx];
    countfreq[t][s1][s2]++;
    s1 = s2;
  }
}

void NAME(adjust_work_replace)(uint64_t *work)
{
  size_t idx, idx2;
  size_t end = compress_state.size;
  T *restrict data = compress_state.data;
  int i, s1, s2, j;

  for (i = 1; i < compress_work_units; i++) {
    idx = work[i];
    if (idx <= work[i - 1]) {
      work[i] = work[i - 1];
      continue;
    }
    if (idx == end) continue;
    s1 = read_symbol(data[idx]);
    idx += symtable[s1].len;
    j = 0;
    while (idx < end) {
      s2 = read_symbol(data[idx]);
      if (newtest[s1][s2]) j = 0;
      else {
	if (j == 1) break;
	j = 1;
	idx2 = idx;
      }
      idx += symtable[s2].len;
      s1 = s2;
    }
    if (idx < end)
      work[i] = idx2;
    else
      work[i] = end;
  }
}

void NAME(replace_pairs)(struct ThreadData *thread)
{
  uint64_t idx = thread->begin;
  uint64_t end = thread->end;
  T *restrict data = compress_state.data;
  int s1, s2, a;
  int t = thread->thread_id;

  if (idx == end) return;

  a = -1;
  s1 = read_symbol(data[idx]);
  idx += symtable[s1].len;
  while (idx < end) {
    s2 = read_symbol(data[idx]);
    idx += symtable[s2].len;
    if (newtest[s1][s2]) {
      struct Symbol *sym = &symtable[newpairs[newtest[s1][s2] - 1].sym];
      write_symbol(data[idx - sym->len], sym);
      if (likely(a >= 0))
        countfirst[t][newtest[s1][s2] - 1][a]++;
      a = newtest[s1][s2] - 1;
      if (unlikely(idx == compress_state.size)) break;
      s1 = read_symbol(data[idx]);
      idx += symtable[s1].len;
      countsecond[t][a][s1]++;
      a = newpairs[a].sym;
    } else {
      a = s1;
      s1 = s2;
    }
  }
}

// might not be a win
static void NAME(remove_dtz_worker)(struct ThreadData *thread)
{
  uint64_t idx, idx2;
  uint64_t end = thread->end;
  T *restrict data = compress_state.data;
  uint64_t size = compress_state.size;
  int s;
  int num = num_vals;
  int max = current_map->high_freq_max;

  idx = thread->begin;
  while (idx < end) {
    for (; idx < end; idx++)
      if (data[idx] == num) break;
    if (idx == end) break;
    for (idx2 = idx + 1; idx2 < end; idx2++)
      if (data[idx2] != num) break;
    if (idx2 - idx >= 32) {
      idx = idx2;
      continue;
    }
    if (idx == 0)
      s = data[idx2];
    else {
      s = data[idx - 1];
      if (idx2 < size && data[idx2] < s)
	s = data[idx2];
    }
    if (s >= max)
      idx = idx2;
    else
      for (; idx < idx2; idx++)
	data[idx] = s;
  }
}

void NAME(adjust_work_dontcares_dtz)(uint64_t *restrict work1,
    uint64_t *restrict work2)
{
  uint64_t idx;
  uint64_t end = work1[compress_work_units];
  T *data = compress_state.data;
  int i;
  int num = num_vals;

  work2[0] = work1[0];
  for (i = 1; i < compress_work_units; i++) {
    idx = work1[i];
    if (idx < work2[i - 1]) {
      work2[i] = work2[i - 1];
      continue;
    }
    while (idx < end && (data[idx - 1] == num || data[idx] == num
					|| data[idx - 1] == data[idx]))
//    while (idx < end && data[idx] == num)
      idx++;
    work2[i] = idx;
  }
  work2[compress_work_units] = work1[compress_work_units];
}

void NAME(calc_block_sizes)(T *restrict data, uint64_t size,
    struct HuffCode *c, int maxsize)
{
  uint64_t idx;
  int i, s, t;
  int64_t block;
  int maxbits, bits, numpos;
  uint32_t avg;

  uint64_t rawsize = calc_size(c);
  printf("calc_size: %"PRIu64"\n", rawsize);

  uint64_t optsize, compsize;

  block = 0;
  compsize = INT64_MAX;
  blocksize = maxsize + 1;
  do {
    optsize = compsize;
    num_blocks = block;
    blocksize--;
    if (((rawsize * ((1 << blocksize) + 2)) >> blocksize) >= optsize) break;
    maxbits = 8 << blocksize;
    bits = 0;
    numpos = 0;
    block = 0;
    for (idx = 0; idx < size;) {
      s = read_symbol(data[idx]);
      t = symtable[s].len;
      if (bits + c->length[s] > maxbits || numpos + t > 65536) {
        block++;
        bits = 0;
        numpos = 0;
      }
      bits += c->length[s];
      numpos += t;
      idx += t;
    }
    if (numpos > 0)
      block++;
    compsize = block << blocksize;
    compsize = (compsize + 0x3f) & ~0x3f;
    avg = size / block;
    idxbits = 0;
    while (avg) {
      idxbits++;
      avg >>= 1;
    }
    idxbits += 4;
    while (idxbits > 1 && (1ULL << (idxbits - 1)) > size)
      idxbits--;

    num_indices = (size + (1ULL << idxbits) - 1) >> idxbits;
    t = ((2 * num_indices - 1) << (idxbits - 1)) - size;
    if (t > 0) block += (t + 65535) >> 16;
    else t = 0;
    compsize += 2 * block + 6 * num_indices;

    printf("bits = %d; blocks = %"PRIi64" (%d); size = %"PRIu64"\n",
        blocksize, block - ((t + 65535) >> 16), (t + 65535) >> 16, compsize);

  } while (compsize < optsize);

  blocksize++;
  maxbits = 8ULL << blocksize;
  sizetable = malloc((num_blocks + 16) * sizeof(uint16_t));

  calc_symbol_tweaks(c);

  bits = 0;
  numpos = 0;
  block = 0;
  for (idx = 0; idx < size;) {
    s = read_symbol(data[idx]);
    t = symtable[s].len;
    if (bits + c->length[s] > maxbits || numpos + t > 65536) {
      if (numpos + t <= 65536) {
        if (bits + c->length[replace[s]] <= maxbits) {
          idx += t;
          sizetable[block++] = numpos + t - 1;
          bits = 0;
          numpos = 0;
          continue;
        }
        if (t > 1) {
          int s1 = symtable[s].pattern[0];
          int s2 = symtable[s].pattern[1];
          if (c->length[s1] != 0 && c->length[s2] != 0) {
            if (bits + c->length[s1] + c->length[replace[s2]] <= maxbits) {
              idx += t;
              sizetable[block++] = numpos + t - 1;
              bits = 0;
              numpos = 0;
              continue;
            }
            if (   c->length[s2] < c->length[s]
                && bits + c->length[replace[s1]] <= maxbits)
            {
              sizetable[block++] = numpos + symtable[s1].len - 1;
              idx += t;
              bits = c->length[s2];
              numpos = symtable[s2].len;
              continue;
            }
          }
        }
      }
      if (numpos + t <= 65536 && bits + c->length[replace[s]] <= maxbits) {
        idx += t;
        sizetable[block++] = numpos + t - 1;
        bits = 0;
        numpos = 0;
        continue;
      }

      sizetable[block++] = numpos - 1;
      bits = 0;
      numpos = 0;
    }
    bits += c->length[s];
    numpos += t;
    idx += t;
  }
  if (numpos > 0)
    sizetable[block++] = numpos - 1;

  real_num_blocks = block;
  num_blocks = block;

  avg = size / num_blocks;
  idxbits = 0;
  while (avg) {
    idxbits++;
    avg >>= 1;
  }
  idxbits += 4;
  while (idxbits > 1 && (1ULL << (idxbits - 1)) > size)
    idxbits--;

  num_indices = (size + (1ULL << idxbits) - 1) >> idxbits;
  indextable = malloc(num_indices * 6);

  uint64_t idx2 = 1ULL << (idxbits-1);
  block = 0;
  idx = 0;
  for (i = 0; i < num_indices; i++) {
    while (block < num_blocks && idx + sizetable[block] < idx2)
      idx += sizetable[block++] + 1;
    if (block == num_blocks) {
      sizetable[num_blocks++] = 65535;
      while (idx + sizetable[block] < idx2) {
        idx += sizetable[block++] + 1;
        sizetable[num_blocks++] = 65535;
      }
    }
    uint16_t diff = idx2 - idx;
    memcpy(indextable + 6 * i, &block, 4);
    memcpy(indextable + 6 * i + 4, &diff, 2);
    idx2 += 1ULL << idxbits;
  }

  printf("real_num_blocks = %"PRIu64"; num_blocks = %"PRIu64"\n", real_num_blocks, num_blocks);
  printf("idxbits = %d\n", idxbits);
  printf("num_indices = %u\n", num_indices);
}

static void NAME(write_ctb_data)(FILE *F, T *restrict data, struct HuffCode *c,
    uint64_t size, int blocksize)
{
  uint64_t idx;
  int s, t, l;
  unsigned bits, numpos;
  unsigned maxbits;

  struct BWState bw;
  init_bits(&bw, F, out_buffer, BUFFER_SIZE);

  maxbits = 8 << blocksize;
  bits = numpos = 0;
  for (idx = 0; idx < size;) {
    s = read_symbol(data[idx]);
    t = symtable[s].len;
    // A block of compressed data is at most maxbits long and represents
    // at most 65536 data values.
    if (bits + c->length[s] > maxbits || numpos + t > 65536) {
      if (numpos + t <= 65536) {
        // The regular Huffman symbol for s does not fit, but we can still
        // try to squeeze in:
        // (1) the shortest symbol representing a string of data values
        //     whose initial part coincides with our data, or
        // (2) two symbols s1 and s2' where s was formed by pairing s1 and
        //     s2, and s2' is the shortest symbol for s2 as in (1), or
        // (3) the shortest symbol s1' for s1, and leave s2 for the next
        //     block if it is shorter than s.
        if (bits + c->length[replace[s]] <= maxbits) {
          s = replace[s];
          l = c->length[s];
          put_bits(&bw, c->base[l] + (c->inv[s] - c->offset[l]), l);
          bits += l;
          idx += t;
          pad_bits(&bw, maxbits - bits);
          bits = 0;
          numpos = 0;
          continue;
        }
        if (t > 1) {
          int s1 = symtable[s].pattern[0];
          int s2 = symtable[s].pattern[1];
          if (c->length[s1] != 0 && c->length[s2] != 0) {
            if (bits + c->length[s1] + c->length[replace[s2]] <= maxbits) {
              l = c->length[s1];
              put_bits(&bw, c->base[l] + (c->inv[s1] - c->offset[l]), l);
              bits += l;
              l = c->length[replace[s2]];
              put_bits(&bw, c->base[l] + (c->inv[replace[s2]] - c->offset[l]),
                  l);
              bits += l;
              idx += t;
              pad_bits(&bw, maxbits - bits);
              bits = 0;
              numpos = 0;
              continue;
            }
            if (   c->length[s2] < c->length[s]
                && bits + c->length[replace[s1]] <= maxbits)
            {
              if (bits + c->length[s1] > maxbits) s1 = replace[s1];
              l = c->length[s1];
              put_bits(&bw, c->base[l] + (c->inv[s1] - c->offset[l]), l);
              bits += l;
              pad_bits(&bw, maxbits - bits);
              l = c->length[s2];
              put_bits(&bw, c->base[l] + (c->inv[s2] - c->offset[l]), l);
              idx += t;
              bits = l;
              numpos = symtable[s2].len;
              continue;
            }
          }
        }
      }

      pad_bits(&bw, maxbits - bits);
      bits = 0;
      numpos = 0;
    }
    l = c->length[s];
    put_bits(&bw, c->base[l] + (c->inv[s] - c->offset[l]), l);
    bits += l;
    numpos += t;
    idx += t;
  }
  if (numpos > 0)
    pad_bits(&bw, maxbits - bits);
  flush_bits(&bw);
}

uint64_t NAME(calc_size_rans)(T *restrict data, uint64_t size,
    struct RansCode *c)
{
  uint64_t csize = 0;
  RansState rans;
  rans_enc_init(&rans, true);
  for (uint64_t idx = 0; idx < size;) {
    int s = read_symbol(data[idx]);
    idx += symtable[s].len;
    int remaining = 16;
    rans_enc_count(&rans, c, c->inv[s], &remaining);
    csize += 16 - remaining;
  }
  csize += sizeof(RansState);
  return csize;
}

void NAME(calc_block_sizes_rans)(T *restrict data, uint64_t size,
    struct RansCode *c, int maxsize)
{
  uint64_t idx;
  int i, s, t;
  int numpos, remaining;
  int64_t block;
  uint32_t avg;

  uint64_t rawsize = NAME(calc_size_rans)(data, size, c);
  printf("calc_size: %"PRIu64"\n", rawsize);

  uint64_t optsize, compsize;

  RansState rans;

  block = 0;
  compsize = INT64_MAX;
  blocksize = maxsize + 1;
  do {
    optsize = compsize;
    num_blocks = block;
    blocksize--;
    if (blocksize < 3) break;
    if (((rawsize * ((1 << blocksize) + 2)) >> blocksize) >= optsize) break;
    block = 1;
    numpos = 0;
    remaining = (1 << blocksize) - sizeof(RansState);
    rans_enc_init(&rans, true);
    for (idx = 0; idx < size;) {
      s = read_symbol(data[idx]);
      t = symtable[s].len;
      if (   numpos + t > 65536
          || !rans_enc_count(&rans, c, c->inv[s], &remaining))
      {
        block++;
        numpos = 0;
        remaining = (1 << blocksize) - sizeof(RansState);
        rans_enc_init(&rans, true);
        continue;
      }
      numpos += t;
      idx += t;
    }
    compsize = block << blocksize;
    compsize = (compsize + 0x3f) & ~0x3f;
    avg = size / block;
    idxbits = 0;
    while (avg) {
      idxbits++;
      avg >>= 1;
    }
    idxbits += 4;
    while (idxbits > 1 && (1ULL << (idxbits - 1)) > size)
      idxbits--;

    num_indices = (size + (1ULL << idxbits) - 1) >> idxbits;
    t = ((2 * num_indices - 1) << (idxbits - 1)) - size;
    if (t > 0) block += (t + 65535) >> 16;
    else t = 0;
    compsize += 2 * block + 6 * num_indices;

    printf("bits = %d; blocks = %"PRIi64" (%d); size = %"PRIu64"\n",
        blocksize, block - ((t + 65535) >> 16), (t + 65535) >> 16, compsize);

  } while (compsize < optsize);

  blocksize++;
  sizetable = malloc((num_blocks + 16) * sizeof(uint16_t));

  // FIXME: try tweaks, and try to use the last bit.
//  calc_symbol_tweaks(c);

  block = 0;
  numpos = 0;
  remaining = (1 << blocksize) - sizeof(RansState);
  rans_enc_init(&rans, true);
  for (idx = 0; idx < size;) {
    s = read_symbol(data[idx]);
    t = symtable[s].len;
    if (numpos + t > 65536 || !rans_enc_count(&rans, c, c->inv[s], &remaining))
    {
      sizetable[block++] = numpos - 1;
      numpos = 0;
      remaining = (1 << blocksize) - sizeof(RansState);
      rans_enc_init(&rans, true);
      continue;
    }
    numpos += t;
    idx += t;
  }
  if (numpos > 0)
    sizetable[block++] = numpos - 1;

  real_num_blocks = block;
  num_blocks = block;

  avg = size / num_blocks;
  idxbits = 0;
  while (avg) {
    idxbits++;
    avg >>= 1;
  }
  idxbits += 4;
  while (idxbits > 1 && (1ULL << (idxbits - 1)) > size)
    idxbits--;

  num_indices = (size + (1ULL << idxbits) - 1) >> idxbits;
  indextable = malloc(num_indices * 6);

  uint64_t idx2 = 1ULL << (idxbits-1);
  block = 0;
  idx = 0;
  for (i = 0; i < num_indices; i++) {
    while (block < num_blocks && idx + sizetable[block] < idx2)
      idx += sizetable[block++] + 1;
    if (block == num_blocks) {
      sizetable[num_blocks++] = 65535;
      while (idx + sizetable[block] < idx2) {
        idx += sizetable[block++] + 1;
        sizetable[num_blocks++] = 65535;
      }
    }
    uint16_t diff = idx2 - idx;
    memcpy(indextable + 6 * i, &block, 4);
    memcpy(indextable + 6 * i + 4, &diff, 2);
    idx2 += 1ULL << idxbits;
  }

  printf("real_num_blocks = %"PRIu64"; num_blocks = %"PRIu64"\n", real_num_blocks, num_blocks);
  printf("idxbits = %d\n", idxbits);
  printf("num_indices = %u\n", num_indices);
}

static void NAME(write_ctb_data_rans)(FILE *F, T *restrict data,
    struct RansCode *c, uint64_t size, int blocksize)
{
  RansState rans;

  uint8_t *block = out_buffer, *p = block + (1 << blocksize);
  int numblocks = BUFFER_SIZE >> blocksize;
  int numpos = 0;
  rans_enc_init(&rans, true);
  for (uint64_t idx = 0; idx < size;) {
    int s = read_symbol(data[idx]);
    int t = symtable[s].len;
    if (numpos + t > 65536 || !rans_enc_put(&rans, c, c->inv[s], &p, block + 8))
    {
      rans_enc_flush(&rans, &p, block);
      block += 1 << blocksize;
      numblocks--;
      if (numblocks == 0) {
        file_write(out_buffer, BUFFER_SIZE, F);
        numblocks = BUFFER_SIZE / (1 << blocksize);
        block = out_buffer;
      }
      p = block + (1 << blocksize);
      numpos = 0;
      rans_enc_init(&rans, true);
      continue;
    }
    numpos += t;
    idx += t;
  }
  rans_enc_flush(&rans, &p, block);
  numblocks--;
  file_write(out_buffer, BUFFER_SIZE - numblocks * (1 << blocksize), F);
}

#undef NAME
