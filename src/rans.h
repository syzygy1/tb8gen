#ifndef RANS_H
#define RANS_H

#include <inttypes.h>

#include "defs.h"
#include "types.h"
#include "util.h"

// rANS encoding and decoding based on https://github.com/rygorous/ryg_rans.

typedef uint64_t RansState;

static constexpr RansState RANS_L = 1ULL << 31;

struct RansEntry {
  uint16_t alias_div;
  uint16_t start;
};

struct RansDecode {
  struct RansEntry entry[4096];
  uint16_t freq[4096];
};

struct RansCode {
  struct RansDecode d;
  uint16_t cum_freq[4096];
  uint16_t alias_remap[1 << 16];
  uint16_t map[4096];
  uint16_t inv[4096];
};

INLINE void rans_enc_init(RansState *r, bool init_zero)
{
  *r = init_zero ? 0 : RANS_L;
}

INLINE bool rans_enc_put(RansState *r, struct RansCode *c, int s, uint8_t **p,
    uint8_t *const end)
{
  uint16_t freq = c->d.freq[s];
  uint64_t x_max = (uint64_t)freq << (63 - 16);
  uint64_t x = *r;

  // Renormalize the state. Allow bit 63 to be used at the end of the block.
  if (x >= x_max && (*p != end || x >= 2 * x_max)) {
    if (*p == end)
      return false;
    *p -= 4;
    write_le_u32(*p, (uint32_t)x);
    x >>= 32;
  }

  // Encode s into the state.
  *r = ((x / freq) << 16) + c->alias_remap[(x % freq) + c->cum_freq[s]];
  return true;
}

// Count the number of output bytes.
INLINE bool rans_enc_count(RansState *r, struct RansCode *c, int s,
    int *remaining)
{
  uint16_t freq = c->d.freq[s];
  uint64_t x_max = (uint64_t)freq << (63 - 16);
  uint64_t x = *r;

  if (x >= x_max && (*remaining > 0 || x >= 2 * x_max)) {
    if (*remaining == 0)
      return false;
    *remaining -= 4;
    x >>= 32;
  }

  *r = ((x / freq) << 16) + c->alias_remap[(x % freq) + c->cum_freq[s]];
  return true;
}

INLINE void rans_enc_flush(RansState *r, uint8_t **p, uint8_t *const start)
{
  *p -= 8;
  write_le_u64(*p, *r);
  while (*p > start) {
    *p -= 4;
    write_le_u32(*p, 0);
  }
}

INLINE void rans_dec_init(RansState *r, const uint8_t **p,
    const uint8_t *const end)
{
  while (read_le_u32(*p) == 0 && *p + 8 < end)
    *p += 4;
  *r = read_le_u64(*p);
  *p += 8;
}

INLINE uint16_t rans_dec_get(RansState *r, struct RansDecode *d)
{
  uint64_t x = *r;

  // Figure out symbol via alias table.
  int s = (x >> 4) & 4095; // bucket id
#if 0
  // branchless seems equal or a bit slower
  uint32_t remainder = x & 15;
  uint32_t alias_div = d->entry[s].alias_div;
  int alias = alias_div >> 4;
  uint32_t mask = -(uint32_t)(remainder >= (alias_div & 15));
  uint32_t start = d->entry[s].start;
  remainder += start & mask;
  remainder &= 0xffff;
  s ^= (s ^ alias) & mask;
#else
  uint16_t remainder = x & 15;
  if (remainder >= (d->entry[s].alias_div & 15)) {
    // switch to alias symbol
    remainder += d->entry[s].start; // can wrap around by design
    s = d->entry[s].alias_div >> 4;
  }
#endif

  // Update state.
  uint16_t freq = d->freq[s];
  *r = freq * (x >> 16) + remainder;

  return s;
}

INLINE void rans_dec_renorm(RansState *r, const uint8_t **p)
{
  uint64_t x = *r;
  // Branchless renormalization is a win.
  int64_t s = x < RANS_L;
  uint64_t tx = x;
  tx = (tx << 32) | (uint64_t)read_le_u32(*p);
  *p += s << 2;
  x = (-s & (tx - x)) + x;
  *r = x;
}

struct RansCode *create_code_rans(int64_t *freq, int num_syms);
void make_alias_table(struct RansDecode *d, struct RansCode *c);
void free_code_rans(struct RansCode *c);
size_t write_freq_table(FILE *F, struct RansCode *c, int num_syms);
const uint8_t *read_freq_table(struct RansDecode *d, int *num_syms,
    const uint8_t *p);

#endif
