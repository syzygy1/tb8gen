#ifndef COMPRESS_H
#define COMPRESS_H

#include <inttypes.h>
#include "types.h"

extern int g_compress_type;

struct HuffCode;
struct RansCode;

static constexpr int MAXSYMB = 4095 + 8;

static constexpr int MAX_VAL = (MAX_STATS / 2 - DRAW_RULE) / 2;

struct DtzMap {
  uint16_t map[4][MAX_VAL];
  uint16_t inv_map[4][MAX_VAL];
//  uint16_t dc_map[2][MAX_VAL];
  uint16_t num[4];
  uint16_t max_num;
  uint16_t high_freq_max;
  uint8_t stm;
  bool wide;
};

struct tb_handle {
  char name[64];
  int num_tables;
  int type;
  bool split;
  uint8_t perm[12][MAX_PIECES];
  int default_blocksize;
  int blocksize[12];
  int idxbits[12];
  uint64_t real_num_blocks[12];
  uint64_t num_blocks[12];
  uint32_t num_indices[12];
  int num_syms[12];
  union {
    struct HuffCode *huff[12];
    struct RansCode *rans[12];
  } code;
  struct Symbol *symtable[12];
  struct DtzMap *map[12];
  int flags[12];
  uint8_t single_val[12];
  FILE *H[12];
};

struct tb_handle *create_tb(const char *tablename, int type, int blocksize);
void compress_tb(struct tb_handle *F, int num, void *data, uint64_t tb_size,
    uint8_t *perm, int minfreq, bool wide);
void merge_tb(struct tb_handle *F);
void compress_init_wdl(bool vals[]);
void compress_init_dtz(struct DtzMap *map, bool wide);
void compress_alloc_wdl(void);
void compress_free_wdl(void);
void compress_alloc_dtz(bool wide);
void compress_free_dtz(void);
int64_t *construct_pairs_wdl(void *data, uint64_t size, int minfreq,
    int maxsymbols, int *nsyms);
int64_t *construct_pairs_dtz(void *data, uint64_t size, int minfreq,
    int maxsymbols, int *nsyms, bool wide);
int64_t *construct_pairs(void *data, uint64_t size, int minfreq,
    int maxsymbols, int *nsyms, bool wide, bool wdl);
void compress_data_462(int s, int stm, int type, void *data, uint64_t tb_size,
    uint8_t *perm, int minfreq, bool wide);

#endif
