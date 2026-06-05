/*
  Copyright (c) 2011-2014, 2017, 2025 Ronald de Man

  This file is distributed under the terms of the GNU GPL, version 2.
*/

#ifndef PROBE_H
#define PROBE_H

#include <inttypes.h>
#include <stdatomic.h>

#define LOOKUP

#include "index.h"
#include "movegen.h"
#include "rans.h"
#include "types.h"
#include "util.h"

enum { WDL = 0, DTM, DTZ };
//enum { WL_BOTH = 0, WL_WTM, WL_BTM, W_ONLY, L_ONLY };
enum {
  LT_PIECE = 0, LT_PAWN_FILE, LT_PAWN_RANK,
  LT_PIECE_K, LT_PIECE_KK,
  LT_PAWN_P, LT_PAWN_PK, LT_PAWN_PP, LT_PAWN_PvP
};

enum {
  TWO_SIDED = 0x01,
  WTM_OR_BTM = 0x02,
  WIN_OR_LOSS = 0x04,
  WIN_ONLY = 0x08,
  WTM_ONLY = 0x10
};

#ifdef LOOKUP
static constexpr int LUBITS = 12;

struct LUEntry {
  uint16_t len;
  uint8_t bits;
  uint8_t cwl;
};

static_assert(sizeof(struct LUEntry) == 4);
#endif

static constexpr int STARTBITS = 8;

struct PairsData {
  const uint8_t *index_table;
  const uint16_t *size_table;
  const uint8_t *data;
  uint8_t *symlen;
  const uint8_t *sympat;
  size_t tb_size; // used for decompressing the whole table
  uint16_t num_syms;
  uint8_t block_size;
  uint8_t idx_bits;
  uint8_t compr_type;
  union {
    // Huffman
    struct {
      const uint16_t *offset;
#ifdef LOOKUP
      struct LUEntry *lookup;
#endif
      uint8_t max_len, min_len;
      uint8_t start[1 << STARTBITS];
    };
    // rANS
    struct RansDecode *rans;
    // constant
    uint8_t const_value[2];
  };
  // for Huffman
  uint64_t base[];
};

struct EncInfo {
  size_t factor[MAX_PIECES];
  uint8_t pieces[MAX_PIECES];
  uint8_t norm[MAX_PIECES];
};

struct DecInfo {
  size_t factor[MAX_PIECES];
  uint32_t fac_iter[MAX_PIECES];
  uint8_t ord_iter[MAX_PIECES];
  uint8_t norm[MAX_PIECES];
  uint8_t order[MAX_PIECES];
};

struct TbTable {
  struct PairsData *precomp;
  struct EncInfo ei;
};

struct WdlTable {
  struct PairsData *precomp;
  struct EncInfo ei;
};

struct DtmTable {
  struct PairsData *precomp;
  struct EncInfo ei;
  const uint16_t *map_dtm;
  uint16_t map_dtm_idx[2];
  uint8_t dtm_flags;
};

struct DtzTable {
  struct PairsData *precomp;
  struct EncInfo ei;
  union {
    const uint8_t *map_dtz;
    const uint16_t *map_dtz16;
  };
  uint16_t map_dtz_idx[4];
  uint8_t dtz_flags;
};

struct TbTable2 {
  struct PairsData *precomp;
  struct RankInfo *ri;
};

struct WdlTable2 {
  struct PairsData *precomp;
  struct RankInfo *ri;
};

struct DtmTable2 {
  struct PairsData *precomp;
  struct RankInfo *ri;
  uint8_t mapped;
  uint8_t dist_format;
  const uint16_t *map_dtm;
  uint16_t map_dtm_idx[2];
};

struct DtzTable2 {
  struct PairsData *precomp;
  struct RankInfo *ri;
  uint8_t mapped;
  uint8_t dist_format;
  uint16_t map_dtz_idx[4];
  union {
    const uint8_t *map_dtz;
    const uint16_t *map_dtz16;
  };
};

struct TbTableConst {
  struct PairsData *precomp;
  uint16_t const_value;
  uint8_t dist_format;
};

struct Tbase {
  const void *data;
  map_t mapping;
  uint8_t pt[MAX_PIECES];
  uint8_t layout;
  uint8_t dist_format;
  uint8_t version;
  bool flipped;
  uint8_t offset;
  void *_Atomic table[];
};

struct TbEntry {
  uint64_t key;
  bool has_dtm, has_dtz, has_pawns, symmetric; // pack in one byte?
  union {
    bool kk_enc;
    uint8_t pawns[2];
  };
  uint8_t num, numsets;
  struct Tbase *_Atomic tbase[3];
};

extern const uint8_t Triangle[64];
extern const uint8_t InvTriangle[10];
extern const int8_t OffDiag[64];
extern const uint8_t PawnMap[48];
extern const uint8_t PawnTwist[2][64];
extern size_t PawnIdx[2][6][24];
extern const uint8_t PawnFlip[2][64];
extern const uint8_t InvPawnFlip[2][24];
extern const uint8_t InvPawnTwist[2][48];
extern int16_t KKPIdx[64][64];

extern const char *suffix[3];
extern uint32_t magic[3], magic2[3];

size_t set_dec_info(struct DecInfo *di, struct TbEntry *e, uint8_t *pcs,
    uint8_t *perm, int order, int order2, int fr, const int enc);

void decode_piece(uint64_t idx, uint8_t *p, struct DecInfo *di,
    struct TbEntry *e);
void decode_pawn_r(uint64_t idx, uint8_t *p, struct DecInfo *di,
    struct TbEntry *e, int rank);

void decode_init(uint32_t *sub, uint64_t idx, struct DecInfo *di);
void decode_piece_iter(uint32_t *sub, uint8_t *p, struct DecInfo *di,
    struct TbEntry *e);
void decode_pawn_r_iter(uint32_t *sub, uint8_t *p, struct DecInfo *di,
    struct TbEntry *e, int rank);

size_t encode_piece(uint8_t *p, struct EncInfo *ei, struct TbEntry *e);
size_t encode_pawn_f(uint8_t *p, struct EncInfo *ei, struct TbEntry *e);
size_t encode_pawn_r(uint8_t *p, struct EncInfo *ei, struct TbEntry *e);

const void *map_tb(const char *name, const int type, map_t *mapping,
    bool use_paths);

void init_tablebases(const char *pathList);
bool init_table(struct TbEntry *e, const char *str, int type, bool use_paths);
void free_tablebases(void);

int probe_capts_wdl(Position *pos, int alpha, int beta);
int probe_capts_dtm(Position *pos, int lower, int upper, const bool won);

int probe_wdl(Position *pos, int alpha, int beta);
int probe_dtm_win(Position *pos, int lower, int upper);
int probe_dtm_loss(Position *pos, int lower, int upper);

#endif
