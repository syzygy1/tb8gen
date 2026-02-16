/*
  Copyright (c) 2011-2014, 2017, 2025 Ronald de Man

  This file is distributed under the terms of the GNU GPL, version 2.
*/

#ifndef PROBE_H
#define PROBE_H

#include <inttypes.h>
#include <stdatomic.h>

#define LOOKUP

#include "movegen.h"
#include "rans.h"
#include "types.h"
#include "util.h"

#if defined(SUICIDE) || defined(GIVEAWAY) || defined(ATOMIC)
#define CONNECTED_KINGS
#endif

enum { WDL = 0, DTM, DTZ };
enum { PIECE_ENC = 0, FILE_ENC, RANK_ENC };

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
  size_t tb_size; // user for decompressing the whole table
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
  struct PairsData *precomp;
  size_t factor[TB_PIECES];
  uint8_t pieces[TB_PIECES];
  uint8_t norm[TB_PIECES];
};

struct DecInfo {
  size_t factor[TB_PIECES];
  uint32_t fac_iter[TB_PIECES];
  uint8_t ord_iter[TB_PIECES];
  uint8_t norm[TB_PIECES];
  uint8_t order[TB_PIECES];
};

struct BaseEntry {
  uint64_t key;
  const uint8_t *data[3];
  map_t mapping[3];
  atomic_bool ready[3];
  uint8_t num;
  bool symmetric, has_pawns, has_dtm, has_dtz;
  union {
    uint8_t kk_enc;
    uint8_t pawns[2];
  };
  bool dtm_losses_only;
};

struct PieceEntry {
  struct BaseEntry be;
  struct EncInfo ei[5]; // 2 + 2 + 1
  uint16_t *map_dtm;
  uint16_t map_dtm_idx[2][2];
  const void *map_dtz;
  uint16_t map_dtz_idx[4];
  uint8_t dtz_flags;
};

struct PawnEntry {
  struct BaseEntry be;
  struct EncInfo ei[24]; // 4 * 2 + 6 * 2 + 4
  uint16_t *map_dtm;
  uint16_t map_dtm_idx[6][2][2];
  const void *map_dtz;
  uint16_t map_dtz_idx[4][4];
  uint8_t dtz_flags[4];
  bool dtm_switched;
};

#ifndef CONNECTED_KINGS
extern const int16_t KKIdx[10][64];
extern uint8_t KKSquare[462][2];
#endif
extern const uint8_t Triangle[64];
extern const uint8_t InvTriangle[10];
extern const int8_t OffDiag[64];
extern const uint8_t FlipDiag[64];
extern size_t Binomial[8][64];
extern const uint8_t PawnMap[48];
extern const uint8_t PawnTwist[2][64];
extern size_t PawnIdx[2][6][24];
extern const uint8_t PawnFlip[2][64];
extern const uint8_t InvPawnFlip[2][24];
extern const uint8_t InvPawnTwist[2][48];

extern const char *suffix[3];
extern uint32_t magic[3];

size_t set_dec_info(struct DecInfo *di, struct BaseEntry *be, uint8_t *pcs,
    uint8_t *perm, int order, int order2, int fr, const int enc);
size_t subfactor(size_t k, size_t n);

void decode_piece(uint64_t idx, uint8_t *p, struct DecInfo *di,
    struct BaseEntry *be);
void decode_pawn_r(uint64_t idx, uint8_t *p, struct DecInfo *di,
    struct BaseEntry *be, int rank);

void decode_init(uint32_t *sub, uint64_t idx, struct DecInfo *di);
void decode_piece_iter(uint32_t *sub, uint8_t *p, struct DecInfo *di,
    struct BaseEntry *be);
void decode_pawn_r_iter(uint32_t *sub, uint8_t *p, struct DecInfo *di,
    struct BaseEntry *be, int rank);

size_t encode_piece(uint8_t *p, struct EncInfo *ei, struct BaseEntry *be);
size_t encode_pawn_f(uint8_t *p, struct EncInfo *ei, struct BaseEntry *be);
size_t encode_pawn_r(uint8_t *p, struct EncInfo *ei, struct BaseEntry *be);

const void *map_tb(const char *name, const int type, map_t *mapping,
    bool use_paths);

void init_tablebases(const char *pathList);
bool init_table(struct BaseEntry *be, const char *str, int type,
    bool use_paths);
void free_tablebases(void);

int probe_capts_wdl(Position *pos, int alpha, int beta);
int probe_capts_dtm(Position *pos, int lower, int upper, const bool won);

int probe_wdl(Position *pos, int alpha, int beta);
int probe_dtm_win(Position *pos, int lower, int upper);
int probe_dtm_loss(Position *pos, int lower, int upper);

#endif
