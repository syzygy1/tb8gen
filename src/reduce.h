#ifndef REDUCE_H
#define REDUCE_H

#include "defs.h"
#include "threads.h"

struct MergeInfo {
  union {
    uint8_t v_u8[MAX_STATS];
    uint16_t v_u16[MAX_STATS];
  };
  union {
    uint16_t v_inv_u8[256];
    uint16_t v_inv_u16[MAX_STATS];
  };
  bool wide;
};

extern struct MergeInfo mi[2];

extern int epoch;
extern int reduce_cnt_win[64];
extern int reduce_cnt_loss[64];

enum {
  RAM_UNRESOLVED = 0,
  RAM_LOSS_IN_0 = 1,
  RAM_CAPT_LOSS = 2,
  RAM_CAPT_BLOSS = 2 + DRAW_RULE,
  RAM_PAWN_DRAW = 126,
  RAM_CAPT_DRAW = 127,
  RAM_CAPT_CWIN = 253 - DRAW_RULE - 1,
  RAM_PAWN_CWIN = 253 - DRAW_RULE - 2,
  RAM_PAWN_WIN = 253,
  RAM_CAPT_WIN = 254,
  RAM_ILLEGAL = 255,

  RAM_REDUCED_LOSS = 1,
  RAM_REDUCED_CAPT_BLOSS = 2,
  RAM_REDUCED_BLOSS = 3,
  RAM_REDUCED_CWIN = 251,
  RAM_REDUCED_CAPT_CWIN = 252,
  RAM_REDUCED_WIN = 253,
};

INLINE uint8_t win_to_byte(int n)
{
  if (epoch == 0)
    return 253 - n - 2 * (n > DRAW_RULE);
  else
    return reduce_cnt_win[epoch - 1] - n;
}

INLINE uint8_t loss_to_byte(int n)
{
  if (epoch == 0)
    return 1 + n + (n > DRAW_RULE);
  else
    return n - reduce_cnt_loss[epoch - 1];
}

INLINE int win_to_stat(int n)
{
  return 2 + n + 2 * (n > DRAW_RULE);
}

INLINE int loss_to_stat(int n)
{
  return MAX_STATS - 1 - n;
}

void reduce_tables(int n);
void unlink_table(int stm);
void unlink_saves(int stm);
void store_table(uint8_t *table, int stm);

void transform_table_u8(struct ThreadData *thread);
void transform_table_u16(struct ThreadData *thread);
void load_table_u8(uint8_t *table, int stm);
void load_table_u16(uint16_t *table, int stm);
//void reconstruct_table_u8(uint8_t *table, int stm, struct DtzMap *map);
//void reconstruct_table_u16(uint16_t *table, int stm, struct DtzMap *map);

#endif
