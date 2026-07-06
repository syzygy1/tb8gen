#ifndef REDUCE_H
#define REDUCE_H

#include "threads.h"

extern int epoch;
extern int reduce_cnt_win[64];
extern int reduce_cnt_loss[64];

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
    return 1 + n;
  else
    return n - reduce_cnt_loss[epoch - 1];
}

void reduce_tables(int n, int epoch);
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
