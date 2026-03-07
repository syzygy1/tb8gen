#ifndef PERMUTE_H
#define PERMUTE_H

#include <inttypes.h>

extern uint64_t tb_size;

void init_permute_piece(uint8_t *pcs, uint8_t *pt);
void permute_piece_wdl(void *tb_table, uint8_t *pcs, uint8_t *pt, void *table,
   uint8_t *best);
void permute_piece_dtz(void *tb_table, uint8_t *pcs, uint8_t *pt, void *table,
   uint8_t *best, bool wide);

void init_permute_pawn(uint8_t *pcs, uint8_t *pt);
void *init_permute_rank(uint8_t *pcs, int rank, void *tb_table,
    bool wide);
void permute_pawn_dtm(void *tb_table, uint8_t *pcs, uint8_t *pt, void *table,
    uint8_t *best, int rank, void *v, bool wide);

#endif
