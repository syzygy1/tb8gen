#ifdef BMI2_FANCY

Bitboard RookMasks[64], RookMasks2[64];
uint16_t *RookAttacks[64];

Bitboard BishopMasks[64], BishopMasks2[64];
uint16_t *BishopAttacks[64];

static uint16_t BishopTable[5248];
static uint16_t RookTable[102400];

static void init_bmi2(uint16_t *table, uint16_t *attacks[], Bitboard *masks,
    Bitboard *masks2, signed char delta[][2])
{
  for (int sq = 0; sq < 64; sq++) {
    attacks[sq] = table;

    // Calculate mask.
    Bitboard mask = 0;
    for (int i = 0; i < 4; i++) {
      if (!valid(sq, delta[i])) continue;
      for (int s = sq + delta[i][0]; valid(s, delta[i]); s += delta[i][0])
        mask |= bit(s); 
    }
    masks[sq] = mask;

    Bitboard mask2 = 0;
    for (int i = 0; i < 4; i++) {
      for (int s = sq; valid(s, delta[i]); s += delta[i][0])
        mask2 |= bit(s + delta[i][0]); 
    }
    masks2[sq] = mask2;

    // Use Carry-Rippler trick.
    Bitboard b = 0;
    do {
      Bitboard atts = 0;
      for (int j = 0; j < 4; j++)
        for (int s = sq; valid(s, delta[j]); s += delta[j][0]) {
          atts |= bit(s + delta[j][0]);
          if (b & bit(s + delta[j][0])) break;
        }
      attacks[sq][_pext_u64(b, mask)] = _pext_u64(atts, mask2);
      b = (b - mask) & mask;
      table++;
    } while (b);
  }
}

static void init_sliding_attacks(void)
{
  init_bmi2(BishopTable, BishopAttacks, BishopMasks, BishopMasks2, BishopDelta);
  init_bmi2(RookTable, RookAttacks, RookMasks, RookMasks2, RookDelta);
}

#endif
