#ifdef BMI2_PLAIN

Bitboard RookMasks[64];
Bitboard *RookAttacks[64];

Bitboard BishopMasks[64];
Bitboard *BishopAttacks[64];

Bitboard BishopTable[5248];
Bitboard RookTable[102400];

static void init_bmi2(Bitboard *table, Bitboard *attacks[], Bitboard *masks,
    signed char delta[][2])
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

    // Use Carry-Rippler trick.
    Bitboard b = 0;
    do {
      Bitboard atts = 0;
      for (int j = 0; j < 4; j++)
        for (int s = sq; valid(s, delta[j]); s += delta[j][0]) {
          atts |= bit(s + delta[j][0]);
          if (b & bit(s + delta[j][0])) break;
        }
      attacks[sq][_pext_u64(b, mask)] = atts;
      b = (b - mask) & mask;
      table++;
    } while (b);
  }
}

static void init_sliding_attacks(void)
{
  init_bmi2(BishopTable, BishopAttacks, BishopMasks, BishopDelta);
  init_bmi2(RookTable, RookAttacks, RookMasks, RookDelta);
}

#endif
