#ifndef TB8GEN_H
#define TB8GEN_H

struct IdxInfo {
  int numsets;   // number of sets of like pieces, excluding kings.
//  int norm[MAX_PIECES];
  uint64_t size;
  uint64_t factor[MAX_SETS];
  uint32_t subfactor[MAX_SETS]; // total number of placements for a set
  int first[MAX_SETS];          // index of first piece of each set
  int mult[MAX_SETS];           // number of like pieces in each set
  int last[MAX_SETS];
};

extern struct IdxInfo ii, capt_ii[MAX_SETS];
extern size_t kslice_size, kslice_sub_size[MAX_SETS];
extern char *g_tablename;
extern Position g_pos;
extern int16_t KKMap[64][64];

#endif
