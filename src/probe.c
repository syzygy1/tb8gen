/*
  Copyright (c) 2011-2017, 2025, 2026 Ronald de Man

  This file is distributed under the terms of the GNU GPL, version 2.
*/

#include <assert.h>
#include <fcntl.h>
#include <limits.h>
#include <stdatomic.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>
#include <x86intrin.h>

#include "defs.h"
#include "index.h"
#include "movegen.h"
#include "probe.h"
#include "tb8gen.h"
#include "threads.h"
#include "types.h"
#include "util.h"

static constexpr int MAX_TABLES = 4031;
static constexpr int TB_HASHBITS = 13;

const char *suffix[] = { ".rtbw", ".rtbm", ".rtbz" };
uint32_t magic[] = { 0x5d23e871, 0x88ac504b, 0xa50c66d7 };
uint32_t magic2[] = { 0x9ca55208, 0xd895e4a4, 0xeaf9b743 };



struct HashEntry {
  uint64_t key;
  struct TbEntry *ptr;
};

static LOCK_T mutex;
static char *path_string = nullptr;
static int num_paths = 0;
static char **paths;

static int num_tbs;
//static int num_wdl, num_dtm;

static struct TbEntry *tb_entry;
static struct HashEntry tb_hash[1 << TB_HASHBITS];

static constexpr uint64_t MaterialPieceKey[16] = {
  [1] = 0x5ced000000000001ULL,
        0xe173000000000010ULL,
        0xd64d000000000100ULL,
        0xab88000000001000ULL,
        0x680b000000010000ULL,
  [9] = 0xf209000000100000ULL,
        0xbb14000001000000ULL,
        0x58df000010000000ULL,
        0xa15f000100000000ULL,
        0x7c94001000000000ULL
};

uint32_t Binomial[8][64];

static uint64_t material_key_from_counts(int white_cnts[8], int black_cnts[8])
{
  uint64_t key = 0;

  for (int i = PAWN; i <= QUEEN ; i++)
    key +=  white_cnts[i] * MaterialPieceKey[i]
          + black_cnts[i] * MaterialPieceKey[i + 8];

  return key;
}

static LOCK_T fail_mutex, dtm_mutex;

void init_indices(void);

static FD open_tb(const char *str, const int type)
{
  char name[256];

  for (int i = 0; i < num_paths; i++) {
    strcpy(name, paths[i]);
    strcat(name, "/");
    strcat(name, str);
    strcat(name, suffix[type]);
    FD fd = open_file(name);
    if (fd != FD_ERR)
      return fd;
  }
  return FD_ERR;
}

static bool test_tb(const char *str, const int type)
{
  FD fd = open_tb(str, type);
  if (fd == FD_ERR) return false;
  close_file(fd);
  return true;
}

const void *map_tb(const char *name, const int type, map_t *mapping,
    bool use_paths)
{
  FD fd;

  if (use_paths)
    fd = open_tb(name, type);
  else {
    char str[256];
    strcpy(str, name);
    strcat(str, suffix[type]);
    fd = open_file(str);
  }
  if (fd == FD_ERR) return NULL;

  const void *data = map_file(fd, true, mapping);
  if (data == NULL) {
    fprintf(stderr, "Could not map %s%s into memory.\n", name, suffix[type]);
    exit(EXIT_FAILURE);
  }

  close_file(fd);

  return data;
}

void add_to_hash(void *ptr, uint64_t key)
{
  int idx;

  idx = key >> (64 - TB_HASHBITS);
  while (tb_hash[idx].ptr)
    idx = (idx + 1) & ((1 << TB_HASHBITS) - 1);

  tb_hash[idx] = (struct HashEntry){ .key = key, .ptr = ptr };
}

static void detect_tb(char *str)
{
  if (!test_tb(str, WDL)) return;

  int pcs[16];
  for (int i = 0; i < 16; i++)
    pcs[i] = 0;
  int color = 0;
  for (char *s = str; *s; s++)
    if (*s == 'v')
      color = 8;
    else
      for (int i = PAWN; i <= KING; i++)
        if (*s == PieceChar[i]) {
          pcs[i | color]++;
          break;
        }
  
  uint64_t key  = material_key_from_counts(pcs, pcs + 8);
  uint64_t key2 = material_key_from_counts(pcs + 8, pcs);

  struct TbEntry *entry = &tb_entry[num_tbs++];
  *entry = (struct TbEntry){ 0 };
  entry->has_pawns = pcs[WPAWN] || pcs[BPAWN];
  entry->key = key;
  entry->symmetric = key == key2;
  for (int i = 0; i < 16; i++) {
    entry->num += pcs[i];
    entry->numsets += pcs[i] != 0;
  }
  entry->numsets -= 2;

//  num_dtm += be->has_dtm = test_tb(str, DTM);

  for (int t = 0; t < 3; t++)
    atomic_init(&entry->tbase[t], NULL);

  if (!entry->has_pawns) {
    int j = 0;
    for (int i = 0; i < 16; i++)
      if (pcs[i] == 1) j++;
    entry->kk_enc = j == 2; // for suicide a bit more may have to be done
  } else {
    entry->pawns[0] = pcs[WPAWN];
    entry->pawns[1] = pcs[BPAWN];
    if (pcs[BPAWN] && (!pcs[WPAWN] || pcs[WPAWN] > pcs[BPAWN]))
      Swap(entry->pawns[0], entry->pawns[1]);
  }

  add_to_hash(entry, key);
  if (key != key2)
    add_to_hash(entry, key2);
}

INLINE int num_tables(struct TbEntry *entry, const int t)
{
  return entry->has_pawns ? t == DTM ? 6 : 4 : 1;
}

static void free_tb_table(struct TbTable *table)
{
  free(table->precomp);
}

static void free_tbase(struct TbEntry *entry, struct Tbase *tb)
{
  unmap_file(tb->data, tb->mapping);
  int num =  tb->layout == LT_PIECE ? 1
           : tb->layout == LT_PIECE_K ? 10
           : tb->layout == LT_PIECE_KK ? 462
           : tb->layout == LT_PAWN_FILE ? 4 : 6;
  if (   !entry->symmetric
      && tb->dist_format != WL_WTM
      && tb->dist_format != WL_BTM)
    num *= 2;
  for (int i = 0; i < num; i++) {
    struct TbTable *table = atomic_load_explicit(&tb->table[i],
        memory_order_relaxed);
    if (table) free_tb_table(table);
  }
  free(tb);
}

static void free_tb_entry(struct TbEntry *entry)
{
  for (int i = 0; i < 3; i++) {
    struct Tbase *tb = atomic_load_explicit(&entry->tbase[i],
        memory_order_relaxed);
    if (tb) free_tbase(entry, tb);
  }
}

void free_tablebases(void)
{
  for (int i = 0; i < num_tbs; i++)
    free_tb_entry(&tb_entry[i]);

  free(tb_entry);
  free(path_string);
  free(paths);
}

void create_piece_string(char *s, int n, int idx)
{
  s[n] = 0;
  if (n == 0) return;
  for (int k = n - 1; k > 0; k--) {
    int l = 0;
    while (idx >= Binomial[k + 1][k + 1 + l])
      l++;
    idx -= Binomial[k + 1][k + l];
    s[n - 1 - k] = PieceChar[l + 1];
  }
  s[n - 1] = PieceChar[idx + 1];
}

void init_tablebases(const char *path_list)
{
  path_string = strdup(path_list ? path_list : ".");
  char *p = path_string;
  for (p = strtok(p, SEP_STR); p; p = strtok(nullptr, SEP_STR))
    num_paths++;
  paths = malloc(num_paths * sizeof(*paths));
  p = path_string;
  for (int i = 0; i < num_paths; i++, p += strlen(p) + 1)
    paths[i] = p;

  LOCK_INIT(mutex);
  LOCK_INIT(dtm_mutex);
  LOCK_INIT(fail_mutex);

  num_tbs = 0;

  tb_entry = malloc(MAX_TABLES * sizeof *tb_entry);

  for (int i = 0; i < (1 << TB_HASHBITS); i++)
    tb_hash[i] = (struct HashEntry){ 0 };

  init_indices();

  char white[16], black[16], name[40];

  for (int p = 1; p <= 5; p++)
    for (int q = 0; q <= min(p, 5 - p); q++)
      for (int k = Binomial[4][p + 4] - 1; k >= 0; k--) {
        create_piece_string(white, p, k);
        for (int l = q < p ? Binomial[4][q + 4] - 1 : k; l >= 0; l--) {
          create_piece_string(black, q, l);
          sprintf(name, "K%svK%s", white, black);
          detect_tb(name);
        }
      }

//  printf("Found %d WDL and %d DTM tablebase files.\n", num_tbs, num_dtz);
  printf("Found %d WDL tablebase files.\n", num_tbs);
}

const int8_t OffDiag[64] = {
  0, -1, -1, -1, -1, -1, -1, -1,
  1,  0, -1, -1, -1, -1, -1, -1,
  1,  1,  0, -1, -1, -1, -1, -1,
  1,  1,  1,  0, -1, -1, -1, -1,
  1,  1,  1,  1,  0, -1, -1, -1,
  1,  1,  1,  1,  1,  0, -1, -1,
  1,  1,  1,  1,  1,  1,  0, -1,
  1,  1,  1,  1,  1,  1,  1,  0
};

const uint8_t Triangle[64] = {
  6, 0, 1, 2, 2, 1, 0, 6,
  0, 7, 3, 4, 4, 3, 7, 0,
  1, 3, 8, 5, 5, 8, 3, 1,
  2, 4, 5, 9, 9, 5, 4, 2,
  2, 4, 5, 9, 9, 5, 4, 2,
  1, 3, 8, 5, 5, 8, 3, 1,
  0, 7, 3, 4, 4, 3, 7, 0,
  6, 0, 1, 2, 2, 1, 0, 6
};

// Why not InvDiag[8]?
static const uint8_t InvDiag[16] = {
  0,  9, 18, 27, 36, 45, 54, 63,
  7, 14, 21, 28, 35, 42, 49, 56
};

const uint8_t InvTriangle[10] = {
  1, 2, 3, 10, 11, 19, 0, 9, 18, 27
};

const uint8_t FlipDiag[64] = {
   0,  8, 16, 24, 32, 40, 48, 56,
   1,  9, 17, 25, 33, 41, 49, 57,
   2, 10, 18, 26, 34, 42, 50, 58,
   3, 11, 19, 27, 35, 43, 51, 59,
   4, 12, 20, 28, 36, 44, 52, 60,
   5, 13, 21, 29, 37, 45, 53, 61,
   6, 14, 22, 30, 38, 46, 54, 62,
   7, 15, 23, 31, 39, 47, 55, 63
};

static const uint8_t Lower[64] = {
  28,  0,  1,  2,  3,  4,  5,  6,
   0, 29,  7,  8,  9, 10, 11, 12,
   1,  7, 30, 13, 14, 15, 16, 17,
   2,  8, 13, 31, 18, 19, 20, 21,
   3,  9, 14, 18, 32, 22, 23, 24,
   4, 10, 15, 19, 22, 33, 25, 26,
   5, 11, 16, 20, 23, 25, 34, 27,
   6, 12, 17, 21, 24, 26, 27, 35
};

static const uint8_t Diag[64] = {
   0,  0,  0,  0,  0,  0,  0,  8,
   0,  1,  0,  0,  0,  0,  9,  0,
   0,  0,  2,  0,  0, 10,  0,  0,
   0,  0,  0,  3, 11,  0,  0,  0,
   0,  0,  0, 12,  4,  0,  0,  0,
   0,  0, 13,  0,  0,  5,  0,  0,
   0, 14,  0,  0,  0,  0,  6,  0,
  15,  0,  0,  0,  0,  0,  0,  7
};

static const uint8_t InvLower[36] = {
   1,  2,  3,  4,  5,  6,  7,
  10, 11, 12, 13, 14, 15,
  19, 20, 21, 22, 23,
  28, 29, 30, 31,
  37, 38, 39,
  46, 47,
  55,
   0,  9, 18, 27, 36, 45, 54, 63
};

const uint8_t PawnFlip[2][64] = {
  { [8] =  0,  6, 12, 18, 18, 12,  6,  0,
           1,  7, 13, 19, 19, 13,  7,  1,
           2,  8, 14, 20, 20, 14,  8,  2,
           3,  9, 15, 21, 21, 15,  9,  3,
           4, 10, 16, 22, 22, 16, 10,  4,
           5, 11, 17, 23, 23, 17, 11,  5  },
  { [8] =  0,  1,  2,  3,  3,  2,  1,  0,
           4,  5,  6,  7,  7,  6,  5,  4,
           8,  9, 10, 11, 11, 10,  9,  8,
          12, 13, 14, 15, 15, 14, 13, 12,
          16, 17, 18, 19, 19, 18, 17, 16,
          20, 21, 22, 23, 23, 22, 21, 20,
           0,  0,  0,  0,  0,  0,  0,  0  }
};

const uint8_t PawnTwist[2][64] = {
  { [8] = 47, 35, 23, 11, 10, 22, 34, 46,
          45, 33, 21,  9,  8, 20, 32, 44,
          43, 31, 19,  7,  6, 18, 30, 42,
          41, 29, 17,  5,  4, 16, 28, 40,
          39, 27, 15,  3,  2, 14, 26, 38,
          37, 25, 13,  1,  0, 12, 24, 36  },
  { [8] = 47, 45, 43, 41, 40, 42, 44, 46,
          39, 37, 35, 33, 32, 34, 36, 38,
          31, 29, 27, 25, 24, 26, 28, 30,
          23, 21, 19, 17, 16, 18, 20, 22,
          15, 13, 11,  9,  8, 10, 12, 14,
           7,  5,  3,  1,  0,  2,  4,  6  }
};

const uint8_t InvPawnFlip[2][24] = {
  {  8, 16, 24, 32, 40, 48,
     9, 17, 25, 33, 41, 49,
    10, 18, 26, 34, 42, 50,
    11, 19, 27, 35, 43, 51  },
  {  8,  9, 10, 11,
    16, 17, 18, 19,
    24, 25, 26, 27,
    32, 33, 34, 35,
    40, 41, 42, 43,
    48, 49, 50, 51  }
};

const uint8_t InvPawnTwist[2][48] = {
  {  52, 51, 44, 43, 36, 35, 28, 27, 20, 19, 12, 11,
     53, 50, 45, 42, 37, 34, 29, 26, 21, 18, 13, 10,
     54, 49, 46, 41, 38, 33, 30, 25, 22, 17, 14,  9,
     55, 48, 47, 40, 39, 32, 31, 24, 23, 16, 15,  8  },
  {  52, 51, 53, 50, 54, 49, 55, 48,
     44, 43, 45, 42, 46, 41, 47, 40,
     36, 35, 37, 34, 38, 33, 39, 32,
     28, 27, 29, 26, 30, 25, 31, 24,
     20, 19, 21, 18, 22, 17, 23, 16,
     12, 11, 13, 10, 14,  9, 15,  8  }
};

const uint8_t PawnMap[48] = {
   8,  9, 10, 11, 12, 13, 14, 15,
  16, 17, 18, 19, 20, 21, 22, 23,
  24, 25, 26, 27, 28, 29, 30, 31,
  32, 33, 34, 35, 36, 37, 38, 39,
  40, 41, 42, 43, 44, 45, 46, 47,
  48, 49, 50, 51, 52, 53, 54, 55
};

const int16_t KKIdx[10][64] = {
  {  -1,  -1,  -1,   0,   1,   2,   3,   4,
     -1,  -1,  -1,   5,   6,   7,   8,   9,
     10,  11,  12,  13,  14,  15,  16,  17,
     18,  19,  20,  21,  22,  23,  24,  25,
     26,  27,  28,  29,  30,  31,  32,  33,
     34,  35,  36,  37,  38,  39,  40,  41,
     42,  43,  44,  45,  46,  47,  48,  49,
     50,  51,  52,  53,  54,  55,  56,  57  },
  {  58,  -1,  -1,  -1,  59,  60,  61,  62,
     63,  -1,  -1,  -1,  64,  65,  66,  67,
     68,  69,  70,  71,  72,  73,  74,  75,
     76,  77,  78,  79,  80,  81,  82,  83,
     84,  85,  86,  87,  88,  89,  90,  91,
     92,  93,  94,  95,  96,  97,  98,  99,
    100, 101, 102, 103, 104, 105, 106, 107,
    108, 109, 110, 111, 112, 113, 114, 115  },
  { 116, 117,  -1,  -1,  -1, 118, 119, 120,
    121, 122,  -1,  -1,  -1, 123, 124, 125,
    126, 127, 128, 129, 130, 131, 132, 133,
    134, 135, 136, 137, 138, 139, 140, 141,
    142, 143, 144, 145, 146, 147, 148, 149,
    150, 151, 152, 153, 154, 155, 156, 157,
    158, 159, 160, 161, 162, 163, 164, 165,
    166, 167, 168, 169, 170, 171, 172, 173  },
  { 174,  -1,  -1,  -1, 175, 176, 177, 178,
    179,  -1,  -1,  -1, 180, 181, 182, 183,
    184,  -1,  -1,  -1, 185, 186, 187, 188,
    189, 190, 191, 192, 193, 194, 195, 196,
    197, 198, 199, 200, 201, 202, 203, 204,
    205, 206, 207, 208, 209, 210, 211, 212,
    213, 214, 215, 216, 217, 218, 219, 220,
    221, 222, 223, 224, 225, 226, 227, 228  },
  { 229, 230,  -1,  -1,  -1, 231, 232, 233,
    234, 235,  -1,  -1,  -1, 236, 237, 238,
    239, 240,  -1,  -1,  -1, 241, 242, 243,
    244, 245, 246, 247, 248, 249, 250, 251,
    252, 253, 254, 255, 256, 257, 258, 259,
    260, 261, 262, 263, 264, 265, 266, 267,
    268, 269, 270, 271, 272, 273, 274, 275,
    276, 277, 278, 279, 280, 281, 282, 283  },
  { 284, 285, 286, 287, 288, 289, 290, 291,
    292, 293,  -1,  -1,  -1, 294, 295, 296,
    297, 298,  -1,  -1,  -1, 299, 300, 301,
    302, 303,  -1,  -1,  -1, 304, 305, 306,
    307, 308, 309, 310, 311, 312, 313, 314,
    315, 316, 317, 318, 319, 320, 321, 322,
    323, 324, 325, 326, 327, 328, 329, 330,
    331, 332, 333, 334, 335, 336, 337, 338  },
  {  -1,  -1, 339, 340, 341, 342, 343, 344,
     -1,  -1, 345, 346, 347, 348, 349, 350,
     -1,  -1, 441, 351, 352, 353, 354, 355,
     -1,  -1,  -1, 442, 356, 357, 358, 359,
     -1,  -1,  -1,  -1, 443, 360, 361, 362,
     -1,  -1,  -1,  -1,  -1, 444, 363, 364,
     -1,  -1,  -1,  -1,  -1,  -1, 445, 365,
     -1,  -1,  -1,  -1,  -1,  -1,  -1, 446  },
  {  -1,  -1,  -1, 366, 367, 368, 369, 370,
     -1,  -1,  -1, 371, 372, 373, 374, 375,
     -1,  -1,  -1, 376, 377, 378, 379, 380,
     -1,  -1,  -1, 447 ,381, 382, 383, 384,
     -1,  -1,  -1,  -1 ,448, 385, 386, 387,
     -1,  -1,  -1,  -1 , -1, 449, 388, 389,
     -1,  -1,  -1,  -1 , -1,  -1, 450, 390,
     -1,  -1,  -1,  -1 , -1,  -1,  -1, 451  },
  { 452, 391, 392, 393 ,394, 395, 396, 397,
     -1,  -1,  -1,  -1 ,398, 399, 400, 401,
     -1,  -1,  -1,  -1 ,402, 403, 404, 405,
     -1,  -1,  -1,  -1 ,406, 407, 408, 409,
     -1,  -1,  -1,  -1 ,453, 410, 411, 412,
     -1,  -1,  -1,  -1 , -1, 454, 413, 414,
     -1,  -1,  -1,  -1 , -1,  -1, 455, 415,
     -1,  -1,  -1,  -1 , -1,  -1,  -1, 456  },
  { 457, 416, 417, 418 ,419, 420, 421, 422,
     -1, 458, 423, 424 ,425, 426, 427, 428,
     -1,  -1,  -1,  -1 , -1, 429, 430, 431,
     -1,  -1,  -1,  -1 , -1, 432, 433, 434,
     -1,  -1,  -1,  -1 , -1, 435, 436, 437,
     -1,  -1,  -1,  -1 , -1, 459, 438, 439,
     -1,  -1,  -1,  -1 , -1,  -1, 460, 440,
     -1,  -1,  -1,  -1 , -1,  -1,  -1, 461  }
};

uint8_t KKSquare[462][2];

static const uint8_t FileToFile[] = { 0, 1, 2, 3, 3, 2, 1, 0 };

size_t PawnIdx[2][6][24];
static size_t PawnFactorFile[6][4];
static size_t PawnFactorRank[6][6];
static uint8_t Off10[10][64];

void init_indices(void)
{
  for (int i = 0; i < 10; i++)
    for (int j = 0; j < 64; j++)
      if (KKIdx[i][j] >= 0) {
	KKSquare[KKIdx[i][j]][0] = InvTriangle[i];
	KKSquare[KKIdx[i][j]][1] = j;
      }

  // Binomial[k][n] = Bin(n, k)
  for (int j = 0; j < 64; j++)
    Binomial[0][j] = 1;
  for (int i = 1; i < 8; i++)
    for (int j = 1; j < 64; j++)
      Binomial[i][j] = Binomial[i - 1][j - 1] + Binomial[i][j - 1];

  for (int i = 0; i < 6; i++) {
    size_t s = 0;
    for (int j = 0; j < 24; j++) {
      PawnIdx[0][i][j] = s;
      s += Binomial[i][PawnTwist[0][InvPawnFlip[0][j]]];
      if ((j + 1) % 6 == 0) {
        PawnFactorFile[i][j / 6] = s;
        s = 0;
      }
    }
  }

  for (int i = 0; i < 6; i++) {
    size_t s = 0;
    for (int j = 0; j < 24; j++) {
      PawnIdx[1][i][j] = s;
      s += Binomial[i][PawnTwist[1][InvPawnFlip[1][j]]];
      if ((j + 1) % 4 == 0) {
        PawnFactorRank[i][j / 4] = s;
        s = 0;
      }
    }
  }

  for (int i = 0; i < 10; i++) {
    int n = 0;
    for (int j = 0; j < 64; j++) {
      if (KKIdx[i][j] < 0) continue;
      Off10[i][j] = n++;
    }
  }
}

INLINE int leading_pawn(uint8_t *p, struct TbEntry *entry, const int lt)
{
  const int enc = lt - LT_PAWN_FILE;
  for (int i = 1; i < entry->pawns[0]; i++)
    if (PawnFlip[enc][p[0]] > PawnFlip[enc][p[i]])
      Swap(p[0], p[i]);

  return lt == LT_PAWN_FILE ? FileToFile[p[0] & 7] : (p[0] - 8) >> 3;
}

INLINE void sort_squares_mapped(int n, uint8_t *p, const uint8_t *map)
{
  for (int i = 0; i < n; i++)
    for (int j = i + 1; j < n; j++)
      if (map[p[i]] > map[p[j]])
        Swap(p[i], p[j]);
}

INLINE size_t encode(uint8_t *p, struct EncInfo *ei, struct TbEntry *entry,
    const int lt)
{
  int n = entry->num;
  size_t idx;
  Bitboard occ = 0;
  int k;

  if (p[0] & 0x04)
    for (int i = 0; i < 8; i++)
      p[i] ^= 0x07;

  if (lt == LT_PIECE) {
    if (p[0] & 0x20)
      for (int i = 0; i < 8; i++)
        p[i] ^= 0x38;

    for (int i = 0; i < n; i++)
      if (OffDiag[p[i]]) {
        if (OffDiag[p[i]] > 0 && i < (entry->kk_enc ? 2 : 3))
          for (int j = 0; j < n; j++)
            p[j] = FlipDiag[p[j]];
        break;
      }

    if (entry->kk_enc) {
      idx = KKIdx[Triangle[p[0]]][p[1]];
      k = 2;
    } else {
      int s1 = (p[1] > p[0]);
      int s2 = (p[2] > p[0]) + (p[2] > p[1]);

      if (OffDiag[p[0]])
        idx = Triangle[p[0]] * 63*62 + (p[1] - s1) * 62 + (p[2] - s2);
      else if (OffDiag[p[1]])
        idx = 6*63*62 + Diag[p[0]] * 28*62 + Lower[p[1]] * 62 + p[2] - s2;
      else if (OffDiag[p[2]])
        idx =  6*63*62 + 4*28*62 + Diag[p[0]] * 7*28
             + (Diag[p[1]] - s1) * 28 + Lower[p[2]];
      else
        idx =  6*63*62 + 4*28*62 + 4*7*28 + Diag[p[0]] * 7*6
             + (Diag[p[1]] - s1) * 6 + (Diag[p[2]] - s2);
      k = 3;
    }
    idx *= ei->factor[0];
    for (int i = 0; i < k; i++)
      occ |= bit(p[i]);
  } else { /* LT_PAWN_FILE/_RANK */
    const int enc = lt - LT_PAWN_FILE;
    for (int i = 1; i < entry->pawns[0]; i++)
      for (int j = i + 1; j < entry->pawns[0]; j++)
        if (PawnTwist[enc][p[i]] < PawnTwist[enc][p[j]])
          Swap(p[i], p[j]);

    k = entry->pawns[0];
    idx = PawnIdx[enc][k-1][PawnFlip[enc][p[0]]];
    for (int i = 1; i < k; i++)
      idx += Binomial[k-i][PawnTwist[enc][p[i]]];
    idx *= ei->factor[0];

    for (int i = 0; i < k; i++)
      occ |= bit(p[i]);

    // Pawns of other color
    if (entry->pawns[1]) {
      int t = k + entry->pawns[1];
      sort_squares(entry->pawns[1], &p[k]);
      size_t s = 0;
      for (int i = k; i < t; i++) {
        int rank = rank_among_free(p[i], occ);
        s += Binomial[i - k + 1][rank - 8];
      }
      idx += s * ei->factor[k];
      for (; k < t; k++)
        occ |= bit(p[k]);
    }
  }

  for (; k < n;) {
    int t = k + ei->norm[k];
    sort_squares(ei->norm[k], &p[k]);
    size_t s = 0;
    for (int i = k; i < t; i++) {
      int rank = rank_among_free(p[i], occ);
      s += Binomial[i - k + 1][rank];
    }
    idx += s * ei->factor[k];
    for (; k < t; k++)
      occ |= bit(p[k]);
  }

  return idx;
}

NOINLINE size_t encode_piece(uint8_t *p, struct EncInfo *ei,
    struct TbEntry *entry)
{
  return encode(p, ei, entry, LT_PIECE);
}

NOINLINE size_t encode_pawn_f(uint8_t *p, struct EncInfo *ei,
    struct TbEntry *entry)
{
  return encode(p, ei, entry, LT_PAWN_FILE);
}

NOINLINE size_t encode_pawn_r(uint8_t *p, struct EncInfo *ei,
    struct TbEntry *entry)
{
  return encode(p, ei, entry, LT_PAWN_RANK);
}

size_t set_enc_info(struct EncInfo *ei, struct TbEntry *entry,
    const uint8_t *tb, int shift, int fr, const int lt)
{
  bool more_pawns = lt != LT_PIECE && entry->pawns[1] > 0;

  for (int i = 0; i < entry->num; i++) {
    ei->pieces[i] = (tb[i + 1 + more_pawns] >> shift) & 0x0f;
    ei->norm[i] = 0;
  }

  int order = (tb[0] >> shift) & 0x0f;
  int order2 = more_pawns ? (tb[1] >> shift) & 0x0f : 0x0f;

  int k = ei->norm[0] =  lt != LT_PIECE ? entry->pawns[0]
                       : entry->kk_enc ? 2 : 3;

  if (more_pawns) {
    ei->norm[k] = entry->pawns[1];
    k += ei->norm[k];
  }

  for (int i = k; i < entry->num; i += ei->norm[i])
    for (int j = i; j < entry->num && ei->pieces[j] == ei->pieces[i]; j++)
      ei->norm[i]++;

  int n = 64 - k;
  size_t f = 1;

  for (int i = 0; k < entry->num || i == order || i == order2; i++) {
    if (i == order) {
      ei->factor[0] = f;
      f *=  lt == LT_PAWN_FILE ? PawnFactorFile[ei->norm[0] - 1][fr]
          : lt == LT_PAWN_RANK ? PawnFactorRank[ei->norm[0] - 1][fr]
          : entry->kk_enc ? 462 : 31332;
    } else if (i == order2) {
      ei->factor[ei->norm[0]] = f;
      f *= Binomial[ei->norm[ei->norm[0]]][48 - ei->norm[0]];
    } else {
      ei->factor[k] = f;
      f *= Binomial[ei->norm[k]][n];
      n -= ei->norm[k];
      k += ei->norm[k];
    }
  }

  return f;
}

size_t set_dec_info(struct DecInfo *di, struct TbEntry *entry, uint8_t *pcs,
    uint8_t *type_perm, int order, int order2, int fr, const int lt)
{
  bool more_pawns = lt != LT_PIECE && entry->pawns[1] > 0;

  for (int i = 0; i < entry->num; i++)
    di->norm[i] = 0;

  int k = di->norm[0] =  lt != LT_PIECE ? entry->pawns[0]
                       : entry->kk_enc ? 2 : 3;

  if (more_pawns) {
    di->norm[k] = entry->pawns[1];
    k += di->norm[k];
  }

  int i, n = 64 - k;
  size_t f = 1;

  for (i = 0; k < entry->num || i == order || i == order2; i++) {
    if (i == order) {
      di->factor[0] = f;
      f *=  lt == LT_PAWN_FILE ? PawnFactorFile[di->norm[0] - 1][fr]
          : lt == LT_PAWN_RANK ? PawnFactorRank[di->norm[0] - 1][fr]
          : entry->kk_enc ? 462 : 31332;
      if (lt == LT_PIECE)
        i += di->norm[0] - 1;
    } else if (i == order2) {
      di->factor[di->norm[0]] = f;
      f *= Binomial[di->norm[di->norm[0]]][48 - di->norm[0]];
    } else {
      di->factor[k] = f;
      di->norm[k] = pcs[type_perm[i]];
      f *= Binomial[di->norm[k]][n];
      n -= di->norm[k];
      k += di->norm[k];
    }
  }

  uint8_t tmp[8];
  k = di->norm[0];
  if (order2 < 0x0f)
    k += di->norm[k];
  for (i = 0; k < entry->num || i == order || i == order2; i++)
    if (i == order)
      tmp[i] = 0;
    else if (i == order2)
      tmp[i] = di->norm[0];
    else {
      tmp[i] = k;
      k += di->norm[k];
    }

  for (k = 0; k < i; k++) {
    di->order[i - k - 1] = tmp[k];
    di->ord_iter[k] = tmp[k];
  }

  for (k = 0; k < i - 1 ; k++)
    di->fac_iter[di->ord_iter[k]] =
      di->factor[di->ord_iter[k + 1]] / di->factor[di->ord_iter[k]];
  di->fac_iter[di->ord_iter[k]] = f / di->factor[di->ord_iter[k]] + 1;

  return f;
}

INLINE Bitboard unrank_binomial_mapped(uint64_t idx, int n, uint8_t *restrict p,
    const uint8_t *restrict const map, Bitboard occ)
{
  if (n == 0)
    return occ;

  Bitboard b = ~occ;
  for (int i = n - 1; i > 0; i--) {
    int r = i;
    while (idx >= Binomial[i + 1][r + 1])
      r++;
    idx -= Binomial[i + 1][r];
//    r = map ? map[r] : r;
    r = map[r];
    Bitboard b1 = _pdep_u64(1ULL << r, b);
    p[i] = lsb(b1);
    occ |= b1;
  }
//  idx = map ? map[idx] : idx;
  idx = map[idx];
  Bitboard b1 = _pdep_u64(1ULL << idx, b);
  p[0] = lsb(b1);
  occ |= b1;

  return occ;
}

INLINE uint8_t place_empty(uint8_t n, Bitboard *b)
{
  *b = 1ULL << n;
  return n;
}

INLINE uint8_t place_no_skip(uint8_t n, Bitboard *b)
{
  *b |= 1ULL << n;
  return n;
}

// set the nth zero bit to one and return its index
INLINE uint8_t place_skip(uint8_t n, Bitboard *b)
{
  Bitboard b1 = _pdep_u64(1ULL << n, ~*b);
  *b |= b1;
  return lsb(b1);
}

INLINE void decode_helper(uint32_t *sub, uint8_t *restrict p,
    struct DecInfo *di, struct TbEntry *entry, const int fr, const int lt)
{
  uint32_t q, r;
  int i;
  Bitboard occ = 0;
  constexpr Bitboard DIAGONAL_BB = 0x8040201008040201ULL;

  q = sub[0];
  if (lt == LT_PIECE) {
    switch (entry->kk_enc) {

    case 0: /* 111 */
      if (q < 6*63*62) {
        r = q / (63*62);
        q -= r * 63*62;
        p[0] = place_empty(InvTriangle[r], &occ);
        r = q / 62;
        q -= r * 62;
        r += (r >= p[0]);
        p[1] = place_no_skip(r, &occ);
        p[2] = place_skip(q, &occ);
      } else if (q < 6*63*62 + 4*28*62) {
        q -= 6*63*62;
        r = q / (28*62);
        q -= r * 28*62;
        p[0] = place_empty(InvDiag[r], &occ);
        r = q / 62;
        q -= r * 62;
        p[1] = place_no_skip(InvLower[r], &occ);
        p[2] = place_skip(q, &occ);
      } else if (q < 6*63*62 + 4*28*62 + 4*7*28) {
        q -= 6*63*62 + 4*28*62;
        r = q / (7*28);
        q -= r * 7*28;
        p[0] = place_empty(InvDiag[r], &occ);
        r = q / 28;
        q -= r * 28;
        r += (InvDiag[r] >= p[0]);
        p[1] = place_no_skip(InvDiag[r], &occ);
        p[2] = place_no_skip(InvLower[q], &occ);
      } else {
        q -= 6*63*62 + 4*28*62 + 4*7*28;
        r = q / (7 * 6);
        q -= r * 7*6;
        p[0] = place_empty(InvDiag[r], &occ);
        r = q / 6;
        q -= r * 6;
        r += (InvDiag[r] >= p[0]);
        p[1] = place_no_skip(InvDiag[r], &occ);
        occ |= ~DIAGONAL_BB;
        p[2] = place_skip(q, &occ); // qth free square on A1-H8
        occ &= DIAGONAL_BB;
      }
      i = 3;
    break;

    default: /* K2 */
      p[0] = place_empty(KKSquare[q][0], &occ);
      p[1] = place_no_skip(KKSquare[q][1], &occ);
      i = 2;
      break;
    }
  }
  else { /* LT_PAWN_FILE, LT_PAWN_RANK */
    const int num = lt == LT_PAWN_RANK ? 6 : 4;
    const int enc = lt - LT_PAWN_FILE;
    int t = entry->pawns[0] - 1;
    for (i = 0; i < num - 1; i++)
      if (PawnIdx[enc][t][num * fr + i + 1] > q) break;
    q -= PawnIdx[enc][t][num * fr + i];
    p[0] = place_empty(InvPawnFlip[enc][num * fr + i], &occ);

    if (t > 0)
      occ = unrank_binomial_mapped(q, t, p + 1, InvPawnTwist[enc], occ);

    i = entry->pawns[0];

    if (entry->pawns[1] > 0) {
      q = sub[entry->pawns[0]];
      occ = unrank_binomial_mapped(q, entry->pawns[1], p + i, PawnMap, occ);
      i += entry->pawns[1];
    }

    for (; i < di->norm[0]; i++)
      p[i] = place_skip(sub[i], &occ);
  }

  for (; i < entry->num;) {
    q = sub[i];
    occ = unrank_binomial(q, di->norm[i], p + i, occ);
    i += di->norm[i];
  }
}

INLINE void decode(uint64_t idx, uint8_t *restrict p, struct DecInfo *di,
    struct TbEntry *entry, const int fr, const int lt)
{
  uint32_t sub[MAX_PIECES];
  int i;

  // TODO: convert into multiplications, e.g. using libdivide.
  // https://github.com/ridiculousfish/libdivide
  for (i = 0; di->factor[di->order[i]] != 1; i++) {
    uint32_t q = idx / di->factor[di->order[i]];
    idx -= q * di->factor[di->order[i]];
    sub[di->order[i]] = q;
  }
  sub[di->order[i]] = idx;

  decode_helper(sub, p, di, entry, fr, lt);
}

NOINLINE void decode_init(uint32_t *sub, uint64_t idx, struct DecInfo *di)
{
  int i;
  for (i = 0; di->factor[di->order[i]] != 1; i++) {
    uint32_t q = idx / di->factor[di->order[i]];
    idx -= q * di->factor[di->order[i]];
    sub[di->order[i]] = q;
  }
  sub[di->order[i]] = idx;
}

INLINE void decode_iter(uint32_t *sub, uint8_t *restrict p, struct DecInfo *di,
    struct TbEntry *entry, const int fr, const int lt)
{
  decode_helper(sub, p, di, entry, fr, lt);

  int i = 0;
  while (true) {
    int j = di->ord_iter[i];
    if (++sub[j] < di->fac_iter[j])
      break;
    sub[j] = 0;
    i++;
  }
}

NOINLINE void decode_piece(uint64_t idx, uint8_t *restrict p,
    struct DecInfo *di, struct TbEntry *entry)
{
  decode(idx, p, di, entry, 0, LT_PIECE);
}

NOINLINE void decode_piece_iter(uint32_t *sub, uint8_t *restrict p,
    struct DecInfo *di, struct TbEntry *entry)
{
  decode_iter(sub, p, di, entry, 0, LT_PIECE);
}

NOINLINE void decode_pawn_r(uint64_t idx, uint8_t *restrict p,
    struct DecInfo *di, struct TbEntry *entry, int rank)
{
  decode(idx, p, di, entry, rank, LT_PAWN_RANK);
}

NOINLINE void decode_pawn_r_iter(uint32_t *sub, uint8_t *restrict p,
    struct DecInfo *di, struct TbEntry *entry, int rank)
{
  decode_iter(sub, p, di, entry, rank, LT_PAWN_RANK);
}


/******* Start of actual probing and decompression code *******/

static void calc_symlen(struct PairsData *d, uint32_t s, bool *tmp)
{
  const uint32_t w = read_le_u32(d->sympat + 3 * s);
  uint32_t s2 = (w >> 12) & 0xfff;
  if (s2 == 0x0fff)
    d->symlen[s] = 0;
  else {
    uint32_t s1 = w & 0xfff;
    if (!tmp[s1]) calc_symlen(d, s1, tmp);
    if (!tmp[s2]) calc_symlen(d, s2, tmp);
    d->symlen[s] = d->symlen[s1] + d->symlen[s2] + 1;
  }
  tmp[s] = true;
}

struct PairsData *setup_huffman(const uint8_t **ptr)
{
  const uint8_t *data = *ptr;
  int max_len = data[8];
  int min_len = data[9];
  int h = max_len - min_len + 1;
  uint32_t num_syms = read_le_u16(&data[10 + 2 * h]);
  uint16_t *offset = (uint16_t *)(&data[10]);

  int hh = h;
#ifdef LOOKUP
  if (max_len < LUBITS)
    hh = LUBITS - min_len + 1;
#endif

  uint64_t tmp_base[32];
  for (int i = h - 1; i < hh; i++)
    tmp_base[i] = 0;
  for (int i = h - 2; i >= 0; i--)
    tmp_base[i] = (tmp_base[i + 1] + offset[i] - offset[i + 1]) / 2;
  for (int i = 0; i < h; i++)
    tmp_base[i] <<= 64 - (min_len + i);

  int num_lu = 0;
#ifdef LOOKUP
  if (min_len <= LUBITS)
    num_lu = (1 << LUBITS) - (tmp_base[LUBITS - min_len] >> (64 - LUBITS));
#endif

  struct PairsData *d =
    malloc(sizeof(*d) + hh * sizeof(uint64_t) + num_lu * 4 + num_syms);
  d->compr_type = 0;
  d->num_syms = num_syms;
  for (int i = 0; i < hh; i++)
    d->base[i] = tmp_base[i];
  d->offset = offset;
#ifdef LOOKUP
  d->lookup = (struct LUEntry *)((uint8_t *)(d + 1) + hh * sizeof(uint64_t));
  d->symlen = (uint8_t *)(d->lookup + num_lu);
#else
  d->symlen = (uint8_t *)(d + 1) + h * sizeof(uint64_t);
#endif
  d->sympat = &data[12 + 2 * h];
  d->max_len = max_len;
  d->min_len = min_len;
  *ptr = &data[12 + 2 * h + 3 * num_syms + (num_syms & 1)];

  bool tmp[4096] = { 0 };
  for (uint32_t s = 0; s < num_syms; s++)
    if (!tmp[s])
      calc_symlen(d, s, tmp);

#ifdef LOOKUP
  for (int i = 0; i < num_lu; i++) {
    uint64_t code = tmp_base[LUBITS - min_len] + ((uint64_t)i << (64 - LUBITS));
    int bits = LUBITS;
    d->lookup[i] = (struct LUEntry){ 0 };
    for (;;) {
      int l = 0;
      while (code < tmp_base[l]) l++;
      if (l + min_len > bits) break;
      int sym = d->offset[l] + ((code - tmp_base[l]) >> (64 - (l + min_len)));
      d->lookup[i].len += d->symlen[sym] + 1;
      d->lookup[i].bits += l + min_len;
      if (d->lookup[i].cwl == 0)
        d->lookup[i].cwl = l + min_len;
      bits -= l + min_len;
      code <<= (l + min_len);
    }
  }
#endif

  for (uint64_t i = 0; i < (1 << STARTBITS); i++) {
    uint64_t code = ((i + 1) << (64 - STARTBITS)) - 1;
    int l = 0;
    while (code < d->base[l]) l++;
    d->start[i] = l + d->min_len;
  }
  
  d->offset -= d->min_len;

  return d;
}

struct PairsData *setup_rans(const uint8_t **ptr)
{
  const uint8_t *data = *ptr;

  struct PairsData *d = malloc(sizeof(struct PairsData));
  d->compr_type = 1;

  int num_syms;
  d->rans = calloc(1, sizeof(struct RansDecode));
  const uint8_t *p = read_freq_table(d->rans, &num_syms, data + 10);
  d->num_syms = num_syms;
  make_alias_table(d->rans, NULL);

  d->sympat = p;
  d->symlen = malloc(num_syms);
  bool tmp[4096] = { 0 };
  for (int s = 0; s < num_syms; s++)
    if (!tmp[s])
      calc_symlen(d, s, tmp);

  *ptr = p + 3 * num_syms + (num_syms & 1);

  return d;
}

static struct PairsData *setup_pairs(const uint8_t **ptr, size_t tb_size,
   size_t *size, uint8_t *flags, int type, bool new)
{
  struct PairsData *d;
  const uint8_t *data = *ptr;

  if (!new) {
    *flags = data[0];
    if (data[0] & 0x80) {
      d = malloc(sizeof *d);
      d->compr_type = 2;
      d->const_value[0] = data[1];
      d->const_value[1] = 0;
      *ptr = data + 2;
      size[0] = size[1] = size[2] = 0;
      d->tb_size = tb_size;
      return d;
    }
    d = (data[0] & 0x40) ? setup_rans(ptr) : setup_huffman(ptr);
  } else {
    d = data[0] == 2 ? setup_rans(ptr) : setup_huffman(ptr);
  }

  uint8_t block_size = data[1];
  uint8_t idx_bits = data[2];
  uint32_t real_num_blocks = read_le_u32(&data[4]);
  uint32_t num_blocks = real_num_blocks + data[3];

  d->block_size = block_size;
  d->idx_bits = idx_bits;
  d->tb_size = tb_size;

  int num_indices = (tb_size + (1ULL << idx_bits) - 1) >> idx_bits;
  size[0] = 6ULL * num_indices;
  size[1] = 2ULL * num_blocks;
  size[2] = (size_t)real_num_blocks << block_size;

  return d;
}

static const size_t TABLE_SIZE[3] = {
  sizeof(struct WdlTable), sizeof(struct DtmTable), sizeof(struct DtzTable)
};

static const size_t TABLE_SIZE2[3] = {
  sizeof(struct WdlTable2), sizeof(struct DtmTable2), sizeof(struct DtzTable2)
};

static const struct TbTableConst const_table[5] = {
  { NULL, 0 }, { NULL, 1 }, { NULL, 2 }, { NULL, 3 }, { NULL, 4 }
};

NOINLINE static struct Tbase *init_old_layout(struct TbEntry *entry,
    struct Tbase *tb, int type, const uint8_t *data, bool new)
{
  bool split;

  if (!new) {
    split = type != DTZ && (data[0] & 0x01);
    tb->flipped = false;
    tb->dist_format = data[0] & 0x04 ? L_ONLY : WL_BOTH;
  } else {
    split = !entry->symmetric;
    if (type != WDL) {
      tb->dist_format = *++data;
      if (tb->dist_format == WL_WTM || tb->dist_format == WL_BTM)
        split = false;
    }
    tb->flipped = false;
  }

  data++;

  size_t tb_size[6][2];
  struct TbTable *table[12];

  int num = num_tables(entry, type);
  for (int t = 0; t < num; t++) {
    table[t] = malloc(TABLE_SIZE[type]);
    atomic_store_explicit(&tb->table[t], table[t], memory_order_relaxed);
  }
  if (split)
    for (int t = 0; t < num; t++) {
      table[t + num] = malloc(TABLE_SIZE[type]);
      atomic_store_explicit(&tb->table[t + num], table[t + num],
          memory_order_relaxed);
    }

  for (int t = 0; t < num; t++) {
    struct EncInfo *ei = &table[t]->ei;
    tb_size[t][0] = set_enc_info(&ei[t], entry, data, 0, t, tb->layout);
    if (split) {
      ei = &table[num + t]->ei;
      tb_size[t][1] = set_enc_info(ei, entry, data, 4, t, tb->layout);
    }
    data += entry->num + 1 + (entry->has_pawns && entry->pawns[1]);
  }
  data += (uintptr_t)data & 1;

  size_t size[6][2][3];
  for (int t = 0; t < num; t++) {
    uint8_t flags;
    table[t]->precomp = setup_pairs(&data, tb_size[t][0], size[t][0], &flags,
        type, false);
    if (type == DTZ) {
      struct DtzTable *dtz = (struct DtzTable *)table[t];
      dtz->dtz_flags = flags;
    }
    if (split)
      table[num + t]->precomp = setup_pairs(&data, tb_size[t][1], size[t][1],
          &flags, type, false);
  }

  // To be revisited later.
  if (type == DTM && tb->dist_format != L_ONLY && tb->dist_format != W_ONLY) {
    const void *map = data;
    for (int t = 0; t < num; t++) {
      struct DtmTable *dtm = (struct DtmTable *)table[t];
      dtm->map_dtm = map;
      for (int i = 0; i < 2; i++) {
        dtm->map_dtm_idx[i] = (uint16_t *)data + 1 - dtm->map_dtm;
        data += 2 + 2 * read_le_u16(data);
      }
    }
  }

  for (int t = 0; t < num; t++) {
    table[t]->precomp->index_table = data;
    data += size[t][0][0];
    if (split) {
      table[num + t]->precomp->index_table = data;
      data += size[t][1][0];
    }
  }

  for (int t = 0; t < num; t++) {
    table[t]->precomp->size_table = (uint16_t *)data;
    data += size[t][0][1];
    if (split) {
      table[num + t]->precomp->size_table = (uint16_t *)data;
      data += size[t][1][1];
    }
  }

  for (int t = 0; t < num; t++) {
    data = (uint8_t *)(((uintptr_t)data + 0x3f) & ~0x3fULL);
    table[t]->precomp->data = data;
    data += size[t][0][2];
    if (split) {
      data = (uint8_t *)(((uintptr_t)data + 0x3f) & ~0x3fULL);
      table[num + t]->precomp->data = data;
      data += size[t][1][2];
    }
  }

  // To be revisited later.
#if 0
  if (type == DTM && entry->has_pawns) {
    int count[16];
    for (int i = 0; i < 16; i++)
      count[i] = 0;
    for (int i = 0; i < entry->num; i++)
      count[ei[0].pieces[i]]++;
    tb->dtm_switched = material_key_from_counts(count, count + 8) != entry->key;
  }
#endif

  return tb;
}

NOINLINE struct Tbase *init_tbase(struct TbEntry *entry, const char *str,
    const int type, bool use_paths)
{
  map_t mapping;
  const uint8_t *restrict data = map_tb(str, type, &mapping, use_paths);
  if (!data) return false;

  if (read_le_u32(data) == magic[type]) {
    int num = entry->has_pawns ? type == DTM ? 6 : 4 : 1;
    if (type != DTZ && !entry->symmetric)
      num *= 2;
    struct Tbase *tb = malloc(sizeof(struct Tbase) + num * sizeof(void *));
    tb->data = data;
    tb->mapping = mapping;
    tb->layout = !entry->has_pawns ? LT_PIECE : LT_PAWN_FILE;
    return init_old_layout(entry, tb, type, data + 4, false);
  }

  if (read_le_u32(data) == magic2[type]) {
    if (data[4] == 0) {
      int num = entry->has_pawns ? type == DTM ? 6 : 4 : 1;
      if (    !entry->symmetric
          && (type == WDL || (data[5] != WL_WTM && data[5] != WL_BTM)))
        num *= 2;
      struct Tbase *tb = malloc(sizeof(struct Tbase) + num * sizeof(void *));
      tb->data = data;
      tb->mapping = mapping;
      tb->layout =  !entry->has_pawns ? LT_PIECE
                  : type == DTM ? LT_PAWN_RANK : LT_PAWN_FILE;
      return init_old_layout(entry, tb, type, data + 5, true);
    }

    if (data[4] != 1)
      return nullptr;

    const uint8_t *p = data + 4 + entry->num;
    bool two_sided = !entry->symmetric;
    int dist_format;
    if (type != WDL) {
      dist_format = *p++;
      two_sided = two_sided && dist_format != WL_WTM && dist_format != WL_BTM;
    }
    int layout = *p++;
    // FIXME: LT_PIECE should probably be left to init_old() ?
    int num = layout == LT_PIECE ? 1 : layout == LT_PIECE_K ? 10 : 462;
    num *= 1 + two_sided;

    struct Tbase *tbase = calloc(1,
        sizeof(struct Tbase) + num * sizeof(void *));
    tbase->data = data;
    tbase->mapping = mapping;
    tbase->layout = layout;
    if (type != WDL)
      tbase->dist_format = dist_format;
    tbase->flipped = false; // Revisit for pawnful tables.
    tbase->offset = ((p - data) + 7) & ~7;
    tbase->pt[0] = 6;
    tbase->pt[1] = 14;
    for (int i = 2; i < entry->num; i++)
      tbase->pt[i] = (data[4 + i] & 0x07) | ((data[4 + i] & 0x80) >> 4);
    return tbase;
  }

  fprintf(stderr, "Invalid tablebase file.\n");
  unmap_file(data, mapping);
  return nullptr;
}

NOINLINE struct TbTable2 *init_new_table(struct TbEntry *entry,
    struct Tbase *tb, int type, int tidx, int tsq)
{
  const uint64_t *offsets = (uint64_t *)((uint8_t *)tb->data + tb->offset);
  const uint8_t *data = (uint8_t *)tb->data + offsets[tidx];

  if (data[0] == 0xff) {
    if (data[1] >= 5) {
      struct TbTableConst *tbl = malloc(sizeof *tbl);
      *tbl = (struct TbTableConst){ .precomp = NULL, .const_value = data[1] };
      return (struct TbTable2 *)tbl;
    }
    return (struct TbTable2 *)&const_table[data[1]];
  }

  struct TbTable2 *table = malloc(TABLE_SIZE2[type]);
  int mapped = *data++;

  uint8_t first[MAX_SETS];
  uint8_t mult[MAX_SETS];
  int k = 0;
  for (int i = 2, l = 0; i < entry->num; i++) {
    if (tb->pt[i] != l) {
      l = tb->pt[i];
      first[k] = i;
      mult[k] = 0;
      k++;
    }
    assume(k >= 1);
    mult[k - 1]++;
  }
  static constexpr uint8_t knum[] = { 58, 58, 58, 55, 55, 55, 33, 30, 30, 30 };
  uint64_t tb_size = tb->layout == LT_PIECE_KK ? 1 : knum[tsq];
  for (int i = 0, n = 62; i < k; i++) {
    table->first[i] = first[data[i]];
    int m = mult[data[i]];
    table->mult[i] = m;
    table->factor[i] = Binomial[m][n];
    tb_size *= table->factor[i];
    n -= m;
  }
  data += k;
  data += (uintptr_t)data & 1;

  size_t size[3];
  table->precomp = setup_pairs(&data, tb_size, size, NULL, type, true);

  if (type == DTM) {
    struct DtmTable2 *dtm = (struct DtmTable2 *)table;
    dtm->mapped = mapped;
    if (mapped) {
      dtm->map_dtm = (uint16_t *)data;
      for (int i = 0; i < 2; i++) {
        dtm->map_dtm_idx[i] = (uint16_t *)data + 1 - dtm->map_dtm;
        data += 2 + 2 * read_le_u16(data);
      }
    }
  }

#if 0
  if (type == DTZ) {
    struct DtzTable2 *dtz = (struct DtzTable2 *)table;
    dtz->mapped = mapped;
    if (mapped == 1) {
      dtz->map_dtz = data;
      for (int i = 0; i < 4; i++) {
        dtz->map_dtz_idx[i] = data + 1 - dtz->map_dtz;
        data += 1 + data[0];
      }
    } else if (mapped == 2) {
      dtz->map_dtz16 = (uint16_t *)data;
      for (int i = 0; i < 4; i++) {
        dtz->dtzMapIdx[i] = (uint16_t *)data + 1 - dtz->map_dtz16;
        data += 2 + 2 * read_le_u16(data);
      }
    }
  }
#endif

  table->precomp->index_table = data;
  data += size[0];
  table->precomp->size_table = (uint16_t *)data;
  data += size[1];
  data = (uint8_t *)(((uintptr_t)data + 0x3f) & ~0x3fULL);
  table->precomp->data = data;

  return table;
}

INLINE const uint8_t *decompress_rans(struct PairsData *d, uint64_t idx)
{
  uint32_t main_idx = idx >> d->idx_bits;
  int litidx =  (int)((uint32_t)idx & ((1u << d->idx_bits) - 1))
              - (int)(1u << (d->idx_bits - 1));
  uint32_t block = read_le_u32(d->index_table + 6 * main_idx);
  litidx += read_le_u16(d->index_table + 6 * main_idx + 4);

  if (litidx < 0)
    while (litidx < 0)
      litidx += d->size_table[--block] + 1;
  else
    while (litidx > d->size_table[block])
      litidx -= d->size_table[block++] + 1;

  // Since the symbols in the block are decoded in reverse order, we
  // need to start counting from the end.
  litidx -= d->size_table[block] + 1;

  const uint8_t *p = d->data + ((size_t)block << d->block_size);
  const uint8_t *end = p + ((size_t)1 << d->block_size);
  int sym;
  RansState rans;
  rans_dec_init(&rans, &p, end);
  for (; litidx < 0 && p < end; litidx += d->symlen[sym] + 1) {
    sym = rans_dec_get(&rans, d->rans);
    rans_dec_renorm(&rans, &p);
  }
  for (; litidx < 0; litidx += d->symlen[sym] + 1)
    sym = rans_dec_get(&rans, d->rans);

  const uint8_t *sympat = d->sympat;
  while (d->symlen[sym] != 0) {
    uint32_t w = read_le_u32(sympat + 3 * sym);
    int s1 = w & 0xfff;
    if (litidx < (int)d->symlen[s1] + 1)
      sym = s1;
    else {
      litidx -= (int)d->symlen[s1] + 1;
      sym = (w >> 12) & 0xfff;
    }
  }

  return &sympat[3 * sym];
}

INLINE const uint8_t *decompress_huff(struct PairsData *d, uint64_t idx)
{
  uint32_t main_idx = idx >> d->idx_bits;
  int litidx =  (int)((uint32_t)idx & ((1u << d->idx_bits) - 1))
              - (int)(1u << (d->idx_bits - 1));
  uint32_t block = read_le_u32(d->index_table + 6 * main_idx);
  litidx += read_le_u16(d->index_table + 6 * main_idx + 4);

  // Add/subtract sizes until 0 <= litidx <= d->size_table[block].
  if (litidx < 0)
    while (litidx < 0)
      litidx += d->size_table[--block] + 1;
  else
    while (litidx > d->size_table[block])
      litidx -= d->size_table[block++] + 1;

  const uint8_t *ptr = d->data + ((size_t)block << d->block_size);

  int l;
  const uint16_t *offset = d->offset;
  uint64_t *base = d->base - d->min_len;
  uint8_t *symlen = d->symlen;
  uint32_t sym, bitcnt;

  uint64_t bitbuf = read_be_u64(ptr), pending = bitbuf;
  ptr += 8;
  bitcnt = 0; // number of consumed bits in bitbuf
  for (;;) {
#ifdef LOOKUP
    if (bitbuf >= base[LUBITS]) {
      int lu = (bitbuf - base[LUBITS]) >> (64 - LUBITS);
      if (litidx < d->lookup[lu].len) {
	for (;;) {
	  l = d->lookup[lu].cwl;
	  sym = from_le_u16(offset[l]) + ((bitbuf - base[l]) >> (64 - l));
	  if (litidx < (int)symlen[sym] + 1) break;
	  litidx -= (int)symlen[sym] + 1;
	  bitbuf <<= l;
          lu = (bitbuf - base[LUBITS]) >> (64 - LUBITS);
	}
	break;
      }
      litidx -= d->lookup[lu].len;
      l = d->lookup[lu].bits;
      goto refill;
    }
#endif
    l = d->start[bitbuf >> (64 - STARTBITS)];
    while (bitbuf < base[l]) l++;
    sym = from_le_u16(offset[l]) + ((bitbuf - base[l]) >> (64 - l));
    if (litidx < (int)symlen[sym] + 1) break;
    litidx -= (int)symlen[sym] + 1;
refill:
    bitbuf = pending << l;
    bitcnt += l;
    pending = bitbuf | read_be_u64(ptr) >> (64 - bitcnt);
    ptr += (bitcnt >> 5) * sizeof(uint32_t);
    bitcnt &= 31;
  }

  const uint8_t *sympat = d->sympat;
  while (symlen[sym] != 0) {
    uint32_t w = read_le_u32(sympat + 3 * sym);
    int s1 = w & 0xfff;
    if (litidx < (int)symlen[s1] + 1)
      sym = s1;
    else {
      litidx -= (int)symlen[s1] + 1;
      sym = (w >> 12) & 0xfff;
    }
  }

  return &sympat[3 * sym];
}

static const uint8_t *decompress_pairs(struct PairsData *d, uint64_t idx)
{
  switch (d->compr_type) {
  case 0:
    return decompress_huff(d, idx);
  case 1:
    return decompress_rans(d, idx);
  default:
    return d->const_value;
  }
}

void create_material_string(Position *pos, char str[16], bool flip)
{
  int cnt[16];

  for (int i = 0; i < 16; i++)
    cnt[i] = 0;
  for (int i = 0; i < pos->num; i++)
    cnt[pos->pt[i] ^ (flip ? 8 : 0)]++;

  int j = 0;
  for (int i = KING; i >= PAWN; i--)
    while (cnt[i]--)
      str[j++] = PieceChar[i];
  str[j++] = 'v';
  flip ^= 8;
  for (int i = KING; i >= PAWN; i--)
    while (cnt[i + 8]--)
      str[j++] = PieceChar[i];
  str[j] = 0;
}

[[noreturn]] NOINLINE static void probe_failed(Position *pos, int type)
{
  char str[16];

  LOCK(fail_mutex);
  create_material_string(pos, str, false);
  fprintf(stderr, "Missing table: %s%s\n", str, suffix[type]);
  exit(EXIT_FAILURE);
}

INLINE void list_squares(Position *pos, const uint8_t *restrict pt, bool flip,
    uint8_t *restrict p)
{
  for (int i = 0; i < pos->num; ) {
    int t = pt[i] ^ (flip << 3);
    for (int j = 0; j < pos->num; j++)
      if (pos->pt[j] == t)
        p[i++] = pos->sq[j] ^ (flip ? 0x38 : 0x00);
  }
}

INLINE int probe_table(Position *pos, int s, const int type)
{
  // Test for KvK
  if (type == WDL && pos->num == 2) return 0;

  uint64_t key = 0;
  for (int i = 0; i < pos->num; i++)
    key += MaterialPieceKey[pos->pt[i]];

  int hash_idx = key >> (64 - TB_HASHBITS);
  while (tb_hash[hash_idx].key && tb_hash[hash_idx].key != key)
    hash_idx = (hash_idx + 1) & ((1 << TB_HASHBITS) - 1);
  if (!tb_hash[hash_idx].ptr)
    probe_failed(pos, type);

  struct TbEntry *entry = tb_hash[hash_idx].ptr;
  if ((type == DTM && !entry->has_dtm) || (type == DTZ && !entry->has_dtz))
    probe_failed(pos, type);

  struct Tbase *tb;

  // Use double-checked locking to reduce locking overhead
  if (!(tb = atomic_load_explicit(&entry->tbase[type], memory_order_acquire))) {
    LOCK(mutex);
    if (!(tb = atomic_load_explicit(&entry->tbase[type], memory_order_relaxed)))
    {
      char str[16];
      create_material_string(pos, str, entry->key != key);
      if (!(tb = init_tbase(entry, str, type, true)))
        probe_failed(pos, type);
      atomic_store_explicit(&entry->tbase[type], tb, memory_order_release);
    }
    UNLOCK(mutex);
  }

  uint8_t p[MAX_PIECES];
  bool flip = !entry->symmetric ? (key != entry->key) != tb->flipped
                                : pos->stm != WHITE;
  int btm_side = (pos->stm == WHITE) == flip;

  if (tb->layout != LT_PIECE_K && tb->layout != LT_PIECE_KK) {

    uint64_t idx;
    int t = 0;
    struct EncInfo *ei;
    struct TbTable *table;

    if (!entry->has_pawns) {
      table = atomic_load_explicit(&tb->table[btm_side], memory_order_relaxed);
      ei = &table->ei;
      list_squares(pos, ei->pieces, flip, p);
      idx = encode_piece(p, ei, entry);
    } else {
      table = atomic_load_explicit(&tb->table[0], memory_order_relaxed);
      list_squares(pos, table->ei.pieces, flip, p);
      t = leading_pawn(p, entry, tb->layout);
      t += btm_side * (type == DTM ? 6 : 4);
      table = atomic_load_explicit(&tb->table[t], memory_order_relaxed);
      ei = &table->ei;
      list_squares(pos, ei->pieces, flip, p);
      // Bring the leading pawn back to the front.
      leading_pawn(p, entry, tb->layout);
      idx = type != DTM ? encode_pawn_f(p, ei, entry)
                        : encode_pawn_r(p, ei, entry);
    }

    const uint8_t *w = decompress_pairs(table->precomp, idx);

    if (type == WDL) return (int)w[0] - 2;

    int v = read_le_u16(w) & 0xfff;

    if (type == DTM) {
      struct DtmTable *dtm = (struct DtmTable *)table;
      v = from_le_u16(dtm->map_dtm[dtm->map_dtm_idx[s] + v]);
    }

    return v;

  } else {

    struct TbTable2 *table;
    list_squares(pos, tb->pt, flip, p);

    if (tb->layout == LT_PIECE_K && btm_side)
      Swap(p[0], p[1]);

    int tsq = tb->layout == LT_PIECE_KK ? KKMap[p[0]][p[1]] : Triangle[p[0]];
    int t = !entry->symmetric ? 2 * tsq + btm_side : tsq;

    if (!(table = atomic_load_explicit(&tb->table[t], memory_order_acquire))) {
      LOCK(mutex);
      if (!(table = atomic_load_explicit(&tb->table[t], memory_order_relaxed)))
      {
        table = init_new_table(entry, tb, type, t, tsq);
        atomic_store_explicit(&tb->table[t], table, memory_order_release);
      }
      UNLOCK(mutex);
    }

    if (!table->precomp)
      return (int)((struct TbTableConst *)table)->const_value
        + (type == WDL ? -2 : 0);

    // Normalize the position.
    uint8_t mask = MirrorMask[p[0]];
    for (int i = 0; i < entry->num; i++)
      p[i] = p[i] ^ mask;

    if (FlipTest[p[0]][p[1]])
      for (int i = 0; i < entry->num; i++)
        p[i] = FlipDiag[p[i]];

    // Calculate index.
    uint64_t idx = tb->layout == LT_PIECE_KK ? 0 : Off10[Triangle[p[0]]][p[1]];
    Bitboard occ = bit(p[0]) | bit(p[1]);
    for (int k = 0; k < entry->numsets; k++) {
      int m = table->first[k];
      sort_squares(table->mult[k], &p[m]);
      size_t s = 0;
      Bitboard occ2 = occ;
      for (int i = 0; i < table->mult[k]; i++, m++) {
        int rank = rank_among_free(p[m], occ);
        occ2 |= bit(p[m]);
        s += Binomial[i + 1][rank];
      }
      idx = idx * table->factor[k] + s;
      occ = occ2;
    }

    const uint8_t *w = decompress_pairs(table->precomp, idx);

    if (type == WDL)
      return (int)w[0] - 2;

    int v = read_le_u16(w) & 0xfff;

    if (type == DTM) {
      struct DtmTable2 *dtm = (struct DtmTable2 *)table;
      if (dtm->mapped)
        v = from_le_u16(dtm->map_dtm[dtm->map_dtm_idx[s] + v]);
    }

    return v;
 
  }
}

// No need to instantiate the DTZ version
NOINLINE static int probe_wdl_table(Position *pos)
{
  return probe_table(pos, 0, WDL);
}

NOINLINE static int probe_dtm_table(Position *pos, bool won)
{
  return probe_table(pos, won, DTM);
}

[[gnu::always_inline]]
inline int probe_capts_wdl(Position *pos, int alpha, int beta)
{
  for (int i = 0; i < pos->num; i++) {
    if ((pos->pt[i] >> 3) != pos->stm)
      continue;
#ifdef HAS_PAWNS
    bool is_pawn = (pos->pt[i] & 7) == PAWN;
#endif
    int from = pos->sq[i];
    Bitboard b = piece_attacks(pos->pt[i], from, pos->occ) & pos->occ;
    while (b) {
      int to = pop_lsb(&b);
      int j = piece_idx(pos, to);
      if (!((pos->pt[i] ^ pos->pt[j]) & 8))
        continue;
      if (do_capture(pos, from, to, i, j)) {
#ifdef HAS_PAWNS
        if (!(is_pawn && rank18(to))) {
          alpha = max(alpha, -probe_wdl(pos, -beta, -alpha));
        } else {
          int l = i == pos->num ? j : i;
          pos->pt[l] += QUEEN - PAWN;
          for (int k = 0; k < 4; k++, pos->pt[l]--)
            if (alpha < beta)
              alpha = max(alpha, -probe_wdl(pos, -beta, -alpha));
        }
#else
        alpha = max(alpha, -probe_wdl(pos, -beta, -alpha));
#endif
        undo_capture(pos, from, to, i, j);
        if (alpha >= beta) return alpha;
      }
    }
  }

  return alpha;
}

int probe_wdl(Position *pos, int alpha, int beta)
{
  alpha = probe_capts_wdl(pos, alpha, beta);
  if (alpha >= beta) return alpha;

  int v = probe_wdl_table(pos);
  return max(alpha, v);
}

#ifdef HAS_PAWNS
// check if a drawn position without ep rights is stalemate
static bool stalemate(Position *pos)
{
  for (int i = 0; i < pos->num; i++) {
    if ((pos->pt[i] >> 3) != pos->stm)
      continue;
    int from = pos->sq[i];
    Bitboard b = piece_attacks(pos->pt[i], from, pos->occ) & pos->occ;
    while (b) {
      int to = pop_lsb(&b);
      int j = piece_idx(pos, to);
      if (!((pos->pt[i] ^ pos->pt[j]) & 8))
        continue;
      if (do_capture(pos, from, to, i, j)) {
        undo_capture(pos, from, to, i, j);
        return false;
      }
    }
    b = piece_moves(pos->pt[i], from, pos->occ);
    while (b) {
      int to = pop_lsb(&b);
      if (do_move(pos, from, to, i)) {
        undo_move(pos, from, to, i);
        return false;
      }
    }
  }
  return true;
}
#endif

int probe_dtm_loss(Position *pos, int lower, int upper);
int probe_dtm_win(Position *pos, int lower, int upper);

[[gnu::always_inline]]
inline int probe_capts_dtm(Position *pos, int alpha, int beta, const bool won)
{
  for (int i = 0; i < pos->num; i++) {
    if ((pos->pt[i] >> 3) != pos->stm)
      continue;
#ifdef HAS_PAWNS
    bool is_pawn = (pos->pt[i] & 7) == PAWN;
#endif
    int from = pos->sq[i];
    Bitboard b = piece_attacks(pos->pt[i], from, pos->occ) & pos->occ;
    while (b) {
      int to = pop_lsb(&b);
      int j = piece_idx(pos, to);
      if (!((pos->pt[i] ^ pos->pt[j]) & 8))
        continue;
      if (do_capture(pos, from, to, i, j)) {
#ifdef HAS_PAWNS
        if (!(is_pawn && rank18(to))) {
          int v =  !won ? -probe_dtm_win(pos, max(1, -beta), -alpha)
                 : probe_wdl(pos, -1, 0) >= 0 ? 10000
                 : 1 - probe_dtm_loss(pos, 1 - beta, 1 - alpha);
          beta = min(beta, v);
        } else {
          int l = i == pos->num ? j : i;
          pos->pt[l] += QUEEN - PAWN;
          for (int k = 0; k < 4; k++, pos->pt[l]--)
            if (alpha < beta) {
              int v =  !won ? -probe_dtm_win(pos, max(1, -beta), -alpha)
                     : probe_wdl(pos, -1, 0) >= 0 ? 10000
                     : 1 - probe_dtm_loss(pos, 1 - beta, 1 - alpha);
              beta = min(beta, v);
            }
        }
#else
        int v =  !won ? -probe_dtm_win(pos, max(1, -beta), -alpha)
               : probe_wdl(pos, -1, 0) >= 0 ? 10000
               : 1 - probe_dtm_loss(pos, 1 - beta, 1 - alpha);
        beta = min(beta, v);
#endif
        undo_capture(pos, from, to, i, j);
        if (beta <= alpha)
          return beta;
      }
    }
  }

  return beta;
}

// alpha < beta <= 0; return value <= 0
int probe_dtm_loss(Position *pos, int alpha, int beta)
{
  beta = min(beta, probe_capts_dtm(pos, alpha, beta, false));
  if (beta <= alpha) return beta;

  return min(beta, -probe_dtm_table(pos, false));
}

// 1 <= alpha < beta; return value >= 1
int probe_dtm_win(Position *pos, int alpha, int beta)
{
  if (beta <= alpha) return beta;

  beta = min(beta, probe_capts_dtm(pos, alpha, beta, true));
  if (beta <= alpha) return beta;

  // Try quiet moves
  for (int i = 0; i < pos->num; i++) {
    if ((pos->pt[i] >> 3) != pos->stm)
      continue;
#ifdef HAS_PAWNS
    bool is_pawn = (pos->pt[i] & 7) == PAWN;
#endif
    int from = pos->sq[i];
    Bitboard b = piece_moves(pos->pt[i], from, pos->occ);
    while (b) {
      int to = pop_lsb(&b);
      if (do_move(pos, from, to, i)) {
#ifdef HAS_PAWNS
        if (is_pawn && rank18(to)) { // pawn promotion
          pos->pt[i] += QUEEN - PAWN;
          for (int k = 0; k < 4; k++, pos->pt[i]--)
            if (alpha < beta) {
              if (probe_wdl(pos, -1, 0) < 0)
                beta = min(beta, 1 - probe_dtm_loss(pos, 1 - beta, 1 - alpha));
            }
        }
        else {
          int best_ep = 0;
          if (is_pawn && (from ^ to) == 16) { // double pawn push
            Bitboard b1 = pawn_attacks(pos->stm ^ 1, to ^ 0x08) & pos->occ;
            while (b1) {
              int s = pop_lsb(&b1);
              int k = piece_idx(pos, s);
              if ((pos->pt[i] ^ pos->pt[k]) != 8) continue;
              if (do_ep_capture(pos, s, to, k, i)) {
                int v1 = probe_wdl(pos, 0, 1);
                if (v1 > 0)
                  best_ep = max(best_ep,
                      probe_dtm_win(pos, max(1, alpha - 1), beta - 1) + 1);
                undo_ep_capture(pos, s, to, k, i);
                if (v1 <= 0)
                  goto skip; // double pawn push is not winning, so skip
              }
            }
          }
          if (best_ep == 0) {
            if (probe_wdl(pos, -1, 0) < 0)
              beta = min(beta, 1 - probe_dtm_loss(pos, 1 - beta, 1 - alpha));
          } else {
            int v1 = probe_wdl(pos, -1, 1);
            if (v1 < 0) {
              beta = min(beta, best_ep);
              if (alpha < beta)
                beta = min(beta, 1 - probe_dtm_loss(pos, 1 - beta, 1 - alpha));
            }
            else if (v1 == 0 && stalemate(pos))
              beta = min(beta, best_ep);
          }
        }
skip:
#else
        if (probe_wdl(pos, -1, 0) < 0)
          beta = min(beta, 1 - probe_dtm_loss(pos, 1 - beta, 1 - alpha));
#endif
        undo_move(pos, from, to, i);
        if (beta <= alpha)
          return beta;
      }
    }
  }

  return beta;
}
