/*
  Copyright (c) 2011-2017, 2025, 2026 Ronald de Man

  This file is distributed under the terms of the GNU GPL, version 2.
*/

#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdint.h>
#include <string.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <x86intrin.h>

#include "defs.h"
#include "movegen.h"
#include "probe.h"
#include "threads.h"
#include "types.h"
#include "util.h"

#define TB_HASHBITS 11
#define MAX_PIECE 650 // MAX_PIECE and MAX_PAWN to be increased for 8-piece TBs
#define MAX_PAWN 861

#define WDLSUFFIX ".rtbw"
#define DTMSUFFIX ".rtbm"

const char *suffix[] = { ".rtbw", ".rtbm", ".rtbz" };
uint32_t magic[] = { 0x5d23e871, 0x88ac504b, 0xa50c66d7 };

struct HashEntry {
  uint64_t key;
  struct BaseEntry *ptr;
};

static LOCK_T mutex;
static char *path_string = nullptr;
static int num_paths = 0;
static char **paths;

static int num_piece, num_pawn;

static int num_wdl, num_dtm;

static struct PieceEntry *piece_entry;
static struct PawnEntry *pawn_entry;
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

size_t Binomial[8][64];

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

  bool has_pawns = pcs[WPAWN] || pcs[BPAWN];

  struct BaseEntry *be = has_pawns ? &pawn_entry[num_pawn++].be
                                   : &piece_entry[num_piece++].be;

  be->has_pawns = has_pawns;
  be->key = key;
  be->symmetric = key == key2;
  be->num = 0;
  for (int i = 0; i < 16; i++)
    be->num += pcs[i];

  num_wdl++;
  num_dtm += be->has_dtm = test_tb(str, DTM);

  for (int t = 0; t < 3; t++)
    be->ready[t] = false;

  if (!be->has_pawns) {
    int j = 0;
    for (int i = 0; i < 16; i++)
      if (pcs[i] == 1) j++;
    be->kk_enc = j == 2; // for suicide a bit more may have to be done
  } else {
    be->pawns[0] = pcs[WPAWN];
    be->pawns[1] = pcs[BPAWN];
    if (pcs[BPAWN] && (!pcs[WPAWN] || pcs[WPAWN] > pcs[BPAWN]))
      Swap(be->pawns[0], be->pawns[1]);
  }

  add_to_hash(be, key);
  if (key != key2)
    add_to_hash(be, key2);
}

#define PCE_E(x) ((struct PieceEntry *)(x))
#define PWN_E(x) ((struct PawnEntry *)(x))

INLINE int num_tables(struct BaseEntry *be, const int t)
{
  return be->has_pawns ? t == DTM ? 6 : 4 : 1;
}

INLINE struct EncInfo *first_ei(struct BaseEntry *be, const int t)
{
  return  be->has_pawns
        ? &PWN_E(be)->ei[t == WDL ? 0 : t == DTM ? 8 : 20]
        : &PCE_E(be)->ei[t == WDL ? 0 : t == DTM ? 2 :  4];
}

static void free_tb_entry(struct BaseEntry *be)
{
  for (int t = 0; t < 3; t++) {
    if (be->ready[t]) {
      unmap_file(be->data[t], be->mapping[t]);
      int num = num_tables(be, t);
      struct EncInfo *ei = first_ei(be, t);
      for (int i = 0; i < num; i++) {
        free(ei[i].precomp);
        if (t != DTZ)
          free(ei[num + i].precomp);
      }
      be->ready[t] = false;
    }
  }
}

void free_tablebases(void)
{
  for (int i = 0; i < num_piece; i++)
    free_tb_entry((struct BaseEntry *)&piece_entry[i]);
  for (int i = 0; i < num_pawn; i++)
    free_tb_entry((struct BaseEntry *)&pawn_entry[i]);
  free(piece_entry);
  free(pawn_entry);
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

  num_piece = num_pawn = 0;

  piece_entry = malloc(MAX_PIECE * sizeof *piece_entry);
  pawn_entry  = malloc(MAX_PAWN  * sizeof *pawn_entry );

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

  printf("Found %d WDL and %d DTM tablebase files.\n", num_wdl, num_dtm);
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
  for (int i = 1; i < 7; i++)
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
}

INLINE int leading_pawn(uint8_t *p, struct BaseEntry *be, const int enc)
{
  for (int i = 1; i < be->pawns[0]; i++)
    if (PawnFlip[enc-1][p[0]] > PawnFlip[enc-1][p[i]])
      Swap(p[0], p[i]);

  return enc == FILE_ENC ? FileToFile[p[0] & 7] : (p[0] - 8) >> 3;
}

INLINE void sort_squares(int n, uint8_t *p)
{
  for (int i = 0; i < n; i++)
    for (int j = i + 1; j < n; j++)
      if (p[i] > p[j])
        Swap(p[i], p[j]);
}

INLINE void sort_squares_mapped(int n, uint8_t *p, const uint8_t *map)
{
  for (int i = 0; i < n; i++)
    for (int j = i + 1; j < n; j++)
      if (map[p[i]] > map[p[j]])
        Swap(p[i], p[j]);
}

INLINE int rank_among_free(uint8_t sq, Bitboard occ)
{
  return sq - popcnt(occ & ((1ULL << sq) - 1));
}

INLINE size_t encode(uint8_t *p, struct EncInfo *ei, struct BaseEntry *be,
    const int enc)
{
  int n = be->num;
  size_t idx;
  Bitboard occ = 0;
  int k;

  if (p[0] & 0x04)
    for (int i = 0; i < 8; i++)
      p[i] ^= 0x07;

  if (enc == PIECE_ENC) {
    if (p[0] & 0x20)
      for (int i = 0; i < 8; i++)
        p[i] ^= 0x38;

    for (int i = 0; i < n; i++)
      if (OffDiag[p[i]]) {
        if (OffDiag[p[i]] > 0 && i < (be->kk_enc ? 2 : 3))
          for (int j = 0; j < n; j++)
            p[j] = FlipDiag[p[j]];
        break;
      }

    if (be->kk_enc) {
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
  } else { /* RANK_ENC */
    for (int i = 1; i < be->pawns[0]; i++)
      for (int j = i + 1; j < be->pawns[0]; j++)
        if (PawnTwist[enc-1][p[i]] < PawnTwist[enc-1][p[j]])
          Swap(p[i], p[j]);

    k = be->pawns[0];
    idx = PawnIdx[enc-1][k-1][PawnFlip[enc-1][p[0]]];
    for (int i = 1; i < k; i++)
      idx += Binomial[k-i][PawnTwist[enc-1][p[i]]];
    idx *= ei->factor[0];

    for (int i = 0; i < k; i++)
      occ |= bit(p[i]);

    // Pawns of other color
    if (be->pawns[1]) {
      int t = k + be->pawns[1];
      sort_squares(be->pawns[1], &p[k]);
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
    struct BaseEntry *be)
{
  return encode(p, ei, be, PIECE_ENC);
}

NOINLINE size_t encode_pawn_f(uint8_t *p, struct EncInfo *ei,
    struct BaseEntry *be)
{
  return encode(p, ei, be, FILE_ENC);
}

NOINLINE size_t encode_pawn_r(uint8_t *p, struct EncInfo *ei,
    struct BaseEntry *be)
{
  return encode(p, ei, be, RANK_ENC);
}

// Count number of placements of k like pieces on n squares
size_t subfactor(size_t k, size_t n)
{
  size_t f = 1;

  for (size_t i = 0; i < k; i++)
    f = (f * (n - i)) / (i + 1);

  return f;
}

size_t set_enc_info(struct EncInfo *ei, struct BaseEntry *be,
    const uint8_t *tb, int shift, int fr, const int enc)
{
  bool more_pawns = enc != PIECE_ENC && be->pawns[1] > 0;

  for (int i = 0; i < be->num; i++) {
    ei->pieces[i] = (tb[i + 1 + more_pawns] >> shift) & 0x0f;
    ei->norm[i] = 0;
  }

  int order = (tb[0] >> shift) & 0x0f;
  int order2 = more_pawns ? (tb[1] >> shift) & 0x0f : 0x0f;

  int k = ei->norm[0] =  enc != PIECE_ENC ? be->pawns[0]
                       : be->kk_enc ? 2 : 3;

  if (more_pawns) {
    ei->norm[k] = be->pawns[1];
    k += ei->norm[k];
  }

  for (int i = k; i < be->num; i += ei->norm[i])
    for (int j = i; j < be->num && ei->pieces[j] == ei->pieces[i]; j++)
      ei->norm[i]++;

  int n = 64 - k;
  size_t f = 1;

  for (int i = 0; k < be->num || i == order || i == order2; i++) {
    if (i == order) {
      ei->factor[0] = f;
      f *=  enc == FILE_ENC ? PawnFactorFile[ei->norm[0] - 1][fr]
          : enc == RANK_ENC ? PawnFactorRank[ei->norm[0] - 1][fr]
          : be->kk_enc ? 462 : 31332;
    } else if (i == order2) {
      ei->factor[ei->norm[0]] = f;
      f *= subfactor(ei->norm[ei->norm[0]], 48 - ei->norm[0]);
    } else {
      ei->factor[k] = f;
      f *= subfactor(ei->norm[k], n);
      n -= ei->norm[k];
      k += ei->norm[k];
    }
  }

  return f;
}

size_t set_dec_info(struct DecInfo *di, struct BaseEntry *be, uint8_t *pcs,
    uint8_t *type_perm, int order, int order2, int fr, const int enc)
{
  bool more_pawns = enc != PIECE_ENC && be->pawns[1] > 0;

  for (int i = 0; i < be->num; i++)
    di->norm[i] = 0;

  int k = di->norm[0] =  enc != PIECE_ENC ? be->pawns[0]
                       : be->kk_enc ? 2 : 3;

  if (more_pawns) {
    di->norm[k] = be->pawns[1];
    k += di->norm[k];
  }

  int i, n = 64 - k;
  size_t f = 1;

  for (i = 0; k < be->num || i == order || i == order2; i++) {
    if (i == order) {
      di->factor[0] = f;
      f *=  enc == FILE_ENC ? PawnFactorFile[di->norm[0] - 1][fr]
          : enc == RANK_ENC ? PawnFactorRank[di->norm[0] - 1][fr]
          : be->kk_enc ? 462 : 31332;
      if (enc == PIECE_ENC)
        i += di->norm[0] - 1;
    } else if (i == order2) {
      di->factor[di->norm[0]] = f;
      f *= subfactor(di->norm[di->norm[0]], 48 - di->norm[0]);
    } else {
      di->factor[k] = f;
      di->norm[k] = pcs[type_perm[i]];
      f *= subfactor(di->norm[k], n);
      n -= di->norm[k];
      k += di->norm[k];
    }
  }

  uint8_t tmp[8];
  k = di->norm[0];
  if (order2 < 0x0f)
    k += di->norm[k];
  for (i = 0; k < be->num || i == order || i == order2; i++)
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

INLINE void unrank_binomial(uint64_t idx, int n, uint8_t *p,
    const uint8_t *const map, Bitboard *occ)
{
  Bitboard b = ~*occ;
  for (int i = n - 1; i > 0; i--) {
    int r = i;
    // FIXME: improve this
    while (true) {
      uint64_t f = Binomial[i][r];
      if (f > idx) break;
      idx -= f;
      r++;
    }
    r = map ? map[r] : r;
    Bitboard b1 = _pdep_u64(1ULL << r, b);
    p[i] = lsb(b1);
    *occ |= b1;
  }
  idx = map ? map[idx] : idx;
  Bitboard b1 = _pdep_u64(1ULL << idx, b);
  p[0] = lsb(b1);
  *occ |= b1;
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
    struct DecInfo *di, struct BaseEntry *be, const int fr, const int enc)
{
  uint32_t q, r;
  int i;
  Bitboard occ = 0;
  constexpr Bitboard DIAGONAL_BB = 0x8040201008040201ULL;

  q = sub[0];
  if (enc == PIECE_ENC) {
    switch (be->kk_enc) {

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
  else { /* FILE_ENC, RANK_ENC */
    const int num = enc == FILE_ENC ? 6 : 4;
    int t = be->pawns[0] - 1;
    for (i = 0; i < num - 1; i++)
      if (PawnIdx[enc-1][t][num * fr + i + 1] > q) break;
    q -= PawnIdx[enc-1][t][num * fr + i];
    p[0] = place_empty(InvPawnFlip[enc-1][num * fr + i], &occ);

    if (t > 0)
      unrank_binomial(q, t, p + 1, InvPawnTwist[enc-1], &occ);

    i = be->pawns[0];

    if (be->pawns[1] > 0) {
      q = sub[be->pawns[0]];
      unrank_binomial(q, be->pawns[1], p + i, PawnMap, &occ);
      i += be->pawns[1];
    }

assert(i >= di->norm[0]);
    for (; i < di->norm[0]; i++)
      p[i] = place_skip(sub[i], &occ);
  }

  for (; i < be->num;) {
    q = sub[i];
    unrank_binomial(q, di->norm[i], p + i, nullptr, &occ);
    i += di->norm[i];
  }
}

INLINE void decode(uint64_t idx, uint8_t *restrict p, struct DecInfo *di,
    struct BaseEntry *be, const int fr, const int enc)
{
  uint32_t sub[TB_PIECES];
  int i;

  // TODO: convert into multiplications, e.g. using libdivide.
  // https://github.com/ridiculousfish/libdivide
  for (i = 0; di->factor[di->order[i]] != 1; i++) {
    uint32_t q = idx / di->factor[di->order[i]];
    idx -= q * di->factor[di->order[i]];
    sub[di->order[i]] = q;
  }
  sub[di->order[i]] = idx;

  decode_helper(sub, p, di, be, fr, enc);
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
    struct BaseEntry *be, const int fr, const int enc)
{
  decode_helper(sub, p, di, be, fr, enc);

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
    struct DecInfo *di, struct BaseEntry *be)
{
  decode(idx, p, di, be, 0, PIECE_ENC);
}

NOINLINE void decode_piece_iter(uint32_t *sub, uint8_t *restrict p,
    struct DecInfo *di, struct BaseEntry *be)
{
  decode_iter(sub, p, di, be, 0, PIECE_ENC);
}

NOINLINE void decode_pawn_r(uint64_t idx, uint8_t *restrict p,
    struct DecInfo *di, struct BaseEntry *be, int rank)
{
  decode(idx, p, di, be, rank, RANK_ENC);
}

NOINLINE void decode_pawn_r_iter(uint32_t *sub, uint8_t *restrict p,
    struct DecInfo *di, struct BaseEntry *be, int rank)
{
  decode_iter(sub, p, di, be, rank, RANK_ENC);
}


/******* start of actual probing and decompression code *******/

static void calc_symlen(struct PairsData *d, uint32_t s, bool *tmp)
{
  const uint8_t *w = d->sympat + 3 * s;
  uint32_t s2 = (w[2] << 4) | (w[1] >> 4);
  if (s2 == 0x0fff)
    d->symlen[s] = 0;
  else {
    uint32_t s1 = ((w[1] & 0xf) << 8) | w[0];
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
   size_t *size, uint8_t *flags, int type)
{
  struct PairsData *d;
  const uint8_t *data = *ptr;

  *flags = data[0];
  if (data[0] & 0x80) {
    d = malloc(sizeof(struct PairsData));
    d->compr_type = 2;
    d->const_value[0] = data[1];
    d->const_value[1] = 0;
    *ptr = data + 2;
    size[0] = size[1] = size[2] = 0;
    d->tb_size = tb_size;
    return d;
  }

  d = (data[0] & 0x40) ? setup_rans(ptr) : setup_huffman(ptr);

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

NOINLINE bool init_table(struct BaseEntry *be, const char *str,
    int type, bool use_paths)
{
  const uint8_t *data = map_tb(str, type, &be->mapping[type], use_paths);
  if (!data) return false;

  if (read_le_u32(data) != magic[type]) {
    fprintf(stderr, "Corrupted table.\n");
    unmap_file(data, be->mapping[type]);
    return false;
  }

  be->data[type] = data;

  bool split = type != DTZ && (data[4] & 0x01);
  if (type == DTM)
    be->dtm_losses_only = data[4] & 0x04;

  data += 5;

  size_t tb_size[6][2];
  int num = num_tables(be, type);
  struct EncInfo *ei = first_ei(be, type);
  int enc = !be->has_pawns ? PIECE_ENC : type != DTM ? FILE_ENC : RANK_ENC;

  for (int t = 0; t < num; t++) {
    tb_size[t][0] = set_enc_info(&ei[t], be, data, 0, t, enc);
    if (split)
      tb_size[t][1] = set_enc_info(&ei[num + t], be, data, 4, t, enc);
    data += be->num + 1 + (be->has_pawns && be->pawns[1]);
  }
  data += (uintptr_t)data & 1;

  size_t size[6][2][3];
  for (int t = 0; t < num; t++) {
    uint8_t flags;
    ei[t].precomp = setup_pairs(&data, tb_size[t][0], size[t][0], &flags, type);
    if (type == DTZ) {
      if (!be->has_pawns)
        PCE_E(be)->dtz_flags = flags;
      else
        PWN_E(be)->dtz_flags[t] = flags;
    }
    if (split)
      ei[num + t].precomp = setup_pairs(&data, tb_size[t][1], size[t][1], &flags, type);
    else if (type != DTZ)
      ei[num + t].precomp = nullptr;
  }

  if (type == DTM && !be->dtm_losses_only) {
    uint16_t *map = (uint16_t *)data;
    *(be->has_pawns ? &PWN_E(be)->map_dtm : &PCE_E(be)->map_dtm) = map;
    uint16_t (*map_idx)[2][2] = be->has_pawns ? &PWN_E(be)->map_dtm_idx[0]
                                              : &PCE_E(be)->map_dtm_idx;
    for (int t = 0; t < num; t++) {
      for (int i = 0; i < 2; i++) {
        map_idx[t][0][i] = (uint16_t *)data + 1 - map;
        data += 2 + 2 * read_le_u16(data);
      }
      if (split) {
        for (int i = 0; i < 2; i++) {
          map_idx[t][1][i] = (uint16_t *)data + 1 - map;
          data += 2 + 2 * read_le_u16(data);
        }
      }
    }
  }

  if (type == DTZ) {
    const void *map = data;
    *(be->has_pawns ? &PWN_E(be)->map_dtz : &PCE_E(be)->map_dtz) = map;
    uint16_t (*map_idx)[4] = be->has_pawns ? &PWN_E(be)->map_dtz_idx[0]
                                           : &PCE_E(be)->map_dtz_idx;
    uint8_t *flags = be->has_pawns ? &PWN_E(be)->dtz_flags[0]
                                   : &PCE_E(be)->dtz_flags;
    for (int t = 0; t < num; t++) {
      if (flags[t] & 2) {
        if (!(flags[t] & 16)) {
          for (int i = 0; i < 4; i++) {
            map_idx[t][i] = data + 1 - (uint8_t *)map;
            data += 1 + data[0];
          }
        } else {
          data += (uintptr_t)data & 0x01;
          for (int i = 0; i < 4; i++) {
            map_idx[t][i] = (uint16_t *)data + 1 - (uint16_t *)map;
            data += 2 + 2 * read_le_u16(data);
          }
        }
      }
    }
    data += (uintptr_t)data & 0x01;
  }

  for (int t = 0; t < num; t++) {
    ei[t].precomp->index_table = data;
    data += size[t][0][0];
    if (split) {
      ei[num + t].precomp->index_table = data;
      data += size[t][1][0];
    }
  }

  for (int t = 0; t < num; t++) {
    ei[t].precomp->size_table = (uint16_t *)data;
    data += size[t][0][1];
    if (split) {
      ei[num + t].precomp->size_table = (uint16_t *)data;
      data += size[t][1][1];
    }
  }

  for (int t = 0; t < num; t++) {
    data = (uint8_t *)(((uintptr_t)data + 0x3f) & ~0x3f);
    ei[t].precomp->data = data;
    data += size[t][0][2];
    if (split) {
      data = (uint8_t *)(((uintptr_t)data + 0x3f) & ~0x3f);
      ei[num + t].precomp->data = data;
      data += size[t][1][2];
    }
  }

  if (type == DTM && be->has_pawns) {
    int count[16];
    for (int i = 0; i < 16; i++)
      count[i] = 0;
    for (int i = 0; i < be->num; i++)
      count[ei[0].pieces[i]]++;
    PWN_E(be)->dtm_switched =
        material_key_from_counts(count, count + 8) != be->key;
  }

  return true;
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

  // TODO: expand Re-Pair dictionary and create pointers from sym to pattern
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

INLINE void list_squares(Position *pos, uint8_t *pt, bool flip, uint8_t *p,
    int n)
{
  for (int i = 0; i < n; ) {
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

  struct BaseEntry *be = tb_hash[hash_idx].ptr;
  if ((type == DTM && !be->has_dtm) || (type == DTZ && !be->has_dtz))
    probe_failed(pos, type);

  // Use double-checked locking to reduce locking overhead
  if (!atomic_load_explicit(&be->ready[type], memory_order_acquire)) {
    LOCK(mutex);
    if (!atomic_load_explicit(&be->ready[type], memory_order_relaxed)) {
      char str[16];
      create_material_string(pos, str, be->key != key);
      if (!init_table(be, str, type, true))
        probe_failed(pos, type);
      atomic_store_explicit(&be->ready[type], true, memory_order_release);
    }
    UNLOCK(mutex);
  }

  bool bside, flip;
  if (!be->symmetric) {
    flip = key != be->key;
    bside = (pos->stm == WHITE) == flip;
    if (type == DTM && be->has_pawns && PWN_E(be)->dtm_switched) {
      flip = !flip;
      bside = !bside;
    }
  } else {
    flip = pos->stm != WHITE;
    bside = false;
  }

  struct EncInfo *ei = first_ei(be, type);
  uint8_t p[8];
  size_t idx;
  int t = 0;
  uint8_t flags;

  if (!be->has_pawns) {
    if (type == DTZ) {
      flags = PCE_E(be)->dtz_flags;
      if ((flags & 1) != bside && !be->symmetric) return INT_MIN;
    }
    ei = type != DTZ ? &ei[bside] : ei;
    list_squares(pos, ei->pieces, flip, p, be->num);
    idx = encode_piece(p, ei, be);
  } else {
    list_squares(pos, ei->pieces, flip, p, be->pawns[0]);
    t = leading_pawn(p, be, type != DTM ? FILE_ENC : RANK_ENC);
    if (type == DTZ) {
      flags = PWN_E(be)->dtz_flags[t];
      if ((flags & 1) != bside && !be->symmetric) return INT_MIN;
    }
    ei =  type == WDL ? &ei[t + 4 * bside]
        : type == DTM ? &ei[t + 6 * bside] : &ei[t];
    list_squares(pos, ei->pieces + be->pawns[0], flip, p + be->pawns[0],
        be->num - be->pawns[0]);
    idx = type != DTM ? encode_pawn_f(p, ei, be) : encode_pawn_r(p, ei, be);
  }

  const uint8_t *w = decompress_pairs(ei->precomp, idx);

  if (type == WDL) return (int)w[0] - 2;

  int v = w[0] + ((w[1] & 0x0f) << 8);

  if (type == DTM) {
    if (!be->dtm_losses_only)
      v =  from_le_u16(be->has_pawns
         ? PWN_E(be)->map_dtm[PWN_E(be)->map_dtm_idx[t][bside][s] + v]
         : PCE_E(be)->map_dtm[PCE_E(be)->map_dtm_idx[bside][s] + v]);
  } else {
#if 0
    if (flags & 2) {
      int m = WdlToMap[s + 2];
      if (!(flags & 16))
        v =  be->has_pawns
           ? ((uint8_t *)PWN_E(be)->map_dtz)[PWN_E(be)->map_dtz_idx[t][m] + v]
           : ((uint8_t *)PCE_E(be)->map_dtz)[PCE_E(be)->map_dtz_idx[m] + v];
      else
        v =  from_le_u16(be->has_pawns
           ? ((uint16_t *)PWN_E(be)->map_dtz)[PWN_E(be)->map_dtz_idx[t][m] + v]
           : ((uint16_t *)PCE_E(be)->map_dtz)[PCE_E(be)->map_dtz_idx[m] + v]);
    }
    if (!(flags & PAFlags[s + 2]) || (s & 1))
      v *= 2;
#endif
  }

  return v;
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
