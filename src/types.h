#ifndef TYPES_H
#define TYPES_H

#include <inttypes.h>
#include <stddef.h>

#include "defs.h"

typedef uint64_t Bitboard;
typedef uint16_t Move;

typedef uint8_t u8;
typedef uint16_t u16;

enum { WHITE = 0, BLACK };

enum { PAWN = 1, KNIGHT, BISHOP, ROOK, QUEEN, KING };

enum {
  WPAWN = 1, WKNIGHT, WBISHOP, WROOK, WQUEEN, WKING,
  BPAWN = 9, BKNIGHT, BBISHOP, BROOK, BQUEEN, BKING
};

#define INLINE [[gnu::always_inline]] static inline
#define NOINLINE [[gnu::noinline]]

#define assume(x) do { if (!(x)) unreachable(); } while (0)
#define likely(x) __builtin_expect(!!(x), 1)
#define unlikely(x) __builtin_expect(!!(x), 0)

#define PASTER(x,y) x##_##y
#define EVALUATOR(x,y) PASTER(x,y)

#undef max
#undef min

#define MAX(T) INLINE T max_##T(T a, T b) { return a > b ? a : b; }
MAX(int)
MAX(int64_t)
MAX(uint8_t)
MAX(uint16_t)
MAX(uint32_t)
MAX(uint64_t)
#undef MAX

#define MIN(T) INLINE T min_##T(T a, T b) { return a < b ? a : b; }
MIN(int)
MIN(int64_t)
MIN(uint8_t)
MIN(uint16_t)
MIN(uint32_t)
MIN(uint64_t)
#undef MIN

#define TEMPLATE(F,a,b) _Generic((a), \
         int: F##_int,             \
     int64_t: F##_int64_t,         \
     uint8_t: F##_uint8_t,         \
    uint16_t: F##_uint16_t,        \
    uint32_t: F##_uint32_t,        \
    uint64_t: F##_uint64_t         \
) (a,b)

#define max(a,b) TEMPLATE(max,a,b)
#define min(a,b) TEMPLATE(min,a,b)

#define SWAP(T) INLINE void Swap_##T(T *a, T *b) { T tmp = *a; *a = *b; *b = tmp; }
SWAP(uint8_t)
SWAP(uint16_t)
SWAP(int)
SWAP(int64_t)
SWAP(uint32_t)
SWAP(uint64_t)
#undef SWAP

#define TEMPLATE2(F,a,b) _Generic(*(a), \
         int: F##_int,               \
     int64_t: F##_int64_t,           \
     uint8_t: F##_uint8_t,           \
    uint16_t: F##_uint16_t,          \
    uint32_t: F##_uint32_t,          \
    uint64_t: F##_uint64_t           \
) (a,b)

#define Swap(a,b) TEMPLATE2(Swap,&(a),&(b))

#endif
