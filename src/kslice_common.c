#if defined(__AVX512F__) || defined(__AVX2__)
#include <immintrin.h>
#endif

static constexpr size_t MIN_STREAMING_STORES_SIZE = (size_t)1024*1024;

static void *work_p, *work_q;

static void set_worker(struct ThreadData *thread)
{
  uint8_t *restrict p = work_p;

  size_t begin = thread->begin << 6, end = thread->end << 6;

#ifdef __AVX512F__

  if (end - begin >= MIN_STREAMING_STORES_SIZE) {
    __m512i ones = _mm512_set1_epi64(-1);
    for (size_t idx = begin; idx < end; idx += 64)
      _mm512_stream_si512((__m512i *)(p + idx), ones);
    return;
  }

#elifdef __AVX2__

  if (end - begin >= MIN_STREAMING_STORES_SIZE) {
    __m256i ones = _mm256_set1_epi32(-1);
    for (idx = begin; idx < end; idx += 64) {
      _mm256_stream_si256((__m256i *)(p + idx), ones);
      _mm256_stream_si256((__m256i *)(p + idx + 32), ones);
    }
    return;
  }

#endif

  memset(p + begin, 0xff, end - begin);
}

static void clear_worker(struct ThreadData *thread)
{
  uint8_t *restrict p = work_p;

  size_t begin = thread->begin << 6, end = thread->end << 6;

#ifdef __AVX512F__

  if (end - begin >= MIN_STREAMING_STORES_SIZE) {
    __m512i z = _mm512_setzero_si512();
    for (size_t idx = begin; idx < end; idx += 64)
      _mm512_stream_si512((__m512i *)(p + idx), z);
    return;
  }

#elifdef __AVX2__

  if (thread->end - thread->begin >= MIN_STREAMING_STORES / 64) {
    __m256i z = _mm256_setzero_si256();
    for (size_t idx = begin; idx < end; idx += 64) {
      _mm256_stream_si256((__m256i *)(p + idx), z);
      _mm256_stream_si256((__m256i *)(p + idx + 32), z);
    }
    return;
  }

#endif

  memset(p + begin, 0x00, end - begin);
}

static void not_worker(struct ThreadData *thread)
{
  uint64_t *restrict p = work_p;

  uint64_t idx = thread->begin << 3;
  uint64_t end = thread->end << 3;

#ifdef __AVX512F__

  __m512i ones = _mm512_set1_epi64(-1);
  for (; idx < end; idx += 8) {
    __m512i a = _mm512_load_si512(p + idx);
    _mm512_store_si512(p + idx, _mm512_xor_si512(a, ones));
  }

#elifdef __AVX2__

  __m256i ones = _mm256_set1_epi32(-1);
  for (; idx < end; idx += 4) {
    __m256i a = _mm256_load_si256(p + idx);
    _mm256_store_si256(p + idx, _mm256_xor_si256(a, ones));
  }

#else

  for (; idx < end; idx++)
    p[idx] = ~p[idx];

#endif
}

static void or_worker(struct ThreadData *thread)
{
  uint64_t *restrict p = work_p;
  uint64_t *restrict q = work_q;

  uint64_t idx = thread->begin << 3;
  uint64_t end = thread->end << 3;

#ifdef __AVX512F__

  for (; idx < end; idx += 8) {
    __m512i a = _mm512_load_si512(p + idx);
    __m512i b = _mm512_load_si512(q + idx);
    _mm512_store_si512(p + idx, _mm512_or_si512(a, b));
  }

#elifdef __AVX2__

  for (; idx < end; idx += 4) {
    __m256i a = _mm256_load_si256(p + idx);
    __m256i b = _mm256_load_si256(q + idx);
    _mm256_store_si256(p + idx, _mm256_or_si256(a, b));
  }

#else

  for (; idx < end; idx++)
    p[idx] |= q[idx];

#endif
}

static void or_not_worker(struct ThreadData *thread)
{
  uint64_t *restrict p = work_p;
  uint64_t *restrict q = work_q;

  uint64_t idx = thread->begin << 3;
  uint64_t end = thread->end << 3;

#ifdef __AVX512F__

  for (; idx < end; idx += 8) {
    __m512i a = _mm512_load_si512(p + idx);
    __m512i b = _mm512_load_si512(q + idx);
    _mm512_store_si512(p + idx, _mm512_ternarylogic_epi64(a, b, b, 0xF3));
  }

#elifdef __AVX2__

  const __m256i ones = _mm256_set1_epi64x(-1);
  for (; idx < end; idx += 4) {
    __m256i a = _mm256_load_si256(p + idx);
    __m256i b = _mm256_xor_si256(_mm256_load_si256(q + idx), ones);
    _mm256_store_si256(p + idx, _mm256_or_si256(a, b));
  }

#else

  for (; idx < end; idx++)
    p[idx] |= ~q[idx];

#endif
}

static void and_worker(struct ThreadData *thread)
{
  uint64_t *restrict p = work_p;
  uint64_t *restrict q = work_q;

  uint64_t idx = thread->begin << 3;
  uint64_t end = thread->end << 3;

#ifdef __AVX512F__

  for (; idx < end; idx += 8) {
    __m512i a = _mm512_load_si512(p + idx);
    __m512i b = _mm512_load_si512(q + idx);
    _mm512_store_si512(p + idx, _mm512_and_si512(a, b));
  }

#elifdef __AVX2__

  for (; idx < end; idx += 4) {
    __m256i a = _mm256_load_si256(p + idx);
    __m256i b = _mm256_load_si256(q + idx);
    _mm256_store_si256(p + idx, _mm256_and_si256(a, b));
  }

#else

  for (; idx < end; idx++)
    p[idx] &= q[idx];

#endif
}

static void and_not_worker(struct ThreadData *thread)
{
  uint64_t *restrict p = work_p;
  uint64_t *restrict q = work_q;

  uint64_t idx = thread->begin << 3;
  uint64_t end = thread->end << 3;

#ifdef __AVX512F__

  for (; idx < end; idx += 8) {
    __m512i a = _mm512_load_si512(p + idx);
    __m512i b = _mm512_load_si512(q + idx);
    _mm512_store_si512(p + idx, _mm512_andnot_epi64(b, a));
  }

#elifdef __AVX2__

  for (; idx < end; idx += 4) {
    __m256i a = _mm256_load_si256(p + idx);
    __m256i b = _mm256_load_si256(q + idx);
    _mm256_store_si256(p + idx, _mm256_andnot_si256(b, a));
  }

#else

  for (; idx < end; idx++)
    p[idx] &= ~q[idx];

#endif
}

static void not_and_worker(struct ThreadData *thread)
{
  uint64_t *restrict p = work_p;
  uint64_t *restrict q = work_q;

  uint64_t idx = thread->begin << 3;
  uint64_t end = thread->end << 3;

#ifdef __AVX512F__

  for (; idx < end; idx += 8) {
    __m512i a = _mm512_load_si512(p + idx);
    __m512i b = _mm512_load_si512(q + idx);
    _mm512_store_si512(p + idx, _mm512_andnot_epi64(a, b));
  }

#elifdef __AVX2__

  for (; idx < end; idx += 4) {
    __m256i a = _mm256_load_si256(p + idx);
    __m256i b = _mm256_load_si256(q + idx);
    _mm256_store_si256(p + idx, _mm256_andnot_si256(a, b));
  }

#else

  for (; idx < end; idx++)
    p[idx] = ~p[idx] & q[idx];

#endif
}

void nor_worker(struct ThreadData *thread)
{
  uint64_t *restrict p = work_p;
  uint64_t *restrict q = work_q;

  uint64_t idx = thread->begin << 3;
  uint64_t end = thread->end << 3;

#ifdef __AVX512F__

  for (; idx < end; idx += 8) {
    __m512i a = _mm512_load_si512(p + idx);
    __m512i b = _mm512_load_si512(q + idx);
    _mm512_store_si512(p + idx, _mm512_ternarylogic_epi64(a, b, b, 0x05));
  }

#elifdef __AVX2__

  const __m256i ones = _mm256_set1_epi64x(-1);
  for (; idx < end; idx += 4) {
    __m256i a = _mm256_load_si256(p + idx);
    __m256i b = _mm256_load_si256(q + idx);
    _mm256_store_si256(p + idx, _mm256_xor_si256(_mm256_or_si256(a, b), ones);
  }

#else

  for (; idx < end; idx++)
    p[idx] = ~(p[idx] | q[idx]);

#endif
}

void split_worker(struct ThreadData *thread)
{
  uint64_t *restrict p = work_p;
  uint64_t *restrict q = work_q;

  uint64_t idx = thread->begin << 3;
  uint64_t end = thread->end << 3;

#ifdef __AVX512F__

  for (; idx < end; idx += 8) {
    __m512i a = _mm512_load_si512(p + idx);
    __m512i b = _mm512_load_si512(q + idx);
    _mm512_store_si512(p + idx, _mm512_andnot_si512(b, a));
    _mm512_store_si512(p + idx, _mm512_and_si512(b, a));
  }

#elifdef __AVX2__

  for (; idx < end; idx += 4) {
    __m256i a = _mm256_load_si256(p + idx);
    __m256i b = _mm256_load_si256(q + idx);
    _mm256_store_si256(p + idx, _mm256_andnot_si256(b, a));
    _mm256_store_si256(p + idx, _mm256_and_si256(b, a));
  }

#else

  for (; idx < end; idx++) {
    uint64_t a = p[idx];
    uint64_t b = q[idx];
    p[idx] = a & ~b;
    q[idx] = a & b;
  }

#endif
}

INLINE void clear_tail(void *p, size_t num_bits, size_t num_words)
{
  uint64_t *restrict q = p;
  size_t idx = num_bits >> 6;
  int r = num_bits & 63;
  if (r)
    q[idx++] &= (1ULL << r) - 1;

  for (; idx < num_words; idx++)
    q[idx] = 0;
}

static void count_worker(struct ThreadData *thread)
{
  uint64_t cnt = 0, *restrict p = work_p;

  uint64_t idx = thread->begin << 3;
  uint64_t end = thread->end << 3;

#ifdef __AVX512VPOPCNTDQ__

  __m512i sum = _mm512_setzero_si512();
  for (; idx < end; idx += 8) {
    __m512i a = _mm512_load_si512(p + idx);
    sum = _mm512_add_epi64(sum, _mm512_popcnt_epi64(a));
  }
  cnt = _mm512_reduce_add_epi64(sum);

#else

  for (; idx < end; idx++)
    cnt += popcnt(p[idx]);

#endif

  thread->cnt += cnt;
}
