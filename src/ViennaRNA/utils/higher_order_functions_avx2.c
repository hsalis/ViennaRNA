#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "ViennaRNA/utils/basic.h"

#if defined(__x86_64__) || defined(_M_AMD64) || defined (_M_X64)
#include <immintrin.h>
#endif


#if defined(__x86_64__) || defined(_M_AMD64) || defined (_M_X64)
static int
horizontal_min_Vec8i(__m256i x);
#endif


PUBLIC int
vrna_fun_zip_add_min_avx2(const int *e1,
                          const int *e2,
                          int       count)
{
#if defined(__x86_64__) || defined(_M_AMD64) || defined (_M_X64)
  int     i       = 0;
  int     decomp  = INF;
  __m256i inf     = _mm256_set1_epi32(INF);

  for (i = 0; i < count - 7; i += 8) {
    __m256i a    = _mm256_loadu_si256((const __m256i *)&e1[i]);
    __m256i b    = _mm256_loadu_si256((const __m256i *)&e2[i]);
    __m256i sum  = _mm256_add_epi32(a, b);
    __m256i mask = _mm256_and_si256(_mm256_cmpgt_epi32(inf, a),
                                    _mm256_cmpgt_epi32(inf, b));
    __m256i kept = _mm256_or_si256(_mm256_and_si256(mask, sum),
                                   _mm256_andnot_si256(mask, inf));

    decomp = MIN2(decomp, horizontal_min_Vec8i(kept));
  }

  for (; i < count; i++) {
    if ((e1[i] != INF) && (e2[i] != INF)) {
      const int en = e1[i] + e2[i];
      decomp = MIN2(decomp, en);
    }
  }

  return decomp;
#else
  int i;
  int decomp = INF;

  for (i = 0; i < count; i++) {
    if ((e1[i] != INF) && (e2[i] != INF)) {
      const int en = e1[i] + e2[i];
      decomp = MIN2(decomp, en);
    }
  }

  return decomp;
#endif
}


#if defined(__x86_64__) || defined(_M_AMD64) || defined (_M_X64)
static int
horizontal_min_Vec8i(__m256i x)
{
  __m128i low   = _mm256_castsi256_si128(x);
  __m128i high  = _mm256_extracti128_si256(x, 1);
  __m128i best  = _mm_min_epi32(low, high);
  __m128i min1  = _mm_shuffle_epi32(best, _MM_SHUFFLE(0, 0, 3, 2));
  __m128i min2  = _mm_min_epi32(best, min1);
  __m128i min3  = _mm_shuffle_epi32(min2, _MM_SHUFFLE(0, 0, 0, 1));
  __m128i min4  = _mm_min_epi32(min2, min3);

  return _mm_cvtsi128_si32(min4);
}
#endif
