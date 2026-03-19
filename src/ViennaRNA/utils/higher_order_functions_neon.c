#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "ViennaRNA/utils/basic.h"

#include <arm_neon.h>


PUBLIC int
vrna_fun_zip_add_min_neon(const int *e1,
                          const int *e2,
                          int       count)
{
  int       i       = 0;
  int       decomp  = INF;
  int32x4_t inf     = vdupq_n_s32(INF);
  int32x4_t best    = inf;

  for (i = 0; i < count - 3; i += 4) {
    int32x4_t a      = vld1q_s32(&e1[i]);
    int32x4_t b      = vld1q_s32(&e2[i]);
    uint32x4_t mask  = vandq_u32(vcltq_s32(a, inf), vcltq_s32(b, inf));
    int32x4_t sum    = vaddq_s32(a, b);
    int32x4_t kept   = vbslq_s32(mask, sum, inf);

    best = vminq_s32(best, kept);
  }

  int32_t lanes[4];
  vst1q_s32(lanes, best);
  decomp = MIN2(MIN2(lanes[0], lanes[1]), MIN2(lanes[2], lanes[3]));

  for (; i < count; i++) {
    if ((e1[i] != INF) && (e2[i] != INF)) {
      const int en = e1[i] + e2[i];
      decomp = MIN2(decomp, en);
    }
  }

  return decomp;
}
