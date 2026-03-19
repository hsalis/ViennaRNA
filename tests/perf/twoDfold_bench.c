#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "ViennaRNA/2Dfold.h"
#include "ViennaRNA/fold_compound.h"
#include "ViennaRNA/mfe/global.h"
#include "ViennaRNA/params/constants.h"


static uint64_t
now_ns(void)
{
  struct timespec ts;

  clock_gettime(CLOCK_MONOTONIC, &ts);

  return ((uint64_t)ts.tv_sec * 1000000000ULL) + (uint64_t)ts.tv_nsec;
}


static char *
synthetic_sequence(size_t        length,
                   unsigned int  seed)
{
  static const char alphabet[] = "ACGU";
  char              *seq;
  unsigned int      state;
  size_t            i;

  seq   = (char *)malloc(length + 1);
  state = seed ^ (unsigned int)length;

  for (i = 0; i < length; i++) {
    state = state * 1664525u + 1013904223u;
    seq[i] = alphabet[(state >> 24) & 0x3u];
  }

  seq[length] = '\0';

  return seq;
}


static void
run_case(size_t        length,
         unsigned int  repeats)
{
  uint64_t      wall_sum;
  double        energy_sum;
  size_t        classes_sum;
  unsigned int  r;

  wall_sum    = 0;
  energy_sum  = 0.;
  classes_sum = 0;

  for (r = 0; r < repeats; r++) {
    char                  *seq, *ref1, *ref2;
    vrna_fold_compound_t  *fc1, *fc2;
    vrna_sol_TwoD_t       *sol;
    uint64_t              t0;
    size_t                i;

    seq = synthetic_sequence(length, 0x2DF0u + r * 41u);

    fc1 = vrna_fold_compound(seq, NULL, VRNA_OPTION_DEFAULT);
    ref1 = (char *)malloc(length + 1);
    (void)vrna_mfe(fc1, ref1);
    vrna_fold_compound_free(fc1);

    ref2 = (char *)malloc(length + 1);
    memset(ref2, '.', length);
    ref2[length] = '\0';

    fc2 = vrna_fold_compound_TwoD(seq,
                                  ref1,
                                  ref2,
                                  NULL,
                                  VRNA_OPTION_MFE);
    t0 = now_ns();
    sol = vrna_mfe_TwoD(fc2, -1, -1);
    wall_sum += now_ns() - t0;

    if (sol) {
      for (i = 0; sol[i].k != INF; i++) {
        energy_sum += sol[i].en;
        free(sol[i].s);
      }
      classes_sum += i;
      free(sol);
    }
    vrna_fold_compound_free(fc2);
    free(ref2);
    free(ref1);
    free(seq);
  }

  printf("%zu,%u,%llu,%.6f,%zu\n",
         length,
         repeats,
         (unsigned long long)(wall_sum / repeats),
         energy_sum / (double)repeats,
         classes_sum / repeats);
}


int
main(int argc,
     char *argv[])
{
  unsigned int  repeats;
  size_t        i;
  const size_t  lengths[] = {
    50, 100
  };

  repeats = 1;
  if (argc > 1)
    repeats = (unsigned int)strtoul(argv[1], NULL, 10);

  printf("length,repeats,wall_ns,mean_energy_sum,mean_classes\n");

  for (i = 0; i < sizeof(lengths) / sizeof(lengths[0]); i++)
    run_case(lengths[i], repeats);

  return 0;
}
