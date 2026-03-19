#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "ViennaRNA/fold_compound.h"
#include "ViennaRNA/mfe/global.h"
#include "ViennaRNA/model.h"
#include "ViennaRNA/params/basic.h"
#include "ViennaRNA/partfunc/global.h"
#include "ViennaRNA/sampling/basic.h"


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
  double        mfe_sum;
  size_t        char_sum;
  unsigned int  r;

  wall_sum = 0;
  mfe_sum  = 0.;
  char_sum = 0;

  for (r = 0; r < repeats; r++) {
    char                  *seq, *mfe_structure;
    char                  **samples;
    vrna_md_t             md;
    vrna_fold_compound_t  *fc;
    double                mfe;
    uint64_t              t0;
    unsigned int          i;

    seq = synthetic_sequence(length, 0x5A11u + r * 31u);
    vrna_md_set_default(&md);
    md.uniq_ML      = 1;
    md.compute_bpp  = 0;

    fc = vrna_fold_compound(seq, &md, VRNA_OPTION_DEFAULT);
    mfe_structure = (char *)malloc(length + 1);

    mfe = vrna_mfe(fc, mfe_structure);
    vrna_exp_params_rescale(fc, &mfe);

    t0 = now_ns();
    (void)vrna_pf(fc, NULL);
    samples = vrna_pbacktrack_num(fc, 20, VRNA_PBACKTRACK_DEFAULT);
    wall_sum += now_ns() - t0;

    mfe_sum += mfe;
    if (samples) {
      for (i = 0; samples[i]; i++) {
        char_sum += strlen(samples[i]);
        free(samples[i]);
      }
      free(samples);
    }

    free(mfe_structure);
    vrna_fold_compound_free(fc);
    free(seq);
  }

  printf("%zu,%u,%llu,%.6f,%zu\n",
         length,
         repeats,
         (unsigned long long)(wall_sum / repeats),
         mfe_sum / (double)repeats,
         char_sum / repeats);
}


int
main(int argc,
     char *argv[])
{
  unsigned int  repeats;
  size_t        i;
  const size_t  lengths[] = {
    100, 200, 500
  };

  repeats = 3;
  if (argc > 1)
    repeats = (unsigned int)strtoul(argv[1], NULL, 10);

  printf("length,repeats,wall_ns,mfe,total_sample_chars\n");

  for (i = 0; i < sizeof(lengths) / sizeof(lengths[0]); i++)
    run_case(lengths[i], repeats);

  return 0;
}
