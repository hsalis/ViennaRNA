#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "ViennaRNA/fold_compound.h"
#include "ViennaRNA/landscape/findpath.h"
#include "ViennaRNA/landscape/paths.h"
#include "ViennaRNA/mfe/global.h"


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
  double        saddle_sum;
  size_t        path_sum;
  unsigned int  r;

  wall_sum    = 0;
  saddle_sum  = 0.;
  path_sum    = 0;

  for (r = 0; r < repeats; r++) {
    char                  *seq, *s1, *s2;
    vrna_fold_compound_t  *fc;
    vrna_path_t           *path, *ptr;
    uint64_t              t0;
    int                   saddle;

    seq = synthetic_sequence(length, 0xF11Du + r * 13u);
    fc = vrna_fold_compound(seq, NULL, VRNA_OPTION_DEFAULT);

    s1 = (char *)malloc(length + 1);
    (void)vrna_mfe(fc, s1);

    s2 = (char *)malloc(length + 1);
    memset(s2, '.', length);
    s2[length] = '\0';

    t0 = now_ns();
    saddle = vrna_path_findpath_saddle(fc, s1, s2, (int)(length / 2));
    path = vrna_path_findpath(fc, s1, s2, (int)(length / 2));
    wall_sum += now_ns() - t0;

    saddle_sum += (double)saddle;
    if (path) {
      for (ptr = path; ptr->s; ptr++)
        path_sum++;
      vrna_path_free(path);
    }

    free(s2);
    free(s1);
    vrna_fold_compound_free(fc);
    free(seq);
  }

  printf("%zu,%u,%llu,%.6f,%zu\n",
         length,
         repeats,
         (unsigned long long)(wall_sum / repeats),
         saddle_sum / (double)repeats,
         path_sum / repeats);
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

  printf("length,repeats,wall_ns,saddle_dcal,mean_path_len\n");

  for (i = 0; i < sizeof(lengths) / sizeof(lengths[0]); i++)
    run_case(lengths[i], repeats);

  return 0;
}
