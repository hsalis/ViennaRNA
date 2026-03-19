#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "ViennaRNA/part_func_up.h"
#include "ViennaRNA/partfunc/global.h"


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
run_case(size_t        len1,
         size_t        len2,
         unsigned int  repeats)
{
  uint64_t      wall_sum;
  double        interaction_sum;
  size_t        hit_sum;
  unsigned int  r;

  wall_sum        = 0;
  interaction_sum = 0.;
  hit_sum         = 0;

  for (r = 0; r < repeats; r++) {
    char        *s1, *s2;
    pu_contrib  *p1, *p2;
    interact    *ix;
    uint64_t    t0;
    int         max_w;

    s1 = synthetic_sequence(len1, 0x0A01u + r * 7u);
    s2 = synthetic_sequence(len2, 0x0A02u + r * 17u);
    max_w = (int)((len2 < 30) ? len2 : 30);

    (void)pf_fold(s1, NULL);
    p1 = pf_unstru(s1, max_w);
    free_pf_arrays();

    (void)pf_fold(s2, NULL);
    p2 = pf_unstru(s2, max_w);
    free_pf_arrays();

    t0 = now_ns();
    ix = pf_interact(s1, s2, p1, p2, max_w, NULL, 0, 0);
    wall_sum += now_ns() - t0;

    if (ix) {
      interaction_sum += ix->Gikjl;
      hit_sum += (size_t)((ix->i > ix->k) ? 1 : 0);
      free_interact(ix);
    }

    free_pu_contrib(p1);
    free_pu_contrib(p2);
    free(s1);
    free(s2);
  }

  printf("%zu,%zu,%u,%llu,%.6f,%zu\n",
         len1,
         len2,
         repeats,
         (unsigned long long)(wall_sum / repeats),
         interaction_sum / (double)repeats,
         hit_sum / repeats);
}


int
main(int argc,
     char *argv[])
{
  unsigned int        repeats;
  const size_t        lengths[] = {
    50, 100, 200
  };
  size_t              i, j;

  repeats = 3;
  if (argc > 1)
    repeats = (unsigned int)strtoul(argv[1], NULL, 10);

  printf("len1,len2,repeats,wall_ns,mean_Gikjl,mean_hits\n");

  for (i = 0; i < sizeof(lengths) / sizeof(lengths[0]); i++)
    for (j = 0; j < sizeof(lengths) / sizeof(lengths[0]); j++)
      run_case(lengths[i], lengths[j], repeats);

  return 0;
}
