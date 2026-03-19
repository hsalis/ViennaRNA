#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "ViennaRNA/heat_capacity.h"
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
  double        last_hc_sum;
  size_t        points_sum;
  unsigned int  r;

  wall_sum    = 0;
  last_hc_sum = 0.;
  points_sum  = 0;

  for (r = 0; r < repeats; r++) {
    char                    *seq;
    vrna_heat_capacity_t    *result;
    uint64_t                t0;
    size_t                  n;

    seq = synthetic_sequence(length, 0xEA77u + r * 23u);

    t0 = now_ns();
    result = vrna_heat_capacity_simple(seq, 0., 100., 5., 2);
    wall_sum += now_ns() - t0;

    for (n = 0; result[n].temperature > (-K0 - 0.5f); n++)
      ;

    if (n > 0)
      last_hc_sum += result[n - 1].heat_capacity;
    points_sum += n;

    free(result);
    free(seq);
  }

  printf("%zu,%u,%llu,%.6f,%zu\n",
         length,
         repeats,
         (unsigned long long)(wall_sum / repeats),
         last_hc_sum / (double)repeats,
         points_sum / repeats);
}


int
main(int argc,
     char *argv[])
{
  unsigned int  repeats;
  size_t        i;
  const size_t  lengths[] = {
    50, 100, 200, 500
  };

  repeats = 3;
  if (argc > 1)
    repeats = (unsigned int)strtoul(argv[1], NULL, 10);

  printf("length,repeats,wall_ns,last_heat_capacity,points\n");

  for (i = 0; i < sizeof(lengths) / sizeof(lengths[0]); i++)
    run_case(lengths[i], repeats);

  return 0;
}
