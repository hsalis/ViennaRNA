#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "ViennaRNA/inverse/basic.h"


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


static char *
target_structure(size_t length)
{
  char    *structure;
  size_t  pairs, i;

  structure = (char *)malloc(length + 1);
  memset(structure, '.', length);
  structure[length] = '\0';

  pairs = length / 4;
  for (i = 0; i < pairs; i++) {
    structure[i] = '(';
    structure[length - 1 - i] = ')';
  }

  return structure;
}


static void
run_case(size_t        length,
         unsigned int  repeats)
{
  uint64_t      wall_sum;
  double        score_sum;
  size_t        gc_sum;
  unsigned int  r;

  wall_sum   = 0;
  score_sum  = 0.;
  gc_sum     = 0;

  for (r = 0; r < repeats; r++) {
    char      *start, *target;
    uint64_t  t0;
    float     score;
    size_t    i;

    start = synthetic_sequence(length, 0x1A90u + r * 37u);
    target = target_structure(length);

    t0 = now_ns();
    score = inverse_fold(start, target);
    wall_sum += now_ns() - t0;

    score_sum += score;
    for (i = 0; i < length; i++)
      if ((start[i] == 'G') || (start[i] == 'C'))
        gc_sum++;

    free(target);
    free(start);
  }

  printf("%zu,%u,%llu,%.6f,%zu\n",
         length,
         repeats,
         (unsigned long long)(wall_sum / repeats),
         score_sum / (double)repeats,
         gc_sum / repeats);
}


int
main(int argc,
     char *argv[])
{
  unsigned int  repeats;
  size_t        i;
  const size_t  lengths[] = {
    50, 100, 150
  };

  repeats = 1;
  if (argc > 1)
    repeats = (unsigned int)strtoul(argv[1], NULL, 10);

  printf("length,repeats,wall_ns,mean_score,mean_gc_count\n");

  for (i = 0; i < sizeof(lengths) / sizeof(lengths[0]); i++)
    run_case(lengths[i], repeats);

  return 0;
}
