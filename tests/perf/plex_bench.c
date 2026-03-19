#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <fcntl.h>

#include "ViennaRNA/plex.h"


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
  double        energy_sum;
  size_t        hit_sum;
  unsigned int  r;

  wall_sum    = 0;
  energy_sum  = 0.;
  hit_sum     = 0;

  for (r = 0; r < repeats; r++) {
    char      *s1, *s2;
    duplexT   **hits;
    uint64_t  t0;
    int       count;
    int       saved_stdout, devnull;

    s1 = synthetic_sequence(len1, 0xA1E0u + r * 11u);
    s2 = synthetic_sequence(len2, 0xA1E1u + r * 19u);

    fflush(stdout);
    saved_stdout = dup(STDOUT_FILENO);
    devnull = open("/dev/null", O_WRONLY);
    if ((saved_stdout != -1) && (devnull != -1))
      (void)dup2(devnull, STDOUT_FILENO);

    t0 = now_ns();
    hits = Lduplexfold(s1, s2, -5, 0, 5, 0, 1, 6, 6, 3, 3);
    wall_sum += now_ns() - t0;

    fflush(stdout);
    if (saved_stdout != -1) {
      (void)dup2(saved_stdout, STDOUT_FILENO);
      close(saved_stdout);
    }
    if (devnull != -1)
      close(devnull);

    count = (hits) ? arraySize(hits) : 0;
    hit_sum += (size_t)count;
    if ((hits) && (count > 0))
      energy_sum += hits[0][0].energy;

    if (hits)
      freeDuplexT(hits);

    free(s1);
    free(s2);
  }

  printf("%zu,%zu,%u,%llu,%.6f,%zu\n",
         len1,
         len2,
         repeats,
         (unsigned long long)(wall_sum / repeats),
         energy_sum / (double)repeats,
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

  printf("len1,len2,repeats,wall_ns,first_hit_energy,mean_hits\n");

  for (i = 0; i < sizeof(lengths) / sizeof(lengths[0]); i++)
    for (j = 0; j < sizeof(lengths) / sizeof(lengths[0]); j++)
      run_case(lengths[i], lengths[j], repeats);

  return 0;
}
