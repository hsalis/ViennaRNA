#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ViennaRNA/fold_compound.h"
#include "ViennaRNA/model.h"
#include "ViennaRNA/subopt/wuchty.h"


static char *
synthetic_sequence(size_t length)
{
  static const char alphabet[] = "ACGU";
  unsigned int      state = 0xFACE1234u;
  char              *seq = (char *)malloc(length + 1);

  for (size_t i = 0; i < length; i++) {
    state = state * 1664525u + 1013904223u;
    seq[i] = alphabet[(state >> 24) & 0x3u];
  }

  seq[length] = '\0';
  return seq;
}


static void
free_subopt(vrna_subopt_solution_t *sol)
{
  if (!sol)
    return;

  for (vrna_subopt_solution_t *s = sol; s->structure; s++)
    free(s->structure);

  free(sol);
}


static void
run_case(size_t        length,
         int           delta,
         unsigned int  repeats)
{
  vrna_md_t md;
  char      *seq = synthetic_sequence(length);
  uint64_t  wall_sum = 0;
  size_t    n_sol = 0;

  vrna_md_set_default(&md);
  md.uniq_ML = 1;

  for (unsigned int r = 0; r < repeats; r++) {
    vrna_fold_compound_t    *fc;
    vrna_subopt_solution_t  *sol;
    struct timespec         ts0, ts1;

    fc = vrna_fold_compound(seq, &md, VRNA_OPTION_DEFAULT);
    clock_gettime(CLOCK_MONOTONIC, &ts0);
    sol = vrna_subopt(fc, delta, 0, NULL);
    clock_gettime(CLOCK_MONOTONIC, &ts1);
    wall_sum += ((uint64_t)(ts1.tv_sec - ts0.tv_sec) * 1000000000ULL) +
                (uint64_t)(ts1.tv_nsec - ts0.tv_nsec);

    n_sol = 0;
    for (vrna_subopt_solution_t *s = sol; s && s->structure; s++)
      n_sol++;

    free_subopt(sol);
    vrna_fold_compound_free(fc);
  }

  printf("%zu,%d,%u,%zu,%llu\n",
         length,
         delta,
         repeats,
         n_sol,
         (unsigned long long)(wall_sum / repeats));

  free(seq);
}


int
main(int argc,
     char *argv[])
{
  unsigned int repeats = 1;
  static const size_t lengths[] = { 50, 100, 150, 200 };
  static const int    deltas[] = { 5, 10 };

  if (argc > 1)
    repeats = (unsigned int)strtoul(argv[1], NULL, 10);

  printf("length,delta_dcal,repeats,solutions,wall_ns\n");

  for (size_t i = 0; i < sizeof(lengths) / sizeof(lengths[0]); i++)
    for (size_t j = 0; j < sizeof(deltas) / sizeof(deltas[0]); j++)
      run_case(lengths[i], deltas[j], repeats);

  return 0;
}
