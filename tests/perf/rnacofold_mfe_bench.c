#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ViennaRNA/fold_compound.h"
#include "ViennaRNA/mfe/global.h"
#include "ViennaRNA/model.h"
#include "ViennaRNA/intern/mfe_profile.h"


static char *
read_first_dimer(const char *path)
{
  FILE *fp;
  char line[4096];
  const char *paths[3] = { path, "../tests/data/rnacofold.small.seq", NULL };

  for (size_t p = 0; paths[p]; p++) {
    fp = fopen(paths[p], "r");
    if (!fp)
      continue;

    while (fgets(line, sizeof(line), fp)) {
      size_t len = strlen(line);

      while ((len > 0) &&
             ((line[len - 1] == '\n') || (line[len - 1] == '\r')))
        line[--len] = '\0';

      if ((len > 0) && (strspn(line, "ACGUTacgut&") == len) && strchr(line, '&')) {
        char *seq = (char *)malloc(len + 1);
        memcpy(seq, line, len + 1);
        fclose(fp);
        return seq;
      }
    }

    fclose(fp);
  }
  return NULL;
}


static char *
synthetic_dimer(size_t total_length,
                size_t split)
{
  static const char alphabet[] = "ACGU";
  unsigned int      state = 0xC0F0D123u;
  char              *seq = (char *)malloc(total_length + 2);
  size_t            pos = 0;

  for (size_t i = 0; i < total_length; i++) {
    if (i == split)
      seq[pos++] = '&';

    state = state * 1664525u + 1013904223u;
    seq[pos++] = alphabet[(state >> 24) & 0x3u];
  }

  seq[pos] = '\0';

  return seq;
}


static void
run_case(const char    *dataset,
         const char    *variant,
         const char    *sequence,
         unsigned int  repeats)
{
  vrna_md_t          md;
  uint64_t           wall_sum = 0;
  vrna_mfe_profile_t profile_sum = { 0 };
  float              mfe = 0.f;

  vrna_md_set_default(&md);
  md.min_loop_size = 0;

  for (unsigned int r = 0; r < repeats; r++) {
    vrna_fold_compound_t      *fc;
    char                      *structure;
    uint64_t                  t0;
    const vrna_mfe_profile_t  *profile;

    fc = vrna_fold_compound(sequence, &md, VRNA_OPTION_DEFAULT);
    structure = (char *)malloc(strlen(sequence) + 1);

    vrna_mfe_profile_enable(1);
    vrna_mfe_profile_reset();
    t0 = vrna_mfe_profile_now();
    mfe = vrna_mfe_dimer(fc, structure);
    wall_sum += (vrna_mfe_profile_now() - t0);
    profile = vrna_mfe_profile_get();

    profile_sum.fill_arrays_ns += profile->fill_arrays_ns;
    profile_sum.decompose_pair_ns += profile->decompose_pair_ns;
    profile_sum.internal_ns += profile->internal_ns;
    profile_sum.multibranch_ns += profile->multibranch_ns;
    profile_sum.exterior_ns += profile->exterior_ns;
    profile_sum.backtrack_ns += profile->backtrack_ns;
    profile_sum.allocations += profile->allocations;

    free(structure);
    vrna_fold_compound_free(fc);
  }

  printf("%s,%s,%zu,%u,%.2f,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu\n",
         dataset,
         variant,
         strlen(sequence) - 1,
         repeats,
         mfe,
         (unsigned long long)(wall_sum / repeats),
         (unsigned long long)(profile_sum.fill_arrays_ns / repeats),
         (unsigned long long)(profile_sum.decompose_pair_ns / repeats),
         (unsigned long long)(profile_sum.internal_ns / repeats),
         (unsigned long long)(profile_sum.multibranch_ns / repeats),
         (unsigned long long)(profile_sum.exterior_ns / repeats),
         (unsigned long long)(profile_sum.backtrack_ns / repeats),
         (unsigned long long)(profile_sum.allocations / repeats));
}


int
main(int argc,
     char *argv[])
{
  unsigned int repeats = 3;
  char         *small;

  if (argc > 1)
    repeats = (unsigned int)strtoul(argv[1], NULL, 10);

  printf("dataset,variant,total_length,repeats,mfe,wall_ns,fill_arrays_ns,decompose_pair_ns,internal_ns,multibranch_ns,exterior_ns,backtrack_ns,allocations\n");

  small = read_first_dimer("tests/data/rnacofold.small.seq");
  if (small) {
    run_case("rnacofold.small.seq", "default", small, repeats);
    free(small);
  }

  for (size_t len_i = 0; len_i < 6; len_i++) {
    static const size_t lengths[] = { 50, 100, 150, 200, 500, 1000 };
    char                label[64];
    char                *equal_split;
    char                *ratio_split;
    size_t              total = lengths[len_i];

    equal_split = synthetic_dimer(total, total / 2);
    snprintf(label, sizeof(label), "synthetic.%zu", total);
    run_case(label, "split_1_1", equal_split, repeats);
    free(equal_split);

    ratio_split = synthetic_dimer(total, total / 4);
    run_case(label, "split_1_3", ratio_split, repeats);
    free(ratio_split);
  }

  return 0;
}
