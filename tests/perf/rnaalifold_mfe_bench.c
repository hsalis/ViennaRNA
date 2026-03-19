#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "ViennaRNA/fold_compound.h"
#include "ViennaRNA/io/file_formats_msa.h"
#include "ViennaRNA/mfe/global.h"
#include "ViennaRNA/model.h"
#include "ViennaRNA/intern/mfe_profile.h"


static void
free_alignment(char **names,
               char **aln,
               char *id,
               char *structure)
{
  if (names) {
    for (size_t i = 0; names[i]; i++)
      free(names[i]);
    free(names);
  }

  if (aln) {
    for (size_t i = 0; aln[i]; i++)
      free(aln[i]);
    free(aln);
  }

  free(id);
  free(structure);
}


static char **
synthetic_alignment(size_t length,
                    unsigned int n_seq)
{
  static const char alphabet[] = "ACGU";
  unsigned int      state = 0x1234ABCDu;
  char              **aln = (char **)calloc(n_seq + 1, sizeof(char *));

  for (unsigned int s = 0; s < n_seq; s++) {
    aln[s] = (char *)malloc(length + 1);
    for (size_t i = 0; i < length; i++) {
      state = state * 1664525u + 1013904223u + s;
      aln[s][i] = alphabet[(state >> 24) & 0x3u];
    }
    aln[s][length] = '\0';
  }

  return aln;
}


static void
run_case(const char    *dataset,
         const char    **alignment,
         unsigned int  repeats)
{
  vrna_md_t          md;
  uint64_t           wall_sum = 0;
  vrna_mfe_profile_t profile_sum = { 0 };
  float              mfe = 0.f;
  size_t             len = strlen(alignment[0]);

  vrna_md_set_default(&md);

  for (unsigned int r = 0; r < repeats; r++) {
    vrna_fold_compound_t      *fc;
    char                      *structure;
    uint64_t                  t0;
    const vrna_mfe_profile_t  *profile;

    fc = vrna_fold_compound_comparative(alignment, &md, VRNA_OPTION_DEFAULT);
    structure = (char *)malloc(len + 1);

    vrna_mfe_profile_enable(1);
    vrna_mfe_profile_reset();
    t0 = vrna_mfe_profile_now();
    mfe = vrna_mfe(fc, structure);
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

  printf("%s,%zu,%u,%.2f,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu\n",
         dataset,
         len,
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

  if (argc > 1)
    repeats = (unsigned int)strtoul(argv[1], NULL, 10);

  printf("dataset,length,repeats,mfe,wall_ns,fill_arrays_ns,decompose_pair_ns,internal_ns,multibranch_ns,exterior_ns,backtrack_ns,allocations\n");

  {
    static const char *paths[] = {
      "tests/data/alignment_fasta.fa",
      "../tests/data/alignment_fasta.fa",
      "../../tests/data/alignment_fasta.fa",
      "../../../tests/data/alignment_fasta.fa",
      "tests/data/alignment_clustal.aln",
      "../tests/data/alignment_clustal.aln",
      "../../tests/data/alignment_clustal.aln",
      "../../../tests/data/alignment_clustal.aln",
      "tests/data/rfam_seed_selected.stk",
      "../tests/data/rfam_seed_selected.stk",
      "../../tests/data/rfam_seed_selected.stk",
      "../../../tests/data/rfam_seed_selected.stk",
      NULL
    };

    for (size_t p = 0; paths[p]; p++) {
      char **names = NULL, **aln = NULL, *id = NULL, *structure = NULL;
      if (access(paths[p], R_OK) != 0)
        continue;

      int  n_seq = vrna_file_msa_read(paths[p],
                                      &names,
                                      &aln,
                                      &id,
                                      &structure,
                                      VRNA_FILE_FORMAT_MSA_DEFAULT);

      if (n_seq > 0)
        run_case(paths[p], (const char **)aln, repeats);

      free_alignment(names, aln, id, structure);
    }
  }

  for (size_t len_i = 0; len_i < 3; len_i++) {
    static const size_t lengths[] = { 100, 300, 500 };
    char                label[64];
    char                **aln = synthetic_alignment(lengths[len_i], 4);

    snprintf(label, sizeof(label), "synthetic.%zu", lengths[len_i]);
    run_case(label, (const char **)aln, repeats);
    free_alignment(NULL, aln, NULL, NULL);
  }

  return 0;
}
