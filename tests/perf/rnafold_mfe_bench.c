#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ViennaRNA/fold_compound.h"
#include "ViennaRNA/mfe/global.h"
#include "ViennaRNA/model.h"
#include "ViennaRNA/constraints/basic.h"
#include "ViennaRNA/constraints/soft.h"
#include "ViennaRNA/intern/mfe_profile.h"


typedef struct {
  int           noLP;
  int           dangles;
  int           max_bp_span;
  int           circ;
  int           gquad;
  unsigned int  add_hard_constraints;
  unsigned int  add_soft_constraints;
} bench_options_t;


static char *
read_first_sequence(const char *path)
{
  FILE          *fp;
  char          line[4096];
  size_t        cap, len, first_line_len;
  char          *seq;
  int           have_header, multi_line_record;

  fp = fopen(path, "r");
  if (!fp)
    return NULL;

  cap = 1024;
  len = 0;
  first_line_len = 0;
  have_header = 0;
  multi_line_record = 0;
  seq = (char *)malloc(cap);
  seq[0] = '\0';

  while (fgets(line, sizeof(line), fp)) {
    size_t line_len = strlen(line);

    if (line[0] == '>') {
      have_header = 1;
      if (len > 0)
        break;
      continue;
    }

    while ((line_len > 0) &&
           ((line[line_len - 1] == '\n') || (line[line_len - 1] == '\r')))
      line[--line_len] = '\0';

    if (line_len == 0)
      break;

    if (strspn(line, "ACGUTacgut") != line_len) {
      if (len > 0)
        break;
      continue;
    }

    if (len > 0) {
      if (!have_header && !multi_line_record) {
        if (line_len != first_line_len)
          multi_line_record = 1;
        else
          break;
      }
    } else {
      first_line_len = line_len;
    }

    if (len + line_len + 1 > cap) {
      while (len + line_len + 1 > cap)
        cap *= 2;
      seq = (char *)realloc(seq, cap);
    }

    memcpy(seq + len, line, line_len);
    len += line_len;
    seq[len] = '\0';
  }

  fclose(fp);

  if (len == 0) {
    free(seq);
    return NULL;
  }

  return seq;
}


static char *
synthetic_sequence(size_t length)
{
  static const char alphabet[] = "ACGU";
  unsigned int      state = 0xC0DEFACEu;
  char              *seq = (char *)malloc(length + 1);

  for (size_t i = 0; i < length; i++) {
    state = state * 1664525u + 1013904223u;
    seq[i] = alphabet[(state >> 24) & 0x3u];
  }

  seq[length] = '\0';

  return seq;
}


static void
apply_hard_constraints(vrna_fold_compound_t *fc)
{
  unsigned int  n;
  char          *constraint;

  n = fc->length;
  constraint = (char *)malloc(n + 1);
  memset(constraint, '.', n);
  constraint[n] = '\0';

  if (n > 10) {
    constraint[1]      = 'x';
    constraint[4]      = '(';
    constraint[n - 6]  = ')';
  }

  vrna_constraints_add(fc, constraint, VRNA_CONSTRAINT_DB_DEFAULT);
  free(constraint);
}


static void
apply_soft_constraints(vrna_fold_compound_t *fc)
{
  unsigned int n = fc->length;

  if (n > 6)
    vrna_sc_add_up(fc, 3, -1.0, VRNA_OPTION_DEFAULT);

  if (n > 12)
    vrna_sc_add_bp(fc, 5, n - 4, -0.5, VRNA_OPTION_DEFAULT);
}


static void
run_case(const char              *dataset,
         const char              *variant,
         const char              *sequence,
         const bench_options_t   *opt,
         unsigned int            repeats)
{
  vrna_md_t               md;
  vrna_fold_compound_t    *fc;
  vrna_mfe_profile_t      profile_sum = {
    0
  };
  uint64_t                wall_sum = 0;
  double                  mfe = 0.;
  size_t                  len = strlen(sequence);

  vrna_md_set_default(&md);
  md.noLP        = opt->noLP;
  md.dangles     = opt->dangles;
  md.max_bp_span = opt->max_bp_span;
  md.circ        = opt->circ;
  md.gquad       = opt->gquad;

  for (unsigned int r = 0; r < repeats; r++) {
    char                    *structure;
    uint64_t                t0;
    const vrna_mfe_profile_t *profile;

    fc = vrna_fold_compound(sequence, &md, VRNA_OPTION_DEFAULT);
    if (!fc) {
      fprintf(stderr, "failed to create fold compound for %s/%s\n", dataset, variant);
      return;
    }

    if (opt->add_hard_constraints)
      apply_hard_constraints(fc);

    if (opt->add_soft_constraints)
      apply_soft_constraints(fc);

    structure = (char *)malloc(len + 1);

    vrna_mfe_profile_enable(1);
    vrna_mfe_profile_reset();
    t0 = vrna_mfe_profile_now();
    mfe = vrna_mfe(fc, structure);
    wall_sum += (vrna_mfe_profile_now() - t0);
    profile = vrna_mfe_profile_get();

    profile_sum.fill_arrays_ns    += profile->fill_arrays_ns;
    profile_sum.decompose_pair_ns += profile->decompose_pair_ns;
    profile_sum.internal_ns       += profile->internal_ns;
    profile_sum.multibranch_ns    += profile->multibranch_ns;
    profile_sum.exterior_ns       += profile->exterior_ns;
    profile_sum.backtrack_ns      += profile->backtrack_ns;
    profile_sum.allocations       += profile->allocations;

    free(structure);
    vrna_fold_compound_free(fc);
  }

  printf("%s,%s,%zu,%u,%.2f,%llu,%.3f,%llu,%llu,%llu,%llu,%llu,%llu,%llu\n",
         dataset,
         variant,
         len,
         repeats,
         mfe,
         (unsigned long long)(wall_sum / repeats),
         (double)(wall_sum / repeats) / ((double)len * (double)len * (double)len),
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
  unsigned int          repeats = 3;
  char                  *small = NULL, *large = NULL, *gquad = NULL;
  const bench_options_t default_opt = {
    .dangles = 2,
    .max_bp_span = -1
  };

  if (argc > 1)
    repeats = (unsigned int)strtoul(argv[1], NULL, 10);

  small = read_first_sequence("tests/data/rnafold.small.seq");
  large = read_first_sequence("tests/data/rnafold.seq");
  gquad = read_first_sequence("tests/data/rnafold.gquad.fa");

  printf("dataset,variant,length,repeats,mfe,wall_ns,ns_per_n3,fill_arrays_ns,decompose_pair_ns,internal_ns,multibranch_ns,exterior_ns,backtrack_ns,allocations\n");

  if (small) {
    bench_options_t opt = default_opt;
    run_case("rnafold.small.seq", "default", small, &opt, repeats);

    opt.noLP = 1;
    run_case("rnafold.small.seq", "noLP", small, &opt, repeats);

    opt = default_opt;
    opt.dangles = 0;
    run_case("rnafold.small.seq", "d0", small, &opt, repeats);

    opt.dangles = 2;
    run_case("rnafold.small.seq", "d2", small, &opt, repeats);

    opt = default_opt;
    opt.max_bp_span = 30;
    run_case("rnafold.small.seq", "maxBPspan30", small, &opt, repeats);

    opt = default_opt;
    opt.circ = 1;
    run_case("rnafold.small.seq", "circ", small, &opt, repeats);

    opt = default_opt;
    opt.add_hard_constraints = 1;
    run_case("rnafold.small.seq", "hard", small, &opt, repeats);

    opt = default_opt;
    opt.add_soft_constraints = 1;
    run_case("rnafold.small.seq", "soft", small, &opt, repeats);
  }

  if (large) {
    run_case("rnafold.seq", "default", large, &default_opt, repeats);
  }

  if (gquad) {
    bench_options_t opt = default_opt;
    opt.gquad = 1;
    run_case("rnafold.gquad.fa", "gquad", gquad, &opt, repeats);
  }

  for (size_t idx = 0; idx < 6; idx++) {
    static const size_t lengths[] = {
      50, 100, 150, 200, 500, 1000
    };
    char *seq = synthetic_sequence(lengths[idx]);

    run_case("synthetic", "default", seq, &default_opt, repeats);
    free(seq);
  }

  free(small);
  free(large);
  free(gquad);

  return 0;
}
