#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "ViennaRNA/fold_compound.h"
#include "ViennaRNA/mfe/local.h"
#include "ViennaRNA/model.h"
#include "ViennaRNA/partfunc/local.h"


static char *
synthetic_sequence(size_t length)
{
  static const char alphabet[] = "ACGU";
  unsigned int      state = 0xA5A5A5A5u;
  char              *seq = (char *)malloc(length + 1);

  for (size_t i = 0; i < length; i++) {
    state = state * 1103515245u + 12345u;
    seq[i] = alphabet[(state >> 16) & 0x3u];
  }

  seq[length] = '\0';
  return seq;
}


static void
noop_mfe_cb(unsigned int start,
            unsigned int end,
            const char   *structure,
            float        en,
            void         *data)
{
  size_t *count = (size_t *)data;
  (void)start;
  (void)end;
  (void)structure;
  (void)en;
  (*count)++;
}


static void
noop_pf_cb(FLT_OR_DBL   *pr,
           int          pr_size,
           int          i,
           int          max,
           unsigned int type,
           void         *data)
{
  size_t *count = (size_t *)data;
  (void)pr;
  (void)pr_size;
  (void)i;
  (void)max;
  (void)type;
  (*count)++;
}


static uint64_t
elapsed_ns(struct timespec a,
           struct timespec b)
{
  return ((uint64_t)(b.tv_sec - a.tv_sec) * 1000000000ULL) +
         (uint64_t)(b.tv_nsec - a.tv_nsec);
}


static void
run_case(size_t        length,
         unsigned int  repeats)
{
  vrna_md_t  md;
  char       *seq = synthetic_sequence(length);
  uint64_t   mfe_sum = 0;
  uint64_t   pf_sum = 0;
  size_t     mfe_hits = 0;
  size_t     pf_hits = 0;

  vrna_md_set_default(&md);
  md.window_size = 200;
  md.max_bp_span = 150;

  for (unsigned int r = 0; r < repeats; r++) {
    vrna_fold_compound_t  *fc_mfe;
    vrna_fold_compound_t  *fc_pf;
    struct timespec       t0, t1;
    size_t                local_mfe_hits = 0;
    size_t                local_pf_hits = 0;

    fc_mfe = vrna_fold_compound(seq, &md, VRNA_OPTION_MFE | VRNA_OPTION_WINDOW);
    clock_gettime(CLOCK_MONOTONIC, &t0);
    (void)vrna_mfe_window_cb(fc_mfe, noop_mfe_cb, &local_mfe_hits);
    clock_gettime(CLOCK_MONOTONIC, &t1);
    mfe_sum += elapsed_ns(t0, t1);
    mfe_hits = local_mfe_hits;
    vrna_fold_compound_free(fc_mfe);

    fc_pf = vrna_fold_compound(seq, &md, VRNA_OPTION_PF | VRNA_OPTION_WINDOW);
    clock_gettime(CLOCK_MONOTONIC, &t0);
    (void)vrna_probs_window(fc_pf,
                            md.window_size,
                            VRNA_PROBS_WINDOW_BPP,
                            noop_pf_cb,
                            &local_pf_hits);
    clock_gettime(CLOCK_MONOTONIC, &t1);
    pf_sum += elapsed_ns(t0, t1);
    pf_hits = local_pf_hits;
    vrna_fold_compound_free(fc_pf);
  }

  printf("%zu,%u,%zu,%llu,%zu,%llu\n",
         length,
         repeats,
         mfe_hits,
         (unsigned long long)(mfe_sum / repeats),
         pf_hits,
         (unsigned long long)(pf_sum / repeats));

  free(seq);
}


int
main(int argc,
     char *argv[])
{
  unsigned int repeats = 1;
  static const size_t lengths[] = { 500, 1000, 2000 };

  if (argc > 1)
    repeats = (unsigned int)strtoul(argv[1], NULL, 10);

  printf("length,repeats,mfe_hits,mfe_wall_ns,pf_hits,pf_wall_ns\n");

  for (size_t i = 0; i < sizeof(lengths) / sizeof(lengths[0]); i++)
    run_case(lengths[i], repeats);

  return 0;
}
