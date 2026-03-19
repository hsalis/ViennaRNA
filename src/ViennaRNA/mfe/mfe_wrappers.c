/*
 * Various wrappers to create simplified interfaces for Minimum Free Energy
 * prediction
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ViennaRNA/fold_compound.h"
#include "ViennaRNA/model.h"
#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/backtrack/global.h"
#include "ViennaRNA/constraints/soft.h"
#include "ViennaRNA/datastructures/sparse_mx.h"
#include "ViennaRNA/mfe/global.h"


#ifndef INLINE
# ifdef __GNUC__
#   define INLINE inline
# else
#   define INLINE
# endif
#endif


PRIVATE INLINE int
get_stored_bp_contributions(vrna_sc_bp_storage_t  *container,
                            unsigned int          j);


PRIVATE int
backtrack_second_monomer(vrna_fold_compound_t *vc,
                         unsigned int         l1,
                         unsigned int         l2,
                         char                 *structure,
                         float                *mfe);


/* wrappers for single sequences */
PUBLIC float
vrna_fold(const char  *string,
          char        *structure)
{
  float                 mfe;
  vrna_fold_compound_t  *vc;
  vrna_md_t             md;

  vrna_md_set_default(&md);
  vc  = vrna_fold_compound(string, &md, 0);
  mfe = vrna_mfe(vc, structure);

  vrna_fold_compound_free(vc);

  return mfe;
}


PUBLIC float
vrna_circfold(const char  *string,
              char        *structure)
{
  float                 mfe;
  vrna_fold_compound_t  *vc;
  vrna_md_t             md;

  vrna_md_set_default(&md);
  md.circ = 1;
  vc      = vrna_fold_compound(string, &md, 0);
  mfe     = vrna_mfe(vc, structure);

  vrna_fold_compound_free(vc);

  return mfe;
}


/* wrappers for multiple sequence alignments */

PUBLIC float
vrna_alifold(const char **strings,
             char       *structure)
{
  float                 mfe;
  vrna_fold_compound_t  *vc;
  vrna_md_t             md;

  vrna_md_set_default(&md);

  vc  = vrna_fold_compound_comparative(strings, &md, VRNA_OPTION_DEFAULT);
  mfe = vrna_mfe(vc, structure);

  vrna_fold_compound_free(vc);

  return mfe;
}


PUBLIC float
vrna_circalifold(const char **sequences,
                 char       *structure)
{
  float                 mfe;
  vrna_fold_compound_t  *vc;
  vrna_md_t             md;

  vrna_md_set_default(&md);
  md.circ = 1;

  vc  = vrna_fold_compound_comparative(sequences, &md, VRNA_OPTION_DEFAULT);
  mfe = vrna_mfe(vc, structure);

  vrna_fold_compound_free(vc);

  return mfe;
}


/* wrappers for RNA-RNA cofolding interaction */
PUBLIC float
vrna_cofold(const char  *seq,
            char        *structure)
{
  float                 mfe;
  vrna_fold_compound_t  *vc;
  vrna_md_t             md;

  vrna_md_set_default(&md);
  md.min_loop_size = 0;  /* set min loop length to 0 */

  /* get compound structure */
  vc = vrna_fold_compound(seq, &md, 0);

  mfe = vrna_mfe_dimer(vc, structure);

  vrna_fold_compound_free(vc);

  return mfe;
}


PUBLIC float
vrna_mfe_dimer(vrna_fold_compound_t *vc,
               char                 *structure)
{
  char          *s2, *ss1, *ss2;
  unsigned int  l1, l2;
  float         mfe, mfe1, mfe2;

  mfe = vrna_mfe(vc, structure);

  /*
   *  for backward compatibility reasons, we alson need
   *  to see whether the unconnected structure is better
   */
  if (vc->strands > 1) {
    l1  = vc->nucleotides[0].length;
    l2  = vc->nucleotides[1].length;
    s2  = vc->nucleotides[1].string;
    ss1 = (char *)vrna_alloc(sizeof(char) * (l1 + 1));
    ss2 = (char *)vrna_alloc(sizeof(char) * (l2 + 1));

    mfe1 = vrna_backtrack5(vc, l1, ss1);

    if (!backtrack_second_monomer(vc, l1, l2, ss2, &mfe2)) {
      vrna_fold_compound_t *fc2 = vrna_fold_compound(s2,
                                                     &(vc->params->model_details),
                                                     VRNA_OPTION_DEFAULT);

      /* extract hard constraints for second sequence, if any non-defaults */
      if (vc->hc) {
        vrna_smx_csr(vrna_uchar) *hc_nondefaults = vrna_hc_nondefaults(vc);

        if (hc_nondefaults) {
          /* single nucleotide constraints first */
          for (unsigned int i = l1 + 1; i <= l1 + l2; ++i) {
#ifndef VRNA_DISABLE_C11_FEATURES
            unsigned char constraint = vrna_smx_csr_get(hc_nondefaults,
                                                        i,
                                                        i,
                                                        VRNA_CONSTRAINT_CONTEXT_NO_REMOVE);
#else
            unsigned char constraint = vrna_smx_csr_vrna_uchar_get(hc_nondefaults,
                                                        i,
                                                        i,
                                                        VRNA_CONSTRAINT_CONTEXT_NO_REMOVE);
#endif
            if (constraint != VRNA_CONSTRAINT_CONTEXT_NO_REMOVE) {
              vrna_hc_add_up(fc2, i - l1, constraint);
              vrna_hc_add_up(fc2, l2 + i - l1, constraint);
            }
          }

          /* base pair constraints next */
          for (unsigned int i = l1 + 1; i <= l1 + l2; ++i) {
            for (unsigned int k = 1; k < l2; k++) {
              unsigned int j = i + k;

              if (j > l1 + l2)
                break;

#ifndef VRNA_DISABLE_C11_FEATURES
              unsigned char constraint = vrna_smx_csr_get(hc_nondefaults,
                                                          i,
                                                          j,
                                                          VRNA_CONSTRAINT_CONTEXT_NO_REMOVE);
#else
              unsigned char constraint = vrna_smx_csr_vrna_uchar_get(hc_nondefaults,
                                                          i,
                                                          j,
                                                          VRNA_CONSTRAINT_CONTEXT_NO_REMOVE);
#endif
              if (constraint != VRNA_CONSTRAINT_CONTEXT_NO_REMOVE) {
                vrna_hc_add_bp(fc2, i - l1, j - l1, constraint | VRNA_CONSTRAINT_CONTEXT_NO_REMOVE);
                vrna_hc_add_bp(fc2, l2 + i - l1, l2 + j - l1, constraint | VRNA_CONSTRAINT_CONTEXT_NO_REMOVE);
              }
            }
          }

#ifndef VRNA_DISABLE_C11_FEATURES
          vrna_smx_csr_free(hc_nondefaults);
#else
          vrna_smx_csr_vrna_uchar_free(hc_nondefaults);
#endif
        }
      }

      /* extract soft constraints for second sequence, if any */
      if (vc->sc) {
        if (vc->sc->up_storage) {
          for (unsigned int i = l1 + 1; i <= l1 + l2; ++i)
            vrna_sc_add_up(fc2, i - l1, (FLT_OR_DBL)vc->sc->up_storage[i] / 100., VRNA_OPTION_DEFAULT);
        }

        if (vc->sc->bp_storage) {
          for (unsigned int i = l1 + 1; i <= l1 + l2; ++i) {
            if (vc->sc->bp_storage[i]) {
              for (unsigned int k = 1; k < l2; k++) {
                unsigned int j = i + k;

                if (j > l1 + l2)
                  break;

                int e = get_stored_bp_contributions(vc->sc->bp_storage[i], j);

                if (e != 0)
                  vrna_sc_add_bp(fc2, i - l1, j - l1, (FLT_OR_DBL)e / 100., VRNA_OPTION_DEFAULT);
              }
            }
          }
        }

        if (vc->sc->energy_stack) {
          for (unsigned int i = l1 + 1; i <= l1 + l2; ++i) {
            if (vc->sc->energy_stack[i] != 0)
              vrna_sc_add_stack(fc2, i - l1, (FLT_OR_DBL)vc->sc->energy_stack[i] / 100., VRNA_OPTION_DEFAULT);
          }
        }
      }

      mfe2 = vrna_mfe(fc2, ss2);
      vrna_fold_compound_free(fc2);
    }

    if (mfe1 + mfe2 < mfe) {
      mfe = mfe1 + mfe2;
      memcpy(structure, ss1, sizeof(char) * l1);
      memcpy(structure + l1, ss2, sizeof(char) * l2);
      structure[l1 + l2] = '\0';
    }

    free(ss1);
    free(ss2);
  }

  return mfe;
}


PRIVATE INLINE int
get_stored_bp_contributions(vrna_sc_bp_storage_t  *container,
                            unsigned int          j)
{
  unsigned int  cnt;
  int           e;

  e = 0;

  /* go through list of constraints for current position i */
  for (cnt = 0; container[cnt].interval_start != 0; cnt++) {
    if (container[cnt].interval_start > j)
      break; /* only constraints for pairs (i,q) with q > j left */

    if (container[cnt].interval_end < j)
      continue; /* constraint for pairs (i,q) with q < j */

    /* constraint has interval [p,q] with p <= j <= q */
    e += container[cnt].e;
  }

  return e;
}


PRIVATE int
backtrack_second_monomer(vrna_fold_compound_t *vc,
                         unsigned int         l1,
                         unsigned int         l2,
                         char                 *structure,
                         float                *mfe)
{
  int               ret;
  sect              bt_stack[1];
  vrna_bp_stack_t   *bp_stack;

  if ((!vc) ||
      (!structure) ||
      (!mfe) ||
      (!vc->matrices) ||
      (!vc->matrices->fms5) ||
      (!vc->matrices->fms5[1]))
    return 0;

  if (vc->matrices->fms5[1][l1 + 1] >= INF)
    return 0;

  memset(structure, '.', sizeof(char) * l2);
  structure[l2]  = '\0';
  bt_stack[0]    = (sect){
    .i = (int)(l1 + 1),
    .j = 1,
    .ml = VRNA_MX_FLAG_MS5
  };

  bp_stack       = (vrna_bp_stack_t *)vrna_alloc(sizeof(vrna_bp_stack_t) * (4 * (1 + l2 / 2)));
  bp_stack[0].i  = 0;
  ret            = vrna_backtrack_from_intervals(vc, bp_stack, bt_stack, 1);

  if (ret) {
    unsigned int n_pairs = (unsigned int)bp_stack[0].i;

    for (unsigned int n = 1; n <= n_pairs; n++) {
      unsigned int i = (unsigned int)bp_stack[n].i;
      unsigned int j = (unsigned int)bp_stack[n].j;

      if ((i <= l1) ||
          (j <= l1) ||
          (i > l1 + l2) ||
          (j > l1 + l2)) {
        ret = 0;
        break;
      }

      structure[i - l1 - 1] = '(';
      structure[j - l1 - 1] = ')';
    }
  }

  free(bp_stack);

  if (ret)
    *mfe = (float)vc->matrices->fms5[1][l1 + 1] / 100.;

  return ret;
}

