#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <limits.h>

#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/params/default.h"
#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/sequences/alphabet.h"
#include "ViennaRNA/constraints/hard.h"
#include "ViennaRNA/constraints/soft.h"
#include "ViennaRNA/eval/gquad.h"
#include "ViennaRNA/structured_domains.h"
#include "ViennaRNA/unstructured_domains.h"
#include "ViennaRNA/eval/structures.h"
#include "ViennaRNA/eval/exterior.h"
#include "ViennaRNA/utils/higher_order_functions.h"
#include "ViennaRNA/intern/mfe_profile.h"
#include "ViennaRNA/intern/mfe_scratch.h"

#include "ViennaRNA/intern/grammar_dat.h"

#ifdef __GNUC__
# define INLINE inline
#else
# define INLINE
#endif

#include "ViennaRNA/constraints/exterior_sc.inc"

#include "ViennaRNA/mfe/exterior.h"

typedef struct {
  int           *stems;
  short         *tmp1;
  short         *tmp2;
  unsigned int  stems_size;
  unsigned int  tmp_size;
  unsigned char owns_tmp;
} f5_helper_data_t;


PRIVATE INLINE unsigned char
use_comparative_sc_fallback(vrna_fold_compound_t *fc);

/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */
PRIVATE INLINE int
reduce_f5_up(vrna_fold_compound_t   *fc,
             unsigned int           j,
             struct sc_f5_dat       *sc_wrapper);


PRIVATE INLINE int *
get_stem_contributions_d0(vrna_fold_compound_t  *fc,
                          unsigned int          j,
                          struct sc_f5_dat      *sc_wrapper,
                          f5_helper_data_t      *cache);


PRIVATE INLINE int *
get_stem_contributions_d2(vrna_fold_compound_t  *fc,
                          unsigned int          j,
                          struct sc_f5_dat      *sc_wrapper,
                          f5_helper_data_t      *cache);


PRIVATE INLINE int *
f5_get_stem_contributions_d5(vrna_fold_compound_t   *fc,
                             unsigned int           j,
                             struct sc_f5_dat       *sc_wrapper,
                             f5_helper_data_t       *cache);


PRIVATE INLINE int *
f5_get_stem_contributions_d3(vrna_fold_compound_t   *fc,
                             unsigned int           j,
                             struct sc_f5_dat       *sc_wrapper,
                             f5_helper_data_t       *cache);


PRIVATE INLINE int *
f5_get_stem_contributions_d53(vrna_fold_compound_t  *fc,
                              unsigned int          j,
                              struct sc_f5_dat      *sc_wrapper,
                              f5_helper_data_t      *cache);


PRIVATE INLINE int
decompose_f5_ext_stem(vrna_fold_compound_t  *fc,
                      unsigned int          j,
                      int                   *stems);


PRIVATE INLINE int
decompose_f5_ext_stem_d0(vrna_fold_compound_t   *fc,
                         unsigned int           j,
                         struct sc_f5_dat       *sc_wrapper,
                         f5_helper_data_t       *cache);


PRIVATE INLINE int
decompose_f5_ext_stem_d2(vrna_fold_compound_t   *fc,
                         unsigned int           j,
                         struct sc_f5_dat       *sc_wrapper,
                         f5_helper_data_t       *cache);


PRIVATE INLINE int
decompose_f5_ext_stem_d1(vrna_fold_compound_t   *fc,
                         unsigned int           j,
                         struct sc_f5_dat       *sc_wrapper,
                         f5_helper_data_t       *cache);


PRIVATE INLINE int
add_f5_gquad(vrna_fold_compound_t   *fc,
             unsigned int           j,
             struct sc_f5_dat       *sc_wrapper);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PRIVATE INLINE unsigned char
use_comparative_sc_fallback(vrna_fold_compound_t *fc)
{
  unsigned int s;

  if ((!fc) ||
      (fc->type != VRNA_FC_TYPE_COMPARATIVE) ||
      (!fc->scs))
    return 0;

  for (s = 0; s < fc->n_seq; s++) {
    vrna_sc_t *sc = fc->scs[s];

    if ((!sc))
      continue;

    if ((sc->f) ||
        (sc->exp_f) ||
        (sc->bt) ||
        (sc->energy_up) ||
        (sc->energy_bp) ||
        (sc->energy_bp_local) ||
        (sc->energy_stack))
      return 1;
  }

  return 0;
}


PUBLIC int
vrna_mfe_exterior_f5(vrna_fold_compound_t *fc)
{
  if (fc) {
    unsigned int          j, length, dangle_model, with_gquad;
    int                   en, *f5;
    vrna_param_t          *P;
    vrna_mfe_scratch_t    *scratch;
    struct sc_f5_dat      sc_wrapper;
    vrna_gr_aux_t         grammar;
    f5_helper_data_t      cache = {
      0
    };
    uint64_t              t0 = 0;

    length        = fc->length;
    f5            = fc->matrices->f5;
    P             = fc->params;
    dangle_model  = P->model_details.dangles;
    with_gquad    = P->model_details.gquad;
    grammar       = fc->aux_grammar;
    scratch       = (vrna_mfe_scratch_t *)fc->matrices->aux_mfe;

    cache.stems       = scratch->ext_stems;
    cache.stems_size  = scratch->ext_stems_size;

    if (fc->type == VRNA_FC_TYPE_COMPARATIVE) {
      if (use_comparative_sc_fallback(fc)) {
        cache.tmp1      = (short *)vrna_alloc(sizeof(short) * fc->n_seq);
        cache.tmp2      = (short *)vrna_alloc(sizeof(short) * fc->n_seq);
        cache.tmp_size  = fc->n_seq;
        cache.owns_tmp  = 1;
      } else {
        cache.tmp1      = scratch->ext_tmp1;
        cache.tmp2      = scratch->ext_tmp2;
        cache.tmp_size  = scratch->ext_tmp_size;
      }
    }

    if (vrna_mfe_profile_enabled())
      t0 = vrna_mfe_profile_now();

    if ((fc->type == VRNA_FC_TYPE_SINGLE) &&
        (fc->strands == 1) &&
        (!P->model_details.circ) &&
        (!fc->domains_up) &&
        (!fc->aux_grammar) &&
        (!fc->hc->f) &&
        (!with_gquad) &&
        ((!fc->sc) ||
         ((!fc->sc->f) &&
          (!fc->sc->exp_f) &&
          (!fc->sc->bt) &&
          (!fc->sc->energy_up) &&
          (!fc->sc->energy_bp) &&
          (!fc->sc->energy_bp_local) &&
          (!fc->sc->energy_stack)))) {
      vrna_hc_t *hc = fc->hc;
      char      *ptype = fc->ptype;
      short     *S = fc->sequence_encoding;
      int       *c = fc->matrices->c;
      int       *indx = fc->jindx;

      f5[0] = 0;
      f5[1] = ((hc->up_ext[1] >= 1) ? 0 : INF);

      for (j = 2; j <= length; j++) {
        unsigned int  i;
        int           ee;

        f5[j] = ((f5[j - 1] != INF) && (hc->up_ext[j] >= 1)) ? f5[j - 1] : INF;

        switch (dangle_model) {
          case 0:
            for (i = j - 1; i > 1; i--) {
              int          ij;
              unsigned int tt;

              cache.stems[i] = INF;
              ij             = indx[j] + i;

              if ((c[ij] == INF) ||
                  !(hc->mx[length * j + i] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP))
                continue;

              tt              = (unsigned int)ptype[ij];
              cache.stems[i]  = c[ij] + vrna_E_exterior_stem((tt == 0) ? 7U : tt, -1, -1, P);
            }

            cache.stems[1] = INF;
            if ((c[indx[j] + 1] != INF) &&
                (hc->mx[length + j] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP))
              cache.stems[1] = c[indx[j] + 1] +
                               vrna_E_exterior_stem((((unsigned int)ptype[indx[j] + 1]) == 0) ?
                                                    7U :
                                                    (unsigned int)ptype[indx[j] + 1],
                                                    -1,
                                                    -1,
                                                    P);

            en    = decompose_f5_ext_stem(fc, j, cache.stems);
            en    = MIN2(en, cache.stems[1]);
            f5[j] = MIN2(f5[j], en);
            break;

          case 2:
            for (i = j - 1; i > 1; i--) {
              int          ij;
              unsigned int tt;

              cache.stems[i] = INF;
              ij             = indx[j] + i;

              if ((c[ij] == INF) ||
                  !(hc->mx[length * j + i] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP))
                continue;

              tt              = (unsigned int)ptype[ij];
              cache.stems[i]  = c[ij] + vrna_E_exterior_stem((tt == 0) ? 7U : tt, S[i - 1], ((j < length) ? S[j + 1] : -1), P);
            }

            cache.stems[1] = INF;
            if ((c[indx[j] + 1] != INF) &&
                (hc->mx[length + j] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP))
              cache.stems[1] = c[indx[j] + 1] +
                               vrna_E_exterior_stem((((unsigned int)ptype[indx[j] + 1]) == 0) ?
                                                    7U :
                                                    (unsigned int)ptype[indx[j] + 1],
                                                    -1,
                                                    ((j < length) ? S[j + 1] : -1),
                                                    P);

            en    = decompose_f5_ext_stem(fc, j, cache.stems);
            en    = MIN2(en, cache.stems[1]);
            f5[j] = MIN2(f5[j], en);
            break;

          default:
            en = INF;

            for (i = j - 1; i > 1; i--) {
              int          ij;
              unsigned int tt;

              /* no dangles */
              cache.stems[i] = INF;
              ij             = indx[j] + i;
              if ((c[ij] != INF) &&
                  (hc->mx[length * j + i] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP)) {
                tt              = (unsigned int)ptype[ij];
                cache.stems[i]  = c[ij] + vrna_E_exterior_stem((tt == 0) ? 7U : tt, -1, -1, P);
              }
            }
            cache.stems[1] = INF;
            if ((c[indx[j] + 1] != INF) &&
                (hc->mx[length + j] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP))
              cache.stems[1] = c[indx[j] + 1] +
                               vrna_E_exterior_stem((((unsigned int)ptype[indx[j] + 1]) == 0) ?
                                                    7U :
                                                    (unsigned int)ptype[indx[j] + 1],
                                                    -1,
                                                    -1,
                                                    P);
            ee  = decompose_f5_ext_stem(fc, j, cache.stems);
            ee  = MIN2(ee, cache.stems[1]);
            en  = MIN2(en, ee);

            /* 5' dangle */
            for (i = j - 1; i > 1; i--) {
              int          ij;
              unsigned int tt;

              cache.stems[i] = INF;
              ij             = indx[j] + i + 1;
              if ((i + 1 < j) &&
                  (c[ij] != INF) &&
                  (hc->mx[length * j + (i + 1)] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP) &&
                  (hc->up_ext[i] >= 1)) {
                tt              = (unsigned int)ptype[ij];
                cache.stems[i]  = c[ij] + vrna_E_exterior_stem((tt == 0) ? 7U : tt, S[i], -1, P);
              }
            }
            cache.stems[1] = INF;
            if ((2 < j) &&
                (c[indx[j] + 2] != INF) &&
                (hc->mx[length * 2 + j] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP) &&
                (hc->up_ext[1] >= 1))
              cache.stems[1] = c[indx[j] + 2] +
                               vrna_E_exterior_stem((((unsigned int)ptype[indx[j] + 2]) == 0) ?
                                                    7U :
                                                    (unsigned int)ptype[indx[j] + 2],
                                                    S[1],
                                                    -1,
                                                    P);
            ee  = decompose_f5_ext_stem(fc, j, cache.stems);
            ee  = MIN2(ee, cache.stems[1]);
            en  = MIN2(en, ee);

            /* 3' dangle */
            for (i = j - 1; i > 1; i--) {
              int          ij;
              unsigned int tt;

              cache.stems[i] = INF;
              ij             = indx[j - 1] + i;
              if ((i + 1 < j) &&
                  (c[ij] != INF) &&
                  (hc->mx[length * (j - 1) + i] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP) &&
                  (hc->up_ext[j] >= 1)) {
                tt              = (unsigned int)ptype[ij];
                cache.stems[i]  = c[ij] + vrna_E_exterior_stem((tt == 0) ? 7U : tt, -1, S[j], P);
              }
            }
            cache.stems[1] = INF;
            if ((2 < j) &&
                (c[indx[j - 1] + 1] != INF) &&
                (hc->mx[length + (j - 1)] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP) &&
                (hc->up_ext[j] >= 1))
              cache.stems[1] = c[indx[j - 1] + 1] +
                               vrna_E_exterior_stem((((unsigned int)ptype[indx[j - 1] + 1]) == 0) ?
                                                    7U :
                                                    (unsigned int)ptype[indx[j - 1] + 1],
                                                    -1,
                                                    S[j],
                                                    P);
            ee  = decompose_f5_ext_stem(fc, j, cache.stems);
            ee  = MIN2(ee, cache.stems[1]);
            en  = MIN2(en, ee);

            /* double dangles */
            for (i = j - 1; i > 1; i--) {
              int          ij;
              unsigned int tt;

              cache.stems[i] = INF;
              ij             = indx[j - 1] + i + 1;
              if ((i + 2 < j) &&
                  (c[ij] != INF) &&
                  (hc->mx[length * (j - 1) + (i + 1)] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP) &&
                  (hc->up_ext[i] >= 1) &&
                  (hc->up_ext[j] >= 1)) {
                tt              = (unsigned int)ptype[ij];
                cache.stems[i]  = c[ij] + vrna_E_exterior_stem((tt == 0) ? 7U : tt, S[i], S[j], P);
              }
            }
            cache.stems[1] = INF;
            if ((3 < j) &&
                (c[indx[j - 1] + 2] != INF) &&
                (hc->mx[length * 2 + (j - 1)] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP) &&
                (hc->up_ext[1] >= 1) &&
                (hc->up_ext[j] >= 1))
              cache.stems[1] = c[indx[j - 1] + 2] +
                               vrna_E_exterior_stem((((unsigned int)ptype[indx[j - 1] + 2]) == 0) ?
                                                    7U :
                                                    (unsigned int)ptype[indx[j - 1] + 2],
                                                    S[1],
                                                    S[j],
                                                    P);
            ee  = decompose_f5_ext_stem(fc, j, cache.stems);
            ee  = MIN2(ee, cache.stems[1]);
            en  = MIN2(en, ee);

            f5[j] = MIN2(f5[j], en);
            break;
        }
      }

      if (vrna_mfe_profile_enabled())
        vrna_mfe_profile_add_exterior(vrna_mfe_profile_now() - t0);

      return f5[length];
    }

    init_sc_f5(fc, &sc_wrapper);

    f5[0] = 0;
    f5[1] = reduce_f5_up(fc, 1, &sc_wrapper);

    if (grammar) {
      for (size_t c = 0; c < vrna_array_size(grammar->f); c++) {
        if (grammar->f[c].cb) {
          en    = grammar->f[c].cb(fc, 1, 1, grammar->f[c].data);
          f5[1] = MIN2(f5[1], en);
        }
      }
    }

    /*
     *  duplicated code may be faster than conditions inside loop or even
     *  using a function pointer ;)
     */
    switch (dangle_model) {
      case 2:
        for (j = 2; j <= length; j++) {
          /* extend previous solution(s) by adding an unpaired region */
          f5[j] = reduce_f5_up(fc, j, &sc_wrapper);

          /* decompose into exterior loop part followed by a stem */
          en    = decompose_f5_ext_stem_d2(fc, j, &sc_wrapper, &cache);
          f5[j] = MIN2(f5[j], en);

          if (with_gquad) {
            en    = add_f5_gquad(fc, j, &sc_wrapper);
            f5[j] = MIN2(f5[j], en);
          }

          if (grammar) {
            for (size_t c = 0; c < vrna_array_size(grammar->f); c++) {
              if (grammar->f[c].cb) {
                en    = grammar->f[c].cb(fc, 1, j, grammar->f[c].data);
                f5[j] = MIN2(f5[j], en);
              }
            }
          }
        }
        break;

      case 0:
        for (j = 2; j <= length; j++) {
          /* extend previous solution(s) by adding an unpaired region */
          f5[j] = reduce_f5_up(fc, j, &sc_wrapper);

          /* decompose into exterior loop part followed by a stem */
          en    = decompose_f5_ext_stem_d0(fc, j, &sc_wrapper, &cache);
          f5[j] = MIN2(f5[j], en);

          if (with_gquad) {
            en    = add_f5_gquad(fc, j, &sc_wrapper);
            f5[j] = MIN2(f5[j], en);
          }

          if (grammar) {
            for (size_t c = 0; c < vrna_array_size(grammar->f); c++) {
              if (grammar->f[c].cb) {
                en    = grammar->f[c].cb(fc, 1, j, grammar->f[c].data);
                f5[j] = MIN2(f5[j], en);
              }
            }
          }
        }
        break;

      default:
        for (j = 2; j <= length; j++) {
          /* extend previous solution(s) by adding an unpaired region */
          f5[j] = reduce_f5_up(fc, j, &sc_wrapper);

          en    = decompose_f5_ext_stem_d1(fc, j, &sc_wrapper, &cache);
          f5[j] = MIN2(f5[j], en);

          if (with_gquad) {
            en    = add_f5_gquad(fc, j, &sc_wrapper);
            f5[j] = MIN2(f5[j], en);
          }

          if (grammar) {
            for (size_t c = 0; c < vrna_array_size(grammar->f); c++) {
              if (grammar->f[c].cb) {
                en    = grammar->f[c].cb(fc, 1, j, grammar->f[c].data);
                f5[j] = MIN2(f5[j], en);
              }
            }
          }
        }
        break;
    }

    free_sc_f5(&sc_wrapper);

    if (cache.owns_tmp) {
      free(cache.tmp1);
      free(cache.tmp2);
    }

    if (vrna_mfe_profile_enabled())
      vrna_mfe_profile_add_exterior(vrna_mfe_profile_now() - t0);

    return f5[length];
  }

  return INF;
}


/*
 #####################################
 # BEGIN OF STATIC HELPER FUNCTIONS  #
 #####################################
 */
PRIVATE INLINE int
decompose_f5_ext_stem_d0(vrna_fold_compound_t   *fc,
                         unsigned int           j,
                         struct sc_f5_dat       *sc_wrapper,
                         f5_helper_data_t       *cache)
{
  int e, *stems;

  stems = get_stem_contributions_d0(fc, j, sc_wrapper, cache);

  /* 1st case, actual decompostion */
  e = decompose_f5_ext_stem(fc, j, stems);

  /* 2nd case, reduce to single stem */
  e = MIN2(e, stems[1]);

  return e;
}


PRIVATE INLINE int
decompose_f5_ext_stem_d2(vrna_fold_compound_t   *fc,
                         unsigned int           j,
                         struct sc_f5_dat       *sc_wrapper,
                         f5_helper_data_t       *cache)
{
  int e, *stems;

  stems = get_stem_contributions_d2(fc, j, sc_wrapper, cache);

  /* 1st case, actual decompostion */
  e = decompose_f5_ext_stem(fc, j, stems);

  /* 2nd case, reduce to single stem */
  e = MIN2(e, stems[1]);

  return e;
}


PRIVATE INLINE int
decompose_f5_ext_stem_d1(vrna_fold_compound_t   *fc,
                         unsigned int           j,
                         struct sc_f5_dat       *sc_wrapper,
                         f5_helper_data_t       *cache)
{
  int e, ee, *stems;

  e = INF;

  /* A) without dangling end contributions */

  /* 1st case, actual decompostion */
  stems = get_stem_contributions_d0(fc, j, sc_wrapper, cache);

  ee = decompose_f5_ext_stem(fc, j, stems);

  /* 2nd case, reduce to single stem */
  ee = MIN2(ee, stems[1]);

  e = MIN2(e, ee);

  /* B) with dangling end contribution on 5' side of stem */
  stems = f5_get_stem_contributions_d5(fc, j, sc_wrapper, cache);

  /* 1st case, actual decompostion */
  ee = decompose_f5_ext_stem(fc, j, stems);

  /* 2nd case, reduce to single stem */
  ee = MIN2(ee, stems[1]);

  e = MIN2(e, ee);

  /* C) with dangling end contribution on 3' side of stem */
  stems = f5_get_stem_contributions_d3(fc, j, sc_wrapper, cache);

  /* 1st case, actual decompostion */
  ee = decompose_f5_ext_stem(fc, j, stems);

  /* 2nd case, reduce to single stem */
  ee = MIN2(ee, stems[1]);

  e = MIN2(e, ee);

  /* D) with dangling end contribution on both sides of stem */
  stems = f5_get_stem_contributions_d53(fc, j, sc_wrapper, cache);

  /* 1st case, actual decompostion */
  ee = decompose_f5_ext_stem(fc, j, stems);

  /* 2nd case, reduce to single stem */
  ee = MIN2(ee, stems[1]);

  e = MIN2(e, ee);

  return e;
}


/*
 *  extend f5 by adding an unpaired nucleotide or an unstructured domain
 *  to the 3' end
 */
PRIVATE INLINE int
reduce_f5_up(vrna_fold_compound_t   *fc,
             unsigned int           j,
             struct sc_f5_dat       *sc_wrapper)
{
  unsigned int          u, k;
  int                   e, en, *f5;
  vrna_ud_t             *domains_up;
  vrna_hc_t             *hc;
  vrna_hc_eval_loop_f   evaluate;
  sc_f5_cb              sc_red_ext;

  f5          = fc->matrices->f5;
  domains_up  = fc->domains_up;
  hc          = fc->hc;
  evaluate    = fc->hc->eval_ext;
  sc_red_ext  = sc_wrapper->red_ext5;
  e           = INF;


  /* check for 3' extension with one unpaired nucleotide */
  if (f5[j - 1] != INF) {
    if (evaluate(1, j, 1, j - 1, VRNA_DECOMP_EXT_EXT, hc)) {
      e = f5[j - 1];

      if (sc_red_ext)
        e += sc_red_ext(j, 1, j - 1, sc_wrapper);
    }
  }

  if ((domains_up) && (domains_up->energy_cb)) {
    for (k = 0; k < (unsigned int)domains_up->uniq_motif_count; k++) {
      u = domains_up->uniq_motif_size[k];
      if ((j >= u) &&
          (f5[j - u] != INF)) {
        if (evaluate(1, j, 1, j - u, VRNA_DECOMP_EXT_EXT, hc)) {
          en = f5[j - u] +
               domains_up->energy_cb(fc,
                                     j - u + 1,
                                     j,
                                     VRNA_UNSTRUCTURED_DOMAIN_EXT_LOOP | VRNA_UNSTRUCTURED_DOMAIN_MOTIF,
                                     domains_up->data);

          if (sc_red_ext)
            en += sc_red_ext(j, 1, j - u, sc_wrapper);

          e = MIN2(e, en);
        }
      }
    }
  }

  return e;
}


PRIVATE INLINE int *
get_stem_contributions_d0(vrna_fold_compound_t  *fc,
                          unsigned int          j,
                          struct sc_f5_dat      *sc_wrapper,
                          f5_helper_data_t      *cache)
{
  char          *ptype;
  short         **S;
  unsigned int  i, s, n_seq, type;
  int           ij, *indx, *c, *stems;
  vrna_param_t  *P;
  vrna_md_t     *md;
  vrna_hc_t             *hc;
  vrna_hc_eval_loop_f   evaluate;

  sc_f5_cb      sc_spl_stem;
  sc_f5_cb      sc_red_stem;

  stems = cache->stems;

  P     = fc->params;
  md    = &(P->model_details);
  indx  = fc->jindx;
  c     = fc->matrices->c;
  ij    = indx[j] + j - 1;
  ptype = (fc->type == VRNA_FC_TYPE_SINGLE) ? fc->ptype : NULL;
  n_seq = (fc->type == VRNA_FC_TYPE_SINGLE) ? 1 : fc->n_seq;
  S     = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S;

  hc          = fc->hc;
  evaluate    = hc->eval_ext;
  sc_spl_stem = sc_wrapper->decomp_stem5;
  sc_red_stem = sc_wrapper->red_stem5;

  switch (fc->type) {
    case VRNA_FC_TYPE_SINGLE:
      for (i = j - 1; i > 1; i--, ij--) {
        stems[i] = INF;

        if ((c[ij] != INF) &&
            (evaluate(1, j, i - 1, i, VRNA_DECOMP_EXT_EXT_STEM, hc))) {
          stems[i]  = c[ij];
          type      = vrna_get_ptype(ij, ptype);
          stems[i]  += vrna_E_exterior_stem(type, -1, -1, P);
        }
      }
      break;

    case VRNA_FC_TYPE_COMPARATIVE:
      for (i = j - 1; i > 1; i--, ij--) {
        stems[i] = INF;
        if ((c[ij] != INF) &&
            (evaluate(1, j, i - 1, i, VRNA_DECOMP_EXT_EXT_STEM, hc))) {
          stems[i] = c[ij];

          for (s = 0; s < n_seq; s++) {
            type      = vrna_get_ptype_md(S[s][i], S[s][j], md);
            stems[i]  += vrna_E_exterior_stem(type, -1, -1, P);
          }
        }
      }
      break;
  }

  if (sc_spl_stem)
    for (i = j - 1; i > 1; i--)
      if (stems[i] != INF)
        stems[i] += sc_spl_stem(j, i - 1, i, sc_wrapper);

  stems[1]  = INF;
  ij        = indx[j] + 1;

  if ((c[ij] != INF) &&
      (evaluate(1, j, 1, j, VRNA_DECOMP_EXT_STEM, hc))) {
    stems[1] = c[ij];

    switch (fc->type) {
      case VRNA_FC_TYPE_SINGLE:
        type      = vrna_get_ptype(ij, ptype);
        stems[1]  += vrna_E_exterior_stem(type, -1, -1, P);
        break;

      case VRNA_FC_TYPE_COMPARATIVE:
        for (s = 0; s < n_seq; s++) {
          type      = vrna_get_ptype_md(S[s][1], S[s][j], md);
          stems[1]  += vrna_E_exterior_stem(type, -1, -1, P);
        }
        break;
    }

    if (sc_red_stem)
      stems[1] += sc_red_stem(j, 1, j, sc_wrapper);
  }

  return stems;
}


PRIVATE INLINE int *
get_stem_contributions_d2(vrna_fold_compound_t  *fc,
                          unsigned int          j,
                          struct sc_f5_dat      *sc_wrapper,
                          f5_helper_data_t      *cache)
{
  char          *ptype;
  short         *S, sj1, *si1, **SS, **S5, **S3, *s3j, *sj;
  unsigned int  n, i, s, n_seq, **a2s, type, *sn;
  int           ij, *indx, *c, *stems, mm5;
  vrna_param_t  *P;
  vrna_md_t     *md;
  vrna_hc_t     *hc;
  vrna_hc_eval_loop_f evaluate;

  sc_f5_cb      sc_spl_stem;
  sc_f5_cb      sc_red_stem;

  stems = cache->stems;

  n     = fc->length;
  sn    = fc->strand_number;
  P     = fc->params;
  md    = &(P->model_details);
  indx  = fc->jindx;
  c     = fc->matrices->c;
  ij    = indx[j] + j - 1;

  hc          = fc->hc;
  evaluate    = hc->eval_ext;
  sc_spl_stem = sc_wrapper->decomp_stem5;
  sc_red_stem = sc_wrapper->red_stem5;

  switch (fc->type) {
    case VRNA_FC_TYPE_SINGLE:
      S     = fc->sequence_encoding;
      ptype = fc->ptype;
      si1   = S + j - 2;
      sj1   = ((j < n) && (sn[j] == sn[j + 1])) ? S[j + 1] : -1;

      for (i = j - 1; i > 1; i--, ij--, si1--) {
        stems[i] = INF;
        if ((c[ij] != INF) &&
            (evaluate(1, j, i - 1, i, VRNA_DECOMP_EXT_EXT_STEM, hc))) {
          type      = vrna_get_ptype(ij, ptype);
          stems[i]  = c[ij] +
                      vrna_E_exterior_stem(type, *si1, sj1, P);
        }
      }

      if (sc_spl_stem)
        for (i = j - 1; i > 1; i--)
          if (stems[i] != INF)
            stems[i] += sc_spl_stem(j, i - 1, i, sc_wrapper);

      stems[1]  = INF;
      ij        = indx[j] + 1;

      if ((c[ij] != INF) && (evaluate(1, j, 1, j, VRNA_DECOMP_EXT_STEM, hc))) {
        type      = vrna_get_ptype(ij, ptype);
        stems[1]  = c[ij] +
                    vrna_E_exterior_stem(type, -1, sj1, P);

        if (sc_red_stem)
          stems[1] += sc_red_stem(j, 1, j, sc_wrapper);
      }

      break;

    case VRNA_FC_TYPE_COMPARATIVE:
      n_seq = fc->n_seq;
      SS    = fc->S;
      S5    = fc->S5;
      S3    = fc->S3;
      a2s   = fc->a2s;

      /* pre-compute S3[s][j - 1] */
      s3j = cache->tmp1;
      sj  = cache->tmp2;
      for (s = 0; s < n_seq; s++) {
        s3j[s]  = (a2s[s][j] < a2s[s][n]) ? S3[s][j] : -1;
        sj[s]   = SS[s][j];
      }

      for (i = j - 1; i > 1; i--, ij--) {
        stems[i] = INF;
        if ((c[ij] != INF) &&
            (evaluate(1, j, i - 1, i, VRNA_DECOMP_EXT_EXT_STEM, hc))) {
          stems[i] = c[ij];
          for (s = 0; s < n_seq; s++) {
            type      = vrna_get_ptype_md(SS[s][i], sj[s], md);
            mm5       = (a2s[s][i] > 1) ? S5[s][i] : -1;
            stems[i]  += vrna_E_exterior_stem(type, mm5, s3j[s], P);
          }
        }
      }

      if (sc_spl_stem)
        for (i = j - 1; i > 1; i--)
          if (stems[i] != INF)
            stems[i] += sc_spl_stem(j, i - 1, i, sc_wrapper);

      stems[1]  = INF;
      ij        = indx[j] + 1;

      if ((c[ij] != INF) && (evaluate(1, j, 1, j, VRNA_DECOMP_EXT_STEM, hc))) {
        stems[1] = c[ij];

        for (s = 0; s < n_seq; s++) {
          type      = vrna_get_ptype_md(SS[s][1], sj[s], md);
          stems[1]  += vrna_E_exterior_stem(type, -1, s3j[s], P);
        }

        if (sc_red_stem)
          stems[1] += sc_red_stem(j, 1, j, sc_wrapper);
      }
      break;
  }

  return stems;
}


PRIVATE INLINE int *
f5_get_stem_contributions_d5(vrna_fold_compound_t   *fc,
                             unsigned int           j,
                             struct sc_f5_dat       *sc_wrapper,
                             f5_helper_data_t       *cache)
{
  char          *ptype;
  short         *S, *si1, **SS, **S5, *sj;
  unsigned int  i, s, n_seq, **a2s, type;
  int           ij, *indx, *c, *stems, mm5;
  vrna_param_t  *P;
  vrna_md_t     *md;

  vrna_hc_t     *hc;
  vrna_hc_eval_loop_f evaluate;
  sc_f5_cb      sc_spl_stem;
  sc_f5_cb      sc_red_stem;

  stems = cache->stems;

  P     = fc->params;
  md    = &(P->model_details);
  indx  = fc->jindx;
  c     = fc->matrices->c;
  ij    = indx[j] + j;

  hc          = fc->hc;
  evaluate    = hc->eval_ext;
  sc_spl_stem = sc_wrapper->decomp_stem5;
  sc_red_stem = sc_wrapper->red_stem5;

  switch (fc->type) {
    case VRNA_FC_TYPE_SINGLE:
      S     = fc->sequence_encoding;
      ptype = fc->ptype;
      si1   = S + j - 1;

      for (i = j - 1; i > 1; i--, ij--, si1--) {
        stems[i] = INF;
        if ((c[ij] != INF) &&
            (evaluate(1, j, i - 1, i + 1, VRNA_DECOMP_EXT_EXT_STEM, hc))) {
          type      = vrna_get_ptype(ij, ptype);
          stems[i]  = c[ij] +
                      vrna_E_exterior_stem(type, *si1, -1, P);
        }
      }

      if (sc_spl_stem)
        for (i = j - 1; i > 1; i--)
          if (stems[i] != INF)
            stems[i] += sc_spl_stem(j, i - 1, i + 1, sc_wrapper);

      stems[1] = INF;
      if (2 < j) {
        ij = indx[j] + 2;

        if ((c[ij] != INF) && (evaluate(1, j, 2, j, VRNA_DECOMP_EXT_STEM, hc))) {
          type      = vrna_get_ptype(ij, ptype);
          stems[1]  = c[ij] +
                      vrna_E_exterior_stem(type, S[1], -1, P);

          if (sc_red_stem)
            stems[1] += sc_red_stem(j, 2, j, sc_wrapper);
        }
      }

      break;

    case VRNA_FC_TYPE_COMPARATIVE:
      n_seq = fc->n_seq;
      SS    = fc->S;
      S5    = fc->S5;
      a2s   = fc->a2s;

      sj = cache->tmp1;
      for (s = 0; s < n_seq; s++)
        sj[s] = SS[s][j];

      for (i = j - 1; i > 1; i--, ij--) {
        stems[i] = INF;
        if ((c[ij] != INF) &&
            (evaluate(1, j, i - 1, i + 1, VRNA_DECOMP_EXT_EXT_STEM, hc))) {
          stems[i] = c[ij];
          for (s = 0; s < n_seq; s++) {
            type      = vrna_get_ptype_md(SS[s][i + 1], sj[s], md);
            mm5       = (a2s[s][i + 1] > 1) ? S5[s][i + 1] : -1;
            stems[i]  = vrna_E_exterior_stem(type, mm5, -1, P);
          }
        }
      }

      if (sc_spl_stem)
        for (i = j - 1; i > 1; i--)
          if (stems[i] != INF)
            stems[i] += sc_spl_stem(j, i - 1, i + 1, sc_wrapper);

      stems[1] = INF;

      if (2 < j) {
        ij = indx[j] + 2;

        if ((c[ij] != INF) && (evaluate(1, j, 2, j, VRNA_DECOMP_EXT_STEM, hc))) {
          stems[1] = c[ij];
          for (s = 0; s < n_seq; s++) {
            type      = vrna_get_ptype_md(SS[s][2], sj[s], md);
            mm5       = (a2s[s][2] > 1) ? S5[s][2] : -1;
            stems[i]  = vrna_E_exterior_stem(type, mm5, -1, P);
          }

          if (sc_red_stem)
            stems[1] += sc_red_stem(j, 2, j, sc_wrapper);
        }
      }

      break;
  }

  return stems;
}


PRIVATE INLINE int *
f5_get_stem_contributions_d3(vrna_fold_compound_t   *fc,
                             unsigned int           j,
                             struct sc_f5_dat       *sc_wrapper,
                             f5_helper_data_t       *cache)
{
  char          *ptype;
  short         *S, sj1, **SS, **S3, *s3j1, *ssj1;
  unsigned int  i, n, s, n_seq, **a2s, type;
  int           ij, *indx, *c, *stems;
  vrna_param_t  *P;
  vrna_md_t     *md;

  vrna_hc_t     *hc;
  vrna_hc_eval_loop_f evaluate;
  sc_f5_cb      sc_spl_stem;
  sc_f5_cb      sc_red_stem;

  stems = cache->stems;

  n     = fc->length;
  P     = fc->params;
  md    = &(P->model_details);
  indx  = fc->jindx;
  c     = fc->matrices->c;
  ij    = indx[j - 1] + j - 1;

  hc          = fc->hc;
  evaluate    = hc->eval_ext;
  sc_spl_stem = sc_wrapper->decomp_stem51;
  sc_red_stem = sc_wrapper->red_stem5;

  switch (fc->type) {
    case VRNA_FC_TYPE_SINGLE:
      S     = fc->sequence_encoding;
      ptype = fc->ptype;
      sj1   = S[j];

      for (i = j - 1; i > 1; i--, ij--) {
        stems[i] = INF;
        if ((i + 1 < j) &&
            (c[ij] != INF) &&
            (evaluate(1, j, i - 1, i, VRNA_DECOMP_EXT_EXT_STEM1, hc))) {
          type      = vrna_get_ptype(ij, ptype);
          stems[i]  = c[ij] +
                      vrna_E_exterior_stem(type, -1, sj1, P);
        }
      }

      if (sc_spl_stem)
        for (i = j - 1; i > 1; i--)
          if (stems[i] != INF)
            stems[i] += sc_spl_stem(j, i - 1, i, sc_wrapper);

      stems[1] = INF;

      if (2 < j) {
        ij = indx[j - 1] + 1;

        if ((c[ij] != INF) && (evaluate(1, j, 1, j - 1, VRNA_DECOMP_EXT_STEM, hc))) {
          type      = vrna_get_ptype(ij, ptype);
          stems[1]  = c[ij] +
                      vrna_E_exterior_stem(type, -1, sj1, P);

          if (sc_red_stem)
            stems[1] += sc_red_stem(j, 1, j - 1, sc_wrapper);
        }
      }

      break;

    case VRNA_FC_TYPE_COMPARATIVE:
      n_seq = fc->n_seq;
      SS    = fc->S;
      S3    = fc->S3;
      a2s   = fc->a2s;

      /* pre-compute S3[s][j - 1] */
      s3j1  = cache->tmp1;
      ssj1  = cache->tmp2;
      for (s = 0; s < n_seq; s++) {
        s3j1[s] = (a2s[s][j - 1] < a2s[s][n]) ? S3[s][j - 1] : -1;
        ssj1[s] = SS[s][j - 1];
      }

      for (i = j - 1; i > 1; i--, ij--) {
        stems[i] = INF;
        if ((i + 1 < j) &&
            (c[ij] != INF) &&
            (evaluate(1, j, i - 1, i, VRNA_DECOMP_EXT_EXT_STEM1, hc))) {
          stems[i] = c[ij];
          for (s = 0; s < n_seq; s++) {
            type      = vrna_get_ptype_md(SS[s][i], ssj1[s], md);
            stems[i]  += vrna_E_exterior_stem(type, -1, s3j1[s], P);
          }
        }
      }

      if (sc_spl_stem)
        for (i = j - 1; i > 1; i--)
          if (stems[i] != INF)
            stems[i] += sc_spl_stem(j, i - 1, i, sc_wrapper);

      stems[1] = INF;

      if (2 < j) {
        ij = indx[j - 1] + 1;

        if ((c[ij] != INF) && (evaluate(1, j, 1, j - 1, VRNA_DECOMP_EXT_STEM, hc))) {
          stems[1] = c[ij];

          for (s = 0; s < n_seq; s++) {
            type      = vrna_get_ptype_md(SS[s][1], ssj1[s], md);
            stems[1]  += vrna_E_exterior_stem(type, -1, s3j1[s], P);
          }

          if (sc_red_stem)
            stems[1] += sc_red_stem(j, 1, j - 1, sc_wrapper);
        }
      }
      break;
  }

  return stems;
}


PRIVATE INLINE int *
f5_get_stem_contributions_d53(vrna_fold_compound_t  *fc,
                              unsigned int          j,
                              struct sc_f5_dat      *sc_wrapper,
                              f5_helper_data_t      *cache)
{
  char          *ptype;
  short         *S, *si1, sj1, **SS, **S5, **S3, *s3j1, *ssj1;
  unsigned int  i, n, s, n_seq, **a2s, type;
  int           ij, *indx, *c, *stems;
  vrna_param_t  *P;
  vrna_md_t     *md;

  vrna_hc_t     *hc;
  vrna_hc_eval_loop_f evaluate;
  sc_f5_cb      sc_spl_stem;
  sc_f5_cb      sc_red_stem;

  stems = cache->stems;

  n     = fc->length;
  P     = fc->params;
  md    = &(P->model_details);
  indx  = fc->jindx;
  c     = fc->matrices->c;
  ij    = indx[j - 1] + j;

  hc          = fc->hc;
  evaluate    = hc->eval_ext;
  sc_spl_stem = sc_wrapper->decomp_stem51;
  sc_red_stem = sc_wrapper->red_stem5;

  switch (fc->type) {
    case VRNA_FC_TYPE_SINGLE:
      S     = fc->sequence_encoding;
      ptype = fc->ptype;
      sj1   = S[j];
      si1   = S + j - 1;

      for (i = j - 1; i > 1; i--, ij--, si1--) {
        stems[i] = INF;
        if ((i + 2 < j) &&
            (c[ij] != INF) &&
            (evaluate(1, j, i - 1, i + 1, VRNA_DECOMP_EXT_EXT_STEM1, hc))) {
          type      = vrna_get_ptype(ij, ptype);
          stems[i]  = c[ij] +
                      vrna_E_exterior_stem(type, *si1, sj1, P);
        }
      }

      if (sc_spl_stem)
        for (i = j - 1; i > 1; i--)
          if (stems[i] != INF)
            stems[i] += sc_spl_stem(j, i - 1, i + 1, sc_wrapper);

      stems[1] = INF;

      if (3 < j) {
        ij = indx[j - 1] + 2;

        if ((c[ij] != INF) && (evaluate(1, j, 2, j - 1, VRNA_DECOMP_EXT_STEM, hc))) {
          type      = vrna_get_ptype(ij, ptype);
          stems[1]  = c[ij] +
                      vrna_E_exterior_stem(type, S[1], sj1, P);

          if (sc_red_stem)
            stems[1] += sc_red_stem(j, 2, j - 1, sc_wrapper);
        }
      }

      break;

    case VRNA_FC_TYPE_COMPARATIVE:
      n_seq = fc->n_seq;
      SS    = fc->S;
      S5    = fc->S5;
      S3    = fc->S3;
      a2s   = fc->a2s;

      /* pre-compute S3[s][j - 1] */
      s3j1  = cache->tmp1;
      ssj1  = cache->tmp2;
      for (s = 0; s < n_seq; s++) {
        s3j1[s] = (a2s[s][j - 1] < a2s[s][n]) ? S3[s][j - 1] : -1;
        ssj1[s] = SS[s][j - 1];
      }

      for (i = j - 1; i > 1; i--, ij--) {
        stems[i] = INF;
        if ((i + 1 < j) &&
            (c[ij] != INF) &&
            (evaluate(1, j, i - 1, i + 1, VRNA_DECOMP_EXT_EXT_STEM1, hc))) {
          stems[i] = c[ij];
          for (s = 0; s < n_seq; s++) {
            type      = vrna_get_ptype_md(SS[s][i + 1], ssj1[s], md);
            stems[i]  += vrna_E_exterior_stem(type,
                                              (a2s[s][i + 1] > 1) ? S5[s][i + 1] : -1,
                                              s3j1[s],
                                              P);
          }
        }
      }

      if (sc_spl_stem)
        for (i = j - 1; i > 1; i--)
          if (stems[i] != INF)
            stems[i] += sc_spl_stem(j, i - 1, i + 1, sc_wrapper);

      stems[1] = INF;
      if (3 < j) {
        ij = indx[j - 1] + 2;

        if ((c[ij] != INF) && (evaluate(1, j, 2, j - 1, VRNA_DECOMP_EXT_STEM, hc))) {
          stems[1] = c[ij];
          for (s = 0; s < n_seq; s++) {
            type      = vrna_get_ptype_md(SS[s][2], ssj1[s], md);
            stems[1]  += vrna_E_exterior_stem(type, (a2s[s][2] > 1) ? S5[s][2] : -1, s3j1[s], P);
          }

          if (sc_red_stem)
            stems[1] += sc_red_stem(j, 2, j - 1, sc_wrapper);
        }
      }
      break;
  }

  return stems;
}


PRIVATE INLINE int
add_f5_gquad(vrna_fold_compound_t   *fc,
             unsigned int           j,
             struct sc_f5_dat       *sc_wrapper VRNA_UNUSED)
{
  unsigned int      i;
  int               e, e_gq, *f5;
  vrna_smx_csr(int) *c_gq;

  f5    = fc->matrices->f5;
  c_gq  = fc->matrices->c_gq;
  e     = INF;

  for (i = j - 1; i > 1; i--) {
#ifndef VRNA_DISABLE_C11_FEATURES
    e_gq = vrna_smx_csr_get(c_gq, i, j, INF);
#else
    e_gq = vrna_smx_csr_int_get(c_gq, i, j, INF);
#endif
    if ((f5[i - 1] != INF) &&
        (e_gq != INF))
      e = MIN2(e, f5[i - 1] + e_gq);
  }

#ifndef VRNA_DISABLE_C11_FEATURES
  e_gq = vrna_smx_csr_get(c_gq, 1, j, INF);
#else
  e_gq = vrna_smx_csr_int_get(c_gq, 1, j, INF);
#endif
  if (e_gq != INF)
    e = MIN2(e, e_gq);

  return e;
}


PRIVATE INLINE int
decompose_f5_ext_stem(vrna_fold_compound_t  *fc,
                      unsigned int          j,
                      int                   *stems)
{
  int       e, *f5;

  f5  = fc->matrices->f5;
  e   = INF;

  const int count = j;

  e = vrna_fun_zip_add_min(f5 + 1, stems + 2, count - 2);

  return e;
}
