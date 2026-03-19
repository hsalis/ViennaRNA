#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/sequences/alphabet.h"
#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/constraints/hard.h"
#include "ViennaRNA/constraints/soft.h"
#include "ViennaRNA/eval/exterior.h"
#include "ViennaRNA/eval/gquad.h"
#include "ViennaRNA/mfe/gquad.h"
#include "ViennaRNA/structured_domains.h"
#include "ViennaRNA/unstructured_domains.h"
#include "ViennaRNA/eval/internal.h"
#include "ViennaRNA/params/salt.h"
#include "ViennaRNA/intern/mfe_profile.h"
#include "ViennaRNA/intern/mfe_scratch.h"


#ifdef __GNUC__
# define INLINE inline
#else
# define INLINE
#endif

#include "ViennaRNA/constraints/internal_sc.inc"

#include "ViennaRNA/mfe/internal.h"


typedef struct {
  struct sc_int_dat     sc_wrapper;

  unsigned int          *tt;
  unsigned int          tt_size;
  unsigned char         owns_tt;
} helper_data_t;


/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */

PRIVATE int
mfe_internal_loop(vrna_fold_compound_t  *fc,
                  unsigned int          i,
                  unsigned int          j);


PRIVATE int
mfe_internal_loop_single_fast(vrna_fold_compound_t  *fc,
                              unsigned int          i,
                              unsigned int          j);


PRIVATE int
mfe_internal_loop_ext(vrna_fold_compound_t  *fc,
                      unsigned int          i,
                      unsigned int          j,
                      unsigned int          *ip,
                      unsigned int          *iq);


PRIVATE void
prepare_intloop_helpers(helper_data_t        *h,
                        vrna_fold_compound_t  *fc,
                        unsigned int          i,
                        unsigned int          j);


PRIVATE void
free_intloop_helpers(helper_data_t *h);


PRIVATE int
mfe_stacks(vrna_fold_compound_t *fc,
           unsigned int         i,
           unsigned int         j,
           helper_data_t        *helpers);


PRIVATE int
mfe_bulges(vrna_fold_compound_t *fc,
           unsigned int         i,
           unsigned int         j,
           helper_data_t        *helpers);


PRIVATE INLINE unsigned int
ptype_fast(char          *ptype,
           unsigned int  ij);


PRIVATE INLINE unsigned char
use_internal_single_fast_path(vrna_fold_compound_t *fc);


PRIVATE INLINE unsigned char
use_comparative_sc_fallback(vrna_fold_compound_t *fc);


PRIVATE INLINE unsigned char
hc_int_linear_fast(vrna_hc_t      *hc,
                   unsigned int   n,
                   unsigned int   i,
                   unsigned int   j,
                   unsigned int   k,
                   unsigned int   l);


PRIVATE INLINE int
mfe_E_internal_single_fast(unsigned int  n1,
                           unsigned int  n2,
                           unsigned int  type,
                           unsigned int  type_2,
                           int           si1,
                           int           sj1,
                           int           sp1,
                           int           sq1,
                           vrna_param_t  *P);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC int
vrna_mfe_internal(vrna_fold_compound_t  *fc,
                  unsigned int          i,
                  unsigned int          j)
{
  if (fc) {
    uint64_t  t0  = 0;
    int       ret;

    if (vrna_mfe_profile_enabled())
      t0 = vrna_mfe_profile_now();

    if (use_internal_single_fast_path(fc))
      ret = mfe_internal_loop_single_fast(fc, i, j);
    else
      ret = mfe_internal_loop(fc, i, j);

    if (vrna_mfe_profile_enabled())
      vrna_mfe_profile_add_internal(vrna_mfe_profile_now() - t0);

    return ret;
  }

  return INF;
}


PUBLIC int
vrna_mfe_internal_ext(vrna_fold_compound_t  *fc,
                      unsigned int          i,
                      unsigned int          j,
                      unsigned int          *ip,
                      unsigned int          *iq)
{
  if (fc)
    return mfe_internal_loop_ext(fc, i, j, ip, iq);

  return INF;
}


/*
 #####################################
 # BEGIN OF STATIC HELPER FUNCTIONS  #
 #####################################
 */
PRIVATE INLINE unsigned int
ptype_fast(char          *ptype,
           unsigned int  ij)
{
  unsigned int tt = (unsigned int)ptype[ij];

  return (tt == 0) ? 7U : tt;
}


PRIVATE INLINE unsigned char
use_internal_single_fast_path(vrna_fold_compound_t *fc)
{
  vrna_md_t *md;
  vrna_sc_t *sc;

  if ((!fc) ||
      (fc->type != VRNA_FC_TYPE_SINGLE) ||
      (fc->strands != 1) ||
      (!fc->params) ||
      (!fc->hc) ||
      (!fc->matrices) ||
      (fc->hc->type == VRNA_HC_WINDOW) ||
      (fc->domains_up))
    return 0;

  md = &(fc->params->model_details);
  sc = fc->sc;

  if ((md->circ) ||
      (md->gquad) ||
      (fc->hc->f))
    return 0;

  if (!sc)
    return 1;

  if ((sc->f) ||
      (sc->exp_f) ||
      (sc->bt) ||
      (sc->energy_up) ||
      (sc->energy_bp) ||
      (sc->energy_bp_local) ||
      (sc->energy_stack))
    return 0;

  return 1;
}


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


PRIVATE INLINE unsigned char
hc_int_linear_fast(vrna_hc_t      *hc,
                   unsigned int   n,
                   unsigned int   i,
                   unsigned int   j,
                   unsigned int   k,
                   unsigned int   l)
{
  unsigned char pij, pkl;
  unsigned int  u1, u2;

  pij = hc->mx[n * i + j];
  pkl = hc->mx[n * k + l];

  if (!(pij & VRNA_CONSTRAINT_CONTEXT_INT_LOOP))
    return 0;

  if (!(pkl & VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC))
    return 0;

  u1 = k - i - 1;
  u2 = j - l - 1;

  if ((u1 != 0) && (hc->up_int[i + 1] < u1))
    return 0;

  if ((u2 != 0) && (hc->up_int[l + 1] < u2))
    return 0;

  return 1;
}


PRIVATE INLINE int
mfe_E_internal_single_fast(unsigned int  n1,
                           unsigned int  n2,
                           unsigned int  type,
                           unsigned int  type_2,
                           int           si1,
                           int           sj1,
                           int           sp1,
                           int           sq1,
                           vrna_param_t  *P)
{
  unsigned int  nl, ns, u, backbones, no_close;
  int           energy, salt_stack_correction, salt_loop_correction;

  no_close              = 0;
  salt_stack_correction = P->SaltStack;
  salt_loop_correction  = 0;

  if ((P->model_details.noGUclosure) &&
      ((type_2 == 3) || (type_2 == 4) || (type == 3) || (type == 4)))
    no_close = 1;

  if (n1 > n2) {
    nl  = n1;
    ns  = n2;
  } else {
    nl  = n2;
    ns  = n1;
  }

  if (nl == 0)
    return P->stack[type][type_2] + salt_stack_correction;

  if (no_close)
    return INF;

  if (P->model_details.salt != VRNA_MODEL_DEFAULT_SALT) {
    backbones = nl + ns + 2;
    if (backbones <= MAXLOOP + 1)
      salt_loop_correction = P->SaltLoop[backbones];
    else
      salt_loop_correction = vrna_salt_loop_int((int)backbones,
                                                P->model_details.salt,
                                                P->temperature + K0,
                                                P->model_details.backbone_length);
  }

  energy = 0;

  switch (ns) {
    case 0:
      energy = (nl <= MAXLOOP) ?
               P->bulge[nl] :
               (P->bulge[30] + (int)(P->lxc * log(nl / 30.)));
      if (nl == 1) {
        energy += P->stack[type][type_2];
      } else {
        if (type > 2)
          energy += P->TerminalAU;

        if (type_2 > 2)
          energy += P->TerminalAU;
      }
      break;

    case 1:
      if (nl == 1) {
        energy = P->int11[type][type_2][si1][sj1];
      } else if (nl == 2) {
        if (n1 == 1)
          energy = P->int21[type][type_2][si1][sq1][sj1];
        else
          energy = P->int21[type_2][type][sq1][si1][sp1];
      } else {
        energy = (nl + 1 <= MAXLOOP) ?
                 P->internal_loop[nl + 1] :
                 (P->internal_loop[30] + (int)(P->lxc * log((nl + 1) / 30.)));
        energy += MIN2(MAX_NINIO, (int)(nl - ns) * P->ninio[2]);
        energy += P->mismatch1nI[type][si1][sj1] +
                  P->mismatch1nI[type_2][sq1][sp1];
      }
      break;

    case 2:
      if (nl == 2) {
        energy = P->int22[type][type_2][si1][sp1][sq1][sj1];
        break;
      } else if (nl == 3) {
        energy = P->internal_loop[5] +
                 P->ninio[2];
        energy += P->mismatch23I[type][si1][sj1] +
                  P->mismatch23I[type_2][sq1][sp1];
        break;
      }
      /* fall through */

    default:
      u       = nl + ns;
      energy  = (u <= MAXLOOP) ?
                P->internal_loop[u] :
                (P->internal_loop[30] + (int)(P->lxc * log(u / 30.)));
      energy += MIN2(MAX_NINIO, (int)(nl - ns) * P->ninio[2]);
      energy += P->mismatchI[type][si1][sj1] +
                P->mismatchI[type_2][sq1][sp1];
      break;
  }

  return energy + salt_loop_correction;
}


PRIVATE void
prepare_intloop_helpers(helper_data_t        *h,
                        vrna_fold_compound_t  *fc,
                        unsigned int          i,
                        unsigned int          j)
{
  /* init soft constraints wrapper */
  init_sc_int(fc, &(h->sc_wrapper));
  h->owns_tt = 0;


  if (fc->type == VRNA_FC_TYPE_COMPARATIVE) {
    vrna_mfe_scratch_t  *scratch = (vrna_mfe_scratch_t *)fc->matrices->aux_mfe;
    vrna_md_t           *md = &(fc->params->model_details);
    short               **S = fc->S;

    if (use_comparative_sc_fallback(fc)) {
      h->tt      = vrna_alloc(sizeof(unsigned int) * fc->n_seq);
      h->tt_size = fc->n_seq;
      h->owns_tt = 1;
    } else {
      h->tt      = scratch->int_tt;
      h->tt_size = scratch->int_tt_size;
    }

    for (unsigned int s = 0; s < fc->n_seq; s++)
      h->tt[s] = vrna_get_ptype_md(S[s][i], S[s][j], md);
  } else {
    h->tt = NULL;
  }
}


PRIVATE void
free_intloop_helpers(helper_data_t *h)
{
  free_sc_int(&(h->sc_wrapper));
  if (h->owns_tt)
    free(h->tt);
  h->tt        = NULL;
  h->tt_size   = 0;
  h->owns_tt   = 0;
}


PRIVATE int
mfe_stacks(vrna_fold_compound_t *fc,
           unsigned int         i,
           unsigned int         j,
           helper_data_t        *helpers)
{
  unsigned char sliding_window;
  char          *ptype, **ptype_local;
  short         *S, **SS, **S5, **S3;
  unsigned int  k, l, *sn, n_seq, s, n, type, type2;
  int           e, eee, *idx, ij, *c, *rtype, **c_local, kl;
  vrna_param_t  *P;
  vrna_md_t     *md;
  vrna_hc_t     *hc;

  e = INF;

  n               = fc->length;
  sliding_window  = (fc->hc->type == VRNA_HC_WINDOW) ? 1 : 0;
  sn              = fc->strand_number;
  n_seq           = (fc->type == VRNA_FC_TYPE_SINGLE) ? 1 : fc->n_seq;
  idx             = fc->jindx;
  ij              = (sliding_window) ? 0 : idx[j] + i;
  ptype           = (fc->type == VRNA_FC_TYPE_SINGLE) ? (sliding_window ? NULL : fc->ptype) : NULL;
  ptype_local     =
    (fc->type == VRNA_FC_TYPE_SINGLE) ? (sliding_window ? fc->ptype_local : NULL) : NULL;
  S       = (fc->type == VRNA_FC_TYPE_SINGLE) ? fc->sequence_encoding : NULL;
  SS      = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S;
  S5      = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S5;
  S3      = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S3;
  c       = (sliding_window) ? NULL : fc->matrices->c;
  c_local = (sliding_window) ? fc->matrices->c_local : NULL;
  P       = fc->params;
  md      = &(P->model_details);
  rtype   = &(md->rtype[0]);
  hc      = fc->hc;

  k = i + 1;
  l = j - 1;

  if (k < l) {
    if (hc->eval_int(i, j, k, l, hc)) {
      type = 0;

      if (fc->type == VRNA_FC_TYPE_SINGLE)
        type = sliding_window ?
               vrna_get_ptype_window(i, j, ptype_local) :
               vrna_get_ptype(ij, ptype);

      kl  = (sliding_window) ? 0 : idx[l] + k;
      eee = (sliding_window) ? c_local[k][l - k] : c[kl];

      if (eee != INF) {
        switch (fc->type) {
          case VRNA_FC_TYPE_SINGLE:
            type2 = sliding_window ?
                    rtype[vrna_get_ptype_window(k, l, ptype_local)] :
                    rtype[vrna_get_ptype(kl, ptype)];

            eee += vrna_E_internal(0, 0, type, type2, S[i + 1], S[j - 1], S[i], S[j], P);
            break;

          case VRNA_FC_TYPE_COMPARATIVE:
            for (s = 0; s < n_seq; s++) {
              type2 = vrna_get_ptype_md(SS[s][l], SS[s][k], md);
              eee   += vrna_E_internal(0,
                                       0,
                                       helpers->tt[s],
                                       type2,
                                       S3[s][i],
                                       S5[s][j],
                                       S5[s][k],
                                       S3[s][l],
                                       P);
            }

            break;
        }
        if (helpers->sc_wrapper.pair)
          eee += helpers->sc_wrapper.pair(i, j, k, l, &(helpers->sc_wrapper));

        e = MIN2(e, eee);
      }
    }
  }

  return e;
}


PRIVATE int
mfe_bulges(vrna_fold_compound_t *fc,
           unsigned int         i,
           unsigned int         j,
           helper_data_t        *helpers)
{
  unsigned char sliding_window, hc_decompose, *hc_mx, **hc_mx_local;
  char          *ptype, **ptype_local;
  short         *S, **SS, **S5, **S3;
  unsigned int  *sn, **a2s, n_seq, s, n, *hc_up, with_ud, noclose,
                type, type2, has_nick, k, l, last_k, first_l, u1, u2,
                noGUclosure, u1_local;
  int           e, eee, *idx, ij, kl, *c, *rtype, **c_local;
  vrna_param_t  *P;
  vrna_md_t     *md;
  vrna_ud_t     *domains_up;
  vrna_hc_t     *hc;

  e = INF;

  n               = fc->length;
  sliding_window  = (fc->hc->type == VRNA_HC_WINDOW) ? 1 : 0;
  sn              = fc->strand_number;
  n_seq           = (fc->type == VRNA_FC_TYPE_SINGLE) ? 1 : fc->n_seq;
  idx             = fc->jindx;
  ij              = (sliding_window) ? 0 : idx[j] + i;
  hc              = fc->hc;
  hc_mx           = (sliding_window) ? NULL : hc->mx;
  hc_mx_local     = (sliding_window) ? hc->matrix_local : NULL;
  hc_up           = fc->hc->up_int;
  ptype           = (fc->type == VRNA_FC_TYPE_SINGLE) ? (sliding_window ? NULL : fc->ptype) : NULL;
  ptype_local     =
    (fc->type == VRNA_FC_TYPE_SINGLE) ? (sliding_window ? fc->ptype_local : NULL) : NULL;
  S           = (fc->type == VRNA_FC_TYPE_SINGLE) ? fc->sequence_encoding : NULL;
  SS          = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S;
  S5          = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S5;
  S3          = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S3;
  a2s         = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->a2s;
  c           = (sliding_window) ? NULL : fc->matrices->c;
  c_local     = (sliding_window) ? fc->matrices->c_local : NULL;
  P           = fc->params;
  md          = &(P->model_details);
  rtype       = &(md->rtype[0]);
  domains_up  = fc->domains_up;
  with_ud     = ((domains_up) && (domains_up->energy_cb)) ? 1 : 0;

  hc_decompose = (sliding_window) ? hc_mx_local[i][j - i] : hc_mx[n * i + j];

  if (hc_decompose & VRNA_CONSTRAINT_CONTEXT_INT_LOOP) {
    has_nick    = sn[i] != sn[j] ? 1 : 0;
    noGUclosure = md->noGUclosure;
    type        = 0;

    if (fc->type == VRNA_FC_TYPE_SINGLE)
      type = sliding_window ?
             vrna_get_ptype_window(i, j, ptype_local) :
             vrna_get_ptype(ij, ptype);

    noclose = ((noGUclosure) && (type == 3 || type == 4)) ? 1 : 0;

    if (!noclose) {
      /* only proceed if the enclosing pair is allowed */

      /* handle bulges in 5' side */
      l = j - 1;
      if (l > i + 2) {
        last_k = l - 1;

        if (last_k > i + 1 + MAXLOOP)
          last_k = i + 1 + MAXLOOP;

        if (last_k > i + 1 + hc_up[i + 1])
          last_k = i + 1 + hc_up[i + 1];

        u1 = 1;

        k   = i + 2;
        kl  = (sliding_window) ? 0 : idx[l] + k;

        hc_mx += n * l;

        for (; k <= last_k; k++, u1++, kl++) {
          if (hc->eval_int(i, j, k, l, hc)) {
            eee = (sliding_window) ? c_local[k][l - k] : c[kl];

            if (eee < INF) {
              switch (fc->type) {
                case VRNA_FC_TYPE_SINGLE:
                  type2 = sliding_window ?
                          rtype[vrna_get_ptype_window(k, l, ptype_local)] :
                          rtype[vrna_get_ptype(kl, ptype)];

                  if ((noGUclosure) && (type2 == 3 || type2 == 4))
                    continue;

                  if ((has_nick) && ((sn[i] != sn[k]) || (sn[j - 1] != sn[j]))) {
                    eee = INF;
                  } else {
                    eee += vrna_E_internal(u1,
                                           0,
                                           type,
                                           type2,
                                           S[i + 1],
                                           S[j - 1],
                                           S[k - 1],
                                           S[l + 1],
                                           P);
                  }

                  break;

                case VRNA_FC_TYPE_COMPARATIVE:
                  for (s = 0; s < n_seq; s++) {
                    u1_local  = a2s[s][k - 1] - a2s[s][i];
                    type2     = vrna_get_ptype_md(SS[s][l], SS[s][k], md);
                    eee       +=
                      vrna_E_internal(u1_local, 0, helpers->tt[s], type2, S3[s][i], S5[s][j],
                                      S5[s][k], S3[s][l],
                                      P);
                  }

                  break;
              }

              if (helpers->sc_wrapper.pair)
                eee += helpers->sc_wrapper.pair(i, j, k, l, &(helpers->sc_wrapper));

              e = MIN2(e, eee);

              if (with_ud) {
                eee += domains_up->energy_cb(fc,
                                             i + 1, k - 1,
                                             VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP,
                                             domains_up->data);
                e = MIN2(e, eee);
              }
            }
          }
        }

        hc_mx -= n * l;
      }

      /* handle bulges in 3' side */
      k = i + 1;
      if (k + 2 < j) {
        first_l = k + 1;
        if (first_l + MAXLOOP + 1 < j)
          first_l = j - 1 - MAXLOOP;

        u2    = 1;
        hc_mx += n * k;

        for (l = j - 2; l >= first_l; l--, u2++) {
          if (u2 > hc_up[l + 1])
            break;

          kl            = (sliding_window) ? 0 : idx[l] + k;
          if (hc->eval_int(i, j, k, l, hc)) {
            eee = (sliding_window) ? c_local[k][l - k] : c[kl];

            if (eee < INF) {
              switch (fc->type) {
                case VRNA_FC_TYPE_SINGLE:
                  type2 = sliding_window ?
                          rtype[vrna_get_ptype_window(k, l, ptype_local)] :
                          rtype[vrna_get_ptype(kl, ptype)];

                  if ((noGUclosure) && (type2 == 3 || type2 == 4))
                    continue;

                  if ((has_nick) && ((sn[i] != sn[i + 1]) || (sn[j] != sn[l]))) {
                    eee = INF;
                  } else {
                    eee += vrna_E_internal(0,
                                           u2,
                                           type,
                                           type2,
                                           S[i + 1],
                                           S[j - 1],
                                           S[k - 1],
                                           S[l + 1],
                                           P);
                  }

                  break;

                case VRNA_FC_TYPE_COMPARATIVE:
                  for (s = 0; s < n_seq; s++) {
                    unsigned int u2_local = a2s[s][j - 1] - a2s[s][l];
                    type2 = vrna_get_ptype_md(SS[s][l], SS[s][k], md);
                    eee   +=
                      vrna_E_internal(0, u2_local, helpers->tt[s], type2, S3[s][i], S5[s][j],
                                      S5[s][k], S3[s][l],
                                      P);
                  }

                  break;
              }

              if (helpers->sc_wrapper.pair)
                eee += helpers->sc_wrapper.pair(i, j, k, l, &(helpers->sc_wrapper));

              e = MIN2(e, eee);

              if (with_ud) {
                eee += domains_up->energy_cb(fc,
                                             l + 1, j - 1,
                                             VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP,
                                             domains_up->data);
                e = MIN2(e, eee);
              }
            }
          }
        }

        hc_mx -= n * k;
      }
    }
  }

  return e;
}


PRIVATE int
mfe_internal_loop_single_fast(vrna_fold_compound_t  *fc,
                              unsigned int          i,
                              unsigned int          j)
{
  char          *ptype;
  short         *S;
  int           *idx, *rtype, *c, kl, e, eee;
  unsigned int  *hc_up, n, ij, k, l, last_k, first_l, u1, u2, type, type2,
                noGUclosure, hc_decompose;
  vrna_param_t  *P;
  vrna_md_t     *md;
  vrna_hc_t     *hc;

  n             = fc->length;
  idx           = fc->jindx;
  ij            = idx[j] + i;
  hc            = fc->hc;
  hc_up         = hc->up_int;
  ptype         = fc->ptype;
  S             = fc->sequence_encoding;
  c             = fc->matrices->c;
  P             = fc->params;
  md            = &(P->model_details);
  rtype         = &(md->rtype[0]);
  noGUclosure   = md->noGUclosure;
  hc_decompose  = hc->mx[n * i + j];
  e             = INF;

  if (!(hc_decompose & VRNA_CONSTRAINT_CONTEXT_INT_LOOP))
    return INF;

  type = ptype_fast(ptype, ij);

  /* stacks */
  k = i + 1;
  l = j - 1;
  if ((k < l) && hc_int_linear_fast(hc, n, i, j, k, l)) {
    kl  = idx[l] + k;
      eee = c[kl];
      if (eee != INF) {
        type2 = rtype[ptype_fast(ptype, kl)];
        eee  += mfe_E_internal_single_fast(0,
                                           0,
                                           type,
                                           type2,
                                           S[i + 1],
                                           S[j - 1],
                                           S[i],
                                           S[j],
                                           P);
        e     = MIN2(e, eee);
      }
  }

  /* bulges */
  if (!((noGUclosure) && ((type == 3) || (type == 4)))) {
    l = j - 1;
    if (l > i + 2) {
      last_k = l - 1;

      if (last_k > i + 1 + MAXLOOP)
        last_k = i + 1 + MAXLOOP;

      if (last_k > i + 1 + hc_up[i + 1])
        last_k = i + 1 + hc_up[i + 1];

      u1 = 1;
      k  = i + 2;
      kl = idx[l] + k;

      for (; k <= last_k; k++, u1++, kl++) {
        if (!hc_int_linear_fast(hc, n, i, j, k, l))
          continue;

        eee = c[kl];
        if (eee == INF)
          continue;

        type2 = rtype[ptype_fast(ptype, kl)];
        if ((noGUclosure) && ((type2 == 3) || (type2 == 4)))
          continue;

        eee  += mfe_E_internal_single_fast(u1,
                                           0,
                                           type,
                                           type2,
                                           S[i + 1],
                                           S[j - 1],
                                           S[k - 1],
                                           S[l + 1],
                                           P);
        e     = MIN2(e, eee);
      }
    }

    k = i + 1;
    if (k + 2 < j) {
      first_l = k + 1;
      if (first_l + MAXLOOP + 1 < j)
        first_l = j - 1 - MAXLOOP;

      u2 = 1;
      for (l = j - 2; l >= first_l; l--, u2++) {
        if (u2 > hc_up[l + 1])
          break;

        kl = idx[l] + k;
        if (!hc_int_linear_fast(hc, n, i, j, k, l))
          continue;

        eee = c[kl];
        if (eee == INF)
          continue;

        type2 = rtype[ptype_fast(ptype, kl)];
        if ((noGUclosure) && ((type2 == 3) || (type2 == 4)))
          continue;

        eee  += mfe_E_internal_single_fast(0,
                                           u2,
                                           type,
                                           type2,
                                           S[i + 1],
                                           S[j - 1],
                                           S[k - 1],
                                           S[l + 1],
                                           P);
        e     = MIN2(e, eee);
      }
    }

    /* all other internal loops */
    first_l = i + 3;
    if (first_l + MAXLOOP + 1 < j)
      first_l = j - 1 - MAXLOOP;

    u2 = 1;
    for (l = j - 2; l >= first_l; l--, u2++) {
      if (u2 > hc_up[l + 1])
        break;

      last_k = l - 1;
      if (last_k + u2 > i + 1 + MAXLOOP)
        last_k = i + 1 + MAXLOOP - u2;

      if (last_k > i + 1 + hc_up[i + 1])
        last_k = i + 1 + hc_up[i + 1];

      u1 = 1;
      k  = i + 2;
      kl = idx[l] + k;

      for (; k <= last_k; k++, u1++, kl++) {
        if (!hc_int_linear_fast(hc, n, i, j, k, l))
          continue;

        eee = c[kl];
        if (eee == INF)
          continue;

        type2 = rtype[ptype_fast(ptype, kl)];
        if ((noGUclosure) && ((type2 == 3) || (type2 == 4)))
          continue;

        eee  += mfe_E_internal_single_fast(u1,
                                           u2,
                                           type,
                                           type2,
                                           S[i + 1],
                                           S[j - 1],
                                           S[k - 1],
                                           S[l + 1],
                                           P);
        e     = MIN2(e, eee);
      }
    }
  }

  return e;
}


PRIVATE int
mfe_internal_loop(vrna_fold_compound_t  *fc,
                  unsigned int          i,
                  unsigned int          j)
{
  unsigned char sliding_window, hc_decompose, *hc_mx, **hc_mx_local;
  char          *ptype, **ptype_local;
  short         *S, **SS, **S5, **S3;
  unsigned int  *sn, **a2s, n_seq, s, n, *hc_up, type, type2, has_nick,
                k, l, last_k, first_l, u1, u2, noGUclosure, with_ud,
                with_gquad, noclose, u1_local, u2_local;
  int           e, eee, e3, e5, *idx, ij, kl, *c, *rtype, **c_local;
  vrna_param_t  *P;
  vrna_md_t     *md;
  vrna_ud_t     *domains_up;
  vrna_hc_t     *hc;
  helper_data_t helpers = {
    0
  };

  e = INF;

  n               = fc->length;
  sliding_window  = (fc->hc->type == VRNA_HC_WINDOW) ? 1 : 0;
  sn              = fc->strand_number;
  n_seq           = (fc->type == VRNA_FC_TYPE_SINGLE) ? 1 : fc->n_seq;
  idx             = fc->jindx;
  ij              = (sliding_window) ? 0 : idx[j] + i;
  hc              = fc->hc;
  hc_mx           = (sliding_window) ? NULL : hc->mx;
  hc_mx_local     = (sliding_window) ? hc->matrix_local : NULL;
  hc_up           = fc->hc->up_int;
  ptype           = (fc->type == VRNA_FC_TYPE_SINGLE) ? (sliding_window ? NULL : fc->ptype) : NULL;
  ptype_local     =
    (fc->type == VRNA_FC_TYPE_SINGLE) ? (sliding_window ? fc->ptype_local : NULL) : NULL;
  S           = (fc->type == VRNA_FC_TYPE_SINGLE) ? fc->sequence_encoding : NULL;
  SS          = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S;
  S5          = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S5;
  S3          = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S3;
  a2s         = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->a2s;
  c           = (sliding_window) ? NULL : fc->matrices->c;
  c_local     = (sliding_window) ? fc->matrices->c_local : NULL;
  P           = fc->params;
  md          = &(P->model_details);
  rtype       = &(md->rtype[0]);
  domains_up  = fc->domains_up;
  with_ud     = ((domains_up) && (domains_up->energy_cb)) ? 1 : 0;
  with_gquad  = md->gquad;

  prepare_intloop_helpers(&helpers, fc, i, j);

  hc_decompose = (sliding_window) ? hc_mx_local[i][j - i] : hc_mx[n * i + j];

  if (hc_decompose & VRNA_CONSTRAINT_CONTEXT_INT_LOOP) {
    has_nick    = sn[i] != sn[j] ? 1 : 0;
    noGUclosure = md->noGUclosure;
    type        = 0;

    if (fc->type == VRNA_FC_TYPE_SINGLE)
      type = sliding_window ?
             vrna_get_ptype_window(i, j, ptype_local) :
             vrna_get_ptype(ij, ptype);

    noclose = ((noGUclosure) && (type == 3 || type == 4)) ? 1 : 0;

    e = MIN2(e, mfe_stacks(fc, i, j, &helpers));

    e = MIN2(e, mfe_bulges(fc, i, j, &helpers));

    if (!noclose) {
      /* only proceed if the enclosing pair is allowed */

      /* last but not least, all other internal loops */
      first_l = i + 2 + 1;
      if (first_l + MAXLOOP + 1 < j)
        first_l = j - 1 - MAXLOOP;

      u2 = 1;
      for (l = j - 2; l >= first_l; l--, u2++) {
        if (u2 > hc_up[l + 1])
          break;

        last_k = l - 1;

        if (last_k + u2 > i + 1 + MAXLOOP)
          last_k = i + 1 + MAXLOOP - u2;

        if (last_k > i + 1 + hc_up[i + 1])
          last_k = i + 1 + hc_up[i + 1];

        u1  = 1;
        k   = i + 2;
        kl  = (sliding_window) ? 0 : idx[l] + k;

        hc_mx += n * l;

        for (; k <= last_k; k++, u1++, kl++) {
          hc_decompose = (sliding_window) ? hc_mx_local[k][l - k] : hc_mx[k];

          if (hc->eval_int(i, j, k, l, hc)) {
            eee = (sliding_window) ? c_local[k][l - k] : c[kl];

            if (eee < INF) {
              switch (fc->type) {
                case VRNA_FC_TYPE_SINGLE:
                  type2 = sliding_window ?
                          rtype[vrna_get_ptype_window(k, l, ptype_local)] :
                          rtype[vrna_get_ptype(kl, ptype)];

                  if ((noGUclosure) && (type2 == 3 || type2 == 4))
                    continue;

                  if ((has_nick) && ((sn[i] != sn[k]) || (sn[j] != sn[l]))) {
                    eee = INF;
                  } else {
                    eee +=
                      vrna_E_internal(u1, u2, type, type2, S[i + 1], S[j - 1], S[k - 1], S[l + 1],
                                      P);
                  }

                  break;

                case VRNA_FC_TYPE_COMPARATIVE:
                  for (s = 0; s < n_seq; s++) {
                    u1_local  = a2s[s][k - 1] - a2s[s][i];
                    u2_local  = a2s[s][j - 1] - a2s[s][l];
                    type2     = vrna_get_ptype_md(SS[s][l], SS[s][k], md);
                    eee       += vrna_E_internal(u1_local,
                                                 u2_local,
                                                 helpers.tt[s],
                                                 type2,
                                                 S3[s][i],
                                                 S5[s][j],
                                                 S5[s][k],
                                                 S3[s][l],
                                                 P);
                  }

                  break;
              }

              if (helpers.sc_wrapper.pair)
                eee += helpers.sc_wrapper.pair(i, j, k, l, &(helpers.sc_wrapper));

              e = MIN2(e, eee);

              if (with_ud) {
                e5 = domains_up->energy_cb(fc,
                                           i + 1, k - 1,
                                           VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP,
                                           domains_up->data);
                e3 = domains_up->energy_cb(fc,
                                           l + 1, j - 1,
                                           VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP,
                                           domains_up->data);

                e = MIN2(e, eee + e5);
                e = MIN2(e, eee + e3);
                e = MIN2(e, eee + e5 + e3);
              }
            }
          }
        }

        hc_mx -= n * l;
      }

      if (with_gquad) {
        /* include all cases where a g-quadruplex may be enclosed by base pair (i,j) */
        eee = vrna_mfe_gquad_internal_loop(fc, i, j);
        e   = MIN2(e, eee);
      }
    }
  }

  free_intloop_helpers(&helpers);

  return e;
}


PRIVATE int
mfe_internal_loop_ext(vrna_fold_compound_t  *fc,
                      unsigned int          i,
                      unsigned int          j,
                      unsigned int          *ip,
                      unsigned int          *iq)
{
  int                   energy, e, *indx, *c;
  unsigned char         *hc_mx;
  unsigned int          n, q, p, s, n_seq, u1, u2, qmin, *tt, *hc_up;
  short                 **SS;
  vrna_md_t             *md;
  vrna_param_t          *P;
  vrna_hc_t             *hc;
  vrna_mfe_scratch_t    *scratch;

  n     = fc->length;
  n_seq = (fc->type == VRNA_FC_TYPE_SINGLE) ? 1 : fc->n_seq;
  SS    = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S;
  indx  = fc->jindx;
  c     = fc->matrices->c;
  hc    = fc->hc;
  hc_mx = hc->mx;
  hc_up = hc->up_int;
  P     = fc->params;
  md    = &(P->model_details);
  tt    = NULL;
  scratch = (vrna_mfe_scratch_t *)fc->matrices->aux_mfe;

  e = INF;

  /* CONSTRAINED internal LOOP start */
  if (hc_mx[n * i + j] & (VRNA_CONSTRAINT_CONTEXT_INT_LOOP | VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC)) {
    /* prepare necessary variables */
    if (fc->type == VRNA_FC_TYPE_COMPARATIVE) {
      if (use_comparative_sc_fallback(fc))
        tt = vrna_alloc(sizeof(unsigned int) * n_seq);
      else
        tt = scratch->int_ext_tt;

      for (s = 0; s < n_seq; s++)
        tt[s] = vrna_get_ptype_md(SS[s][j], SS[s][i], md);
    }

    for (p = j + 1; p < n; p++) {
      u1 = p - j - 1;
      if (u1 + i - 1 > MAXLOOP)
        break;

      if (hc_up[j + 1] < u1)
        break;

      qmin = p + 1;

      if (u1 + i + n > p + MAXLOOP + 2)
        qmin = u1 + i + n - MAXLOOP - 1;

      for (q = n; q >= qmin; q--) {
        u2 = i - 1 + n - q;
        if (hc_up[q + 1] < u2)
          break;

        if (u1 + u2 > MAXLOOP)
          continue;

        if (hc->eval_int(i, j, p, q, hc)) {
          int pq = indx[q] + p;
          energy = c[pq];
          if (energy < INF) {
            energy += vrna_eval_internal(fc, i, j, p, q, VRNA_EVAL_LOOP_NO_HC);

            if (energy < e) {
              e = energy;
              if ((ip != NULL) && (iq != NULL)) {
                *ip = p;
                *iq = q;
              }
            }
          }
        }
      }
    }
  }

  if ((fc->type == VRNA_FC_TYPE_COMPARATIVE) &&
      use_comparative_sc_fallback(fc))
    free(tt);

  return e;
}


/*
 * ###########################################
 * # deprecated functions below              #
 *###########################################
 */

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

PUBLIC int
vrna_E_int_loop(vrna_fold_compound_t  *fc,
                int                   i,
                int                   j)
{
  return vrna_mfe_internal(fc,
                           (unsigned int)i,
                           (unsigned int)j);
}


PUBLIC int
vrna_E_ext_int_loop(vrna_fold_compound_t  *fc,
                    int                   i,
                    int                   j,
                    int                   *ip,
                    int                   *iq)
{
  int e = INF;

  if (fc) {
    unsigned int p, q;
    e = vrna_mfe_internal_ext(fc,
                              (unsigned int)i,
                              (unsigned int)j,
                              &p,
                              &q);
    if (ip)
      *ip = (int)p;

    if (iq)
      *iq = (int)q;
  }

  return e;
}


#endif
