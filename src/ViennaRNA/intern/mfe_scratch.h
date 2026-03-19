#ifndef VIENNA_RNA_PACKAGE_INTERN_MFE_SCRATCH_H
#define VIENNA_RNA_PACKAGE_INTERN_MFE_SCRATCH_H

#include "ViennaRNA/datastructures/basic.h"
#include "ViennaRNA/mfe/multibranch.h"

typedef struct vrna_mfe_scratch_s {
  int                   *cc;
  int                   *cc1;
  vrna_mx_mfe_aux_ml_t  ml_helpers;

  int                   *ext_stems;
  short                 *ext_tmp1;
  short                 *ext_tmp2;
  unsigned int          ext_stems_size;
  unsigned int          ext_tmp_size;

  unsigned int          *int_tt;
  unsigned int          int_tt_size;
  unsigned int          *int_ext_tt;
  unsigned int          int_ext_tt_size;

  vrna_bts_t            bt_stack;
  vrna_bps_t            bp_stack;
} vrna_mfe_scratch_t;

#endif
