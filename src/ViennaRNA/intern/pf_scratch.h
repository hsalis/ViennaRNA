#ifndef VIENNA_RNA_PACKAGE_INTERN_PF_SCRATCH_H
#define VIENNA_RNA_PACKAGE_INTERN_PF_SCRATCH_H

#include "ViennaRNA/datastructures/basic.h"
#include "ViennaRNA/partfunc/exterior.h"
#include "ViennaRNA/partfunc/multibranch.h"

typedef struct {
  vrna_mx_pf_aux_el_t ext_helpers;
  vrna_mx_pf_aux_ml_t ml_helpers;
  FLT_OR_DBL          *qm1_tmp;
  unsigned int        *int_tt;
  unsigned int        *ext_int_tt;
} vrna_pf_scratch_t;

#endif
