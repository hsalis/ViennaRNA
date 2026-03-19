#ifndef VIENNA_RNA_PACKAGE_INTERN_MFE_PROFILE_H
#define VIENNA_RNA_PACKAGE_INTERN_MFE_PROFILE_H

#include <stdint.h>

typedef struct {
  uint64_t  fill_arrays_ns;
  uint64_t  decompose_pair_ns;
  uint64_t  internal_ns;
  uint64_t  multibranch_ns;
  uint64_t  exterior_ns;
  uint64_t  backtrack_ns;
  uint64_t  allocations;
} vrna_mfe_profile_t;

uint64_t
vrna_mfe_profile_now(void);


void
vrna_mfe_profile_enable(int enabled);


int
vrna_mfe_profile_enabled(void);


void
vrna_mfe_profile_reset(void);


const vrna_mfe_profile_t *
vrna_mfe_profile_get(void);


void
vrna_mfe_profile_add_fill_arrays(uint64_t delta_ns);


void
vrna_mfe_profile_add_decompose_pair(uint64_t delta_ns);


void
vrna_mfe_profile_add_internal(uint64_t delta_ns);


void
vrna_mfe_profile_add_multibranch(uint64_t delta_ns);


void
vrna_mfe_profile_add_exterior(uint64_t delta_ns);


void
vrna_mfe_profile_add_backtrack(uint64_t delta_ns);


void
vrna_mfe_profile_count_alloc(uint64_t count);


#endif
