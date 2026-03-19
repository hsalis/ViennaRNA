#ifndef VIENNA_RNA_PACKAGE_INTERN_PF_PROFILE_H
#define VIENNA_RNA_PACKAGE_INTERN_PF_PROFILE_H

#include <stdint.h>

typedef struct {
  uint64_t  fill_arrays_ns;
  uint64_t  decompose_pair_ns;
  uint64_t  internal_ns;
  uint64_t  multibranch_ns;
  uint64_t  exterior_ns;
  uint64_t  circular_ns;
  uint64_t  bpp_ns;
  uint64_t  allocations;
} vrna_pf_profile_t;

uint64_t
vrna_pf_profile_now(void);


void
vrna_pf_profile_enable(int enabled);


int
vrna_pf_profile_enabled(void);


void
vrna_pf_profile_reset(void);


const vrna_pf_profile_t *
vrna_pf_profile_get(void);


void
vrna_pf_profile_add_fill_arrays(uint64_t delta_ns);


void
vrna_pf_profile_add_decompose_pair(uint64_t delta_ns);


void
vrna_pf_profile_add_internal(uint64_t delta_ns);


void
vrna_pf_profile_add_multibranch(uint64_t delta_ns);


void
vrna_pf_profile_add_exterior(uint64_t delta_ns);


void
vrna_pf_profile_add_circular(uint64_t delta_ns);


void
vrna_pf_profile_add_bpp(uint64_t delta_ns);


void
vrna_pf_profile_count_alloc(uint64_t count);


#endif
