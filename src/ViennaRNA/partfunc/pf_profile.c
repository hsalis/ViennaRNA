#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <string.h>
#include <time.h>

#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/intern/pf_profile.h"


static vrna_pf_profile_t profile_data;
static int               profile_enabled = 0;


PUBLIC uint64_t
vrna_pf_profile_now(void)
{
  struct timespec ts;

  clock_gettime(CLOCK_MONOTONIC, &ts);

  return ((uint64_t)ts.tv_sec * 1000000000ULL) + (uint64_t)ts.tv_nsec;
}


PUBLIC void
vrna_pf_profile_enable(int enabled)
{
  profile_enabled = enabled ? 1 : 0;
}


PUBLIC int
vrna_pf_profile_enabled(void)
{
  return profile_enabled;
}


PUBLIC void
vrna_pf_profile_reset(void)
{
  memset(&profile_data, 0, sizeof(profile_data));
}


PUBLIC const vrna_pf_profile_t *
vrna_pf_profile_get(void)
{
  return &profile_data;
}


PUBLIC void
vrna_pf_profile_add_fill_arrays(uint64_t delta_ns)
{
  if (profile_enabled)
    profile_data.fill_arrays_ns += delta_ns;
}


PUBLIC void
vrna_pf_profile_add_decompose_pair(uint64_t delta_ns)
{
  if (profile_enabled)
    profile_data.decompose_pair_ns += delta_ns;
}


PUBLIC void
vrna_pf_profile_add_internal(uint64_t delta_ns)
{
  if (profile_enabled)
    profile_data.internal_ns += delta_ns;
}


PUBLIC void
vrna_pf_profile_add_multibranch(uint64_t delta_ns)
{
  if (profile_enabled)
    profile_data.multibranch_ns += delta_ns;
}


PUBLIC void
vrna_pf_profile_add_exterior(uint64_t delta_ns)
{
  if (profile_enabled)
    profile_data.exterior_ns += delta_ns;
}


PUBLIC void
vrna_pf_profile_add_circular(uint64_t delta_ns)
{
  if (profile_enabled)
    profile_data.circular_ns += delta_ns;
}


PUBLIC void
vrna_pf_profile_add_bpp(uint64_t delta_ns)
{
  if (profile_enabled)
    profile_data.bpp_ns += delta_ns;
}


PUBLIC void
vrna_pf_profile_count_alloc(uint64_t count)
{
  if (profile_enabled)
    profile_data.allocations += count;
}
