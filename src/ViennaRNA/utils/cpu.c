/*
 * A collection of useful functions to detect various CPU features
 *
 * (c) 2018 - Ronny Lorenz - ViennaRNA Package
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#if defined(_WIN32) && (defined(__x86_64__) || defined(_M_AMD64) || defined (_M_X64))
#include <immintrin.h>
#endif

#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/utils/cpu.h"

#ifdef __GNUC__
# define INLINE inline
#else
# define INLINE
#endif

/*
 *  cpuid register flags for (64bit) x86 CPUs
 *  See also https://en.wikipedia.org/wiki/CPUID
 */

#define bit_SSE2      (1 << 26) /* stored in EDX after cpuid with EAX=1 */
#define bit_SSE3      (1 << 0)  /* stored in ECX after cpuid with EAX=1 */
#define bit_SSE41     (1 << 19) /* stored in ECX after cpuid with EAX=1 */
#define bit_SSE42     (1 << 20) /* stored in ECX after cpuid with EAX=1 */
#define bit_AVX       (1 << 28) /* stored in ECX after cpuid with EAX=1 */
#define bit_AVX2      (1 << 5)  /* stored in EBX after cpuid with EAX=7, ECX=0 */
#define bit_AVX512F   (1 << 16) /* stored in EBX after cpuid with EAX=7, ECX=0 */

/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */
PRIVATE INLINE int
execute_cpuid(uint32_t *regs);


PRIVATE INLINE uint64_t
xgetbv(uint32_t idx);


PRIVATE unsigned int
cpu_feature_bits(void);


PRIVATE unsigned int
cpu_extended_feature_bits(void);


PRIVATE int
os_supports_avx_state(void);


PRIVATE int
os_supports_avx512_state(void);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC char *
vrna_cpu_vendor_string(void)
{
  static char name[13] = {
    0
  };
  uint32_t    regs[4] = {
    0
  };

  if (execute_cpuid(&regs[0])) {
    memcpy(name + 0, &regs[1], 4);
    memcpy(name + 4, &regs[3], 4);
    memcpy(name + 8, &regs[2], 4);
    name[12] = '\0';
  }

  return name;
}


PUBLIC unsigned int
vrna_cpu_simd_capabilities(void)
{
  unsigned int capabilities = VRNA_CPU_SIMD_NONE;

#if defined(__ARM_NEON) || defined(__ARM_NEON__) || defined(__aarch64__)
  capabilities |= VRNA_CPU_SIMD_NEON;
#endif
  capabilities  |= cpu_feature_bits();
  capabilities  |= cpu_extended_feature_bits();

  return capabilities;
}


/*
 #################################
 # STATIC helper functions below #
 #################################
 */

/*
 * execute 'cpuid' command with values stored in registers regs[0]..regs[3]
 * that directly correspond to EAX, EBX, ECX, and EDX. The result of the
 * 'cpuid' command is then returned in the same register array
 */
PRIVATE INLINE int
execute_cpuid(uint32_t *regs)
{
#if defined(__x86_64__) || defined(_M_AMD64) || defined (_M_X64)
# ifdef _WIN32
#   ifndef __MINGW32__
  __cpuidex(regs, regs[0], regs[2]);
#   else
  __asm__ __volatile__ ("cpuid"
                        : "=a" (regs[0]),
                        "=b" (regs[1]),
                        "=c" (regs[2]),
                        "=d" (regs[3])
                        : "0" (regs[0]),
                        "2" (regs[2]));
#   endif
# else
  __asm__ __volatile__ ("cpuid"
                        : "=a" (regs[0]),
                        "=b" (regs[1]),
                        "=c" (regs[2]),
                        "=d" (regs[3])
                        : "0" (regs[0]),
                        "2" (regs[2]));
# endif
  return 1;
#else
  return 0;
#endif
}


PRIVATE INLINE uint64_t
xgetbv(uint32_t idx)
{
#if defined(__x86_64__) || defined(_M_AMD64) || defined (_M_X64)
# ifdef _WIN32
  return _xgetbv(idx);
# else
  uint32_t eax;
  uint32_t edx;

  __asm__ __volatile__(".byte 0x0f, 0x01, 0xd0"
                       : "=a" (eax), "=d" (edx)
                       : "c" (idx));

  return ((uint64_t)edx << 32) | eax;
# endif
#else
  (void)idx;
  return 0;
#endif
}


PRIVATE unsigned int
cpu_feature_bits(void)
{
  unsigned int  features = VRNA_CPU_SIMD_NONE;

  uint32_t      regs[4] = {
    1, 0, 0, 0
  };

  if (execute_cpuid(&regs[0])) {
    if (regs[3] & bit_SSE2)
      features |= VRNA_CPU_SIMD_SSE2;

    if (regs[2] & bit_SSE3)
      features |= VRNA_CPU_SIMD_SSE3;

    if (regs[2] & bit_SSE41)
      features |= VRNA_CPU_SIMD_SSE41;

    if (regs[2] & bit_SSE42)
      features |= VRNA_CPU_SIMD_SSE42;

    if ((regs[2] & bit_AVX) && os_supports_avx_state())
      features |= VRNA_CPU_SIMD_AVX;
  }

  return features;
}


PRIVATE unsigned int
cpu_extended_feature_bits(void)
{
  unsigned int  features = VRNA_CPU_SIMD_NONE;

  uint32_t      regs[4] = {
    7, 0, 0, 0
  };

  if (execute_cpuid(&regs[0]) && os_supports_avx_state()) {
    if (regs[1] & bit_AVX2)
      features |= VRNA_CPU_SIMD_AVX2;

    if ((regs[1] & bit_AVX512F) && os_supports_avx512_state())
      features |= VRNA_CPU_SIMD_AVX512F;

    if ((regs[1] & (1u << 30)) && os_supports_avx512_state())
      features |= VRNA_CPU_SIMD_AVX512BW;

    if ((regs[1] & (1u << 31)) && os_supports_avx512_state())
      features |= VRNA_CPU_SIMD_AVX512VL;
  }

  return features;
}


PRIVATE int
os_supports_avx_state(void)
{
#if defined(__x86_64__) || defined(_M_AMD64) || defined (_M_X64)
  uint32_t regs[4] = { 1, 0, 0, 0 };

  if (!execute_cpuid(&regs[0]))
    return 0;

  if (!(regs[2] & (1u << 27)))
    return 0;

  return ((xgetbv(0) & 0x6) == 0x6);
#else
  return 0;
#endif
}


PRIVATE int
os_supports_avx512_state(void)
{
#if defined(__x86_64__) || defined(_M_AMD64) || defined (_M_X64)
  if (!os_supports_avx_state())
    return 0;

  return ((xgetbv(0) & 0xe6) == 0xe6);
#else
  return 0;
#endif
}
