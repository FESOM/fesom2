/*
 * fesom_ftz.c -- enable flush-to-zero (FTZ) and denormals-are-zero (DAZ) on the
 * calling thread's SSE control/status register (MXCSR).
 *
 * Motivation: in single-precision EVP sea-ice dynamics, tiny strain-rate products
 * (~1e-40) are normal in double precision but SUBNORMAL in single precision. Every
 * FP op on a subnormal operand triggers a slow x86 microcode assist (~100+ cycles),
 * inflating the fixed EVP subcycle loop in SP builds. FESOM's precision-preserving
 * compiler flags (Intel -fp-model precise, GNU -ffloat-store, neither with FTZ) keep
 * denormals, so the penalty is live. Flushing subnormals to zero removes it. The
 * flushed values (< 1.2e-38) are far below any physical strain (>~1e-9) and below the
 * delta_min=1e-11 rheology clamp, so the physical effect is negligible.
 *
 * Called once per rank at startup; MXCSR is per-thread, so with OpenMP each worker
 * thread must call it -- do so from inside a parallel region if ENABLE_OPENMP is used.
 * No-op (and safe) on non-x86 architectures.
 */

#if defined(__x86_64__) || defined(__i386__)
#include <xmmintrin.h>   /* _MM_SET_FLUSH_ZERO_MODE   (FTZ, all SSE) */
#include <pmmintrin.h>   /* _MM_SET_DENORMALS_ZERO_MODE (DAZ, SSE3)  */

void fesom_set_ftz_daz(void)
{
    _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
    _MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON);
}
#else
void fesom_set_ftz_daz(void)
{
    /* no-op on non-x86; other ISAs handle subnormals in hardware without penalty */
}
#endif
