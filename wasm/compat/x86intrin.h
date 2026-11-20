/* <x86intrin.h> for Emscripten, which has none of its own. Without this,
 * clang's x86 header is picked up and pulls in instructions WebAssembly
 * cannot express. It maps to Emscripten's SSE/AVX emulation on top of WASM
 * SIMD (<immintrin.h>, needs -msimd128 and e.g. -msse4.1), so that
 * gaussian.cpp and the VCL library compile unchanged. See
 * https://emscripten.org/docs/porting/simd.html.
 *
 * Makefile.wasm force-includes this file into gaussian.cpp, because VCL
 * includes <smmintrin.h> directly and would not see the definitions below. */
#ifndef EVMLR_COMPAT_X86INTRIN_H
#define EVMLR_COMPAT_X86INTRIN_H

#include <immintrin.h>

/* WebAssembly has no MXCSR: it cannot change the rounding mode or flush
 * denormals, and Emscripten does not declare _mm_setcsr. VCL calls it only
 * from set_control_word/no_subnormals and from its rounding fallbacks below
 * SSE4.1; the sampler reaches none of them. Trap rather than silently do
 * nothing, so a future caller fails loudly. */
#ifndef _mm_setcsr
static inline void evmlr_compat_mm_setcsr(unsigned int w) {
    (void) w;
    __builtin_trap();
}
#define _mm_setcsr evmlr_compat_mm_setcsr
#endif

#endif /* EVMLR_COMPAT_X86INTRIN_H */
