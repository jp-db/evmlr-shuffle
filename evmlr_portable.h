#ifndef EVMLR_PORTABLE_H
#define EVMLR_PORTABLE_H

/* EVMLR_PORTABLE picks the portable C versions of the code that otherwise
 * needs x86: the AES-NI stream in fastrandombytes.c, the SSE4.1/VCL Gaussian
 * in gaussian.cpp, and rdtsc in cpucycles.c. The portable versions live in
 * fastrandombytes_portable.c, gaussian_portable.c and cpucycles.c.
 *
 * It is switched on for every Emscripten (WebAssembly) build, which has
 * neither AES-NI nor SSE4.1. A native build can opt in with -DEVMLR_PORTABLE
 * to compare the two versions; it is never on by default, so the native
 * build compiles the same code as before.
 *
 * The portable AES-256-CTR stream gives the same bytes as the AES-NI one for
 * the same seed, so challenges derived from a transcript match on both
 * platforms and a proof made in WebAssembly verifies natively. */
#if defined(__EMSCRIPTEN__) && !defined(EVMLR_PORTABLE)
#define EVMLR_PORTABLE 1
#endif

#endif /* EVMLR_PORTABLE_H */
