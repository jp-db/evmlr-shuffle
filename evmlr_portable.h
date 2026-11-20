#ifndef EVMLR_PORTABLE_H
#define EVMLR_PORTABLE_H

/* EVMLR_PORTABLE picks the portable AES-256-CTR stream
 * (fastrandombytes_portable.c) over the AES-NI one in fastrandombytes.c.
 *
 * It is switched on for every Emscripten (WebAssembly) build: WebAssembly
 * has no AES instruction, and Emscripten's x86 intrinsic emulation offers
 * none. A native build can opt in with -DEVMLR_PORTABLE to compare the two;
 * it is never on by default, so the native build compiles the same code as
 * before.
 *
 * The portable stream gives the same bytes as the AES-NI one for the same
 * seed, so challenges derived from a transcript match on both platforms and
 * a proof made in WebAssembly verifies natively.
 *
 * The Gaussian sampler needs no portable version: gaussian.cpp compiles
 * unchanged for WebAssembly through Emscripten's SSE4.1 emulation (see
 * wasm/compat/x86intrin.h and Makefile.wasm). */
#if defined(__EMSCRIPTEN__) && !defined(EVMLR_PORTABLE)
#define EVMLR_PORTABLE 1
#endif

#endif /* EVMLR_PORTABLE_H */
