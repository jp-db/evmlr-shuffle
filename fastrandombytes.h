/* ****************************** *
 * Implemented by Raymond K. ZHAO *
 *                                *
 * PRG                            *
 * ****************************** */
 
#ifndef FASTRANDOMBYTES_H
#define FASTRANDOMBYTES_H

#include "evmlr_portable.h"

#ifdef EVMLR_PORTABLE
#include <stdint.h>
#else
#include <x86intrin.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

/* State of an AES-256-CTR stream. */
#ifdef EVMLR_PORTABLE
/* Bitsliced (constant-time) AES, see fastrandombytes_portable.c: each of the
 * 15 round keys is held as 8 bitsliced 32-bit words. */
typedef struct {
	uint32_t round_key[15 * 8];
	uint8_t iv[16];
} prg_t;
#else
typedef struct {
	__m128i round_key[15];
	__m128i iv;
} prg_t;
#endif

/* Key a stream from a 32-byte seed and reset its counter. */
void prg_setseed(prg_t *prg, const unsigned char *seed);
/* Read rlen bytes from a stream. */
void prg_bytes(prg_t *prg, unsigned char *r, unsigned long long rlen);

/* The global stream that feeds the Gaussian sampler and the rejection coins.
 * It holds the secret masks, so it seeds itself from the operating system on
 * first use and must never be keyed with public data: anything derived from a
 * hash, such as a challenge, belongs in a local prg_t instead. */
void fastrandombytes(unsigned char *r, unsigned long long rlen);
/* Replace the global seed, for reproducible sampler tests only. */
void fastrandombytes_setseed(const unsigned char *randomness);

#ifdef __cplusplus
}
#endif

#endif
