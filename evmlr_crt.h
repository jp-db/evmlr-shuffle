#ifndef EVMLR_SHUFFLE_EVMLR_CRT_H
#define EVMLR_SHUFFLE_EVMLR_CRT_H

#include "evmlr_params.h"

/**
 * @file
 *
 * Multiplication in R_q = Z_q[x]/(x^n + 1) via the CRT decomposition.
 *
 * For q = 5 (mod 8) -- which is what the parameter set uses -- Lemma 1 of the
 * paper gives x^n + 1 = (x^(n/2) - r)(x^(n/2) + r) mod q, where r^2 = -1, and
 * both factors are irreducible. So
 *
 *     R_q = Z_q[x]/(x^(n/2) - r)  x  Z_q[x]/(x^(n/2) + r).
 *
 * Multiplying in the decomposition replaces one product of degree-n operands
 * plus a division by the modulus with two products of degree-n/2 operands and
 * no division at all.
 *
 * These routines are drop-in replacements for their nmod_poly counterparts and
 * fall back to them whenever the decomposition does not apply (a different
 * modulus, a different cyclotomic, or no square root of -1), so they are safe
 * to use unconditionally.
 */

/** res = a * b mod (cyclo), i.e. a drop-in nmod_poly_mulmod. */
void evmlr_crt_mulmod(nmod_poly_t res, const nmod_poly_t a, const nmod_poly_t b,
                      const nmod_poly_t cyclo);

/**
 * res = a mod (x^n + 1), for any coefficient modulus.
 *
 * Reducing modulo x^n + 1 needs no division: x^n = -1, so the coefficients
 * above the degree bound just fold back down with alternating signs. This is
 * linear in the degree, where nmod_poly_rem is not, and unlike the CRT routine
 * above it does not care what the coefficient modulus is.
 *
 * Falls back to nmod_poly_rem if cyclo is not of the form x^n + 1.
 */
void evmlr_crt_negacyclic_rem(nmod_poly_t res, const nmod_poly_t a, const nmod_poly_t cyclo);

/** Whether the CRT path is actually in use for the current parameters. */
int evmlr_crt_available(void);

#endif // EVMLR_SHUFFLE_EVMLR_CRT_H
