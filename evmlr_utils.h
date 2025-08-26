#ifndef EVMLR_SHUFFLE_EVMLR_UTILS_H
#define EVMLR_SHUFFLE_EVMLR_UTILS_H
#include "evmlr_params.h"
#include "flint/fmpz_poly.h" // for signed coefficients

// 2 Preliminaries
// 2.2 Polynomial Ring

int evmlr_utils_is_poly_bin(const nmod_poly_t poly);
// R_q^* -> {True, False}. Takes a vector over Rq and evaluates to True
// if all the coefficients of all its elements are in {0, 1}, and
// to False otherwise.
/**
 * Check if a all coefficients from all polynomials in the vector are in {0, 1}.
 * @param poly Pointer to an array of nmod_poly_t polynomials.
 * @param len Length of the array.
 * @return 1 if all coefficients are in {0, 1}, 0 otherwise.
 */
int evmlr_utils_is_bin(const nmod_poly_t* poly, slong len);

// {0, ..., 2^{n} - 1} -> R_q. Maps an integer in the domain to the ring element
/**
 * Maps an integer to a polynomial whose coefficients are the binary representation of the integer.
 * The least significant bit is in the constant coefficient.
 * @param poly The polynomial to be filled with the binary representation.
 * @param n the integer to be converted.
 */
void evmlr_utils_int_to_bin(nmod_poly_t poly, ulong n);

// R_q^N \times \N -> (R_q^N)^*. Bina7ry decomposition of (a vector of) ring elements
/**
 * Decomposes a vector of polynomials into b binary polynomial vectors.
 * @param len Length of the input array.
 * @param bits Number of bits for the binary decomposition.
 * @param bin_vec Output array of (nmod_poly_t *) vectors, size b.
 * @param ring_vec Input array of polynomials to be decomposed.
 * @param mod Modulus for the polynomials.
 */
void evmlr_utils_ring_to_bin(slong len, int bits, nmod_poly_t bin_vec[bits][len], const fmpz_poly_t ring_vec[len], ulong mod);

/**
 * Stacks the binary vectors into a single vector.
 * @param len length of each binary vector
 * @param bits number of binary vectors
 * @param stack output stacked vector of length bits * len
 * @param bin_vec input binary vectors of size bits, each of length len
 * @param mod modulus for the polynomials
 */
void evmlr_utils_stack(slong len, int bits, nmod_poly_t stack[bits * len], const nmod_poly_t bin_vec[bits][len], ulong mod);

// TODO
void evmlr_utils_gadget_matrix();


// B_\eta centered binomial distribution for some positive integer \eta
/**
 * Sample from a centered binomial distribution with a given center.
 * @param center  The center of the binomial distribution, which is also the number of trials.
 * @return An integer sampled from the centered binomial distribution.
 */
int evmlr_utils_binom_sample(int center);

// we are sampling an element of Rq by sampling each of its coefficients independently according to B^*_η
void evmlr_utils_binom_sample_ring(nmod_poly_t poly, int center);

void evmlr_utils_sample_binary_poly(nmod_poly_t poly, slong degree, flint_rand_t state);

/**
 * Generates a new random permutation of length `len` and stores it in the `perm` array.
 * @param perm Pointer to an array where the permutation will be stored.
 * @param len Length of the permutation array.
 */
void evmlr_new_perm(ulong* perm, size_t len);


#endif // EVMLR_SHUFFLE_EVMLR_UTILS_H