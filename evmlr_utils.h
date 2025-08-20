#ifndef EVMLR_SHUFFLE_EVMLR_UTILS_H
#define EVMLR_SHUFFLE_EVMLR_UTILS_H
#include "evmlr_params.h"

// 2 Preliminaries
// 2.2 Polynomial Ring

// R_q^* -> {True, False}. Takes a vector over Rq and evaluates to True
// if all the coefficients of all its elements are in {0, 1}, and
// to False otherwise.
/**
 * Check if a all coefficients from all polynomials in the vector are in {0, 1}.
 * @param poly Pointer to an array of nmod_poly_t polynomials.
 * @param len Length of the array.
 * @return 1 if all coefficients are in {0, 1}, 0 otherwise.
 */
int evmlr_utils_is_bin(nmod_poly_t* poly, slong len);

// {0, ..., 2^{n} - 1} -> R_q. Maps an integer in the domain to the ring element
/**
 * Maps an integer to a polynomial whose coefficients are the binary representation of the integer.
 * The least significant bit is in the constant coefficient.
 * @param poly The polynomial to be filled with the binary representation.
 * @param n the integer to be converted.
 */
void evmlr_utils_int_to_bin(nmod_poly_t poly, ulong n);

// R_q^N \times \N -> (R_q^N)^*. Binary decomposition of (a vector of) ring elements
// TODO
void evmlr_utils_ring_to_bin(nmod_poly_t* poly, ulong n, slong len);

// TODO
void evmlr_utils_stack();

// TODO
void evmlr_utils_gadget_matrix();


// B_\eta centered binomial distribution for some positive integer \eta
/**
 * Sample from a centered binomial distribution with a given center.
 * @param center  The center of the binomial distribution, which is also the number of trials.
 * @return An integer sampled from the centered binomial distribution.
 */
int evmlr_utils_binom_sample(int center);

#endif // EVMLR_SHUFFLE_EVMLR_UTILS_H