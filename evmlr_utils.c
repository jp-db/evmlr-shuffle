#include "evmlr_utils.h"
#include "sys/random.h"
#include "stdint.h"

int evmlr_utils_is_poly_bin(const nmod_poly_t poly) {
    for (slong i = 0; i < poly->length; i++) {
        ulong coeff = nmod_poly_get_coeff_ui(poly, i);
        if (coeff != 0 && coeff != 1) {
            return 0; // Return false if any coefficient is not binary
        }
    }
    return 1; // All coefficients are binary
}

int evmlr_utils_is_bin(const nmod_poly_t* poly, slong len) {
    for (slong i = 0; i < len; i++) {
        if (!evmlr_utils_is_poly_bin(poly[i])) {
            return 0; // Return false if any polynomial is not binary
        }
    }
    return 1; // All polynomials are binary
}

void evmlr_utils_int_to_bin(nmod_poly_t poly, ulong n) {
    slong degree = 0;
    ulong temp = n;

    // Calculate the degree based on the binary representation of n
    while (temp > 0) {
        temp >>= 1;
        degree++;
    }

    // Initialize the polynomial with the calculated degree
    nmod_poly_fit_length(poly, degree);

    // Fill the polynomial coefficients with the binary representation
    for (slong i = 0; i < degree; i++) {
        nmod_poly_set_coeff_ui(poly, i, n & 1); // Get the least significant bit
        n >>= 1; // Shift right to process the next bit
    }
}

int evmlr_utils_binom_sample(int center) {
    int sample = 0;

    // Sample from a centered binomial distribution
    for (int i = 0; i < center; i++) {
        // Generate a random bit (0 or 1)
        u_int8_t buff[1];
        getrandom(buff, sizeof(u_int8_t), 0);
        // get first and second bits from buff
        int a = (buff[0] & 1); // Random bit a
        int b = (buff[0] >> 1) & 1; // Random bit b

        sample += a - b; // Either 1, 0, or -1
    }

    return sample;
}

void evmlr_utils_binom_sample_ring(nmod_poly_t poly, int center) {
    nmod_poly_fit_length(poly, DEGREE_N);
    for (slong i = 0; i < DEGREE_N; i++) {
        int sample = evmlr_utils_binom_sample(center);
        // Ensure the coefficient is non-negative modulo MOD_Q
        ulong coeff = (sample + MOD_Q) % MOD_Q;
        nmod_poly_set_coeff_ui(poly, i, coeff);
    }
}

void evmlr_new_perm(ulong* perm, size_t len) {
    // Initialize the permutation array with the identity permutation
    for (size_t i = 0; i < len; i++) {
        perm[i] = i;
    }

    // Shuffle the array using Fisher-Yates algorithm
    for (size_t i = len - 1; i > 0; i--) {
        u_int8_t buff[4];
        getrandom(buff, sizeof(u_int8_t) * 4, 0);
        size_t j = *((uint32_t*)buff) % (i + 1); // Random index from 0 to i

        // Swap perm[i] and perm[j]
        ulong temp = perm[i];
        perm[i] = perm[j];
        perm[j] = temp;
    }
}
