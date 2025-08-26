#include "evmlr_utils.h"
#include "sys/random.h"
#include "stdint.h"
#include "flint/fmpz.h"

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

void evmlr_utils_ring_to_bin(slong len, int bits, nmod_poly_t bin_vec[bits][len], const fmpz_poly_t ring_vec[len], ulong mod) {
    fmpz_t coeff, two_pow_b_minus_1, two_pow_b;
    fmpz_init(coeff);
    fmpz_init(two_pow_b_minus_1);
    fmpz_init(two_pow_b);

    fmpz_set_ui(two_pow_b, 1);
    fmpz_mul_2exp(two_pow_b, two_pow_b, bits); // 2^b

    fmpz_set_ui(two_pow_b_minus_1, 1);
    fmpz_mul_2exp(two_pow_b_minus_1, two_pow_b_minus_1, bits - 1); // 2^(b-1)

    for (int b = 0; b < bits; b++) {
        for (slong i = 0; i < len; i++) {
            nmod_poly_init(bin_vec[b][i], mod);
            nmod_poly_zero(bin_vec[b][i]);
        }
    }

    for (slong i = 0; i < len; i++) {
        slong deg = fmpz_poly_degree(ring_vec[i]);
        for (int j = 0; j < bits; j++) {
            nmod_poly_fit_length(bin_vec[j][i], deg + 1);
        }
        for (slong k = 0; k < deg; k++) {
            fmpz_poly_get_coeff_fmpz(coeff, ring_vec[i], k);
            // The range is [-2^(b-1), 2^(b-1)). If coeff is negative,
            // its two's complement representation is 2^b + coeff.
            if (fmpz_sgn(coeff) < 0) {
                fmpz_add(coeff, coeff, two_pow_b); // Make coeff non-negative
            }
            for (int j = 0; j < bits; j++) {
                if (fmpz_tstbit(coeff, j)) {
                    nmod_poly_set_coeff_ui(bin_vec[j][i], k, 1);
                }
            }
        }
    }

    fmpz_clear(coeff);
    fmpz_clear(two_pow_b_minus_1);
    fmpz_clear(two_pow_b);
}

void evmlr_utils_stack(slong len, int bits, nmod_poly_t stack[bits * len], const nmod_poly_t bin_vec[bits][len], ulong mod) {
    for (int b = 0; b < bits; b++) {
        for (slong i = 0; i < len; i++) {
            nmod_poly_init(stack[b * len + i], mod);
            nmod_poly_set(stack[b * len + i], bin_vec[b][i]);
        }
    }
}

void evmlr_utils_sample_binary_poly(nmod_poly_t poly, slong degree, flint_rand_t state) {
    nmod_poly_fit_length(poly, degree);
    for (slong i = 0; i < degree; i++) {
        ulong rand_bit = n_randint(state, 2); // Random bit 0 or 1
        nmod_poly_set_coeff_ui(poly, i, rand_bit);
    }
}
