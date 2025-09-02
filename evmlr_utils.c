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
    if (center == 0) return 0;

    int sample = 0;
    // Calculate how many bytes we need. Each byte provides 4 samples (8 bits / 2 bits per sample).
    ulong bytes_needed = (center + 3) / 4;
    uint8_t* buff = (uint8_t*) malloc(bytes_needed);

    // Single system call to get all the random data we need.
    getrandom(buff, bytes_needed, 0);

    int i = 0;
    for (int j = 0; j < bytes_needed; ++j) {
        uint8_t random_byte = buff[j];

        // Process 4 samples from a single byte
        for (int k = 0; k < 4 && i < center; ++k, ++i) {
            int a = (random_byte >> (k * 2)) & 1;     // Bit a
            int b = (random_byte >> (k * 2 + 1)) & 1; // Bit b
            sample += a - b;
        }
    }

    free(buff);
    return sample;
}

void evmlr_utils_binom_sample_ring(nmod_poly_t poly, int center) {
    if (center == 0) {
        // Handle zero-center case (results in a zero polynomial)
        nmod_poly_zero(poly);
        return;
    }

    nmod_poly_fit_length(poly, DEGREE_N);

    // 1. Calculate total random data needed and fetch it in one call.
    // Each sample needs 'center' pairs of bits (2 * center bits).
    // Each 64-bit word provides 32 pairs of bits.
    size_t total_bit_pairs = (size_t)DEGREE_N * center;
    size_t words_needed = (total_bit_pairs + 31) / 32;
    uint64_t* buff = (uint64_t*) malloc(words_needed * sizeof(uint64_t));

    getrandom(buff, words_needed * sizeof(uint64_t), 0);

    size_t word_idx = 0;
    int bit_pair_idx = 0;

    for (slong i = 0; i < DEGREE_N; i++) {
        int sample = 0;

        for (int j = 0; j < center; j++) {
            uint64_t current_word = buff[word_idx];

            uint8_t a = (current_word >> (bit_pair_idx * 2)) & 1;
            uint8_t b = (current_word >> (bit_pair_idx * 2 + 1)) & 1;
            sample += a - b;

            // Advance our position in the buffer.
            bit_pair_idx++;
            if (bit_pair_idx == 32) {
                bit_pair_idx = 0;
                word_idx++;
            }
        }

        ulong coeff = (sample + MOD_Q) % MOD_Q;
        nmod_poly_set_coeff_ui(poly, i, coeff);
    }

    free(buff);
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

void evmlr_utils_ring_to_bin(size_t len, int bits, nmod_poly_t bin_vec[bits][len], const fmpz_poly_t ring_vec[len], ulong mod) {
    fmpz_t coeff, two_pow_b_minus_1, two_pow_b;
    fmpz_init(coeff);
    fmpz_init(two_pow_b_minus_1);
    fmpz_init(two_pow_b);

    fmpz_set_ui(two_pow_b, 1);
    fmpz_mul_2exp(two_pow_b, two_pow_b, bits); // 2^b

    fmpz_set_ui(two_pow_b_minus_1, 1);
    fmpz_mul_2exp(two_pow_b_minus_1, two_pow_b_minus_1, bits - 1); // 2^(b-1)

    for (int b = 0; b < bits; b++) {
        for (size_t i = 0; i < len; i++) {
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

void evmlr_utils_stack(size_t len, int bits, nmod_poly_t stack[bits * len], const nmod_poly_t bin_vec[bits][len], ulong mod) {
    for (int b = 0; b < bits; b++) {
        for (size_t i = 0; i < len; i++) {
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
