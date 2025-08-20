#include "evmlr_utils.h"
#include "sys/random.h"

int _is_poly_bin(nmod_poly_t poly) {
    for (slong i = 0; i < poly->length; i++) {
        if (poly->coeffs[i] != 0 && poly->coeffs[i] != 1) {
            return 0; // Return false if any coefficient is not binary
        }
    }
    return 1; // All coefficients are binary
}

int evmlr_utils_is_bin(nmod_poly_t* poly, slong len) {
    for (slong i = 0; i < len; i++) {
        if (!_is_poly_bin(poly[i])) {
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
        poly->coeffs[i] = (n & 1); // Get the least significant bit
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