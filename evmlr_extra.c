#include "evmlr_params.h"
#include "flint/fmpz.h"
#include "flint/fmpz_mod.h"
#include "flint/fmpz_mod_poly.h"

// If we have lists of elements (a_1, ..., a_N ) and (b_1, ..., b_N ), where all the a_i, b_i
// are in some field F, then the two lists are a permutation of each other if and only
// if over F[X], we have:
// \perm{i=1}^{N} (a_i - X) = \perm{i=1}^{N} (b_i - X)

int is_perm_over_f(const fmpz_t* a, const fmpz_t* b, size_t len, fmpz_mod_ctx_t ctx) {
    if (len == 0) {
        return 1; // Empty lists are permutations of each other
    }

    fmpz_mod_poly_t prod_a, prod_b; // Polynomials to hold the products
    fmpz_mod_poly_t tmp;            // Temporary polynomial
    fmpz_mod_poly_t X_poly;         // Polynomial representing X
    fmpz_t coeff;                   // Temporary coefficient

    fmpz_mod_poly_init(prod_a, ctx);
    fmpz_mod_poly_init(prod_b, ctx);
    fmpz_mod_poly_init(tmp, ctx);
    fmpz_mod_poly_init(X_poly, ctx);
    fmpz_init(coeff);

    // X_poly = 1 * X^1 + 0 * X^0 = X
    fmpz_mod_poly_set_coeff_ui(X_poly, 1, 1, ctx); // Coefficient of X^1 is 1

    // \perm{i=1}^{N} (a_i - X)
    fmpz_mod_poly_one(prod_a, ctx);
    for (size_t i = 0; i < len; i++) {
        // tmp = a_i - X
        fmpz_mod_poly_set_fmpz(tmp, a[i], ctx);
        fmpz_mod_poly_sub(tmp, tmp, X_poly, ctx);
        // prod_a *= tmp
        fmpz_mod_poly_mul(prod_a, prod_a, tmp, ctx);
    }

    // \perm{i=1}^{N} (b_i - X)
    fmpz_mod_poly_one(prod_b, ctx);
    for (size_t i = 0; i < len; i++) {
        // tmp = b_i - X
        fmpz_mod_poly_set_fmpz(tmp, b[i], ctx);
        fmpz_mod_poly_sub(tmp, tmp, X_poly, ctx);
        // prod_b *= tmp
        fmpz_mod_poly_mul(prod_b, prod_b, tmp, ctx);
    }

    // Check if prod_a == prod_b
}