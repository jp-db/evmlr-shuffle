#ifndef EVMLR_SHUFFLE_EVMLR_BIN_PROOF_H
#define EVMLR_SHUFFLE_EVMLR_BIN_PROOF_H

#include "flint/flint.h"
#include "flint/nmod_poly.h"
#include "flint/nmod_poly_mat.h"
#include "evmlr_params.h"
#include "evmlr_commit.h"

// Context for the binary proof
typedef struct {
    slong k;   // K_SIS
    slong m;   // 2 * K_SIS
    slong L;   // message length (number of polynomials)
    nmod_poly_t cyclo_poly;
} evmlr_bin_proof_ctx_struct;
typedef evmlr_bin_proof_ctx_struct evmlr_bin_proof_ctx_t[1];

typedef struct {
    // commitments for garbage/quadratic helpers
    evmlr_commit_t c_g;
    evmlr_commit_t c_g1;
    
    // masking commitments/values
    nmod_poly_mat_t w;
    nmod_poly_mat_t w_g;
    nmod_poly_mat_t w_g1;
    nmod_poly_t h;
    nmod_poly_t v;
    
    // responses
    nmod_poly_mat_t z_m;
    nmod_poly_mat_t z_r;
    nmod_poly_t z_g;
    nmod_poly_mat_t z_rg;
    nmod_poly_t z_g1;
    nmod_poly_mat_t z_rg1;
} evmlr_bin_proof_struct;
typedef evmlr_bin_proof_struct evmlr_bin_proof_t[1];

void evmlr_bin_proof_ctx_init(evmlr_bin_proof_ctx_t ctx, slong k, slong m, slong L);
void evmlr_bin_proof_ctx_clear(evmlr_bin_proof_ctx_t ctx);

void evmlr_bin_proof_init(evmlr_bin_proof_t proof, const evmlr_bin_proof_ctx_t ctx);
void evmlr_bin_proof_clear(evmlr_bin_proof_t proof);

/**
 * Generate a proof for the relation A_1 * m + A_2 * r = c AND that m has binary coefficients.
 */
int evmlr_bin_prove(evmlr_bin_proof_t proof,
                    const nmod_poly_mat_t A_1,
                    const nmod_poly_mat_t A_2,
                    const nmod_poly_mat_t m,
                    const nmod_poly_mat_t r,
                    const nmod_poly_mat_t c,
                    const evmlr_bin_proof_ctx_t ctx,
                    flint_rand_t state);

/**
 * Verify a proof for the relation A_1 * m + A_2 * r = c AND that m has binary coefficients.
 */
int evmlr_bin_verify(const evmlr_bin_proof_t proof,
                     const nmod_poly_mat_t A_1,
                     const nmod_poly_mat_t A_2,
                     const nmod_poly_mat_t c,
                     const evmlr_bin_proof_ctx_t ctx);

#endif // EVMLR_SHUFFLE_EVMLR_BIN_PROOF_H
