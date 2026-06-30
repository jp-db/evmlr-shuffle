#ifndef EVMLR_SHUFFLE_EVMLR_LIN_PROOF_H
#define EVMLR_SHUFFLE_EVMLR_LIN_PROOF_H

#include "flint/flint.h"
#include "flint/nmod_poly.h"
#include "flint/nmod_poly_mat.h"
#include "evmlr_params.h"

// Context for the linear proof: A s_1 + s_2 = t
struct evmlr_lin_proof_ctx_struct;

typedef void (*mask_sampler_t)(nmod_poly_mat_t y, const struct evmlr_lin_proof_ctx_struct* ctx);
typedef int (*rejection_check_t)(const nmod_poly_mat_t z_1, const nmod_poly_mat_t z_2, const struct evmlr_lin_proof_ctx_struct* ctx);
typedef int (*norm_check_t)(const nmod_poly_mat_t z_1, const nmod_poly_mat_t z_2, const struct evmlr_lin_proof_ctx_struct* ctx);

typedef struct evmlr_lin_proof_ctx_struct {
    slong k;
    slong m;
    nmod_poly_t cyclo_poly;
    slong eta;
    slong beta;
    mask_sampler_t sample_mask;
    rejection_check_t rejection_check;
    norm_check_t norm_check;
} evmlr_lin_proof_ctx_struct;

typedef evmlr_lin_proof_ctx_struct evmlr_lin_proof_ctx_t[1];

typedef struct {
    nmod_poly_mat_t w;   // k x 1
    nmod_poly_mat_t z_1; // m x 1
    nmod_poly_mat_t z_2; // k x 1
} evmlr_lin_proof_struct;

typedef evmlr_lin_proof_struct evmlr_lin_proof_t[1];

void evmlr_lin_proof_ctx_init(evmlr_lin_proof_ctx_t ctx, slong k, slong m);
void evmlr_lin_proof_ctx_clear(evmlr_lin_proof_ctx_t ctx);

void evmlr_lin_proof_init(evmlr_lin_proof_t proof, const evmlr_lin_proof_ctx_t ctx);
void evmlr_lin_proof_clear(evmlr_lin_proof_t proof);

/**
 * Generate a proof for the relation A * s_1 + s_2 = t.
 * Note: Uses binomial sampling for masks y_1, y_2 with parameter ETA.
 * 
 * @param proof Output proof structure containing responses and challenge
 * @param A Public matrix of size k x m
 * @param s_1 Secret vector of size m x 1
 * @param s_2 Secret vector of size k x 1
 * @param t Public vector of size k x 1
 * @param ctx Proof context
 * @param state FLINT RNG state for challenge generation
 * @return 1 on success
 */
int evmlr_lin_prove(evmlr_lin_proof_t proof, 
                    const nmod_poly_mat_t A, 
                    const nmod_poly_mat_t s_1, 
                    const nmod_poly_mat_t s_2, 
                    const nmod_poly_mat_t t,
                    const evmlr_lin_proof_ctx_t ctx);

/**
 * Verify a proof for the relation A * s_1 + s_2 = t.
 * 
 * @param proof Proof to verify
 * @param A Public matrix of size k x m
 * @param t Public vector of size k x 1
 * @param ctx Proof context
 * @return 1 if proof is valid, 0 otherwise
 */
int evmlr_lin_verify(const evmlr_lin_proof_t proof, 
                     const nmod_poly_mat_t A, 
                     const nmod_poly_mat_t t, 
                     const evmlr_lin_proof_ctx_t ctx);

void evmlr_lin_proof_ctx_set_linear(evmlr_lin_proof_ctx_t ctx, slong beta);

#endif // EVMLR_SHUFFLE_EVMLR_LIN_PROOF_H
