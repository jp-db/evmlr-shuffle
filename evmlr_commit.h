#ifndef EVMLR_SHUFFLE_EVMLR_COMMIT_H
#define EVMLR_SHUFFLE_EVMLR_COMMIT_H
#include "evmlr_params.h"

// B := \sqrt{nL + n 2 K_SIS ETA^2}
// public_parameters pp := (K_SIS, M_LEN, ETA, B, A_1, A_2)
typedef struct {
    ulong b_sqr; // B^2
    slong N; // message length
    // A_1 \sample R_q^{K_SIS x N}, A_2 \sample R_q^{K_SIS x 2K_SIS}
    nmod_poly_mat_t A_1;
    nmod_poly_mat_t A_2; // K_SIS x 2K_SIS matrix
    nmod_poly_t cyclo_poly; // x^N + 1
} evmlr_commit_ctx_struct;
typedef evmlr_commit_ctx_struct evmlr_commit_ctx_t[1];

typedef struct {
    nmod_poly_mat_t c; // vector of K_SIS polynomials
} evmlr_commit_struct;
typedef evmlr_commit_struct evmlr_commit_t[1];

void evmlr_commit_ctx_init(evmlr_commit_ctx_t ctx, slong N, flint_rand_t state);

void evmlr_commit_ctx_clear(evmlr_commit_ctx_t ctx);

// sample from binom r \sample B^{2K_SIS}_eta
void evmlr_commit_sample_r(nmod_poly_mat_t r);

// c = A_1 * msg + A_2 * r
void evmlr_commit(evmlr_commit_t com, const nmod_poly_mat_t msg, const nmod_poly_mat_t r, const evmlr_commit_ctx_t ctx);

void evmlr_commit_clear(evmlr_commit_t com);

// Verify that ||\binom{m}{r}|| <= B ^ c = A_1 * msg + A_2 * r
int evmlr_commit_verify(const evmlr_commit_t com, const nmod_poly_mat_t msg, const nmod_poly_mat_t r, const evmlr_commit_ctx_t ctx);

#endif