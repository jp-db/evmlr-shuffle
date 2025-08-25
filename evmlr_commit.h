#ifndef EVMLR_SHUFFLE_EVMLR_COMMIT_H
#define EVMLR_SHUFFLE_EVMLR_COMMIT_H
#include "evmlr_params.h"
#define L 1

// B := \sqrt{nL + n 2 K_SIS ETA^2}
// public_parameters pp := (K_SIS, L, ETA, B, A_1, A_2)
typedef struct {
    double b;
    // A_1 \sample R_q^{K_SIS x L}, A_2 \sample R_q^{K_SIS x 2K_SIS}
    nmod_poly_t A_1[K_SIS][L];
    nmod_poly_t A_2[K_SIS][2 * K_SIS];
} evmlr_commit_scheme_struct;
typedef evmlr_commit_scheme_struct evmlr_commit_scheme_t[1];

typedef struct {
 nmod_poly_t c[K_SIS];
} evmlr_commit_struct;
typedef evmlr_commit_struct evmlr_commit_t[1];

void evmlr_commit_scheme_init(evmlr_commit_scheme_t pp, flint_rand_t state);

void evmlr_commit_scheme_clear(evmlr_commit_scheme_t pp);

// sample from binom r \sample B^{2K_SIS}_eta
void evmlr_sample_r(nmod_poly_t r[2*K_SIS]);

// c = A_1 * msg + A_2 * r
void evmlr_commit(evmlr_commit_t com, const nmod_poly_t msg[L], const nmod_poly_t r[2*K_SIS], const evmlr_commit_scheme_t pp);

void evmlr_commit_clear(evmlr_commit_t com, const evmlr_commit_scheme_t pp);

// Verify that ||\binom{m}{r}|| <= B ^ c = A_1 * msg + A_2 * r
int evmlr_commit_verify(const evmlr_commit_t com, const nmod_poly_t msg[L], const nmod_poly_t r[2*K_SIS], const evmlr_commit_scheme_t pp);

#endif