#ifndef EVMLR_SHUFFLE_EVMLR_SHUFFLE_H
#define EVMLR_SHUFFLE_EVMLR_SHUFFLE_H

#include "evmlr_params.h"
#include "evmlr_commit.h"
#include "evmlr_hpke.h"

typedef struct {
    evmlr_commit_ctx_t com_ctx;
    evmlr_hpke_ctx_t hpke_ctx;
    size_t N; // number of elements
    size_t L; // vector length
} evmlr_shuffle_ctx_struct;
typedef evmlr_shuffle_ctx_struct evmlr_shuffle_ctx_t[1];

typedef struct {
    size_t* pi; // permutation array (size N)
    nmod_poly_t** d_dagger; // (2*ETA*(K_LWE + M_LEN)) x N matrix
    nmod_poly_t r_D[2 * K_SIS];
} evmlr_shuffle_sp_struct; // Secret parameters
typedef evmlr_shuffle_sp_struct evmlr_shuffle_sp_t[1];

typedef struct {
    evmlr_commit_t D; // commitment to d_dagger
    nmod_poly_t** c_hat; // shuffled ciphertexts c_hat_i = c_star_pi[i] - H'B(η)d_dagger_i
    evmlr_otse_ciphertext_t* c_star; // symmetric ciphertexts
} evmlr_shuffle_pp_struct; // Public parameters
typedef evmlr_shuffle_pp_struct evmlr_shuffle_pp_t[1];

void evmlr_shuffle_ctx_init(evmlr_shuffle_ctx_t ctx, size_t N,size_t L, flint_rand_t state);

void evmlr_shuffle_ctx_clear(evmlr_shuffle_ctx_t ctx);

void evmlr_shuffle_sample_challenge(nmod_poly_t chal, flint_rand_t state);

#endif //EVMLR_SHUFFLE_EVMLR_SHUFFLE_H
