#ifndef EVMLR_SHUFFLE_EVMLR_SHUFFLE_H
#define EVMLR_SHUFFLE_EVMLR_SHUFFLE_H

#include "evmlr_params.h"
#include "evmlr_commit.h"
#include "evmlr_hpke.h"

typedef struct {
    evmlr_commit_ctx_t com_ctx;
    evmlr_hpke_ctx_t hpke_ctx;
} evmlr_shuffle_ctx_struct;

typedef struct {
    ulong N; // number of elements
    ulong* pi; // permutation array (size N)
    nmod_poly_t* d_dagger[(2*ETA) * (K_LWE + M_LEN)]; // (2*ETA*(K_LWE + M_LEN)) x N matrix
    nmod_poly_t* r_D[2 * K_SIS]; // 2K_SIS x N matrix
} evmlr_shuffle_sp_struct; // Secret parameters
typedef evmlr_shuffle_sp_struct evmlr_shuffle_sp_t[1];

typedef evmlr_shuffle_ctx_struct evmlr_shuffle_ctx_t[1];
typedef struct {
    ulong N; // number of elements
    evmlr_commit_t* D; // commitments to messages (size N)
    evmlr_hpke_cipher_t* C; // encryptions of messages (size N)
} evmlr_shuffle_pp_struct; // Public parameters

typedef evmlr_shuffle_pp_struct evmlr_shuffle_pp_t[1];

void evmlr_shuffle_ctx_init(evmlr_shuffle_ctx_t ctx, ulong L, flint_rand_t state);

void evmlr_shuffle_ctx_clear(evmlr_shuffle_ctx_t ctx);

void evmlr_shuffle_sample_challenge(nmod_poly_t chal, flint_rand_t state);

#endif //EVMLR_SHUFFLE_EVMLR_SHUFFLE_H
