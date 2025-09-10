#ifndef EVMLR_SHUFFLE_EVMLR_OTSE_H
#define EVMLR_SHUFFLE_EVMLR_OTSE_H
#include "evmlr_params.h"

typedef struct {
    slong L;
    nmod_poly_mat_t H; // H[K_LWE + M_LEN][K_LWR]
    nmod_poly_mat_t H_prime; // [M_LEN][K_LWE + M_LEN]
    nmod_poly_t cyclo_poly; // x^n + 1
} evmlr_otse_ctx_struct;
typedef evmlr_otse_ctx_struct evmlr_otse_ctx_t[1];

typedef struct {
    nmod_poly_mat_t s; // s \sample {0,1}^{K_LWR}
} evmlr_otse_key_struct;
typedef evmlr_otse_key_struct evmlr_otse_key_t[1];

typedef struct {
    nmod_poly_mat_t c; // ciphertext c \in R_q^{L}
} evmlr_otse_ciphertext_struct;
typedef evmlr_otse_ciphertext_struct evmlr_otse_ciphertext_t[1];

void evmlr_otse_ctx_init(evmlr_otse_ctx_t ctx, slong L, flint_rand_t state);

void evmlr_otse_ctx_clear(evmlr_otse_ctx_t ctx);

// s \sample {0,1}^{K_LWR}
void evmlr_otse_keygen(evmlr_otse_key_t key, flint_rand_t state);

void evmlr_otse_keyclear(evmlr_otse_key_t key);

// a = H' * B^(n) d_dagger
void evmlr_calc_a(nmod_poly_mat_t a, const nmod_poly_mat_t d_stack, const evmlr_otse_ctx_t ctx);

void evmlr_otse_encrypt(evmlr_otse_ciphertext_t ct, nmod_poly_mat_t d_dagger, const nmod_poly_mat_t m,
                        const evmlr_otse_key_t key, const evmlr_otse_ctx_t ctx);

void evmlr_otse_decrypt(nmod_poly_mat_t m, const evmlr_otse_ciphertext_t ct, const evmlr_otse_key_t key, const evmlr_otse_ctx_t ctx);

void evmlr_otse_ciphertext_clear(evmlr_otse_ciphertext_t ct);

#endif //EVMLR_SHUFFLE_EVMLR_OTSE_H
