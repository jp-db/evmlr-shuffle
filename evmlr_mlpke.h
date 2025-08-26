#ifndef EVMLR_SHUFFLE_EVMLR_MLPKE_H
#define EVMLR_SHUFFLE_EVMLR_MLPKE_H
#include "evmlr_params.h"

typedef struct {
    nmod_poly_t cyclo_poly; // x^N + 1
} evmlr_mlpke_ctx_struct;
typedef evmlr_mlpke_ctx_struct evmlr_mlpke_ctx_t[1];

typedef struct {
    nmod_poly_t s[K_LWE];
} evmlr_mlpke_sk_struct;
typedef evmlr_mlpke_sk_struct evmlr_mlpke_sk_t[1];

typedef struct {
    nmod_poly_t A[K_LWE][K_LWE];
    // t = A*s + e
    nmod_poly_t t[K_LWE];
} evmlr_mlpke_pk_struct;
typedef evmlr_mlpke_pk_struct evmlr_mlpke_pk_t[1];

typedef struct {
    evmlr_mlpke_sk_t sk;
    evmlr_mlpke_pk_t pk;
} evmlr_mlpke_keypair_struct;
typedef evmlr_mlpke_keypair_struct evmlr_mlpke_keypair_t[1];

typedef struct {
    // r, e_2 \sample B^{K_LWE}_\eta
    // u^T = r^T A + e_2^T
    nmod_poly_t uT[K_LWE];
    // e_3 \sample B_\eta
    // v = r^T t + e_3 + \lfloor q/2 \rfloor * msg
    nmod_poly_t v;
} evmlr_mlpke_cipher_struct;
typedef evmlr_mlpke_cipher_struct evmlr_mlpke_cipher_t[1];

void evmlr_mlpke_ctx_init(evmlr_mlpke_ctx_t ctx);

void evmlr_mlpke_ctx_clear(evmlr_mlpke_ctx_t ctx);

void evmlr_mlpke_keypair_gen(evmlr_mlpke_keypair_t keypair, flint_rand_t state, const evmlr_mlpke_ctx_t ctx);

void evmlr_mlpke_keypair_clear(evmlr_mlpke_keypair_t keypair);

void evmlr_mlpke_enc(evmlr_mlpke_cipher_t cipher, const nmod_poly_t msg, const evmlr_mlpke_pk_t pk, const evmlr_mlpke_ctx_t ctx);

void evmlr_mlpke_dec(nmod_poly_t msg, const evmlr_mlpke_cipher_t cipher, const evmlr_mlpke_sk_t sk, const evmlr_mlpke_ctx_t ctx);

void evmlr_mlpke_cipher_clear(evmlr_mlpke_cipher_t cipher);

#endif