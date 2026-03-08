#ifndef EVMLR_SHUFFLE_EVMLR_MLPKE_H
#define EVMLR_SHUFFLE_EVMLR_MLPKE_H
#include "evmlr_params.h"
#include "flint/nmod_poly_mat.h"
#include "flint/flint.h"
#include "evmlr_crt.h"

typedef struct {
    nmod_poly_t cyclo_poly; // x^N + 1
    evmlr_crt_ctx_t crt_ctx;
} evmlr_mlpke_ctx_struct;
typedef evmlr_mlpke_ctx_struct evmlr_mlpke_ctx_t[1];

typedef struct {
    nmod_poly_mat_t s; // secret key s \in R_q^{K_LWE}
} evmlr_mlpke_sk_struct;
typedef evmlr_mlpke_sk_struct evmlr_mlpke_sk_t[1];

typedef struct {
    nmod_poly_mat_t A; // public matrix A \in R_q^{K_LWE x K_LWE}
    // t = A*s + e
    nmod_poly_mat_t t; // public vector t \in R_q^{K_LWE}
} evmlr_mlpke_pk_struct;
typedef evmlr_mlpke_pk_struct evmlr_mlpke_pk_t[1];

typedef struct {
    evmlr_mlpke_sk_t sk;
    evmlr_mlpke_pk_t pk;
} evmlr_mlpke_keypair_struct;
typedef evmlr_mlpke_keypair_struct evmlr_mlpke_keypair_t[1];

typedef struct {
    nmod_poly_mat_t uT; // vector of K_LWE polynomials
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