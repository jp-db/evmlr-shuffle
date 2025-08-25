#ifndef EVMLR_SHUFFLE_EVMLR_COMMIT_H
#define EVMLR_SHUFFLE_EVMLR_COMMIT_H
#include "evmlr_params.h"

typedef struct {
    nmod_poly_t s[K_LWE];
} evmlr_enc_sk_struct;
typedef evmlr_enc_sk_struct evmlr_enc_sk_t[1];

typedef struct {
    nmod_poly_t A[K_LWE][K_LWE];
    // t = A*s + e
    nmod_poly_t t[K_LWE];
} evmlr_enc_pk_struct;
typedef evmlr_enc_pk_struct evmlr_enc_pk_t[1];

typedef struct {
    evmlr_enc_sk_t sk;
    evmlr_enc_pk_t pk;
} evmlr_enc_keypair_struct;
typedef evmlr_enc_keypair_struct evmlr_enc_keypair_t[1];

typedef struct {
    // r, e_2 \sample B^{K_LWE}_\eta
    // u^T = r^T A + e_2^T
    nmod_poly_t uT[K_LWE];
    // e_3 \sample B_\eta
    // v = r^T t + e_3 + \lfloor q/2 \rfloor * msg
    nmod_poly_t v;
} evmlr_enc_cipher_struct;
typedef evmlr_enc_cipher_struct evmlr_enc_cipher_t[1];

void evmlr_enc_keypair_gen(evmlr_enc_keypair_t keypair, flint_rand_t state);

void evmlr_enc_clear_keypair(evmlr_enc_keypair_t keypair);

void evmlr_enc(evmlr_enc_cipher_t cipher, const nmod_poly_t msg, const evmlr_enc_pk_t pk
);

void evmlr_dec(nmod_poly_t msg, const evmlr_enc_cipher_t cipher, const evmlr_enc_sk_t sk);

void evmlr_enc_clear_cipher(evmlr_enc_cipher_t cipher);

#endif