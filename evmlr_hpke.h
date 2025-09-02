#ifndef EVMLR_SHUFFLE_EVMLR_HPKE_H
#define EVMLR_SHUFFLE_EVMLR_HPKE_H

#include "evmlr_mlpke.h"
#include "evmlr_otse.h"

typedef struct {
    evmlr_mlpke_ctx_t enc_ctx;
    evmlr_otse_ctx_t otse_ctx;
} evmlr_hpke_ctx_struct;
typedef evmlr_hpke_ctx_struct evmlr_hpke_ctx_t[1];

typedef struct {
    evmlr_mlpke_keypair_t enc_keypair;
} evmlr_hpke_keypair_struct;
typedef evmlr_hpke_keypair_struct evmlr_hpke_keypair_t[1];

typedef struct {
    evmlr_mlpke_cipher_t enc_cipher[K_LWR];
    evmlr_otse_ciphertext_t otse_cipher;
} evmlr_hpke_cipher_struct;
typedef evmlr_hpke_cipher_struct evmlr_hpke_cipher_t[1];

void evmlr_hpke_ctx_init(evmlr_hpke_ctx_t ctx, size_t L, flint_rand_t state);

void evmlr_hpke_ctx_clear(evmlr_hpke_ctx_t ctx);

void evmlr_hpke_keypair_gen(evmlr_hpke_keypair_t keypair, flint_rand_t state, const evmlr_hpke_ctx_t ctx);

void evmlr_hpke_keypair_clear(evmlr_hpke_keypair_t keypair);

void evmlr_hpke_encrypt(evmlr_hpke_cipher_t cipher, const nmod_poly_t msg[], const evmlr_mlpke_pk_t pk, const evmlr_hpke_ctx_t ctx, flint_rand_t state);

void evmlr_hpke_decrypt(nmod_poly_t msg[], const evmlr_hpke_cipher_t cipher, const evmlr_mlpke_sk_t sk, const evmlr_hpke_ctx_t ctx);

void evmlr_hpke_cipher_clear(evmlr_hpke_cipher_t cipher, const evmlr_hpke_ctx_t ctx);

#endif //EVMLR_SHUFFLE_EVMLR_HPKE_H
