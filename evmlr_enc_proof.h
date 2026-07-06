#ifndef EVMLR_SHUFFLE_EVMLR_ENC_PROOF_H
#define EVMLR_SHUFFLE_EVMLR_ENC_PROOF_H

#include "evmlr_lin_proof.h"
#include "evmlr_mlpke.h"
#include "evmlr_otse.h"
#include "evmlr_hpke.h"

typedef struct {
    evmlr_lin_proof_struct proof_hpke;
    evmlr_lin_proof_struct proof_otse;
} evmlr_enc_proof_struct;
typedef evmlr_enc_proof_struct evmlr_enc_proof_t[1];

void evmlr_enc_proof_init(evmlr_enc_proof_t proof, slong L_curr);
void evmlr_enc_proof_clear(evmlr_enc_proof_t proof);

void evmlr_enc_proof_prove(evmlr_enc_proof_t proof,
                           const evmlr_mlpke_pk_t pk,
                           const evmlr_otse_ctx_t otse_ctx,
                           const nmod_poly_mat_t* r_mats,
                           const nmod_poly_mat_t* e2_mats,
                           const nmod_poly_t* e3_polys,
                           const evmlr_otse_key_t key,
                           const evmlr_hpke_cipher_t cipher,
                           const nmod_poly_mat_t d_dagger,
                           const nmod_poly_mat_t plaintext,
                           slong L_in,
                           slong L_max);

int evmlr_enc_proof_verify(const evmlr_enc_proof_t proof,
                           const evmlr_mlpke_pk_t pk,
                           const evmlr_otse_ctx_t otse_ctx,
                           const evmlr_hpke_cipher_t cipher,
                           const nmod_poly_mat_t decrypted_plaintext,
                           slong L_max);

void evmlr_enc_proof_copy(evmlr_enc_proof_t dest, const evmlr_enc_proof_t src);

#endif // EVMLR_SHUFFLE_EVMLR_ENC_PROOF_H
