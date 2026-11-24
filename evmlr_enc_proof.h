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

/* label (label_len bytes, may be empty) is hashed into both Fiat-Shamir
 * challenges, binding the proof to a context such as an election and a
 * voter: a proof made under one label does not verify under another.
 *
 * The ML-PKE half is also bound to the OTSE half and to the OTSE ciphertext:
 * evmlr_enc_proof_verify_hpke, which needs no secret, fails if either was
 * changed, so the whole of a proof and its ciphertext can be checked for
 * tampering by anyone. */
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
                           slong L_max,
                           const uint8_t* label,
                           size_t label_len);

int evmlr_enc_proof_verify(const evmlr_enc_proof_t proof,
                           const evmlr_mlpke_pk_t pk,
                           const evmlr_otse_ctx_t otse_ctx,
                           const evmlr_hpke_cipher_t cipher,
                           const nmod_poly_mat_t decrypted_plaintext,
                           slong L_max,
                           const uint8_t* label,
                           size_t label_len);

/* Only the first half of evmlr_enc_proof_verify: that the ML-PKE ciphertexts
 * encrypt a binary OTSE key under pk, and that neither the OTSE half nor the
 * OTSE ciphertext was changed since. It needs neither the plaintext nor a
 * secret key, so anyone can run it on a ballot as soon as it is cast. */
int evmlr_enc_proof_verify_hpke(const evmlr_enc_proof_t proof,
                                const evmlr_mlpke_pk_t pk,
                                const evmlr_hpke_cipher_t cipher,
                                const uint8_t* label,
                                size_t label_len);

/* The second half: that the OTSE ciphertext encrypts decrypted_plaintext
 * (whose row count is the message length). Needs the decryption. */
int evmlr_enc_proof_verify_otse(const evmlr_enc_proof_t proof,
                                const evmlr_otse_ctx_t otse_ctx,
                                const evmlr_hpke_cipher_t cipher,
                                const nmod_poly_mat_t decrypted_plaintext,
                                const uint8_t* label,
                                size_t label_len);

/* A hash of the whole proof (both halves), for binding it into another
 * proof's label. */
void evmlr_enc_proof_digest(uint8_t out[32], const evmlr_enc_proof_t proof);

void evmlr_enc_proof_copy(evmlr_enc_proof_t dest, const evmlr_enc_proof_t src);

#endif // EVMLR_SHUFFLE_EVMLR_ENC_PROOF_H
