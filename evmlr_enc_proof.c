#include "evmlr_enc_proof.h"
#include "evmlr_utils.h"
#include <stdlib.h>

static void build_A_hpke(nmod_poly_mat_t A_hpke, const evmlr_mlpke_pk_t pk) {
    nmod_poly_mat_zero(A_hpke);
    ulong shift = (MOD_Q + 1) / 2;
    for (int i = 0; i < K_LWR; i++) {
        slong offset = i * (K_LWE + 1);
        // Copy A^T
        for (int r = 0; r < K_LWE; r++) {
            for (int c = 0; c < K_LWE; c++) {
                nmod_poly_set(nmod_poly_mat_entry(A_hpke, offset + r, offset + c),
                              nmod_poly_mat_entry(pk->A, c, r));
            }
        }
        // Copy t^T
        for (int c = 0; c < K_LWE; c++) {
            nmod_poly_set(nmod_poly_mat_entry(A_hpke, offset + K_LWE, offset + c),
                          nmod_poly_mat_entry(pk->t, c, 0));
        }
        // Set shift in the bottom-right of the block
        nmod_poly_set_coeff_ui(nmod_poly_mat_entry(A_hpke, offset + K_LWE, offset + K_LWE), 0, shift);
    }
}

static void build_A_otse(nmod_poly_mat_t A_otse, const nmod_poly_mat_t H_prime, slong L_in) {
    slong cols = 2 * ETA * (K_LWE + L_in);
    nmod_poly_mat_zero(A_otse);
    for (slong r = 0; r < L_in; r++) {
        for (slong c = 0; c < cols; c++) {
            slong b = c / (K_LWE + L_in);
            slong idx = c % (K_LWE + L_in);
            nmod_poly_struct* entry = nmod_poly_mat_entry(A_otse, r, c);
            nmod_poly_struct* h_prime_entry = nmod_poly_mat_entry(H_prime, r, idx);
            if (b < ETA) {
                nmod_poly_set(entry, h_prime_entry);
            } else {
                nmod_poly_neg(entry, h_prime_entry);
            }
        }
    }
}

static void truncate_d_dagger(nmod_poly_mat_t d_dagger_trunc, const nmod_poly_mat_t d_dagger, slong L_in, slong L_max) {
    nmod_poly_mat_zero(d_dagger_trunc);
    for (int b = 0; b < 2 * ETA; b++) {
        slong src_offset = b * (K_LWE + L_max);
        slong dst_offset = b * (K_LWE + L_in);
        for (int i = 0; i < K_LWE + L_in; i++) {
            nmod_poly_set(nmod_poly_mat_entry(d_dagger_trunc, dst_offset + i, 0),
                          nmod_poly_mat_entry(d_dagger, src_offset + i, 0));
        }
    }
}

void evmlr_enc_proof_init(evmlr_enc_proof_t proof, slong L_curr) {
    evmlr_lin_proof_ctx_t ctx_hpke;
    evmlr_lin_proof_ctx_init(ctx_hpke, K_LWR * (K_LWE + 1), K_LWR * (K_LWE + 1));
    evmlr_lin_proof_init(&proof->proof_hpke, ctx_hpke);
    evmlr_lin_proof_ctx_clear(ctx_hpke);

    evmlr_lin_proof_ctx_t ctx_otse;
    evmlr_lin_proof_ctx_init(ctx_otse, L_curr, 2 * ETA * (K_LWE + L_curr));
    evmlr_lin_proof_init(&proof->proof_otse, ctx_otse);
    evmlr_lin_proof_ctx_clear(ctx_otse);
}

void evmlr_enc_proof_clear(evmlr_enc_proof_t proof) {
    evmlr_lin_proof_clear(&proof->proof_hpke);
    evmlr_lin_proof_clear(&proof->proof_otse);
}

void evmlr_enc_proof_copy(evmlr_enc_proof_t dest, const evmlr_enc_proof_t src) {
    nmod_poly_mat_set(dest->proof_hpke.w, src->proof_hpke.w);
    nmod_poly_mat_set(dest->proof_hpke.z_1, src->proof_hpke.z_1);
    nmod_poly_mat_set(dest->proof_hpke.z_2, src->proof_hpke.z_2);

    nmod_poly_mat_set(dest->proof_otse.w, src->proof_otse.w);
    nmod_poly_mat_set(dest->proof_otse.z_1, src->proof_otse.z_1);
    nmod_poly_mat_set(dest->proof_otse.z_2, src->proof_otse.z_2);
}

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
                           slong L_max) {
    // 1. HPKE proof of correct encryption
    slong hpke_dim = K_LWR * (K_LWE + 1);
    nmod_poly_mat_t A_hpke, s1_hpke, s2_hpke, t_hpke;
    nmod_poly_mat_init(A_hpke, hpke_dim, hpke_dim, MOD_Q);
    nmod_poly_mat_init(s1_hpke, hpke_dim, 1, MOD_Q);
    nmod_poly_mat_init(s2_hpke, hpke_dim, 1, MOD_Q);
    nmod_poly_mat_init(t_hpke, hpke_dim, 1, MOD_Q);

    build_A_hpke(A_hpke, pk);

    for (int i = 0; i < K_LWR; i++) {
        slong offset = i * (K_LWE + 1);
        for (int r = 0; r < K_LWE; r++) {
            nmod_poly_set(nmod_poly_mat_entry(s1_hpke, offset + r, 0),
                          nmod_poly_mat_entry(r_mats[i], 0, r)); // r_mats is 1 x K_LWE
            nmod_poly_set(nmod_poly_mat_entry(s2_hpke, offset + r, 0),
                          nmod_poly_mat_entry(e2_mats[i], 0, r)); // e2_mats is 1 x K_LWE
            nmod_poly_set(nmod_poly_mat_entry(t_hpke, offset + r, 0),
                          nmod_poly_mat_entry(cipher->enc_cipher[i]->uT, 0, r)); // uT is 1 x K_LWE
        }
        nmod_poly_set(nmod_poly_mat_entry(s1_hpke, offset + K_LWE, 0),
                      nmod_poly_mat_entry(key->s, i, 0));
        nmod_poly_set(nmod_poly_mat_entry(s2_hpke, offset + K_LWE, 0),
                      e3_polys[i]);
        nmod_poly_set(nmod_poly_mat_entry(t_hpke, offset + K_LWE, 0),
                      cipher->enc_cipher[i]->v);
    }

    evmlr_lin_proof_ctx_t ctx_hpke;
    evmlr_lin_proof_ctx_init(ctx_hpke, hpke_dim, hpke_dim);
    evmlr_lin_prove(&proof->proof_hpke, A_hpke, s1_hpke, s2_hpke, t_hpke, ctx_hpke);
    evmlr_lin_proof_ctx_clear(ctx_hpke);

    nmod_poly_mat_clear(A_hpke);
    nmod_poly_mat_clear(s1_hpke);
    nmod_poly_mat_clear(s2_hpke);
    nmod_poly_mat_clear(t_hpke);

    // 2. OTSE proof of correct encryption
    slong otse_m = 2 * ETA * (K_LWE + L_in);
    nmod_poly_mat_t A_otse, s1_otse, s2_otse, t_otse;
    nmod_poly_mat_init(A_otse, L_in, otse_m, MOD_Q);
    nmod_poly_mat_init(s1_otse, otse_m, 1, MOD_Q);
    nmod_poly_mat_init(s2_otse, L_in, 1, MOD_Q);
    nmod_poly_mat_init(t_otse, L_in, 1, MOD_Q);

    build_A_otse(A_otse, otse_ctx->H_prime, L_in);
    truncate_d_dagger(s1_otse, d_dagger, L_in, L_max);
    nmod_poly_mat_zero(s2_otse);

    for (slong i = 0; i < L_in; i++) {
        nmod_poly_sub(nmod_poly_mat_entry(t_otse, i, 0),
                      nmod_poly_mat_entry(cipher->otse_cipher->c, i, 0),
                      nmod_poly_mat_entry(plaintext, i, 0));
    }

    evmlr_lin_proof_ctx_t ctx_otse;
    evmlr_lin_proof_ctx_init(ctx_otse, L_in, otse_m);
    evmlr_lin_prove(&proof->proof_otse, A_otse, s1_otse, s2_otse, t_otse, ctx_otse);
    evmlr_lin_proof_ctx_clear(ctx_otse);

    nmod_poly_mat_clear(A_otse);
    nmod_poly_mat_clear(s1_otse);
    nmod_poly_mat_clear(s2_otse);
    nmod_poly_mat_clear(t_otse);
}

int evmlr_enc_proof_verify(const evmlr_enc_proof_t proof,
                           const evmlr_mlpke_pk_t pk,
                           const evmlr_otse_ctx_t otse_ctx,
                           const evmlr_hpke_cipher_t cipher,
                           const nmod_poly_mat_t decrypted_plaintext,
                           slong L_max) {
    slong L_in = decrypted_plaintext->r;

    // 1. Verify HPKE proof
    slong hpke_dim = K_LWR * (K_LWE + 1);
    nmod_poly_mat_t A_hpke, t_hpke;
    nmod_poly_mat_init(A_hpke, hpke_dim, hpke_dim, MOD_Q);
    nmod_poly_mat_init(t_hpke, hpke_dim, 1, MOD_Q);

    build_A_hpke(A_hpke, pk);

    for (int i = 0; i < K_LWR; i++) {
        slong offset = i * (K_LWE + 1);
        for (int r = 0; r < K_LWE; r++) {
            nmod_poly_set(nmod_poly_mat_entry(t_hpke, offset + r, 0),
                          nmod_poly_mat_entry(cipher->enc_cipher[i]->uT, 0, r));
        }
        nmod_poly_set(nmod_poly_mat_entry(t_hpke, offset + K_LWE, 0),
                      cipher->enc_cipher[i]->v);
    }

    evmlr_lin_proof_ctx_t ctx_hpke;
    evmlr_lin_proof_ctx_init(ctx_hpke, hpke_dim, hpke_dim);
    int hpke_ok = evmlr_lin_verify(&proof->proof_hpke, A_hpke, t_hpke, ctx_hpke);
    evmlr_lin_proof_ctx_clear(ctx_hpke);

    nmod_poly_mat_clear(A_hpke);
    nmod_poly_mat_clear(t_hpke);

    if (!hpke_ok) return 0;

    // 2. Verify OTSE proof
    slong otse_m = 2 * ETA * (K_LWE + L_in);
    nmod_poly_mat_t A_otse, t_otse;
    nmod_poly_mat_init(A_otse, L_in, otse_m, MOD_Q);
    nmod_poly_mat_init(t_otse, L_in, 1, MOD_Q);

    build_A_otse(A_otse, otse_ctx->H_prime, L_in);

    for (slong i = 0; i < L_in; i++) {
        nmod_poly_sub(nmod_poly_mat_entry(t_otse, i, 0),
                      nmod_poly_mat_entry(cipher->otse_cipher->c, i, 0),
                      nmod_poly_mat_entry(decrypted_plaintext, i, 0));
    }

    evmlr_lin_proof_ctx_t ctx_otse;
    evmlr_lin_proof_ctx_init(ctx_otse, L_in, otse_m);
    int otse_ok = evmlr_lin_verify(&proof->proof_otse, A_otse, t_otse, ctx_otse);
    evmlr_lin_proof_ctx_clear(ctx_otse);

    nmod_poly_mat_clear(A_otse);
    nmod_poly_mat_clear(t_otse);

    return otse_ok;
}
