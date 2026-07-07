#include "evmlr_voting.h"
#include "evmlr_utils.h"
#include <string.h>

#ifdef MAIN
#include "test.h"
#include "bench.h"
#include "cpucycles.h"
#include <sys/random.h>
#endif

// Serialization of HPKE ciphertext layers
static void serialize_layer(nmod_poly_mat_t msg, const evmlr_voting_layer_t ct, slong next_layer_len) {
    slong idx = 0;
    for (int i = 0; i < K_LWR; i++) {
        for (int j = 0; j < K_LWE; j++) {
            nmod_poly_set(nmod_poly_mat_entry(msg, idx++, 0), nmod_poly_mat_entry(ct->enc_cipher[i]->uT, 0, j));
        }
        nmod_poly_set(nmod_poly_mat_entry(msg, idx++, 0), ct->enc_cipher[i]->v);
    }
    for (int i = 0; i < next_layer_len; i++) {
        nmod_poly_set(nmod_poly_mat_entry(msg, idx++, 0), nmod_poly_mat_entry(ct->otse_cipher->c, i, 0));
    }
}

static void deserialize_layer(evmlr_voting_layer_t ct, const nmod_poly_mat_t msg, slong next_layer_len) {
    slong idx = 0;
    for (int i = 0; i < K_LWR; i++) {
        nmod_poly_mat_init(ct->enc_cipher[i]->uT, 1, K_LWE, MOD_Q);
        nmod_poly_init(ct->enc_cipher[i]->v, MOD_Q);
        for (int j = 0; j < K_LWE; j++) {
            nmod_poly_set(nmod_poly_mat_entry(ct->enc_cipher[i]->uT, 0, j), nmod_poly_mat_entry(msg, idx++, 0));
        }
        nmod_poly_set(ct->enc_cipher[i]->v, nmod_poly_mat_entry(msg, idx++, 0));
    }
    nmod_poly_mat_init(ct->otse_cipher->c, next_layer_len, 1, MOD_Q);
    for (int i = 0; i < next_layer_len; i++) {
        nmod_poly_set(nmod_poly_mat_entry(ct->otse_cipher->c, i, 0), nmod_poly_mat_entry(msg, idx++, 0));
    }
}

static void voting_encrypt_and_prove_layer(evmlr_voting_layer_t ct, evmlr_enc_proof_t proof, const nmod_poly_mat_t msg, const evmlr_mlpke_pk_t pk, const evmlr_voting_pp_t pp, flint_rand_t state) {
    slong L_in = msg->r;

    // Pad msg to L_max
    nmod_poly_mat_t msg_padded;
    nmod_poly_mat_init(msg_padded, pp->L_max, 1, MOD_Q);
    nmod_poly_mat_zero(msg_padded);
    for (slong i = 0; i < L_in; i++) {
        nmod_poly_set(nmod_poly_mat_entry(msg_padded, i, 0), nmod_poly_mat_entry(msg, i, 0));
    }

    // Manually perform HPKE encryption to capture the secrets
    // Sample key/seed
    evmlr_otse_key_t key;
    evmlr_otse_keygen(key, state);

    // Allocate ML-PKE secrets
    nmod_poly_mat_t r_mats[K_LWR];
    nmod_poly_mat_t e2_mats[K_LWR];
    nmod_poly_t e3_polys[K_LWR];

    evmlr_hpke_cipher_t cipher;
    for (int i = 0; i < K_LWR; i++) {
        nmod_poly_mat_init(r_mats[i], 1, K_LWE, MOD_Q);
        nmod_poly_mat_init(e2_mats[i], 1, K_LWE, MOD_Q);
        evmlr_utils_binom_sample_mat_ring(r_mats[i], ETA);
        evmlr_utils_binom_sample_mat_ring(e2_mats[i], ETA);

        nmod_poly_init(e3_polys[i], MOD_Q);
        evmlr_utils_binom_sample_ring(e3_polys[i], ETA);

        evmlr_mlpke_enc_with_secrets(cipher->enc_cipher[i], r_mats[i], e2_mats[i], e3_polys[i],
                                     nmod_poly_mat_entry(key->s, i, 0), pk, pp->hpke_ctx->enc_ctx);
    }

    nmod_poly_mat_t d_dagger;
    nmod_poly_mat_init(d_dagger, 2 * ETA * (K_LWE + pp->L_max), 1, MOD_Q);
    evmlr_otse_encrypt(cipher->otse_cipher, d_dagger, msg_padded, key, pp->hpke_ctx->otse_ctx);

    // Copy cipher to ct
    evmlr_voting_layer_clear(ct);
    evmlr_voting_layer_init(ct, L_in);
    for (int i = 0; i < K_LWR; i++) {
        nmod_poly_mat_set(ct->enc_cipher[i]->uT, cipher->enc_cipher[i]->uT);
        nmod_poly_set(ct->enc_cipher[i]->v, cipher->enc_cipher[i]->v);
    }
    for (slong i = 0; i < L_in; i++) {
        nmod_poly_set(nmod_poly_mat_entry(ct->otse_cipher->c, i, 0), nmod_poly_mat_entry(cipher->otse_cipher->c, i, 0));
    }

    // Generate ZK proofs
    evmlr_enc_proof_prove(proof, pk, pp->hpke_ctx->otse_ctx, r_mats, e2_mats, e3_polys, key, cipher, d_dagger, msg, L_in, pp->L_max);

    // Clear secrets
    for (int i = 0; i < K_LWR; i++) {
        nmod_poly_mat_clear(r_mats[i]);
        nmod_poly_mat_clear(e2_mats[i]);
        nmod_poly_clear(e3_polys[i]);
    }
    evmlr_otse_keyclear(key);
    nmod_poly_mat_clear(d_dagger);
    nmod_poly_mat_clear(msg_padded);
    evmlr_hpke_cipher_clear(cipher);
}

extern void calc_a(nmod_poly_mat_t a, nmod_poly_mat_t d_dagger, const evmlr_otse_key_t key, const evmlr_otse_ctx_t ctx);

static void voting_decrypt_layer(nmod_poly_mat_t msg, nmod_poly_mat_t d_dagger, nmod_poly_mat_t a, const evmlr_voting_layer_t ct, const evmlr_mlpke_sk_t sk, const evmlr_voting_pp_t pp) {
    slong L_in = ct->otse_cipher->c->r;

    // Reconstruct key s
    evmlr_otse_key_t key;
    nmod_poly_t key_poly;
    nmod_poly_init(key_poly, MOD_Q);
    nmod_poly_mat_init(key->s, K_LWR, 1, MOD_Q);

    for (int i = 0; i < K_LWR; i++) {
        evmlr_mlpke_dec(key_poly, ct->enc_cipher[i], sk, pp->hpke_ctx->enc_ctx);
        nmod_poly_set(nmod_poly_mat_entry(key->s, i, 0), key_poly);
    }
    nmod_poly_clear(key_poly);

    // Call calc_a to get a and extract d_dagger
    nmod_poly_mat_init(a, pp->hpke_ctx->otse_ctx->L, 1, MOD_Q);
    calc_a(a, d_dagger, key, pp->hpke_ctx->otse_ctx);

    // m = c - a (truncated to L_in)
    nmod_poly_mat_init(msg, L_in, 1, MOD_Q);
    for (slong i = 0; i < L_in; i++) {
        nmod_poly_struct* msg_val = nmod_poly_mat_entry(msg, i, 0);
        nmod_poly_sub(msg_val, nmod_poly_mat_entry(ct->otse_cipher->c, i, 0), nmod_poly_mat_entry(a, i, 0));
    }

    evmlr_otse_keyclear(key);
}

static void sync_hpke_ctx(evmlr_hpke_ctx_t dst, const evmlr_hpke_ctx_t src) {
    nmod_poly_mat_set(dst->otse_ctx->H, src->otse_ctx->H);
    nmod_poly_mat_set(dst->otse_ctx->H_prime, src->otse_ctx->H_prime);
}

// setup(k, state) -> (pp, sk)
void evmlr_voting_setup(evmlr_voting_pp_t pp, evmlr_voting_sk_t sk, slong k, flint_rand_t state) {
    pp->k = k;
    sk->k = k;

    // Outer message input has length k * K_LWR * (K_LWE + 1) + 1.
    // Decrypting layer 1 yields length (k - 1) * K_LWR * (K_LWE + 1) + 1.
    // This is the maximum length of the messages being shuffled.
    pp->L_max = (k - 1) * K_LWR * (K_LWE + 1) + 1;

    // Initialize HPKE and Shuffle contexts.
    evmlr_hpke_ctx_init(pp->hpke_ctx, pp->L_max, state);
    evmlr_shuffle_ctx_init(pp->shuf_ctx, SHUFFLE_N_MSGS, pp->L_max, state);
    sync_hpke_ctx(pp->shuf_ctx->hpke_ctx, pp->hpke_ctx);

    // Allocate keys
    pp->pks = (evmlr_mlpke_pk_t*) malloc(k * sizeof(evmlr_mlpke_pk_t));
    sk->sks = (evmlr_mlpke_sk_t*) malloc(k * sizeof(evmlr_mlpke_sk_t));

    for (slong j = 0; j < k; j++) {
        nmod_poly_mat_init(pp->pks[j]->A, K_LWE, K_LWE, MOD_Q);
        nmod_poly_mat_init(pp->pks[j]->t, K_LWE, 1, MOD_Q);
        nmod_poly_mat_init(sk->sks[j]->s, K_LWE, 1, MOD_Q);

        evmlr_mlpke_keypair_t kp;
        evmlr_mlpke_keypair_gen(kp, state, pp->hpke_ctx->enc_ctx);
        
        nmod_poly_mat_set(pp->pks[j]->A, kp->pk->A);
        nmod_poly_mat_set(pp->pks[j]->t, kp->pk->t);
        nmod_poly_mat_set(sk->sks[j]->s, kp->sk->s);

        evmlr_mlpke_keypair_clear(kp);
    }
}

// casting(submission, voter_id, vote, pp, state)
void evmlr_voting_casting(evmlr_voting_submission_t submission, slong voter_id, const nmod_poly_mat_t vote, const evmlr_voting_pp_t pp, flint_rand_t state) {
    submission->voter_id = voter_id;
    slong k = pp->k;

    nmod_poly_mat_t current_msg;
    nmod_poly_mat_init(current_msg, vote->r, vote->c, MOD_Q);
    nmod_poly_mat_set(current_msg, vote);

    for (slong j = k - 1; j >= 0; j--) {
        if (j > 0) {
            evmlr_voting_layer_t temp_layer;
            slong L_otse = (k - j - 1) * K_LWR * (K_LWE + 1) + 1;
            evmlr_voting_layer_init(temp_layer, L_otse);

            // Encrypt current_msg under pks[j] and generate proof
            voting_encrypt_and_prove_layer(temp_layer, &submission->proofs[j]->proof, current_msg, pp->pks[j], pp, state);

            // Next plaintext is the serialized ciphertext
            nmod_poly_mat_clear(current_msg);
            slong L_curr = (k - j) * K_LWR * (K_LWE + 1) + 1;
            nmod_poly_mat_init(current_msg, L_curr, 1, MOD_Q);
            serialize_layer(current_msg, temp_layer, L_curr - K_LWR * (K_LWE + 1));

            evmlr_voting_layer_clear(temp_layer);
        } else {
            // Encrypt current_msg under pks[0] to the outermost layer and generate proof
            voting_encrypt_and_prove_layer(submission->ballot->layer, &submission->proofs[0]->proof, current_msg, pp->pks[0], pp, state);
        }
    }

    nmod_poly_mat_clear(current_msg);
}

// verify_proof_of_encryption(pp, proof, ciphertext_layer, decrypted_plaintext, layer_idx)
int evmlr_voting_verify_enc(const evmlr_voting_pp_t pp, const evmlr_voting_proof_enc_t proof, const evmlr_voting_layer_t ct, const nmod_poly_mat_t decrypted_plaintext, slong layer_idx) {
    // Reconstruct the hpke cipher representation
    evmlr_hpke_cipher_t cipher;
    for (int i = 0; i < K_LWR; i++) {
        nmod_poly_mat_init(cipher->enc_cipher[i]->uT, 1, K_LWE, MOD_Q);
        nmod_poly_mat_set(cipher->enc_cipher[i]->uT, ct->enc_cipher[i]->uT);
        nmod_poly_init(cipher->enc_cipher[i]->v, MOD_Q);
        nmod_poly_set(cipher->enc_cipher[i]->v, ct->enc_cipher[i]->v);
    }
    // OTSE ciphertext representation in voting is of size L_in, but verify expects L_max
    nmod_poly_mat_init(cipher->otse_cipher->c, pp->L_max, 1, MOD_Q);
    nmod_poly_mat_zero(cipher->otse_cipher->c);
    for (slong i = 0; i < ct->otse_cipher->c->r; i++) {
        nmod_poly_set(nmod_poly_mat_entry(cipher->otse_cipher->c, i, 0),
                      nmod_poly_mat_entry(ct->otse_cipher->c, i, 0));
    }

    int ok = evmlr_enc_proof_verify(&proof->proof, pp->pks[layer_idx], pp->hpke_ctx->otse_ctx, cipher, decrypted_plaintext, pp->L_max);

    evmlr_hpke_cipher_clear(cipher);
    return ok;
}

// counting(...)
slong evmlr_voting_counting(nmod_poly_mat_t* results, evmlr_shuffle_proof_t* shuffle_proofs, 
                          const evmlr_voting_submission_t* submissions, slong num_submissions,
                          const evmlr_voting_pp_t pp, const evmlr_voting_sk_t sk, flint_rand_t state) {
    slong k = pp->k;

    // 1. Apply revote policy to get the list of unique voter submissions.
    evmlr_voting_submission_struct** active_subs = (evmlr_voting_submission_struct**) malloc(num_submissions * sizeof(evmlr_voting_submission_struct*));
    slong N = 0;
    for (slong i = num_submissions - 1; i >= 0; i--) {
        slong voter_id = submissions[i]->voter_id;
        int found = 0;
        for (slong j = 0; j < N; j++) {
            if (active_subs[j]->voter_id == voter_id) {
                found = 1;
                break;
            }
        }
        if (!found) {
            active_subs[N++] = (evmlr_voting_submission_struct*) submissions[i];
        }
    }
    // Reverse active_subs to maintain chronological order
    for (slong i = 0; i < N / 2; i++) {
        evmlr_voting_submission_struct* tmp = active_subs[i];
        active_subs[i] = active_subs[N - 1 - i];
        active_subs[N - 1 - i] = tmp;
    }

    if (N == 0) {
        free(active_subs);
        return 0;
    }

    // Re-initialize shuffle context for N unique submissions
    evmlr_shuffle_ctx_clear((evmlr_shuffle_ctx_struct *)pp->shuf_ctx);
    evmlr_shuffle_ctx_init((evmlr_shuffle_ctx_struct *)pp->shuf_ctx, N, pp->L_max, state);
    sync_hpke_ctx(((evmlr_shuffle_ctx_struct *)pp->shuf_ctx)->hpke_ctx, pp->hpke_ctx);

    // 2. Initialize current ciphertexts and proofs
    evmlr_voting_layer_struct* curr_cts = (evmlr_voting_layer_struct*) malloc(N * sizeof(evmlr_voting_layer_struct));

    // Array of remaining proofs for each voter
    evmlr_voting_proof_enc_t** remaining_proofs = (evmlr_voting_proof_enc_t**) malloc(N * sizeof(evmlr_voting_proof_enc_t*));
    for (slong i = 0; i < N; i++) {
        remaining_proofs[i] = (evmlr_voting_proof_enc_t*) malloc(k * sizeof(evmlr_voting_proof_enc_t));
        for (slong l = 0; l < k; l++) {
            slong L_curr_l = (k - l - 1) * K_LWR * (K_LWE + 1) + 1;
            evmlr_enc_proof_init(&remaining_proofs[i][l]->proof, L_curr_l);
            evmlr_enc_proof_copy(&remaining_proofs[i][l]->proof, &active_subs[i]->proofs[l]->proof);
        }
    }

    slong L_curr = pp->L_max;
    for (slong i = 0; i < N; i++) {
        evmlr_voting_layer_init(&curr_cts[i], L_curr);
        for (int r = 0; r < K_LWR; r++) {
            nmod_poly_mat_set(curr_cts[i].enc_cipher[r]->uT, active_subs[i]->ballot->layer->enc_cipher[r]->uT);
            nmod_poly_set(curr_cts[i].enc_cipher[r]->v, active_subs[i]->ballot->layer->enc_cipher[r]->v);
        }
        nmod_poly_mat_set(curr_cts[i].otse_cipher->c, active_subs[i]->ballot->layer->otse_cipher->c);
    }

    nmod_poly_mat_t* decrypted_msgs = (nmod_poly_mat_t*) malloc(N * sizeof(nmod_poly_mat_t));
    nmod_poly_mat_t* shuffled_msgs = (nmod_poly_mat_t*) malloc(N * sizeof(nmod_poly_mat_t));

    // Process layer by layer
    for (slong j = 0; j < k; j++) {
        slong L_next = (k - j - 1) * K_LWR * (K_LWE + 1) + 1;
        slong L_curr_j = (k - j) * K_LWR * (K_LWE + 1) + 1;

        // Dynamically allocate and initialize curr_proofs for the current layer j
        evmlr_voting_proof_enc_t* curr_proofs = (evmlr_voting_proof_enc_t*) malloc(N * sizeof(evmlr_voting_proof_enc_t));
        for (slong i = 0; i < N; i++) {
            evmlr_enc_proof_init(&curr_proofs[i]->proof, L_curr_j - K_LWR * (K_LWE + 1));
            evmlr_enc_proof_copy(&curr_proofs[i]->proof, &remaining_proofs[i][j]->proof);
        }

        nmod_poly_mat_t* d_dagger = (nmod_poly_mat_t*) malloc(N * sizeof(nmod_poly_mat_t));
        nmod_poly_mat_t* a_mats = (nmod_poly_mat_t*) malloc(N * sizeof(nmod_poly_mat_t));

        // Decrypt and Verify
        for (slong i = 0; i < N; i++) {
            voting_decrypt_layer(decrypted_msgs[i], d_dagger[i], a_mats[i], &curr_cts[i], sk->sks[j], pp);
            int ok = evmlr_voting_verify_enc(pp, curr_proofs[i], &curr_cts[i], decrypted_msgs[i], j);
            if (!ok) {
                printf("Warning: proof of encryption failed for voter %ld at server %ld\n", active_subs[i]->voter_id, j);
            }
        }

        // Shuffle
        size_t* pi = (size_t *) malloc(N * sizeof(size_t));
        evmlr_utils_new_perm(pi, N);

        for (slong i = 0; i < N; i++) {
            nmod_poly_mat_init(shuffled_msgs[i], L_next, 1, MOD_Q);
            nmod_poly_mat_set(shuffled_msgs[i], decrypted_msgs[pi[i]]);
        }

        // Shuffle remaining proofs
        evmlr_voting_proof_enc_t** temp_remaining_proofs = (evmlr_voting_proof_enc_t**) malloc(N * sizeof(evmlr_voting_proof_enc_t*));
        for (slong i = 0; i < N; i++) {
            temp_remaining_proofs[i] = (evmlr_voting_proof_enc_t*) malloc(k * sizeof(evmlr_voting_proof_enc_t));
            for (slong l = 0; l < k; l++) {
                slong L_curr_l = (k - l - 1) * K_LWR * (K_LWE + 1) + 1;
                evmlr_enc_proof_init(&temp_remaining_proofs[i][l]->proof, L_curr_l);
                evmlr_enc_proof_copy(&temp_remaining_proofs[i][l]->proof, &remaining_proofs[pi[i]][l]->proof);
            }
        }
        for (slong i = 0; i < N; i++) {
            for (slong l = 0; l < k; l++) {
                evmlr_enc_proof_clear(&remaining_proofs[i][l]->proof);
            }
            free(remaining_proofs[i]);
            remaining_proofs[i] = temp_remaining_proofs[i];
        }
        free(temp_remaining_proofs);

        // Generate Proof of Shuffle
        evmlr_shuffle_sp_t shuf_sp;
        evmlr_shuffle_pp_t shuf_pp;

        shuf_sp->pi = pi;
        shuf_sp->d_dagger = d_dagger;
        nmod_poly_mat_init(shuf_sp->r_D, 2*K_SIS, 1, MOD_Q);
        evmlr_commit_sample_r(shuf_sp->r_D);

        nmod_poly_mat_t d_flat;
        nmod_poly_mat_init(d_flat, pp->shuf_ctx->com_ctx->N, 1, MOD_Q);
        nmod_poly_mat_zero(d_flat);
        slong row_len = (2*ETA) * (K_LWE + pp->L_max);
        for (slong i = 0; i < N; i++) {
            for (slong r = 0; r < row_len; r++) {
                nmod_poly_struct* d_dag_ir = nmod_poly_mat_entry(shuf_sp->d_dagger[i], r, 0);
                nmod_poly_struct* poly = nmod_poly_mat_entry(d_flat, i * row_len + r, 0);
                nmod_poly_set(poly, d_dag_ir);
            }
        }
        evmlr_commit(shuf_pp->D, d_flat, shuf_sp->r_D, pp->shuf_ctx->com_ctx);
        nmod_poly_mat_clear(d_flat);

        shuf_pp->c_hat = (nmod_poly_mat_t *) malloc(N * sizeof(nmod_poly_mat_t));
        for (slong i = 0; i < N; i++) {
            nmod_poly_mat_init(shuf_pp->c_hat[i], pp->L_max, 1, MOD_Q);
            nmod_poly_mat_zero(shuf_pp->c_hat[i]);
            for (slong r = 0; r < L_next; r++) {
                nmod_poly_set(nmod_poly_mat_entry(shuf_pp->c_hat[i], r, 0), nmod_poly_mat_entry(shuffled_msgs[i], r, 0));
            }
        }
        for (slong i = 0; i < N; i++) {
            for (slong r = L_next; r < pp->L_max; r++) {
                nmod_poly_struct* hat_val = nmod_poly_mat_entry(shuf_pp->c_hat[i], r, 0);
                nmod_poly_struct* a_val = nmod_poly_mat_entry(a_mats[pi[i]], r, 0);
                nmod_poly_neg(hat_val, a_val);
            }
        }

        shuf_pp->c_star = (evmlr_otse_ciphertext_t *) malloc(N * sizeof(evmlr_otse_ciphertext_t));
        for (slong i = 0; i < N; i++) {
            nmod_poly_mat_init(shuf_pp->c_star[i]->c, pp->L_max, 1, MOD_Q);
            nmod_poly_mat_zero(shuf_pp->c_star[i]->c);
            for (slong r = 0; r < curr_cts[i].otse_cipher->c->r; r++) {
                nmod_poly_set(nmod_poly_mat_entry(shuf_pp->c_star[i]->c, r, 0), nmod_poly_mat_entry(curr_cts[i].otse_cipher->c, r, 0));
            }
        }

        evmlr_proof_init(shuffle_proofs[j], pp->shuf_ctx);
        evmlr_shuffle_prove(shuffle_proofs[j], shuf_sp, shuf_pp, pp->shuf_ctx, state);

        int shuffle_ok = evmlr_shuffle_verify(shuffle_proofs[j], pp->shuf_ctx, shuf_pp);
        if (!shuffle_ok) {
            printf("Warning: Shuffle verification failed at server %ld!\n", j);
        }

        // Clean up pp/sp
        for (slong i = 0; i < N; i++) {
            nmod_poly_mat_clear(shuf_pp->c_hat[i]);
            nmod_poly_mat_clear(shuf_pp->c_star[i]->c);
            nmod_poly_mat_clear(a_mats[i]);
        }
        free(shuf_pp->c_hat);
        free(shuf_pp->c_star);
        free(a_mats);
        evmlr_commit_clear(shuf_pp->D);
        nmod_poly_mat_clear(shuf_sp->r_D);
        free(shuf_sp->pi);
        for (slong i = 0; i < N; i++) {
            nmod_poly_mat_clear(d_dagger[i]);
        }
        free(d_dagger);

        if (j == k - 1) {
            for (slong i = 0; i < N; i++) {
                nmod_poly_mat_init(results[i], L_next, 1, MOD_Q);
                nmod_poly_mat_set(results[i], shuffled_msgs[i]);
            }
        } else {
            slong L_next_ct = L_next;
            for (slong i = 0; i < N; i++) {
                evmlr_voting_layer_clear(&curr_cts[i]);
                deserialize_layer(&curr_cts[i], shuffled_msgs[i], L_next_ct - K_LWR * (K_LWE + 1));
            }
        }

        for (slong i = 0; i < N; i++) {
            evmlr_enc_proof_clear(&curr_proofs[i]->proof);
        }
        free(curr_proofs);

        for (slong i = 0; i < N; i++) {
            nmod_poly_mat_clear(decrypted_msgs[i]);
            nmod_poly_mat_clear(shuffled_msgs[i]);
        }
    }

    free(decrypted_msgs);
    free(shuffled_msgs);
    for (slong i = 0; i < N; i++) {
        evmlr_voting_layer_clear(&curr_cts[i]);
        for (slong l = 0; l < k; l++) {
            evmlr_enc_proof_clear(&remaining_proofs[i][l]->proof);
        }
        free(remaining_proofs[i]);
    }
    free(curr_cts);
    free(remaining_proofs);
    free(active_subs);

    return N;
}

void evmlr_voting_layer_init(evmlr_voting_layer_t layer, slong otse_len) {
    for (int i = 0; i < K_LWR; i++) {
        nmod_poly_mat_init(layer->enc_cipher[i]->uT, 1, K_LWE, MOD_Q);
        nmod_poly_init(layer->enc_cipher[i]->v, MOD_Q);
    }
    nmod_poly_mat_init(layer->otse_cipher->c, otse_len, 1, MOD_Q);
}

void evmlr_voting_layer_clear(evmlr_voting_layer_t layer) {
    for (int i = 0; i < K_LWR; i++) {
        nmod_poly_mat_clear(layer->enc_cipher[i]->uT);
        nmod_poly_clear(layer->enc_cipher[i]->v);
    }
    nmod_poly_mat_clear(layer->otse_cipher->c);
}

void evmlr_voting_submission_init(evmlr_voting_submission_t sub, const evmlr_voting_pp_t pp) {
    slong k = pp->k;
    evmlr_voting_layer_init(sub->ballot->layer, pp->L_max);
    sub->proofs = (evmlr_voting_proof_enc_t*) malloc(k * sizeof(evmlr_voting_proof_enc_t));
    for (slong j = 0; j < k; j++) {
        slong L_curr = (k - j - 1) * K_LWR * (K_LWE + 1) + 1;
        evmlr_enc_proof_init(&sub->proofs[j]->proof, L_curr);
    }
}

void evmlr_voting_submission_clear(evmlr_voting_submission_t sub, const evmlr_voting_pp_t pp) {
    slong k = pp->k;
    evmlr_voting_layer_clear(sub->ballot->layer);
    for (slong j = 0; j < k; j++) {
        evmlr_enc_proof_clear(&sub->proofs[j]->proof);
    }
    free(sub->proofs);
}

void evmlr_voting_pp_clear(evmlr_voting_pp_t pp) {
    for (slong j = 0; j < pp->k; j++) {
        nmod_poly_mat_clear(pp->pks[j]->A);
        nmod_poly_mat_clear(pp->pks[j]->t);
    }
    free(pp->pks);
    evmlr_hpke_ctx_clear(pp->hpke_ctx);
    evmlr_shuffle_ctx_clear(pp->shuf_ctx);
}

void evmlr_voting_sk_clear(evmlr_voting_sk_t sk) {
    for (slong j = 0; j < sk->k; j++) {
        nmod_poly_mat_clear(sk->sks[j]->s);
    }
    free(sk->sks);
}


#ifdef MAIN

static void test_voting_protocol(flint_rand_t state) {
    slong k_servers = 3;
    evmlr_voting_pp_t pp;
    evmlr_voting_sk_t sk;

    evmlr_voting_setup(pp, sk, k_servers, state);

    // Prepare 6 submissions for 5 voters (voter 2 submits twice)
    slong num_subs = 6;
    evmlr_voting_submission_t subs[num_subs];
    for (slong i = 0; i < num_subs; i++) {
        evmlr_voting_submission_init(subs[i], pp);
    }

    // Set voter IDs and votes (each vote is a 1x1 matrix with a distinct candidate polynomial)
    slong voter_ids[6] = {1, 2, 3, 2, 4, 5};
    nmod_poly_mat_t expected_votes[5];
    for (int i = 0; i < 5; i++) {
        nmod_poly_mat_init(expected_votes[i], 1, 1, MOD_Q);
    }
    
    // Voter 1: vote candidate 10
    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(expected_votes[0], 0, 0), 0, 10);
    // Voter 2: first vote is candidate 99 (superseded), second vote is candidate 20
    nmod_poly_mat_t superseded_vote;
    nmod_poly_mat_init(superseded_vote, 1, 1, MOD_Q);
    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(superseded_vote, 0, 0), 0, 99);
    
    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(expected_votes[1], 0, 0), 0, 20);
    // Voter 3: vote candidate 30
    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(expected_votes[2], 0, 0), 0, 30);
    // Voter 4: vote candidate 40
    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(expected_votes[3], 0, 0), 0, 40);
    // Voter 5: vote candidate 50
    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(expected_votes[4], 0, 0), 0, 50);

    // Perform casting
    evmlr_voting_casting(subs[0], voter_ids[0], expected_votes[0], pp, state); // Voter 1
    evmlr_voting_casting(subs[1], voter_ids[1], superseded_vote, pp, state);   // Voter 2 (first vote)
    evmlr_voting_casting(subs[2], voter_ids[2], expected_votes[2], pp, state); // Voter 3
    evmlr_voting_casting(subs[3], voter_ids[3], expected_votes[1], pp, state); // Voter 2 (second vote)
    evmlr_voting_casting(subs[4], voter_ids[4], expected_votes[3], pp, state); // Voter 4
    evmlr_voting_casting(subs[5], voter_ids[5], expected_votes[4], pp, state); // Voter 5

    // Run counting
    nmod_poly_mat_t results[5];
    evmlr_shuffle_proof_t shuffle_proofs[3];
    
    slong num_results = evmlr_voting_counting(results, shuffle_proofs, subs, num_subs, pp, sk, state);

    TEST_BEGIN("voting protocol: number of valid results is correct") {
        TEST_ASSERT(num_results == 5, end);
    } TEST_END;

    TEST_BEGIN("voting protocol: decrypted votes match expected votes") {
        // Since the votes are shuffled, we check that each result matches one expected vote, and all are unique
        int matched[5] = {0, 0, 0, 0, 0};
        for (int i = 0; i < 5; i++) {
            int found = 0;
            for (int j = 0; j < 5; j++) {
                if (!matched[j] && nmod_poly_mat_equal(results[i], expected_votes[j]) == 1) {
                    matched[j] = 1;
                    found = 1;
                    break;
                }
            }
            TEST_ASSERT(found == 1, end);
        }
    } TEST_END;

end:
    nmod_poly_mat_clear(superseded_vote);
    for (int i = 0; i < 5; i++) {
        nmod_poly_mat_clear(expected_votes[i]);
        nmod_poly_mat_clear(results[i]);
    }
    for (int i = 0; i < 3; i++) {
        evmlr_proof_clear(shuffle_proofs[i], pp->shuf_ctx);
    }
    for (slong i = 0; i < num_subs; i++) {
        evmlr_voting_submission_clear(subs[i], pp);
    }
    evmlr_voting_pp_clear(pp);
    evmlr_voting_sk_clear(sk);
}

void bench_voting_protocol(flint_rand_t state) {
    slong k_servers = 2;
    slong N = 100; // 10

    printf("\n=== Voting Protocol Benchmark (N = %ld voters, k = %ld servers) ===\n", N, k_servers);

    // 1. Setup Phase
    uint64_t start_setup = cpucycles();
    evmlr_voting_pp_t pp;
    evmlr_voting_sk_t sk;
    evmlr_voting_setup(pp, sk, k_servers, state);
    uint64_t end_setup = cpucycles();
    printf("Setup Phase: %lu cycles\n", (unsigned long)(end_setup - start_setup));

    // Initialize N submissions
    evmlr_voting_submission_t subs[N];
    for (slong i = 0; i < N; i++) {
        evmlr_voting_submission_init(subs[i], pp);
    }

    nmod_poly_mat_t votes[N];
    for (slong i = 0; i < N; i++) {
        nmod_poly_mat_init(votes[i], 1, 1, MOD_Q);
        nmod_poly_set_coeff_ui(nmod_poly_mat_entry(votes[i], 0, 0), 0, i + 1);
    }

    // 2. Casting Phase (Single Voter vs Total)
    uint64_t start_casting = cpucycles();
    for (slong i = 0; i < N; i++) {
        evmlr_voting_casting(subs[i], i + 1, votes[i], pp, state);
    }
    uint64_t end_casting = cpucycles();
    printf("Casting Phase (total for %ld voters): %lu cycles (avg %lu per voter)\n",
           N, (unsigned long)(end_casting - start_casting), (unsigned long)((end_casting - start_casting) / N));

    // 3. Counting/Decryption/Shuffling/Verification Phases
    // We will measure each sub-phase inside a custom counting loop to isolate the cycles
    nmod_poly_mat_t results[N];
    evmlr_shuffle_proof_t shuffle_proofs[k_servers];

    // Re-initialize shuffle context for N unique submissions
    evmlr_shuffle_ctx_clear((evmlr_shuffle_ctx_struct *)pp->shuf_ctx);
    evmlr_shuffle_ctx_init((evmlr_shuffle_ctx_struct *)pp->shuf_ctx, N, pp->L_max, state);
    sync_hpke_ctx(((evmlr_shuffle_ctx_struct *)pp->shuf_ctx)->hpke_ctx, pp->hpke_ctx);

    // Initialize current ciphertexts and proofs
    evmlr_voting_layer_struct* curr_cts = (evmlr_voting_layer_struct*) malloc(N * sizeof(evmlr_voting_layer_struct));

    // Array of remaining proofs for each voter
    evmlr_voting_proof_enc_t** remaining_proofs = (evmlr_voting_proof_enc_t**) malloc(N * sizeof(evmlr_voting_proof_enc_t*));
    for (slong i = 0; i < N; i++) {
        remaining_proofs[i] = (evmlr_voting_proof_enc_t*) malloc(k_servers * sizeof(evmlr_voting_proof_enc_t));
        for (slong l = 0; l < k_servers; l++) {
            slong L_curr_l = (k_servers - l - 1) * K_LWR * (K_LWE + 1) + 1;
            evmlr_enc_proof_init(&remaining_proofs[i][l]->proof, L_curr_l);
            evmlr_enc_proof_copy(&remaining_proofs[i][l]->proof, &subs[i]->proofs[l]->proof);
        }
    }

    slong L_curr = pp->L_max;
    for (slong i = 0; i < N; i++) {
        evmlr_voting_layer_init(&curr_cts[i], L_curr);
        for (int r = 0; r < K_LWR; r++) {
            nmod_poly_mat_set(curr_cts[i].enc_cipher[r]->uT, subs[i]->ballot->layer->enc_cipher[r]->uT);
            nmod_poly_set(curr_cts[i].enc_cipher[r]->v, subs[i]->ballot->layer->enc_cipher[r]->v);
        }
        nmod_poly_mat_set(curr_cts[i].otse_cipher->c, subs[i]->ballot->layer->otse_cipher->c);
    }

    nmod_poly_mat_t* decrypted_msgs = (nmod_poly_mat_t*) malloc(N * sizeof(nmod_poly_mat_t));
    nmod_poly_mat_t* shuffled_msgs = (nmod_poly_mat_t*) malloc(N * sizeof(nmod_poly_mat_t));

    uint64_t total_dec_verify_enc = 0;
    uint64_t total_shuffle_prove = 0;
    uint64_t total_shuffle_verify = 0;

    for (slong j = 0; j < k_servers; j++) {
        slong L_next = (k_servers - j - 1) * K_LWR * (K_LWE + 1) + 1;
        slong L_curr_j = (k_servers - j) * K_LWR * (K_LWE + 1) + 1;

        // Dynamically allocate and initialize curr_proofs for the current layer j
        evmlr_voting_proof_enc_t* curr_proofs = (evmlr_voting_proof_enc_t*) malloc(N * sizeof(evmlr_voting_proof_enc_t));
        for (slong i = 0; i < N; i++) {
            evmlr_enc_proof_init(&curr_proofs[i]->proof, L_curr_j - K_LWR * (K_LWE + 1));
            evmlr_enc_proof_copy(&curr_proofs[i]->proof, &remaining_proofs[i][j]->proof);
        }

        nmod_poly_mat_t* d_dagger = (nmod_poly_mat_t*) malloc(N * sizeof(nmod_poly_mat_t));
        nmod_poly_mat_t* a_mats = (nmod_poly_mat_t*) malloc(N * sizeof(nmod_poly_mat_t));

        // Decrypt and Verify Encryption Proofs
        uint64_t start_dec_verify = cpucycles();
        for (slong i = 0; i < N; i++) {
            voting_decrypt_layer(decrypted_msgs[i], d_dagger[i], a_mats[i], &curr_cts[i], sk->sks[j], pp);
            int ok = evmlr_voting_verify_enc(pp, curr_proofs[i], &curr_cts[i], decrypted_msgs[i], j);
            if (!ok) {
                printf("Warning: proof of encryption failed for voter %ld at server %ld\n", i + 1, j);
            }
        }
        uint64_t end_dec_verify = cpucycles();
        total_dec_verify_enc += (end_dec_verify - start_dec_verify);

        // Shuffle
        size_t* pi = (size_t *) malloc(N * sizeof(size_t));
        evmlr_utils_new_perm(pi, N);

        for (slong i = 0; i < N; i++) {
            nmod_poly_mat_init(shuffled_msgs[i], L_next, 1, MOD_Q);
            nmod_poly_mat_set(shuffled_msgs[i], decrypted_msgs[pi[i]]);
        }

        // Shuffle remaining proofs
        evmlr_voting_proof_enc_t** temp_remaining_proofs = (evmlr_voting_proof_enc_t**) malloc(N * sizeof(evmlr_voting_proof_enc_t*));
        for (slong i = 0; i < N; i++) {
            temp_remaining_proofs[i] = (evmlr_voting_proof_enc_t*) malloc(k_servers * sizeof(evmlr_voting_proof_enc_t));
            for (slong l = 0; l < k_servers; l++) {
                slong L_curr_l = (k_servers - l - 1) * K_LWR * (K_LWE + 1) + 1;
                evmlr_enc_proof_init(&temp_remaining_proofs[i][l]->proof, L_curr_l);
                evmlr_enc_proof_copy(&temp_remaining_proofs[i][l]->proof, &remaining_proofs[pi[i]][l]->proof);
            }
        }
        for (slong i = 0; i < N; i++) {
            for (slong l = 0; l < k_servers; l++) {
                evmlr_enc_proof_clear(&remaining_proofs[i][l]->proof);
            }
            free(remaining_proofs[i]);
            remaining_proofs[i] = temp_remaining_proofs[i];
        }
        free(temp_remaining_proofs);

        // Generate Proof of Shuffle
        uint64_t start_shuf_prove = cpucycles();
        evmlr_shuffle_sp_t shuf_sp;
        evmlr_shuffle_pp_t shuf_pp;

        shuf_sp->pi = pi;
        shuf_sp->d_dagger = d_dagger;
        nmod_poly_mat_init(shuf_sp->r_D, 2*K_SIS, 1, MOD_Q);
        evmlr_commit_sample_r(shuf_sp->r_D);

        nmod_poly_mat_t d_flat;
        nmod_poly_mat_init(d_flat, pp->shuf_ctx->com_ctx->N, 1, MOD_Q);
        nmod_poly_mat_zero(d_flat);
        slong row_len = (2*ETA) * (K_LWE + pp->L_max);
        for (slong i = 0; i < N; i++) {
            for (slong r = 0; r < row_len; r++) {
                nmod_poly_struct* d_dag_ir = nmod_poly_mat_entry(shuf_sp->d_dagger[i], r, 0);
                nmod_poly_struct* poly = nmod_poly_mat_entry(d_flat, i * row_len + r, 0);
                nmod_poly_set(poly, d_dag_ir);
            }
        }
        evmlr_commit(shuf_pp->D, d_flat, shuf_sp->r_D, pp->shuf_ctx->com_ctx);
        nmod_poly_mat_clear(d_flat);

        shuf_pp->c_hat = (nmod_poly_mat_t *) malloc(N * sizeof(nmod_poly_mat_t));
        for (slong i = 0; i < N; i++) {
            nmod_poly_mat_init(shuf_pp->c_hat[i], pp->L_max, 1, MOD_Q);
            nmod_poly_mat_zero(shuf_pp->c_hat[i]);
            for (slong r = 0; r < L_next; r++) {
                nmod_poly_set(nmod_poly_mat_entry(shuf_pp->c_hat[i], r, 0), nmod_poly_mat_entry(shuffled_msgs[i], r, 0));
            }
        }
        for (slong i = 0; i < N; i++) {
            for (slong r = L_next; r < pp->L_max; r++) {
                nmod_poly_struct* hat_val = nmod_poly_mat_entry(shuf_pp->c_hat[i], r, 0);
                nmod_poly_struct* a_val = nmod_poly_mat_entry(a_mats[pi[i]], r, 0);
                nmod_poly_neg(hat_val, a_val);
            }
        }

        shuf_pp->c_star = (evmlr_otse_ciphertext_t *) malloc(N * sizeof(evmlr_otse_ciphertext_t));
        for (slong i = 0; i < N; i++) {
            nmod_poly_mat_init(shuf_pp->c_star[i]->c, pp->L_max, 1, MOD_Q);
            nmod_poly_mat_zero(shuf_pp->c_star[i]->c);
            for (slong r = 0; r < curr_cts[i].otse_cipher->c->r; r++) {
                nmod_poly_set(nmod_poly_mat_entry(shuf_pp->c_star[i]->c, r, 0), nmod_poly_mat_entry(curr_cts[i].otse_cipher->c, r, 0));
            }
        }

        evmlr_proof_init(shuffle_proofs[j], pp->shuf_ctx);
        evmlr_shuffle_prove(shuffle_proofs[j], shuf_sp, shuf_pp, pp->shuf_ctx, state);
        uint64_t end_shuf_prove = cpucycles();
        total_shuffle_prove += (end_shuf_prove - start_shuf_prove);

        // Verify Shuffle Proof
        uint64_t start_shuf_verify = cpucycles();
        int shuffle_ok = evmlr_shuffle_verify(shuffle_proofs[j], pp->shuf_ctx, shuf_pp);
        uint64_t end_shuf_verify = cpucycles();
        if (!shuffle_ok) {
            printf("Warning: Shuffle verification failed at server %ld!\n", j);
        }
        total_shuffle_verify += (end_shuf_verify - start_shuf_verify);

        // Clean up pp/sp
        for (slong i = 0; i < N; i++) {
            nmod_poly_mat_clear(shuf_pp->c_hat[i]);
            nmod_poly_mat_clear(shuf_pp->c_star[i]->c);
            nmod_poly_mat_clear(a_mats[i]);
        }
        free(shuf_pp->c_hat);
        free(shuf_pp->c_star);
        free(a_mats);
        evmlr_commit_clear(shuf_pp->D);
        nmod_poly_mat_clear(shuf_sp->r_D);
        free(shuf_sp->pi);
        for (slong i = 0; i < N; i++) {
            nmod_poly_mat_clear(d_dagger[i]);
        }
        free(d_dagger);

        if (j == k_servers - 1) {
            for (slong i = 0; i < N; i++) {
                nmod_poly_mat_init(results[i], L_next, 1, MOD_Q);
                nmod_poly_mat_set(results[i], shuffled_msgs[i]);
            }
        } else {
            slong L_next_ct = L_next;
            for (slong i = 0; i < N; i++) {
                evmlr_voting_layer_clear(&curr_cts[i]);
                deserialize_layer(&curr_cts[i], shuffled_msgs[i], L_next_ct - K_LWR * (K_LWE + 1));
            }
        }

        for (slong i = 0; i < N; i++) {
            evmlr_enc_proof_clear(&curr_proofs[i]->proof);
        }
        free(curr_proofs);

        for (slong i = 0; i < N; i++) {
            nmod_poly_mat_clear(decrypted_msgs[i]);
            nmod_poly_mat_clear(shuffled_msgs[i]);
        }
    }

    printf("Peel & Verify Enc Proofs Phase (total for %ld servers): %lu cycles (avg %lu per server)\n",
           k_servers, (unsigned long)total_dec_verify_enc, (unsigned long)(total_dec_verify_enc / k_servers));
    printf("Shuffle Prove Phase (total for %ld servers): %lu cycles (avg %lu per server)\n",
           k_servers, (unsigned long)total_shuffle_prove, (unsigned long)(total_shuffle_prove / k_servers));
    printf("Shuffle Verify Phase (total for %ld servers): %lu cycles (avg %lu per server)\n",
           k_servers, (unsigned long)total_shuffle_verify, (unsigned long)(total_shuffle_verify / k_servers));

    // Calculate total cost and the percentage impact of each phase
    uint64_t total_cycles = (end_setup - start_setup) + (end_casting - start_casting) +
                            total_dec_verify_enc + total_shuffle_prove + total_shuffle_verify;
    printf("\n--- Breakdown of Phases (%ld voters, %ld servers) ---\n", N, k_servers);
    printf("Total protocol execution cost: %lu cycles\n", (unsigned long)total_cycles);
    printf("  1. Setup Phase:                %5.2f%% (%lu cycles)\n",
           100.0 * (double)(end_setup - start_setup) / (double)total_cycles, (unsigned long)(end_setup - start_setup));
    printf("  2. Casting Phase:              %5.2f%% (%lu cycles)\n",
           100.0 * (double)(end_casting - start_casting) / (double)total_cycles, (unsigned long)(end_casting - start_casting));
    printf("  3. Decrypt & Verify ZK Enc:    %5.2f%% (%lu cycles)\n",
           100.0 * (double)total_dec_verify_enc / (double)total_cycles, (unsigned long)total_dec_verify_enc);
    printf("  4. Shuffle Proving:            %5.2f%% (%lu cycles)\n",
           100.0 * (double)total_shuffle_prove / (double)total_cycles, (unsigned long)total_shuffle_prove);
    printf("  5. Shuffle Verification:       %5.2f%% (%lu cycles)\n",
           100.0 * (double)total_shuffle_verify / (double)total_cycles, (unsigned long)total_shuffle_verify);
    printf("=======================================================================\n\n");

    // Clean up
    free(decrypted_msgs);
    free(shuffled_msgs);
    for (slong i = 0; i < N; i++) {
        evmlr_voting_layer_clear(&curr_cts[i]);
        for (slong l = 0; l < k_servers; l++) {
            evmlr_enc_proof_clear(&remaining_proofs[i][l]->proof);
        }
        free(remaining_proofs[i]);
        nmod_poly_mat_clear(votes[i]);
        nmod_poly_mat_clear(results[i]);
        evmlr_voting_submission_clear(subs[i], pp);
    }
    free(curr_cts);
    free(remaining_proofs);
    for (int i = 0; i < k_servers; i++) {
        evmlr_proof_clear(shuffle_proofs[i], pp->shuf_ctx);
    }
    evmlr_voting_pp_clear(pp);
    evmlr_voting_sk_clear(sk);
}

int main() {
    flint_rand_t state;
    flint_rand_init(state);
    ulong seed[2];
    getrandom(seed, sizeof(ulong)*2, 0);
    flint_rand_set_seed(state, seed[0], seed[1]);

    test_voting_protocol(state);
    bench_voting_protocol(state);

    flint_rand_clear(state);
    return 0;
}

#endif
