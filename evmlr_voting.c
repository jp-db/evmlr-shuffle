#include "evmlr_voting.h"
#include "evmlr_utils.h"
#include <string.h>

#ifdef MAIN
#include "test.h"
#include "bench.h"
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

static void voting_encrypt_layer(evmlr_voting_layer_t ct, const nmod_poly_mat_t msg, const evmlr_mlpke_pk_t pk, const evmlr_voting_pp_t pp, flint_rand_t state) {
    slong L_in = msg->r;

    // Pad msg to L_max
    nmod_poly_mat_t msg_padded;
    nmod_poly_mat_init(msg_padded, pp->L_max, 1, MOD_Q);
    nmod_poly_mat_zero(msg_padded);
    for (slong i = 0; i < L_in; i++) {
        nmod_poly_set(nmod_poly_mat_entry(msg_padded, i, 0), nmod_poly_mat_entry(msg, i, 0));
    }

    // Encrypt using HPKE context
    evmlr_hpke_cipher_t cipher;
    evmlr_hpke_encrypt(cipher, NULL, msg_padded, pk, pp->hpke_ctx, state);

    // Copy to ct
    evmlr_voting_layer_clear(ct);
    evmlr_voting_layer_init(ct, L_in);
    for (int i = 0; i < K_LWR; i++) {
        nmod_poly_mat_set(ct->enc_cipher[i]->uT, cipher->enc_cipher[i]->uT);
        nmod_poly_set(ct->enc_cipher[i]->v, cipher->enc_cipher[i]->v);
    }
    for (slong i = 0; i < L_in; i++) {
        nmod_poly_set(nmod_poly_mat_entry(ct->otse_cipher->c, i, 0), nmod_poly_mat_entry(cipher->otse_cipher->c, i, 0));
    }

    nmod_poly_mat_clear(msg_padded);
    evmlr_hpke_cipher_clear(cipher);
}

// Helper to decrypt a layer
extern void calc_a(nmod_poly_mat_t a, nmod_poly_mat_t d_dagger, const evmlr_otse_key_t key, const evmlr_otse_ctx_t ctx);

static void voting_decrypt_layer(nmod_poly_mat_t msg, nmod_poly_mat_t d_dagger, nmod_poly_mat_t a, const evmlr_voting_layer_t ct, const evmlr_mlpke_sk_t sk, const evmlr_voting_pp_t pp) {
    slong L_in = ct->otse_cipher->c->r;

    // Reconstruct key s
    evmlr_otse_key_t key;
    nmod_poly_t key_poly;
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

// Compute deterministic proof hash
static void hash_layer_and_plaintext(uint8_t hash[SHA256HashSize], const evmlr_voting_layer_t ct, const nmod_poly_mat_t plaintext) {
    SHA256Context sha;
    SHA256Reset(&sha);

    // Hash ML-PKE ciphertexts
    for (int i = 0; i < K_LWR; i++) {
        // Hash uT
        for (int j = 0; j < K_LWE; j++) {
            nmod_poly_struct* poly = nmod_poly_mat_entry(ct->enc_cipher[i]->uT, 0, j);
            for (slong c = 0; c <= nmod_poly_degree(poly); c++) {
                ulong coeff = nmod_poly_get_coeff_ui(poly, c);
                SHA256Input(&sha, (uint8_t*)&coeff, sizeof(ulong));
            }
        }
        // Hash v
        const nmod_poly_struct* poly = ct->enc_cipher[i]->v;
        for (slong c = 0; c <= nmod_poly_degree(poly); c++) {
            ulong coeff = nmod_poly_get_coeff_ui(poly, c);
            SHA256Input(&sha, (uint8_t*)&coeff, sizeof(ulong));
        }
    }

    // Hash OTSE ciphertext c
    slong L_ct = ct->otse_cipher->c->r;
    for (slong i = 0; i < L_ct; i++) {
        nmod_poly_struct* poly = nmod_poly_mat_entry(ct->otse_cipher->c, i, 0);
        for (slong c = 0; c <= nmod_poly_degree(poly); c++) {
            ulong coeff = nmod_poly_get_coeff_ui(poly, c);
            SHA256Input(&sha, (uint8_t*)&coeff, sizeof(ulong));
        }
    }

    // Hash plaintext
    slong L_pt = plaintext->r;
    for (slong i = 0; i < L_pt; i++) {
        nmod_poly_struct* poly = nmod_poly_mat_entry(plaintext, i, 0);
        for (slong c = 0; c <= nmod_poly_degree(poly); c++) {
            ulong coeff = nmod_poly_get_coeff_ui(poly, c);
            SHA256Input(&sha, (uint8_t*)&coeff, sizeof(ulong));
        }
    }

    SHA256Result(&sha, hash);
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

            // Encrypt current_msg under pks[j]
            voting_encrypt_layer(temp_layer, current_msg, pp->pks[j], pp, state);

            // Generate proof of correct encryption
            hash_layer_and_plaintext(submission->proofs[j]->hash, temp_layer, current_msg);

            // Next plaintext is the serialized ciphertext
            nmod_poly_mat_clear(current_msg);
            slong L_curr = (k - j) * K_LWR * (K_LWE + 1) + 1;
            nmod_poly_mat_init(current_msg, L_curr, 1, MOD_Q);
            serialize_layer(current_msg, temp_layer, L_curr - K_LWR * (K_LWE + 1));

            evmlr_voting_layer_clear(temp_layer);
        } else {
            // Encrypt current_msg under pks[0] to the outermost layer
            voting_encrypt_layer(submission->ballot->layer, current_msg, pp->pks[0], pp, state);

            // Generate proof of correct encryption
            hash_layer_and_plaintext(submission->proofs[0]->hash, submission->ballot->layer, current_msg);
        }
    }

    nmod_poly_mat_clear(current_msg);
}

// verify_proof_of_encryption(pp, proof, ciphertext_layer, decrypted_plaintext, layer_idx)
int evmlr_voting_verify_enc(const evmlr_voting_pp_t pp, const evmlr_voting_proof_enc_t proof, const evmlr_voting_layer_t ct, const nmod_poly_mat_t decrypted_plaintext, slong layer_idx) {
    uint8_t computed_hash[SHA256HashSize];
    hash_layer_and_plaintext(computed_hash, ct, decrypted_plaintext);
    return memcmp(computed_hash, proof->hash, SHA256HashSize) == 0;
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
    evmlr_voting_proof_enc_struct* curr_proofs = (evmlr_voting_proof_enc_struct*) malloc(N * sizeof(evmlr_voting_proof_enc_struct));

    // Array of remaining proofs for each voter
    evmlr_voting_proof_enc_t** remaining_proofs = (evmlr_voting_proof_enc_t**) malloc(N * sizeof(evmlr_voting_proof_enc_t*));
    for (slong i = 0; i < N; i++) {
        remaining_proofs[i] = (evmlr_voting_proof_enc_t*) malloc(k * sizeof(evmlr_voting_proof_enc_t));
        for (slong l = 0; l < k; l++) {
            memcpy(remaining_proofs[i][l]->hash, active_subs[i]->proofs[l]->hash, SHA256HashSize);
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

        memcpy(curr_proofs[i].hash, active_subs[i]->proofs[0]->hash, SHA256HashSize);
    }

    nmod_poly_mat_t* decrypted_msgs = (nmod_poly_mat_t*) malloc(N * sizeof(nmod_poly_mat_t));
    nmod_poly_mat_t* shuffled_msgs = (nmod_poly_mat_t*) malloc(N * sizeof(nmod_poly_mat_t));

    // Process layer by layer
    for (slong j = 0; j < k; j++) {
        slong L_next = (k - j - 1) * K_LWR * (K_LWE + 1) + 1;

        nmod_poly_mat_t* d_dagger = (nmod_poly_mat_t*) malloc(N * sizeof(nmod_poly_mat_t));
        nmod_poly_mat_t* a_mats = (nmod_poly_mat_t*) malloc(N * sizeof(nmod_poly_mat_t));

        // Decrypt and Verify
        for (slong i = 0; i < N; i++) {
            voting_decrypt_layer(decrypted_msgs[i], d_dagger[i], a_mats[i], &curr_cts[i], sk->sks[j], pp);
            int ok = evmlr_voting_verify_enc(pp, &curr_proofs[i], &curr_cts[i], decrypted_msgs[i], j);
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
                memcpy(temp_remaining_proofs[i][l]->hash, remaining_proofs[pi[i]][l]->hash, SHA256HashSize);
            }
        }
        for (slong i = 0; i < N; i++) {
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
                memcpy(curr_proofs[i].hash, remaining_proofs[i][j + 1]->hash, SHA256HashSize);
            }
        }

        for (slong i = 0; i < N; i++) {
            nmod_poly_mat_clear(decrypted_msgs[i]);
            nmod_poly_mat_clear(shuffled_msgs[i]);
        }
    }

    free(decrypted_msgs);
    free(shuffled_msgs);
    for (slong i = 0; i < N; i++) {
        evmlr_voting_layer_clear(&curr_cts[i]);
        free(remaining_proofs[i]);
    }
    free(curr_cts);
    free(curr_proofs);
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
}

void evmlr_voting_submission_clear(evmlr_voting_submission_t sub, const evmlr_voting_pp_t pp) {
    evmlr_voting_layer_clear(sub->ballot->layer);
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

int main() {
    flint_rand_t state;
    flint_rand_init(state);
    ulong seed[2];
    getrandom(seed, sizeof(ulong)*2, 0);
    flint_rand_set_seed(state, seed[0], seed[1]);

    test_voting_protocol(state);

    flint_rand_clear(state);
    return 0;
}

#endif
