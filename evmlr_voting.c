#include "evmlr_voting.h"
#include "evmlr_utils.h"
#include "evmlr_challenge.h"
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

/* Encrypts msg and proves it. The proof claims `claimed` as the plaintext,
 * which is msg except in the tests of a dishonest voter. */
static void voting_encrypt_and_prove_layer(evmlr_voting_layer_t ct, evmlr_enc_proof_t proof, const nmod_poly_mat_t msg,
                                           const nmod_poly_mat_t claimed, const evmlr_mlpke_pk_t pk, const evmlr_voting_pp_t pp,
                                           const uint8_t* label, size_t label_len, flint_rand_t state) {
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
    evmlr_hpke_cipher_init(cipher, pp->hpke_ctx->otse_ctx->L);
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
    evmlr_enc_proof_prove(proof, pk, pp->hpke_ctx->otse_ctx, r_mats, e2_mats, e3_polys, key, cipher, d_dagger, claimed, L_in, pp->L_max,
                          label, label_len);

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

#define HDR (K_LWR * (K_LWE + 1))

slong evmlr_voting_otse_len(slong k, slong layer) {
    return (k - layer - 1) * HDR + 1;
}

slong evmlr_voting_msg_len(slong k, slong layer) {
    return (k - layer) * HDR + 1;
}

slong evmlr_voting_shuffle_len(slong N) {
    return N == 1 ? 2 : N;
}

static void absorb_u64(evmlr_challenge_t chal, uint64_t v) {
    uint8_t b[8];
    for (int i = 0; i < 8; i++) {
        b[i] = (uint8_t) (v >> (8 * i));
    }
    evmlr_challenge_add_bytes(chal, b, sizeof(b));
}

/* The Fiat-Shamir label of the proof of encryption of layer `layer`: the
 * election, the number of servers and the layer, and for the outermost layer
 * the voter. Without it a proof says nothing about who cast it: anyone could
 * resubmit a voter's ballot under their own id and have the vote counted
 * twice, or replay it in another election.
 *
 * Inner layers are not bound to the voter. Copying needs the outer layer
 * (the inner ones are hidden inside it until server 0 decrypts), and the
 * mixing servers would otherwise have to carry each item's voter id past
 * the shuffles.
 *
 * The outermost label also hashes the proofs of the inner layers (`inner`,
 * the k - 1 proofs of layers 1..k-1). With the binding inside every proof of
 * encryption (evmlr_enc_proof.h), the public check of the outermost layer
 * then covers the whole ballot: its voter id, ciphertext and every proof.
 * Inner proofs cannot be checked before the servers decrypt, but they can no
 * longer be changed without the board noticing. */
static void ballot_label(uint8_t out[SHA256HashSize], const evmlr_voting_pp_struct* pp, slong layer, int64_t voter_id,
                         const evmlr_enc_proof_struct* const* inner) {
    static const char domain[] = "evmlr ballot proof v1";
    evmlr_challenge_t chal;

    evmlr_challenge_init(chal);
    evmlr_challenge_add_bytes(chal, (const uint8_t *) domain, sizeof(domain) - 1);
    evmlr_challenge_add_bytes(chal, pp->election_id, EVMLR_ELECTION_ID_LEN);
    absorb_u64(chal, (uint64_t) pp->k);
    absorb_u64(chal, (uint64_t) layer);
    if (layer == 0) {
        absorb_u64(chal, 1);
        absorb_u64(chal, (uint64_t) voter_id);
        for (slong l = 1; l < pp->k; l++) {
            uint8_t d[SHA256HashSize];
            evmlr_enc_proof_digest(d, inner[l - 1]);
            evmlr_challenge_add_bytes(chal, d, sizeof(d));
        }
    } else {
        absorb_u64(chal, 0);
    }
    evmlr_challenge_get_hash(out, chal);
}

/* The inner proofs of a ballot (layers 1..k-1), as ballot_label takes them.
 * Returns NULL if out of memory; free it. */
static const evmlr_enc_proof_struct** inner_of_submission(const evmlr_voting_submission_struct* sub, slong k) {
    const evmlr_enc_proof_struct** p = malloc((k > 1 ? k - 1 : 1) * sizeof(*p));
    for (slong l = 1; p != NULL && l < k; l++) {
        p[l - 1] = &sub->proofs[l]->proof;
    }
    return p;
}

/* The same for item i of a layer-0 batch, into p (room for k - 1). */
static const evmlr_enc_proof_struct** inner_of_batch(const evmlr_enc_proof_struct** p, const evmlr_voting_batch_struct* b,
                                                     slong i) {
    for (slong l = 1; l < b->k; l++) {
        p[l - 1] = &b->proofs[i * b->k + l];
    }
    return p;
}

/* The OTSE matrices H and H'' used to be sampled by whoever ran setup, the
 * same party that generated every key. Every party now expands them from a
 * hash of the election id and k, so they can be recomputed by anyone and
 * chosen by nobody. */
static void derive_hpke_ctx(evmlr_hpke_ctx_t ctx, const uint8_t id[EVMLR_ELECTION_ID_LEN], slong k, slong L) {
    static const char domain[] = "evmlr otse matrices v1";
    uint8_t seed[SHA256HashSize];
    evmlr_challenge_t chal;
    nmod_poly_mat_t Hpp;
    prg_t prg;

    evmlr_challenge_init(chal);
    evmlr_challenge_add_bytes(chal, (const uint8_t *) domain, sizeof(domain) - 1);
    evmlr_challenge_add_bytes(chal, id, EVMLR_ELECTION_ID_LEN);
    absorb_u64(chal, (uint64_t) k);
    absorb_u64(chal, (uint64_t) L);
    evmlr_challenge_get_hash(seed, chal);

    evmlr_hpke_ctx_init_empty(ctx, L);
    prg_setseed(&prg, seed);
    evmlr_expand_uniform_mat(ctx->otse_ctx->H, &prg);
    nmod_poly_mat_init(Hpp, L, K_LWE, MOD_Q);
    evmlr_expand_uniform_mat(Hpp, &prg);
    for (slong r = 0; r < L; r++) {
        for (slong c = 0; c < K_LWE; c++) {
            nmod_poly_set(nmod_poly_mat_entry(ctx->otse_ctx->H_prime, r, c), nmod_poly_mat_entry(Hpp, r, c));
        }
    }
    nmod_poly_mat_clear(Hpp);
    memset(&prg, 0, sizeof(prg));
}

/* The shuffle's commitment key used to be sampled by whoever ran setup or
 * counting. The counting servers are the provers of the shuffle, and a
 * prover who picks the commitment key can pick one with a trapdoor, which
 * breaks binding and with it the soundness of the shuffle. The key is now
 * expanded from a hash of the public parameters and N: nobody chooses it,
 * and every verifier recomputes it. */
void evmlr_voting_shuffle_ctx_init(evmlr_shuffle_ctx_t ctx, const evmlr_voting_pp_t pp, slong N) {
    static const char domain[] = "evmlr-shuffle commitment key v1";
    uint8_t seed[SHA256HashSize];
    evmlr_challenge_t chal;

    evmlr_challenge_init(chal);
    evmlr_challenge_add_bytes(chal, (const uint8_t *) domain, sizeof(domain) - 1);
    absorb_u64(chal, (uint64_t) pp->k);
    absorb_u64(chal, (uint64_t) N);
    evmlr_challenge_add_matrix(chal, pp->hpke_ctx->otse_ctx->H);
    evmlr_challenge_add_matrix(chal, pp->hpke_ctx->otse_ctx->H_prime);
    for (slong j = 0; j < pp->k; j++) {
        evmlr_challenge_add_matrix(chal, pp->pks[j]->A);
        evmlr_challenge_add_matrix(chal, pp->pks[j]->t);
    }
    evmlr_challenge_get_hash(seed, chal);
    evmlr_shuffle_ctx_init_seeded(ctx, N, pp->hpke_ctx, seed);
}

void evmlr_voting_keygen(evmlr_mlpke_keypair_t kp, flint_rand_t state) {
    evmlr_mlpke_ctx_t ctx;

    evmlr_mlpke_ctx_init(ctx);
    evmlr_mlpke_keypair_gen(kp, state, ctx);
    evmlr_mlpke_ctx_clear(ctx);
}

int evmlr_voting_pp_init(evmlr_voting_pp_t pp, const uint8_t election_id[EVMLR_ELECTION_ID_LEN],
                         slong k, const evmlr_mlpke_pk_t* pks) {
    if (k < 1) {
        return -1;
    }
    pp->k = k;
    // Outer message input has length k * K_LWR * (K_LWE + 1) + 1.
    // Decrypting layer 1 yields length (k - 1) * K_LWR * (K_LWE + 1) + 1.
    // This is the maximum length of the messages being shuffled.
    pp->L_max = evmlr_voting_otse_len(k, 0);
    memcpy(pp->election_id, election_id, EVMLR_ELECTION_ID_LEN);
    derive_hpke_ctx(pp->hpke_ctx, election_id, k, pp->L_max);

    pp->pks = (evmlr_mlpke_pk_t*) malloc(k * sizeof(evmlr_mlpke_pk_t));
    for (slong j = 0; j < k; j++) {
        nmod_poly_mat_init(pp->pks[j]->A, K_LWE, K_LWE, MOD_Q);
        nmod_poly_mat_init(pp->pks[j]->t, K_LWE, 1, MOD_Q);
        nmod_poly_mat_set(pp->pks[j]->A, pks[j]->A);
        nmod_poly_mat_set(pp->pks[j]->t, pks[j]->t);
    }
    return 0;
}

// setup(k, state) -> (pp, sk), with a random election id
void evmlr_voting_setup(evmlr_voting_pp_t pp, evmlr_voting_sk_t sk, slong k, flint_rand_t state) {
    uint8_t id[EVMLR_ELECTION_ID_LEN];
    evmlr_mlpke_keypair_t* kps = (evmlr_mlpke_keypair_t*) malloc(k * sizeof(evmlr_mlpke_keypair_t));
    evmlr_mlpke_pk_t* pks = (evmlr_mlpke_pk_t*) malloc(k * sizeof(evmlr_mlpke_pk_t));

    evmlr_utils_random(id, sizeof(id));
    sk->k = k;
    sk->sks = (evmlr_mlpke_sk_t*) malloc(k * sizeof(evmlr_mlpke_sk_t));
    for (slong j = 0; j < k; j++) {
        evmlr_voting_keygen(kps[j], state);
        pks[j][0] = kps[j]->pk[0]; // shallow: pp_init copies
        nmod_poly_mat_init(sk->sks[j]->s, K_LWE, 1, MOD_Q);
        nmod_poly_mat_set(sk->sks[j]->s, kps[j]->sk->s);
    }
    evmlr_voting_pp_init(pp, id, k, (const evmlr_mlpke_pk_t*) pks);
    for (slong j = 0; j < k; j++) {
        evmlr_mlpke_keypair_clear(kps[j]);
    }
    free(pks);
    free(kps);
}

/* A read-only HPKE view of a layer (the two structs have the same fields). */
static void shallow_cipher(evmlr_hpke_cipher_struct* c, const evmlr_voting_layer_struct* ct) {
    for (int i = 0; i < K_LWR; i++) {
        c->enc_cipher[i][0] = ct->enc_cipher[i][0];
    }
    c->otse_cipher[0] = ct->otse_cipher[0];
}

// casting(submission, voter_id, vote, pp, state)
//
// Layers are made from the inside out, so the outermost proof can hash the
// inner ones. Every layer is checked here, against the plaintext only the
// voter knows, before the ballot is handed out: no one else can check the
// inner layers before the servers decrypt them.
//
// cast() takes a description of how to cheat, for the tests of what the
// servers reject: a dishonest voter can make any layer wrong and still bind
// the whole ballot. NULL for an honest ballot.
typedef struct {
    slong wrong_plaintext; // the proof of this layer claims another plaintext (-1: none)
    slong broken_hpke;     // the ML-PKE half of this layer's proof is changed once made (-1: none)
} cheat_t;

static int cast(evmlr_voting_submission_t submission, int64_t voter_id, const nmod_poly_mat_t vote,
                const evmlr_voting_pp_t pp, flint_rand_t state, const cheat_t* cheat) {
    const slong k = pp->k;
    uint8_t label[SHA256HashSize];
    const evmlr_enc_proof_struct** inner;
    evmlr_hpke_cipher_struct c;
    nmod_poly_mat_t current_msg;
    int ok = 1;

    submission->voter_id = voter_id;
    inner = inner_of_submission(submission, k);
    if (inner == NULL) {
        return 0;
    }
    nmod_poly_mat_init(current_msg, vote->r, vote->c, MOD_Q);
    nmod_poly_mat_set(current_msg, vote);

    for (slong j = k - 1; j >= 0; j--) {
        evmlr_enc_proof_struct* prf = &submission->proofs[j]->proof;
        nmod_poly_mat_t claimed;

        nmod_poly_mat_init_set(claimed, current_msg);
        if (cheat != NULL && cheat->wrong_plaintext == j) {
            nmod_poly_struct* m0 = nmod_poly_mat_entry(claimed, 0, 0);
            nmod_poly_set_coeff_ui(m0, 0, (nmod_poly_get_coeff_ui(m0, 0) + 1) % MOD_Q);
        }
        ballot_label(label, pp, j, voter_id, inner);
        if (j > 0) {
            evmlr_voting_layer_t temp_layer;
            slong L_otse = evmlr_voting_otse_len(k, j);
            evmlr_voting_layer_init(temp_layer, L_otse);

            // Encrypt current_msg under pks[j] and generate proof
            voting_encrypt_and_prove_layer(temp_layer, prf, current_msg, claimed, pp->pks[j], pp, label, sizeof(label), state);
            shallow_cipher(&c, temp_layer);
            ok = ok && evmlr_enc_proof_verify(prf, pp->pks[j], pp->hpke_ctx->otse_ctx, &c, current_msg, pp->L_max,
                                              label, sizeof(label));

            // Next plaintext is the serialized ciphertext
            nmod_poly_mat_clear(current_msg);
            slong L_curr = evmlr_voting_msg_len(k, j);
            nmod_poly_mat_init(current_msg, L_curr, 1, MOD_Q);
            serialize_layer(current_msg, temp_layer, L_otse);

            evmlr_voting_layer_clear(temp_layer);
        } else {
            // Encrypt current_msg under pks[0] to the outermost layer and generate proof
            voting_encrypt_and_prove_layer(submission->ballot->layer, prf, current_msg, claimed, pp->pks[0], pp,
                                           label, sizeof(label), state);
            shallow_cipher(&c, submission->ballot->layer);
            ok = ok && evmlr_enc_proof_verify(prf, pp->pks[0], pp->hpke_ctx->otse_ctx, &c, current_msg, pp->L_max,
                                              label, sizeof(label));
        }
        nmod_poly_mat_clear(claimed);
        if (cheat != NULL && cheat->broken_hpke == j) {
            nmod_poly_struct* z = nmod_poly_mat_entry(prf->proof_hpke.z_1, 0, 0);
            nmod_poly_set_coeff_ui(z, 0, (nmod_poly_get_coeff_ui(z, 0) + 1) % MOD_Q);
        }
    }

    nmod_poly_mat_clear(current_msg);
    free(inner);
    return ok;
}

int evmlr_voting_casting(evmlr_voting_submission_t submission, int64_t voter_id, const nmod_poly_mat_t vote,
                         const evmlr_voting_pp_t pp, flint_rand_t state) {
    return cast(submission, voter_id, vote, pp, state, NULL);
}

int evmlr_voting_check_ballot(const evmlr_voting_pp_t pp, const evmlr_voting_submission_t sub) {
    const evmlr_voting_layer_struct* ct = sub->ballot->layer;
    uint8_t label[SHA256HashSize];
    const evmlr_enc_proof_struct** inner;
    evmlr_hpke_cipher_struct c;
    int ok;

    if (ct->otse_cipher->c->r != pp->L_max || (inner = inner_of_submission(sub, pp->k)) == NULL) {
        return 0;
    }
    ballot_label(label, pp, 0, sub->voter_id, inner);
    shallow_cipher(&c, ct);
    ok = evmlr_enc_proof_verify_hpke(&sub->proofs[0]->proof, pp->pks[0], &c, label, sizeof(label));
    free(inner);
    return ok;
}

int evmlr_voting_verify_ballot(const evmlr_voting_pp_t pp, const evmlr_voting_submission_t sub,
                               const nmod_poly_mat_t plaintext) {
    uint8_t label[SHA256HashSize];
    const evmlr_enc_proof_struct** inner;
    evmlr_hpke_cipher_struct c;
    int ok;

    if (sub->ballot->layer->otse_cipher->c->r != pp->L_max || (inner = inner_of_submission(sub, pp->k)) == NULL) {
        return 0;
    }
    ballot_label(label, pp, 0, sub->voter_id, inner);
    shallow_cipher(&c, sub->ballot->layer);
    ok = evmlr_enc_proof_verify(&sub->proofs[0]->proof, pp->pks[0], pp->hpke_ctx->otse_ctx, &c,
                                plaintext, pp->L_max, label, sizeof(label));
    free(inner);
    return ok;
}

// Revote policy: of the ballots that pass evmlr_voting_check_ballot, the
// last one of each voter counts. Fills active (room for num_submissions) in
// chronological order and returns their number.
//
// Ballots that fail the check are skipped before the policy, not after: a
// ballot copied under Alice's id fails it, and must not replace her vote.
static slong active_submissions(evmlr_voting_submission_struct** active_subs,
                                const evmlr_voting_submission_t* submissions, slong num_submissions,
                                const evmlr_voting_pp_t pp) {
    slong N = 0;
    for (slong i = num_submissions - 1; i >= 0; i--) {
        int64_t voter_id = submissions[i]->voter_id;
        int found = 0;
        for (slong j = 0; j < N; j++) {
            if (active_subs[j]->voter_id == voter_id) {
                found = 1;
                break;
            }
        }
        if (!found && evmlr_voting_check_ballot(pp, submissions[i])) {
            active_subs[N++] = (evmlr_voting_submission_struct*) submissions[i];
        }
    }
    // Reverse active_subs to maintain chronological order
    for (slong i = 0; i < N / 2; i++) {
        evmlr_voting_submission_struct* tmp = active_subs[i];
        active_subs[i] = active_subs[N - 1 - i];
        active_subs[N - 1 - i] = tmp;
    }
    return N;
}

/* Whether rows [0, len) of got equal rows [offset, offset + len) of want and
 * the rest of got is zero. */
static int rows_match(const nmod_poly_mat_t got, const nmod_poly_mat_t want, slong offset, slong len) {
    if (nmod_poly_mat_ncols(got) != 1 || nmod_poly_mat_nrows(got) < len
        || nmod_poly_mat_nrows(want) < offset + len) {
        return 0;
    }
    for (slong r = 0; r < nmod_poly_mat_nrows(got); r++) {
        const nmod_poly_struct* g = nmod_poly_mat_entry(got, r, 0);
        if (r < len ? !nmod_poly_equal(g, nmod_poly_mat_entry(want, offset + r, 0)) : !nmod_poly_is_zero(g)) {
            return 0;
        }
    }
    return 1;
}

/* The dimensions of a shuffle of N messages of length L, without the
 * matrices: all the clear functions read. */
static void shape_ctx(evmlr_shuffle_ctx_struct* shape, slong N, slong L) {
    memset(shape, 0, sizeof(*shape));
    shape->N = N;
    shape->L = L;
    shape->com_ctx->N = evmlr_shuffle_com_len(N, L);
}

/* ---------------------------------------------------------------------------
 * Batches and mix steps
 * ------------------------------------------------------------------------- */

int evmlr_voting_batch_init(evmlr_voting_batch_t b, slong k, slong layer, slong n) {
    const slong nprf = k - layer;

    memset(b, 0, sizeof(evmlr_voting_batch_struct));
    b->k = k;
    b->layer = layer;
    if (k < 1 || layer < 0 || layer > k || n < 0) {
        return EVMLR_VOTING_ERR_ARG;
    }
    if (n == 0) {
        return EVMLR_VOTING_OK;
    }
    b->msgs = (nmod_poly_mat_t*) malloc(n * sizeof(nmod_poly_mat_t));
    b->proofs = nprf > 0 ? (evmlr_enc_proof_struct*) malloc(n * nprf * sizeof(evmlr_enc_proof_struct)) : NULL;
    b->voter_ids = layer == 0 ? (int64_t*) calloc(n, sizeof(int64_t)) : NULL;
    if (b->msgs == NULL || (nprf > 0 && b->proofs == NULL) || (layer == 0 && b->voter_ids == NULL)) {
        free(b->msgs);
        free(b->proofs);
        free(b->voter_ids);
        memset(b, 0, sizeof(evmlr_voting_batch_struct));
        return EVMLR_VOTING_ERR_NOMEM;
    }
    b->n = n;
    for (slong i = 0; i < n; i++) {
        nmod_poly_mat_init(b->msgs[i], evmlr_voting_msg_len(k, layer), 1, MOD_Q);
        for (slong l = 0; l < nprf; l++) {
            evmlr_enc_proof_init(&b->proofs[i * nprf + l], evmlr_voting_otse_len(k, layer + l));
        }
    }
    return EVMLR_VOTING_OK;
}

void evmlr_voting_batch_clear(evmlr_voting_batch_t b) {
    const slong nprf = b->k - b->layer;

    for (slong i = 0; i < b->n; i++) {
        nmod_poly_mat_clear(b->msgs[i]);
        for (slong l = 0; l < nprf; l++) {
            evmlr_enc_proof_clear(&b->proofs[i * nprf + l]);
        }
    }
    free(b->msgs);
    free(b->proofs);
    free(b->voter_ids);
    b->msgs = NULL;
    b->proofs = NULL;
    b->voter_ids = NULL;
    b->n = 0;
}

int evmlr_voting_batch_from_submissions(evmlr_voting_batch_t batch, const evmlr_voting_submission_t* submissions,
                                        slong num_submissions, const evmlr_voting_pp_t pp) {
    const slong k = pp->k;
    evmlr_voting_submission_struct** active;
    slong N;
    int st;

    if (num_submissions < 0) {
        return EVMLR_VOTING_ERR_ARG;
    }
    active = (evmlr_voting_submission_struct**) malloc((num_submissions ? num_submissions : 1) * sizeof(*active));
    if (active == NULL) {
        return EVMLR_VOTING_ERR_NOMEM;
    }
    N = active_submissions(active, submissions, num_submissions, pp);
    st = evmlr_voting_batch_init(batch, k, 0, N);
    for (slong i = 0; i < N && st == EVMLR_VOTING_OK; i++) {
        serialize_layer(batch->msgs[i], active[i]->ballot->layer, pp->L_max);
        batch->voter_ids[i] = active[i]->voter_id;
        for (slong l = 0; l < k; l++) {
            evmlr_enc_proof_copy(&batch->proofs[i * k + l], &active[i]->proofs[l]->proof);
        }
    }
    free(active);
    return st;
}

void evmlr_voting_mix_clear(evmlr_voting_mix_t mix) {
    if (mix->N > 0) {
        evmlr_shuffle_ctx_struct shape;

        shape_ctx(&shape, evmlr_voting_shuffle_len(mix->N), evmlr_voting_otse_len(mix->k, 0));
        evmlr_shuffle_pp_clear(mix->statement, &shape);
        evmlr_proof_clear(mix->proof, &shape);
    }
    free(mix->rejected);
    free(mix->reasons);
    mix->rejected = NULL;
    mix->reasons = NULL;
    mix->n_rejected = 0;
    mix->N = 0;
}

int evmlr_voting_mix_step(evmlr_voting_batch_t out, evmlr_voting_mix_t mix, const evmlr_voting_batch_t in,
                          const evmlr_mlpke_keypair_t kp, const evmlr_voting_pp_t pp, flint_rand_t state) {
    const slong k = pp->k, j = in->layer, n = in->n;
    const slong nprf = k - j, L_max = pp->L_max;
    const slong L_ct = evmlr_voting_otse_len(k, j); // OTSE length of the inputs = length of the outputs
    evmlr_voting_layer_struct* cts;
    nmod_poly_mat_t *dec, *dd, *a;
    slong *keep, N = 0, Ns = 0;
    const evmlr_enc_proof_struct** inner; // the inner proofs of an outer layer (j = 0)
    int st;

    if (in->k != k || j < 0 || j >= k || n < 0 || (j == 0 && n > 0 && in->voter_ids == NULL)) {
        return EVMLR_VOTING_ERR_ARG;
    }
    if (!nmod_poly_mat_equal(kp->pk->A, pp->pks[j]->A) || !nmod_poly_mat_equal(kp->pk->t, pp->pks[j]->t)) {
        return EVMLR_VOTING_ERR_KEY;
    }
    for (slong i = 0; i < n; i++) {
        if (nmod_poly_mat_nrows(in->msgs[i]) != evmlr_voting_msg_len(k, j) || nmod_poly_mat_ncols(in->msgs[i]) != 1) {
            return EVMLR_VOTING_ERR_ARG;
        }
    }

    memset(mix, 0, sizeof(evmlr_voting_mix_struct));
    mix->k = k;
    mix->j = j;
    mix->n_in = n;
    if (n == 0) {
        return evmlr_voting_batch_init(out, k, j + 1, 0);
    }

    // Room for one more decryption than inputs: a single one is shuffled twice.
    cts = (evmlr_voting_layer_struct*) malloc(n * sizeof(*cts));
    dec = (nmod_poly_mat_t*) malloc((n + 1) * sizeof(nmod_poly_mat_t));
    dd = (nmod_poly_mat_t*) malloc((n + 1) * sizeof(nmod_poly_mat_t));
    a = (nmod_poly_mat_t*) malloc((n + 1) * sizeof(nmod_poly_mat_t));
    keep = (slong*) malloc((n + 1) * sizeof(slong));
    inner = malloc(k * sizeof(*inner));
    mix->rejected = (slong*) malloc(n * sizeof(slong));
    mix->reasons = (uint8_t*) malloc(n);
    if (!cts || !dec || !dd || !a || !keep || !inner || !mix->rejected || !mix->reasons) {
        free(cts); free(dec); free(dd); free(a); free(keep); free(inner);
        free(mix->rejected); free(mix->reasons);
        return EVMLR_VOTING_ERR_NOMEM;
    }

    // Decrypt and check every proof of encryption; drop what fails.
    for (slong i = 0; i < n; i++) {
        const evmlr_enc_proof_struct* prf = &in->proofs[i * nprf];
        uint8_t label[SHA256HashSize];
        evmlr_hpke_cipher_struct c;

        deserialize_layer(&cts[i], in->msgs[i], L_ct);
        shallow_cipher(&c, &cts[i]);
        ballot_label(label, pp, j, j == 0 ? in->voter_ids[i] : 0, j == 0 ? inner_of_batch(inner, in, i) : NULL);
        if (!evmlr_enc_proof_verify_hpke(prf, pp->pks[j], &c, label, sizeof(label))) {
            mix->rejected[mix->n_rejected] = i;
            mix->reasons[mix->n_rejected++] = EVMLR_VOTING_REJECT_HPKE;
            continue;
        }
        // d_dagger is owned by the caller; decryption only fills it in.
        nmod_poly_mat_init(dd[N], 2 * ETA * (K_LWE + L_max), 1, MOD_Q);
        voting_decrypt_layer(dec[N], dd[N], a[N], &cts[i], kp->sk, pp);
        if (!evmlr_enc_proof_verify_otse(prf, pp->hpke_ctx->otse_ctx, &c, dec[N], label, sizeof(label))) {
            nmod_poly_mat_clear(dec[N]);
            nmod_poly_mat_clear(dd[N]);
            nmod_poly_mat_clear(a[N]);
            mix->rejected[mix->n_rejected] = i;
            mix->reasons[mix->n_rejected++] = EVMLR_VOTING_REJECT_OTSE;
            continue;
        }
        keep[N++] = i;
    }
    mix->N = N;
    if (mix->n_rejected == 0) {
        free(mix->rejected);
        free(mix->reasons);
        mix->rejected = NULL;
        mix->reasons = NULL;
    }

    /* The proof of shuffle needs two ciphertexts. One left is shuffled with
     * a copy of itself: the statement holds it twice (inputs and outputs),
     * the verifier checks that the two are equal, and only one is passed
     * on. That proves what a shuffle of two proves about each, and hides
     * nothing, as no shuffle of one could. */
    Ns = evmlr_voting_shuffle_len(N);
    if (N == 1) {
        nmod_poly_mat_init_set(dec[1], dec[0]);
        nmod_poly_mat_init_set(dd[1], dd[0]);
        nmod_poly_mat_init_set(a[1], a[0]);
        keep[1] = keep[0];
    }
    st = evmlr_voting_batch_init(out, k, j + 1, N);
    if (st == EVMLR_VOTING_OK && N > 0) {
        evmlr_shuffle_ctx_t sctx;
        evmlr_shuffle_sp_t sp;
        evmlr_shuffle_pp_struct* stmt = mix->statement;
        nmod_poly_mat_t d_flat;
        const slong row_len = (2 * ETA) * (K_LWE + L_max);

        evmlr_voting_shuffle_ctx_init(sctx, pp, Ns);

        // The witness: the permutation, the OTSE randomness and the opening of D
        sp->pi = (size_t*) malloc(Ns * sizeof(size_t));
        evmlr_utils_new_perm(sp->pi, Ns);
        sp->d_dagger = dd;
        nmod_poly_mat_init(sp->r_D, 2 * K_SIS, 1, MOD_Q);
        evmlr_commit_sample_r(sp->r_D);

        nmod_poly_mat_init(d_flat, sctx->com_ctx->N, 1, MOD_Q);
        nmod_poly_mat_zero(d_flat);
        for (slong m = 0; m < Ns; m++) {
            for (slong r = 0; r < row_len; r++) {
                nmod_poly_set(nmod_poly_mat_entry(d_flat, m * row_len + r, 0), nmod_poly_mat_entry(dd[m], r, 0));
            }
        }
        evmlr_commit_init(stmt->D);
        evmlr_commit(stmt->D, d_flat, sp->r_D, sctx->com_ctx);
        nmod_poly_mat_clear(d_flat);

        // Outputs: the decryptions in shuffled order, padded with -a
        stmt->c_hat = (nmod_poly_mat_t*) malloc(Ns * sizeof(nmod_poly_mat_t));
        stmt->c_star = (evmlr_otse_ciphertext_t*) malloc(Ns * sizeof(evmlr_otse_ciphertext_t));
        for (slong m = 0; m < Ns; m++) {
            const slong p = (slong) sp->pi[m];

            nmod_poly_mat_init(stmt->c_hat[m], L_max, 1, MOD_Q);
            for (slong r = 0; r < L_max; r++) {
                if (r < L_ct) {
                    nmod_poly_set(nmod_poly_mat_entry(stmt->c_hat[m], r, 0), nmod_poly_mat_entry(dec[p], r, 0));
                } else {
                    nmod_poly_neg(nmod_poly_mat_entry(stmt->c_hat[m], r, 0), nmod_poly_mat_entry(a[p], r, 0));
                }
            }
            // Inputs: the OTSE parts of the ciphertexts kept, in input order
            nmod_poly_mat_init(stmt->c_star[m]->c, L_max, 1, MOD_Q);
            for (slong r = 0; r < L_ct; r++) {
                nmod_poly_set(nmod_poly_mat_entry(stmt->c_star[m]->c, r, 0),
                              nmod_poly_mat_entry(cts[keep[m]].otse_cipher->c, r, 0));
            }
        }

        evmlr_proof_init(mix->proof, sctx);
        evmlr_shuffle_prove(mix->proof, sp, stmt, sctx, state);

        // The next batch: the same outputs, with the proofs of the layers left
        for (slong m = 0; m < N; m++) {
            const slong src = keep[sp->pi[m]];

            nmod_poly_mat_set(out->msgs[m], dec[sp->pi[m]]);
            for (slong l = 0; l + 1 < nprf; l++) {
                evmlr_enc_proof_copy(&out->proofs[m * (nprf - 1) + l], &in->proofs[src * nprf + l + 1]);
            }
        }

        nmod_poly_mat_clear(sp->r_D);
        free(sp->pi);
        evmlr_shuffle_ctx_clear(sctx);
    } else if (st != EVMLR_VOTING_OK) {
        mix->N = 0; // nothing proven to free
    }

    for (slong m = 0; m < Ns; m++) {
        nmod_poly_mat_clear(dec[m]);
        nmod_poly_mat_clear(dd[m]);
        nmod_poly_mat_clear(a[m]);
    }
    for (slong i = 0; i < n; i++) {
        evmlr_voting_layer_clear(&cts[i]);
    }
    free(cts); free(dec); free(dd); free(a); free(keep); free(inner);
    if (st != EVMLR_VOTING_OK) {
        evmlr_voting_mix_clear(mix);
    }
    return st;
}

int evmlr_voting_batch_matches_mix(const evmlr_voting_batch_t batch, const evmlr_voting_mix_t mix,
                                   const evmlr_voting_pp_t pp) {
    if (batch->k != pp->k || mix->k != pp->k || batch->layer != mix->j + 1 || batch->n != mix->N) {
        return 0;
    }
    for (slong m = 0; m < batch->n; m++) {
        if (!rows_match(batch->msgs[m], mix->statement->c_hat[m], 0, evmlr_voting_msg_len(pp->k, batch->layer))) {
            return 0;
        }
    }
    return 1;
}

/* ---------------------------------------------------------------------------
 * Counting and its verification
 * ------------------------------------------------------------------------- */

slong evmlr_voting_counting_transcript(nmod_poly_mat_t* results, evmlr_voting_transcript_t tr,
                                       const evmlr_voting_submission_t* submissions, slong num_submissions,
                                       const evmlr_voting_pp_t pp, const evmlr_voting_sk_t sk, flint_rand_t state) {
    const slong k = pp->k;
    evmlr_voting_batch_struct cur, next;
    slong N;

    tr->k = k;
    tr->mixes = (evmlr_voting_mix_struct*) calloc(k, sizeof(evmlr_voting_mix_struct));
    if (tr->mixes == NULL || sk->k != k
        || evmlr_voting_batch_from_submissions(&cur, submissions, num_submissions, pp) != EVMLR_VOTING_OK) {
        free(tr->mixes);
        tr->mixes = NULL;
        tr->k = 0;
        return -1;
    }
    for (slong j = 0; j < k; j++) {
        evmlr_mlpke_keypair_struct kp; // server j's key, borrowed

        kp.sk[0] = sk->sks[j][0];
        kp.pk[0] = pp->pks[j][0];
        if (evmlr_voting_mix_step(&next, &tr->mixes[j], &cur, &kp, pp, state) != EVMLR_VOTING_OK) {
            evmlr_voting_batch_clear(&cur);
            for (slong l = 0; l < j; l++) {
                evmlr_voting_mix_clear(&tr->mixes[l]);
            }
            free(tr->mixes);
            tr->mixes = NULL;
            tr->k = 0;
            return -1;
        }
        evmlr_voting_batch_clear(&cur);
        cur = next;
    }
    N = cur.n;
    for (slong i = 0; i < N; i++) {
        nmod_poly_mat_init(results[i], 1, 1, MOD_Q);
        nmod_poly_mat_set(results[i], cur.msgs[i]);
    }
    evmlr_voting_batch_clear(&cur);
    return N;
}

void evmlr_voting_transcript_clear(evmlr_voting_transcript_t tr) {
    for (slong j = 0; tr->mixes != NULL && j < tr->k; j++) {
        evmlr_voting_mix_clear(&tr->mixes[j]);
    }
    free(tr->mixes);
    tr->mixes = NULL;
    tr->k = 0;
}

slong evmlr_voting_verify_mixes(const evmlr_voting_mix_struct* mixes, slong m,
                                const evmlr_voting_submission_t* submissions, slong num_submissions,
                                const evmlr_voting_pp_t pp) {
    const slong k = pp->k;
    evmlr_voting_submission_struct** active;
    evmlr_shuffle_ctx_t ctx;
    slong ctx_N = 0, expected;
    int ok = 1;

    if (m < 0 || m > k || (m > 0 && mixes == NULL) || num_submissions < 0) {
        return -1;
    }
    active = (evmlr_voting_submission_struct**) malloc((num_submissions ? num_submissions : 1) * sizeof(*active));
    if (active == NULL) {
        return -1;
    }
    expected = active_submissions(active, submissions, num_submissions, pp);

    for (slong j = 0; j < m && ok; j++) {
        // Server j's inputs: the OTSE parts of the ballots' outer layers, then
        // those of the ciphertexts the previous server output.
        const evmlr_voting_mix_struct* mx = &mixes[j];
        const slong L_in = evmlr_voting_otse_len(k, j);
        slong m = 0, r = 0;

        ok = mx->k == k && mx->j == j && mx->n_in == expected && mx->n_rejected >= 0
          && mx->N == mx->n_in - mx->n_rejected && mx->N >= 0
          && (mx->n_rejected == 0 || (mx->rejected != NULL && mx->reasons != NULL));
        for (slong i = 0; ok && i < mx->n_in; i++) {
            if (r < mx->n_rejected && mx->rejected[r] == i) {
                const uint8_t why = mx->reasons[r++];

                ok = why == EVMLR_VOTING_REJECT_HPKE || why == EVMLR_VOTING_REJECT_OTSE;
                // At server 0 this rejection is checkable: the ballot is public.
                if (ok && j == 0 && why == EVMLR_VOTING_REJECT_HPKE) {
                    ok = !evmlr_voting_check_ballot(pp, active[i]);
                }
                continue;
            }
            if (j == 0) {
                ok = rows_match(mx->statement->c_star[m]->c, active[i]->ballot->layer->otse_cipher->c, 0, L_in);
            } else {
                ok = rows_match(mx->statement->c_star[m]->c, mixes[j - 1].statement->c_hat[i], HDR, L_in);
            }
            m++;
        }
        // Every rejected index was met, so they are increasing and in range.
        ok = ok && r == mx->n_rejected;
        // A single ciphertext is shuffled with its copy (see the mix step).
        if (ok && mx->N == 1) {
            ok = nmod_poly_mat_equal(mx->statement->c_star[1]->c, mx->statement->c_star[0]->c)
              && nmod_poly_mat_equal(mx->statement->c_hat[1], mx->statement->c_hat[0]);
        }
        if (ok && mx->N > 0) {
            const slong Ns = evmlr_voting_shuffle_len(mx->N);
            if (ctx_N != Ns) {
                if (ctx_N > 0) {
                    evmlr_shuffle_ctx_clear(ctx);
                }
                evmlr_voting_shuffle_ctx_init(ctx, pp, Ns);
                ctx_N = Ns;
            }
            ok = evmlr_shuffle_verify(mx->proof, ctx, mx->statement);
        }
        expected = mx->N;
    }
    if (ctx_N > 0) {
        evmlr_shuffle_ctx_clear(ctx);
    }
    free(active);
    return ok ? expected : -1;
}

slong evmlr_voting_verify_counting(nmod_poly_mat_t* results, const evmlr_voting_transcript_t tr,
                                   const evmlr_voting_submission_t* submissions, slong num_submissions,
                                   const evmlr_voting_pp_t pp) {
    const slong k = pp->k;
    slong N;

    if (tr->k != k || tr->mixes == NULL) {
        return -1;
    }
    N = evmlr_voting_verify_mixes(tr->mixes, k, submissions, num_submissions, pp);
    // The last server's output rows 0..0 are the votes.
    if (N > 0 && results != NULL) {
        for (slong i = 0; i < N; i++) {
            nmod_poly_mat_init(results[i], 1, 1, MOD_Q);
            nmod_poly_set(nmod_poly_mat_entry(results[i], 0, 0),
                          nmod_poly_mat_entry(tr->mixes[k - 1].statement->c_hat[i], 0, 0));
        }
    }
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
        evmlr_enc_proof_init(&sub->proofs[j]->proof, evmlr_voting_otse_len(k, j));
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
}

void evmlr_voting_sk_clear(evmlr_voting_sk_t sk) {
    for (slong j = 0; j < sk->k; j++) {
        nmod_poly_mat_clear(sk->sks[j]->s);
    }
    free(sk->sks);
}

#ifdef MAIN
#include "evmlr_main.h"

static void set_vote(nmod_poly_mat_t vote, ulong value) {
    nmod_poly_mat_init(vote, 1, 1, MOD_Q);
    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(vote, 0, 0), 0, value);
}

static void submission_copy(evmlr_voting_submission_t dst, const evmlr_voting_submission_t src, const evmlr_voting_pp_t pp) {
    evmlr_voting_submission_init(dst, pp);
    dst->voter_id = src->voter_id;
    for (int i = 0; i < K_LWR; i++) {
        nmod_poly_mat_set(dst->ballot->layer->enc_cipher[i]->uT, src->ballot->layer->enc_cipher[i]->uT);
        nmod_poly_set(dst->ballot->layer->enc_cipher[i]->v, src->ballot->layer->enc_cipher[i]->v);
    }
    nmod_poly_mat_set(dst->ballot->layer->otse_cipher->c, src->ballot->layer->otse_cipher->c);
    for (slong j = 0; j < pp->k; j++) {
        evmlr_enc_proof_copy(&dst->proofs[j]->proof, &src->proofs[j]->proof);
    }
}

/* Counts how many of results[0..n-1] hold value. */
static int count_value(nmod_poly_mat_t* results, slong n, ulong value) {
    int c = 0;
    for (slong i = 0; i < n; i++) {
        c += nmod_poly_get_coeff_ui(nmod_poly_mat_entry(results[i], 0, 0), 0) == value;
    }
    return c;
}

/* The protocol as separate parties run it: three servers that each generate
 * their own key, voters (one revotes), a ballot copied under another voter
 * id and a ballot whose outer OTSE ciphertext was altered, then the servers
 * mixing one after the other and a verifier with public data only. */
static void test_voting_protocol(flint_rand_t state) {
    enum { K = 3, SUBS = 10 };
    const int64_t voter_ids[SUBS] = { 1, 2, 3, 2, 4, 5, 6, 7, 8, 1 };
    const ulong votes[SUBS]     = { 10, 99, 30, 20, 40, 50, 0, 70, 80, 0 };
    evmlr_mlpke_keypair_t kps[K];
    evmlr_mlpke_pk_t pks[K];
    uint8_t id[EVMLR_ELECTION_ID_LEN], other_id[EVMLR_ELECTION_ID_LEN];
    evmlr_voting_pp_t pp, pp_other;
    evmlr_voting_submission_t subs[SUBS];
    evmlr_voting_batch_struct batch[K + 1];
    evmlr_voting_transcript_t tr = { { 0, NULL } };
    nmod_poly_mat_t results[SUBS];
    evmlr_voting_submission_t tampered[3];
    int self_checks[SUBS] = { 0 };
    slong counted = -1;
    int batches = 0;

    for (int i = 0; i < EVMLR_ELECTION_ID_LEN; i++) {
        id[i] = (uint8_t) i;
        other_id[i] = (uint8_t) (i + 1);
    }
    for (int j = 0; j < K; j++) {
        evmlr_voting_keygen(kps[j], state);
        pks[j][0] = kps[j]->pk[0];
    }
    evmlr_voting_pp_init(pp, id, K, (const evmlr_mlpke_pk_t*) pks);
    evmlr_voting_pp_init(pp_other, other_id, K, (const evmlr_mlpke_pk_t*) pks);

    // Voters 1-5 vote honestly (voter 2 twice; the second ballot counts).
    for (int i = 0; i < 6; i++) {
        nmod_poly_mat_t vote;
        set_vote(vote, votes[i]);
        evmlr_voting_submission_init(subs[i], pp);
        self_checks[i] = evmlr_voting_casting(subs[i], voter_ids[i], vote, pp, state);
        nmod_poly_mat_clear(vote);
    }
    // Voter 6 copies voter 1's ballot under its own id, and someone posts a
    // copy of voter 3's ballot under voter 1's id, after voter 1 voted.
    submission_copy(subs[6], subs[0], pp);
    subs[6]->voter_id = voter_ids[6];
    submission_copy(subs[9], subs[2], pp);
    subs[9]->voter_id = voter_ids[9];
    // Two dishonest voters bind wrong layers into ballots that pass the
    // public check. Voter 7's outer proof claims another plaintext: the
    // ML-PKE half holds, the OTSE half fails once server 0 decrypts. Voter
    // 8's layer-1 proof is broken after it was made and checked: server 1
    // rejects it. The casting check catches voter 7's, which an honest client
    // would act on.
    for (int i = 7; i < 9; i++) {
        const cheat_t cheat = { i == 7 ? 0 : -1, i == 8 ? 1 : -1 };
        nmod_poly_mat_t vote;
        set_vote(vote, votes[i]);
        evmlr_voting_submission_init(subs[i], pp);
        self_checks[i] = cast(subs[i], voter_ids[i], vote, pp, state, &cheat);
        nmod_poly_mat_clear(vote);
    }
    // A ballot changed after casting: in the ciphertext, an inner proof, or
    // the OTSE half of the outer proof.
    for (int t = 0; t < 3; t++) {
        nmod_poly_struct* c;
        submission_copy(tampered[t], subs[4], pp);
        c = t == 0 ? nmod_poly_mat_entry(tampered[t]->ballot->layer->otse_cipher->c, 0, 0)
          : t == 1 ? nmod_poly_mat_entry(tampered[t]->proofs[K - 1]->proof.proof_otse.z_1, 0, 0)
                   : nmod_poly_mat_entry(tampered[t]->proofs[0]->proof.proof_otse.w, 0, 0);
        nmod_poly_set_coeff_ui(c, 0, (nmod_poly_get_coeff_ui(c, 0) + 1) % MOD_Q);
    }

    TEST_ONCE("voting: honest ballots pass the public check") {
        for (int i = 0; i < 6; i++) {
            TEST_ASSERT(evmlr_voting_check_ballot(pp, subs[i]) == 1, end);
        }
        TEST_ASSERT(evmlr_voting_check_ballot(pp, subs[7]) == 1, end);
        TEST_ASSERT(evmlr_voting_check_ballot(pp, subs[8]) == 1, end);
    } TEST_END;

    TEST_ONCE("voting: casting checks every layer, and catches a wrong one") {
        for (int i = 0; i < 6; i++) {
            TEST_ASSERT(self_checks[i] == 1, end);
        }
        // Voter 7's claims a wrong plaintext; voter 8's breaks its proof after the check.
        TEST_ASSERT(self_checks[7] == 0, end);
    } TEST_END;

    TEST_ONCE("voting: a ballot changed anywhere fails the public check") {
        for (int t = 0; t < 3; t++) {
            TEST_ASSERT(evmlr_voting_check_ballot(pp, tampered[t]) == 0, end);
        }
    } TEST_END;

    TEST_ONCE("voting: a ballot copied under another voter id fails it") {
        TEST_ASSERT(evmlr_voting_check_ballot(pp, subs[6]) == 0, end);
        TEST_ASSERT(evmlr_voting_check_ballot(pp, subs[9]) == 0, end);
    } TEST_END;

    TEST_ONCE("voting: a ballot replayed in another election fails it") {
        TEST_ASSERT(evmlr_voting_check_ballot(pp_other, subs[0]) == 0, end);
    } TEST_END;

    TEST_ONCE("voting: each server mixes with its own key only") {
        TEST_ASSERT(evmlr_voting_batch_from_submissions(&batch[0], (const evmlr_voting_submission_t*) subs, SUBS, pp) == EVMLR_VOTING_OK, end);
        batches = 1;
        // 10 submissions: voter 2 twice, and the two copies do not count.
        TEST_ASSERT(batch[0].n == 7, end);
        tr->k = K;
        tr->mixes = (evmlr_voting_mix_struct*) calloc(K, sizeof(evmlr_voting_mix_struct));
        // A key that is not server 0's is refused.
        TEST_ASSERT(evmlr_voting_mix_step(&batch[1], &tr->mixes[0], &batch[0], kps[1], pp, state) == EVMLR_VOTING_ERR_KEY, end);
        for (int j = 0; j < K; j++) {
            TEST_ASSERT(evmlr_voting_mix_step(&batch[j + 1], &tr->mixes[j], &batch[j], kps[j], pp, state) == EVMLR_VOTING_OK, end);
            batches = j + 2;
            TEST_ASSERT(evmlr_voting_batch_matches_mix(&batch[j + 1], &tr->mixes[j], pp) == 1, end);
        }
    } TEST_END;

    TEST_ONCE("voting: each prefix of the chain verifies before the next step") {
        nmod_poly_struct* c;
        ulong c0;
        TEST_ASSERT(evmlr_voting_verify_mixes(tr->mixes, 0, (const evmlr_voting_submission_t*) subs, SUBS, pp) == 7, end);
        for (int j = 0; j < K; j++) {
            TEST_ASSERT(evmlr_voting_verify_mixes(tr->mixes, j + 1, (const evmlr_voting_submission_t*) subs, SUBS, pp)
                        == batch[j + 1].n, end);
        }
        // A bad proof at server 0 is found before server 1 mixes.
        c = nmod_poly_mat_entry(tr->mixes[0].proof->u, 0, 0);
        c0 = nmod_poly_get_coeff_ui(c, 0);
        nmod_poly_set_coeff_ui(c, 0, (c0 + 1) % MOD_Q);
        TEST_ASSERT(evmlr_voting_verify_mixes(tr->mixes, 1, (const evmlr_voting_submission_t*) subs, SUBS, pp) == -1, end);
        nmod_poly_set_coeff_ui(c, 0, c0);
        TEST_ASSERT(evmlr_voting_verify_mixes(tr->mixes, K + 1, (const evmlr_voting_submission_t*) subs, SUBS, pp) == -1, end);
    } TEST_END;

    TEST_ONCE("voting: the servers reject the broken ballots, and say why") {
        const evmlr_voting_mix_struct *m0 = &tr->mixes[0], *m1 = &tr->mixes[1];
        // Server 0's inputs in cast order: voters 1, 3, 2, 4, 5, 7, 8.
        TEST_ASSERT(m0->n_in == 7 && m0->n_rejected == 1 && m0->N == 6, end);
        TEST_ASSERT(m0->rejected[0] == 5 && m0->reasons[0] == EVMLR_VOTING_REJECT_OTSE, end);
        // Server 1 cannot know whose it is; it drops one for its inner proof.
        TEST_ASSERT(m1->n_in == 6 && m1->n_rejected == 1 && m1->N == 5, end);
        TEST_ASSERT(m1->reasons[0] == EVMLR_VOTING_REJECT_HPKE, end);
        TEST_ASSERT(tr->mixes[2].n_in == 5 && tr->mixes[2].n_rejected == 0, end);
    } TEST_END;

    TEST_ONCE("voting: the last batch holds exactly the honest votes") {
        TEST_ASSERT(batch[K].n == 5, end);
        TEST_ASSERT(count_value(batch[K].msgs, 5, 10) == 1 && count_value(batch[K].msgs, 5, 20) == 1
                    && count_value(batch[K].msgs, 5, 30) == 1 && count_value(batch[K].msgs, 5, 40) == 1
                    && count_value(batch[K].msgs, 5, 50) == 1, end);
    } TEST_END;

    TEST_ONCE("voting: the count verifies from public data") {
        counted = evmlr_voting_verify_counting(results, tr, (const evmlr_voting_submission_t*) subs, SUBS, pp);
        TEST_ASSERT(counted == 5, end);
        for (slong i = 0; i < 5; i++) {
            TEST_ASSERT(nmod_poly_mat_equal(results[i], batch[K].msgs[i]), end);
        }
    } TEST_END;

    TEST_ONCE("voting: a tampered count does not") {
        evmlr_voting_mix_struct* m0 = &tr->mixes[0];
        int ok;
        // Voter 7's ballot rejected for its public proof, which it passes.
        m0->reasons[0] = EVMLR_VOTING_REJECT_HPKE;
        ok = evmlr_voting_verify_counting(NULL, tr, (const evmlr_voting_submission_t*) subs, SUBS, pp) == -1;
        m0->reasons[0] = EVMLR_VOTING_REJECT_OTSE;
        TEST_ASSERT(ok, end);
        // A rejected index that is not among the inputs.
        m0->rejected[0] = 7;
        ok = evmlr_voting_verify_counting(NULL, tr, (const evmlr_voting_submission_t*) subs, SUBS, pp) == -1;
        m0->rejected[0] = 5;
        TEST_ASSERT(ok, end);
        TEST_ASSERT(evmlr_voting_verify_counting(NULL, tr, (const evmlr_voting_submission_t*) subs, SUBS, pp) == 5, end);
        // A ballot left out of the input (voter 8's).
        TEST_ASSERT(evmlr_voting_verify_counting(NULL, tr, (const evmlr_voting_submission_t*) subs, SUBS - 2, pp) == -1, end);
        // An output of the second server changed.
        nmod_poly_struct* c = nmod_poly_mat_entry(tr->mixes[1].statement->c_hat[0], 0, 0);
        nmod_poly_set_coeff_ui(c, 0, (nmod_poly_get_coeff_ui(c, 0) + 1) % MOD_Q);
        TEST_ASSERT(evmlr_voting_verify_counting(NULL, tr, (const evmlr_voting_submission_t*) subs, SUBS, pp) == -1, end);
        TEST_ASSERT(evmlr_voting_batch_matches_mix(&batch[2], &tr->mixes[1], pp) == 0, end);
    } TEST_END;

end:
    if (counted > 0) {
        for (slong i = 0; i < counted; i++) {
            nmod_poly_mat_clear(results[i]);
        }
    }
    for (int b = 0; b < batches; b++) {
        evmlr_voting_batch_clear(&batch[b]);
    }
    if (tr->mixes != NULL) {
        evmlr_voting_transcript_clear(tr);
    }
    for (int i = 0; i < SUBS; i++) {
        evmlr_voting_submission_clear(subs[i], pp);
    }
    for (int t = 0; t < 3; t++) {
        evmlr_voting_submission_clear(tampered[t], pp);
    }
    for (int j = 0; j < K; j++) {
        evmlr_mlpke_keypair_clear(kps[j]);
    }
    evmlr_voting_pp_clear(pp_other);
    evmlr_voting_pp_clear(pp);
}

/* The one-process helpers still count, and an empty election counts zero. */
static void test_voting_dealer(flint_rand_t state) {
    evmlr_voting_pp_t pp;
    evmlr_voting_sk_t sk;
    evmlr_voting_submission_t subs[2];
    evmlr_voting_transcript_t tr;
    nmod_poly_mat_t votes[2], results[2];
    slong n;

    evmlr_voting_setup(pp, sk, 2, state);
    for (int i = 0; i < 2; i++) {
        set_vote(votes[i], 7 + i);
        evmlr_voting_submission_init(subs[i], pp);
        evmlr_voting_casting(subs[i], 42 + i, votes[i], pp, state);
    }

    TEST_ONCE("voting: setup and counting in one process") {
        n = evmlr_voting_counting_transcript(results, tr, (const evmlr_voting_submission_t*) subs, 2, pp, sk, state);
        TEST_ASSERT(n == 2, end);
        TEST_ASSERT(count_value(results, 2, 7) == 1 && count_value(results, 2, 8) == 1, end);
        TEST_ASSERT(evmlr_voting_verify_counting(NULL, tr, (const evmlr_voting_submission_t*) subs, 2, pp) == 2, end);
        nmod_poly_mat_clear(results[0]);
        nmod_poly_mat_clear(results[1]);
        evmlr_voting_transcript_clear(tr);
    } TEST_END;

    TEST_ONCE("voting: a single ballot is mixed, and verifies") {
        nmod_poly_struct* c;
        n = evmlr_voting_counting_transcript(results, tr, (const evmlr_voting_submission_t*) subs, 1, pp, sk, state);
        TEST_ASSERT(n == 1 && tr->mixes[0].N == 1 && tr->mixes[1].N == 1, end);
        TEST_ASSERT(count_value(results, 1, 7) == 1, end);
        nmod_poly_mat_clear(results[0]);
        TEST_ASSERT(evmlr_voting_verify_counting(results, tr, (const evmlr_voting_submission_t*) subs, 1, pp) == 1, end);
        TEST_ASSERT(count_value(results, 1, 7) == 1, end);
        nmod_poly_mat_clear(results[0]);
        // The copy it is shuffled with must be the same ciphertext.
        c = nmod_poly_mat_entry(tr->mixes[1].statement->c_hat[1], 0, 0);
        nmod_poly_set_coeff_ui(c, 0, (nmod_poly_get_coeff_ui(c, 0) + 1) % MOD_Q);
        TEST_ASSERT(evmlr_voting_verify_counting(NULL, tr, (const evmlr_voting_submission_t*) subs, 1, pp) == -1, end);
        evmlr_voting_transcript_clear(tr);
    } TEST_END;

    TEST_ONCE("voting: no ballots, no votes") {
        n = evmlr_voting_counting_transcript(results, tr, NULL, 0, pp, sk, state);
        TEST_ASSERT(n == 0 && tr->mixes[0].n_in == 0 && tr->mixes[1].N == 0, end);
        TEST_ASSERT(evmlr_voting_verify_counting(NULL, tr, NULL, 0, pp) == 0, end);
        evmlr_voting_transcript_clear(tr);
    } TEST_END;

end:
    for (int i = 0; i < 2; i++) {
        nmod_poly_mat_clear(votes[i]);
        evmlr_voting_submission_clear(subs[i], pp);
    }
    evmlr_voting_pp_clear(pp);
    evmlr_voting_sk_clear(sk);
}

void bench_voting_protocol(flint_rand_t state) {
    enum { K = 2, N = 100 };
    evmlr_mlpke_keypair_t kps[K];
    evmlr_mlpke_pk_t pks[K];
    uint8_t id[EVMLR_ELECTION_ID_LEN] = { 0 };
    evmlr_voting_pp_t pp;
    evmlr_voting_submission_t* subs = (evmlr_voting_submission_t*) malloc(N * sizeof(evmlr_voting_submission_t));
    evmlr_voting_batch_struct batch[K + 1];
    evmlr_voting_transcript_t tr;
    uint64_t t0, t_setup, t_cast, t_check = 0, t_mix[K], t_verify;

    printf("\n=== Voting Protocol Benchmark (N = %d voters, k = %d servers) ===\n", N, K);

    t0 = cpucycles();
    for (int j = 0; j < K; j++) {
        evmlr_voting_keygen(kps[j], state);
        pks[j][0] = kps[j]->pk[0];
    }
    evmlr_voting_pp_init(pp, id, K, (const evmlr_mlpke_pk_t*) pks);
    t_setup = cpucycles() - t0;

    t0 = cpucycles();
    for (slong i = 0; i < N; i++) {
        nmod_poly_mat_t vote;
        set_vote(vote, (ulong) (i + 1));
        evmlr_voting_submission_init(subs[i], pp);
        evmlr_voting_casting(subs[i], i + 1, vote, pp, state);
        nmod_poly_mat_clear(vote);
    }
    t_cast = cpucycles() - t0;

    for (slong i = 0; i < N; i++) {
        t0 = cpucycles();
        evmlr_voting_check_ballot(pp, subs[i]);
        t_check += cpucycles() - t0;
    }

    evmlr_voting_batch_from_submissions(&batch[0], (const evmlr_voting_submission_t*) subs, N, pp);
    tr->k = K;
    tr->mixes = (evmlr_voting_mix_struct*) calloc(K, sizeof(evmlr_voting_mix_struct));
    for (int j = 0; j < K; j++) {
        t0 = cpucycles();
        evmlr_voting_mix_step(&batch[j + 1], &tr->mixes[j], &batch[j], kps[j], pp, state);
        t_mix[j] = cpucycles() - t0;
    }

    t0 = cpucycles();
    slong verified = evmlr_voting_verify_counting(NULL, tr, (const evmlr_voting_submission_t*) subs, N, pp);
    t_verify = cpucycles() - t0;

    printf("Key generation + parameters:   %llu cycles\n", (unsigned long long) t_setup);
    printf("Casting (total):               %llu cycles (avg %llu per voter)\n",
           (unsigned long long) t_cast, (unsigned long long) (t_cast / N));
    printf("Public ballot check (total):   %llu cycles (avg %llu per ballot)\n",
           (unsigned long long) t_check, (unsigned long long) (t_check / N));
    for (int j = 0; j < K; j++) {
        printf("Mix step, server %d:            %llu cycles (decrypt, check proofs, shuffle, prove)\n",
               j, (unsigned long long) t_mix[j]);
    }
    printf("Verify counting (%ld votes):   %llu cycles\n", (long) verified, (unsigned long long) t_verify);
    printf("=======================================================================\n\n");

    for (int b = 0; b <= K; b++) {
        evmlr_voting_batch_clear(&batch[b]);
    }
    evmlr_voting_transcript_clear(tr);
    for (slong i = 0; i < N; i++) {
        evmlr_voting_submission_clear(subs[i], pp);
    }
    free(subs);
    for (int j = 0; j < K; j++) {
        evmlr_mlpke_keypair_clear(kps[j]);
    }
    evmlr_voting_pp_clear(pp);
}

int main(int argc, char *argv[]) {
    evmlr_mode_t mode = evmlr_mode(argc, argv);

    flint_rand_t state;
    evmlr_rand_init(state);

    if (evmlr_runs_tests(mode)) {
        test_voting_protocol(state);
        test_voting_dealer(state);
    }
    if (evmlr_runs_benches(mode)) bench_voting_protocol(state);

    flint_rand_clear(state);
    return test_status();
}

#endif
