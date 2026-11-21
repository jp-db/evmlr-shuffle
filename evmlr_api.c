#include "evmlr_api.h"

#include <stdlib.h>
#include <string.h>

#include "evmlr_serialize.h"
#include "evmlr_utils.h"

#ifdef MAIN
#include "test.h"
#include "bench.h"
#endif

/* The decoding errors are evmlr_serialize.h's, with the same numbers. */
_Static_assert((int) EVMLR_ERR_TRUNCATED == (int) EVMLR_SER_ERR_TRUNCATED, "status codes");
_Static_assert((int) EVMLR_ERR_ARG == (int) EVMLR_SER_ERR_ARG, "status codes");
_Static_assert((int) EVMLR_ERR_KEY == (int) EVMLR_SER_ERR_KEY, "status codes");
_Static_assert((int) EVMLR_ERR_INVALID == (int) EVMLR_SER_ERR_INVALID, "status codes");
_Static_assert(EVMLR_ELECTION_ID_BYTES == EVMLR_ELECTION_ID_LEN, "election id length");

void evmlr_bytes_free(evmlr_bytes *b) {
    if (b != NULL) {
        free(b->data);
        b->data = NULL;
        b->len = 0;
    }
}

const char *evmlr_strerror(int status) {
    switch (status) {
        case EVMLR_ERR_TOO_FEW: return "only one ciphertext left to mix; a shuffle needs two";
        case EVMLR_ERR_SPACE:   return "output array too small";
        default:                return evmlr_ser_strerror(status);
    }
}

uint32_t evmlr_modulus(void) {
    return MOD_Q;
}

uint32_t evmlr_vote_len(void) {
    return DEGREE_N;
}

/* ---------------------------------------------------------------------------
 * Helpers
 * ------------------------------------------------------------------------- */

static void empty(evmlr_bytes *b) {
    if (b != NULL) {
        b->data = NULL;
        b->len = 0;
    }
}

/* Hands a finished buffer to the caller, or frees it if anything failed. */
static int take(evmlr_bytes *out, evmlr_buf_t *b, int st) {
    if (st == EVMLR_OK && b->failed) {
        st = EVMLR_ERR_NOMEM;
    }
    if (st != EVMLR_OK) {
        evmlr_buf_clear(b);
        return st;
    }
    out->data = b->data;
    out->len = b->len;
    return EVMLR_OK;
}

static int from_voting(int st) {
    switch (st) {
        case EVMLR_VOTING_OK:          return EVMLR_OK;
        case EVMLR_VOTING_ERR_KEY:     return EVMLR_ERR_KEY;
        case EVMLR_VOTING_ERR_NOMEM:   return EVMLR_ERR_NOMEM;
        case EVMLR_VOTING_ERR_TOO_FEW: return EVMLR_ERR_TOO_FEW;
        default:                       return EVMLR_ERR_ARG;
    }
}

static void rand_init(flint_rand_t state) {
    ulong seed[2];

    evmlr_utils_random(seed, sizeof(seed));
    flint_rand_init(state);
    flint_rand_set_seed(state, seed[0], seed[1]);
}

/* Decodes n ballots into a fresh array; on failure nothing is left to free. */
static int decode_ballots(evmlr_voting_submission_struct **out, const evmlr_voting_pp_t pp,
                          const uint8_t *const *ballots, const size_t *lens, size_t n) {
    evmlr_voting_submission_struct *subs;

    *out = NULL;
    if (n > 0 && (ballots == NULL || lens == NULL)) {
        return EVMLR_ERR_ARG;
    }
    if (n > (size_t) WORD_MAX || (subs = (evmlr_voting_submission_struct *) malloc((n ? n : 1) * sizeof(*subs))) == NULL) {
        return EVMLR_ERR_NOMEM;
    }
    for (size_t i = 0; i < n; i++) {
        int st = evmlr_deserialize_submission(&subs[i], pp, ballots[i], lens[i]);
        if (st != EVMLR_OK) {
            for (size_t l = 0; l < i; l++) {
                evmlr_voting_submission_clear(&subs[l], pp);
            }
            free(subs);
            return st;
        }
    }
    *out = subs;
    return EVMLR_OK;
}

static void free_ballots(evmlr_voting_submission_struct *subs, size_t n, const evmlr_voting_pp_t pp) {
    for (size_t i = 0; i < n; i++) {
        evmlr_voting_submission_clear(&subs[i], pp);
    }
    free(subs);
}

/* Writes the n votes (1 x 1 matrices) as DEGREE_N coefficients each. */
static int write_votes(uint16_t *votes, size_t cap, size_t *n_votes, const nmod_poly_mat_struct *v, slong n) {
    if (n_votes != NULL) {
        *n_votes = (size_t) n;
    }
    if ((size_t) n > cap || (n > 0 && votes == NULL)) {
        return EVMLR_ERR_SPACE;
    }
    for (slong i = 0; i < n; i++) {
        const nmod_poly_struct *p = nmod_poly_mat_entry(&v[i], 0, 0);
        for (slong c = 0; c < DEGREE_N; c++) {
            votes[i * DEGREE_N + c] = (uint16_t) nmod_poly_get_coeff_ui(p, c);
        }
    }
    return EVMLR_OK;
}

/* ---------------------------------------------------------------------------
 * Keys and parameters
 * ------------------------------------------------------------------------- */

int evmlr_keygen(evmlr_bytes *pk, evmlr_bytes *keypair) {
    evmlr_mlpke_keypair_t kp;
    evmlr_buf_t b_pk, b_kp;
    flint_rand_t state;
    int st;

    empty(pk);
    empty(keypair);
    if (pk == NULL || keypair == NULL) {
        return EVMLR_ERR_ARG;
    }
    rand_init(state);
    evmlr_voting_keygen(kp, state);
    evmlr_buf_init(&b_pk);
    evmlr_buf_init(&b_kp);
    st = take(pk, &b_pk, evmlr_serialize_pk(&b_pk, kp->pk));
    if (st == EVMLR_OK) {
        st = take(keypair, &b_kp, evmlr_serialize_keypair(&b_kp, kp));
        if (st != EVMLR_OK) {
            evmlr_bytes_free(pk);
        }
    }
    /* The secret is in keypair now (or nowhere); wipe the working copy. */
    nmod_poly_mat_zero(kp->sk->s);
    evmlr_mlpke_keypair_clear(kp);
    flint_rand_clear(state);
    return st;
}

int evmlr_params(evmlr_bytes *pp_out, const uint8_t election_id[EVMLR_ELECTION_ID_BYTES], uint32_t k,
                 const uint8_t *const *pks, const size_t *pk_lens) {
    evmlr_mlpke_pk_t *keys;
    evmlr_voting_pp_t pp;
    evmlr_buf_t b;
    uint32_t done = 0;
    int st = EVMLR_OK;

    empty(pp_out);
    if (pp_out == NULL || election_id == NULL || k < 1 || pks == NULL || pk_lens == NULL) {
        return EVMLR_ERR_ARG;
    }
    if ((keys = (evmlr_mlpke_pk_t *) malloc((size_t) k * sizeof(evmlr_mlpke_pk_t))) == NULL) {
        return EVMLR_ERR_NOMEM;
    }
    for (; done < k && st == EVMLR_OK; done++) {
        st = evmlr_deserialize_pk(keys[done], pks[done], pk_lens[done]);
        if (st != EVMLR_OK) {
            break;
        }
    }
    if (st == EVMLR_OK) {
        evmlr_voting_pp_init(pp, election_id, (slong) k, (const evmlr_mlpke_pk_t *) keys);
        evmlr_buf_init(&b);
        st = take(pp_out, &b, evmlr_serialize_pp(&b, pp));
        evmlr_voting_pp_clear(pp);
    }
    for (uint32_t j = 0; j < done; j++) {
        nmod_poly_mat_clear(keys[j]->A);
        nmod_poly_mat_clear(keys[j]->t);
    }
    free(keys);
    return st;
}

int evmlr_params_info(const uint8_t *pp_in, size_t pp_len, uint32_t *k, uint8_t election_id[EVMLR_ELECTION_ID_BYTES]) {
    evmlr_voting_pp_t pp;
    int st = evmlr_deserialize_pp(pp, pp_in, pp_len);

    if (st != EVMLR_OK) {
        return st;
    }
    if (k != NULL) {
        *k = (uint32_t) pp->k;
    }
    if (election_id != NULL) {
        memcpy(election_id, pp->election_id, EVMLR_ELECTION_ID_BYTES);
    }
    evmlr_voting_pp_clear(pp);
    return EVMLR_OK;
}

/* ---------------------------------------------------------------------------
 * Ballots
 * ------------------------------------------------------------------------- */

int evmlr_cast(evmlr_bytes *ballot, const uint8_t *pp_in, size_t pp_len, int64_t voter_id,
               const uint16_t *vote, size_t vote_len) {
    evmlr_voting_submission_t sub;
    evmlr_voting_pp_t pp;
    nmod_poly_mat_t v;
    flint_rand_t state;
    evmlr_buf_t b;
    int st;

    empty(ballot);
    if (ballot == NULL || vote_len > DEGREE_N || (vote_len > 0 && vote == NULL)
        || (int64_t) (slong) voter_id != voter_id) { /* slong is 32 bits on wasm32 */
        return EVMLR_ERR_ARG;
    }
    for (size_t i = 0; i < vote_len; i++) {
        if (vote[i] >= MOD_Q) {
            return EVMLR_ERR_ARG;
        }
    }
    if ((st = evmlr_deserialize_pp(pp, pp_in, pp_len)) != EVMLR_OK) {
        return st;
    }
    nmod_poly_mat_init(v, 1, 1, MOD_Q);
    for (size_t i = 0; i < vote_len; i++) {
        nmod_poly_set_coeff_ui(nmod_poly_mat_entry(v, 0, 0), (slong) i, vote[i]);
    }
    rand_init(state);
    evmlr_voting_submission_init(sub, pp);
    evmlr_voting_casting(sub, (slong) voter_id, v, pp, state);
    evmlr_buf_init(&b);
    st = take(ballot, &b, evmlr_serialize_submission(&b, sub, pp));

    evmlr_voting_submission_clear(sub, pp);
    nmod_poly_mat_clear(v);
    flint_rand_clear(state);
    evmlr_voting_pp_clear(pp);
    return st;
}

int evmlr_check_ballot(const uint8_t *pp_in, size_t pp_len, const uint8_t *ballot, size_t ballot_len,
                       int64_t *voter_id) {
    evmlr_voting_submission_t sub;
    evmlr_voting_pp_t pp;
    int st;

    if ((st = evmlr_deserialize_pp(pp, pp_in, pp_len)) != EVMLR_OK) {
        return st;
    }
    st = evmlr_deserialize_submission(sub, pp, ballot, ballot_len);
    if (st == EVMLR_OK) {
        if (voter_id != NULL) {
            *voter_id = (int64_t) sub->voter_id;
        }
        st = evmlr_voting_check_ballot(pp, sub) ? EVMLR_OK : EVMLR_ERR_INVALID;
        evmlr_voting_submission_clear(sub, pp);
    }
    evmlr_voting_pp_clear(pp);
    return st;
}

/* ---------------------------------------------------------------------------
 * Mixing
 * ------------------------------------------------------------------------- */

int evmlr_first_batch(evmlr_bytes *batch_out, const uint8_t *pp_in, size_t pp_len,
                      const uint8_t *const *ballots, const size_t *ballot_lens, size_t n) {
    evmlr_voting_submission_struct *subs;
    evmlr_voting_batch_t batch;
    evmlr_voting_pp_t pp;
    evmlr_buf_t b;
    int st;

    empty(batch_out);
    if (batch_out == NULL) {
        return EVMLR_ERR_ARG;
    }
    if ((st = evmlr_deserialize_pp(pp, pp_in, pp_len)) != EVMLR_OK) {
        return st;
    }
    if ((st = decode_ballots(&subs, pp, ballots, ballot_lens, n)) == EVMLR_OK) {
        st = from_voting(evmlr_voting_batch_from_submissions(batch, (const evmlr_voting_submission_t *) subs,
                                                             (slong) n, pp));
        if (st == EVMLR_OK) {
            evmlr_buf_init(&b);
            st = take(batch_out, &b, evmlr_serialize_batch(&b, batch, pp));
            evmlr_voting_batch_clear(batch);
        }
        free_ballots(subs, n, pp);
    }
    evmlr_voting_pp_clear(pp);
    return st;
}

int evmlr_mix(evmlr_bytes *out_batch, evmlr_bytes *mix_out, const uint8_t *pp_in, size_t pp_len,
              const uint8_t *keypair, size_t keypair_len, const uint8_t *in_batch, size_t in_len,
              const uint8_t *prev_mix, size_t prev_len) {
    evmlr_voting_pp_t pp;
    evmlr_mlpke_keypair_t kp;
    evmlr_voting_batch_t in, out;
    evmlr_voting_mix_t mix;
    int st, have_kp = 0, have_in = 0;

    empty(out_batch);
    empty(mix_out);
    if (out_batch == NULL || mix_out == NULL) {
        return EVMLR_ERR_ARG;
    }
    if ((st = evmlr_deserialize_pp(pp, pp_in, pp_len)) != EVMLR_OK) {
        return st;
    }
    if ((st = evmlr_deserialize_keypair(kp, keypair, keypair_len)) == EVMLR_OK) {
        have_kp = 1;
        st = evmlr_deserialize_batch(in, pp, in_batch, in_len);
    }
    if (st == EVMLR_OK) {
        have_in = 1;
        if (in->layer >= pp->k) {
            st = EVMLR_ERR_ARG; /* the votes: nothing left to mix */
        } else if (in->layer == 0) {
            st = prev_mix == NULL ? EVMLR_OK : EVMLR_ERR_ARG;
        } else if (prev_mix == NULL) {
            st = EVMLR_ERR_ARG;
        } else {
            evmlr_voting_mix_t prev;

            st = evmlr_deserialize_mix(prev, pp, prev_mix, prev_len);
            if (st == EVMLR_OK) {
                if (prev->j != in->layer - 1 || !evmlr_voting_batch_matches_mix(in, prev, pp)) {
                    st = EVMLR_ERR_INVALID;
                }
                evmlr_voting_mix_clear(prev);
            }
        }
    }
    if (st == EVMLR_OK) {
        flint_rand_t state;

        rand_init(state);
        st = from_voting(evmlr_voting_mix_step(out, mix, in, kp, pp, state));
        flint_rand_clear(state);
        if (st == EVMLR_OK) {
            evmlr_buf_t b_out, b_mix;

            evmlr_buf_init(&b_out);
            evmlr_buf_init(&b_mix);
            st = take(out_batch, &b_out, evmlr_serialize_batch(&b_out, out, pp));
            if (st == EVMLR_OK) {
                st = take(mix_out, &b_mix, evmlr_serialize_mix(&b_mix, mix, pp));
                if (st != EVMLR_OK) {
                    evmlr_bytes_free(out_batch);
                }
            }
            evmlr_voting_batch_clear(out);
            evmlr_voting_mix_clear(mix);
        }
    }
    if (have_in) {
        evmlr_voting_batch_clear(in);
    }
    if (have_kp) {
        nmod_poly_mat_zero(kp->sk->s);
        evmlr_mlpke_keypair_clear(kp);
    }
    evmlr_voting_pp_clear(pp);
    return st;
}

int evmlr_mix_info(const uint8_t *pp_in, size_t pp_len, const uint8_t *mix_in, size_t mix_len,
                   uint32_t *server, uint32_t *n_in, uint32_t *n_out, uint32_t *n_rejected,
                   uint32_t *rejected, uint8_t *reasons, size_t cap) {
    evmlr_voting_pp_t pp;
    evmlr_voting_mix_t mix;
    int st;

    if ((st = evmlr_deserialize_pp(pp, pp_in, pp_len)) != EVMLR_OK) {
        return st;
    }
    if ((st = evmlr_deserialize_mix(mix, pp, mix_in, mix_len)) == EVMLR_OK) {
        if (server != NULL) *server = (uint32_t) mix->j;
        if (n_in != NULL) *n_in = (uint32_t) mix->n_in;
        if (n_out != NULL) *n_out = (uint32_t) mix->N;
        if (n_rejected != NULL) *n_rejected = (uint32_t) mix->n_rejected;
        if ((size_t) mix->n_rejected > cap || (mix->n_rejected > 0 && (rejected == NULL || reasons == NULL))) {
            st = EVMLR_ERR_SPACE;
        } else {
            for (slong r = 0; r < mix->n_rejected; r++) {
                rejected[r] = (uint32_t) mix->rejected[r];
                reasons[r] = mix->reasons[r];
            }
        }
        evmlr_voting_mix_clear(mix);
    }
    evmlr_voting_pp_clear(pp);
    return st;
}

/* ---------------------------------------------------------------------------
 * The count
 * ------------------------------------------------------------------------- */

int evmlr_transcript(evmlr_bytes *transcript, const uint8_t *pp_in, size_t pp_len,
                     const uint8_t *const *mixes, const size_t *mix_lens, size_t k) {
    evmlr_voting_transcript_t tr;
    evmlr_voting_pp_t pp;
    evmlr_buf_t b;
    int st;

    empty(transcript);
    if (transcript == NULL || mixes == NULL || mix_lens == NULL) {
        return EVMLR_ERR_ARG;
    }
    if ((st = evmlr_deserialize_pp(pp, pp_in, pp_len)) != EVMLR_OK) {
        return st;
    }
    if (k != (size_t) pp->k) {
        evmlr_voting_pp_clear(pp);
        return EVMLR_ERR_ARG;
    }
    tr->k = pp->k;
    tr->mixes = (evmlr_voting_mix_struct *) calloc(k, sizeof(evmlr_voting_mix_struct));
    if (tr->mixes == NULL) {
        evmlr_voting_pp_clear(pp);
        return EVMLR_ERR_NOMEM;
    }
    for (size_t j = 0; j < k && st == EVMLR_OK; j++) {
        st = evmlr_deserialize_mix(&tr->mixes[j], pp, mixes[j], mix_lens[j]);
        if (st == EVMLR_OK && tr->mixes[j].j != (slong) j) {
            st = EVMLR_ERR_ARG; /* out of order */
        }
    }
    if (st == EVMLR_OK) {
        evmlr_buf_init(&b);
        st = take(transcript, &b, evmlr_serialize_transcript(&b, tr, pp));
    }
    evmlr_voting_transcript_clear(tr);
    evmlr_voting_pp_clear(pp);
    return st;
}

int evmlr_verify(uint16_t *votes, size_t cap, size_t *n_votes, const uint8_t *pp_in, size_t pp_len,
                 const uint8_t *const *ballots, const size_t *ballot_lens, size_t n,
                 const uint8_t *transcript, size_t transcript_len) {
    evmlr_voting_submission_struct *subs;
    evmlr_voting_transcript_t tr;
    evmlr_voting_pp_t pp;
    nmod_poly_mat_t *results;
    slong counted;
    int st;

    if (n_votes != NULL) {
        *n_votes = 0;
    }
    if ((st = evmlr_deserialize_pp(pp, pp_in, pp_len)) != EVMLR_OK) {
        return st;
    }
    if ((st = evmlr_deserialize_transcript(tr, pp, transcript, transcript_len)) != EVMLR_OK) {
        evmlr_voting_pp_clear(pp);
        return st;
    }
    if ((st = decode_ballots(&subs, pp, ballots, ballot_lens, n)) == EVMLR_OK) {
        results = (nmod_poly_mat_t *) malloc((n ? n : 1) * sizeof(nmod_poly_mat_t));
        if (results == NULL) {
            st = EVMLR_ERR_NOMEM;
        } else {
            counted = evmlr_voting_verify_counting(results, tr, (const evmlr_voting_submission_t *) subs, (slong) n, pp);
            if (counted < 0) {
                st = EVMLR_ERR_INVALID;
            } else {
                st = write_votes(votes, cap, n_votes, (const nmod_poly_mat_struct *) results, counted);
                for (slong i = 0; i < counted; i++) {
                    nmod_poly_mat_clear(results[i]);
                }
            }
            free(results);
        }
        free_ballots(subs, n, pp);
    }
    evmlr_voting_transcript_clear(tr);
    evmlr_voting_pp_clear(pp);
    return st;
}

int evmlr_results(uint16_t *votes, size_t cap, size_t *n_votes, const uint8_t *pp_in, size_t pp_len,
                  const uint8_t *batch_in, size_t batch_len) {
    evmlr_voting_batch_t batch;
    evmlr_voting_pp_t pp;
    int st;

    if (n_votes != NULL) {
        *n_votes = 0;
    }
    if ((st = evmlr_deserialize_pp(pp, pp_in, pp_len)) != EVMLR_OK) {
        return st;
    }
    if ((st = evmlr_deserialize_batch(batch, pp, batch_in, batch_len)) == EVMLR_OK) {
        st = batch->layer == pp->k
           ? write_votes(votes, cap, n_votes, (const nmod_poly_mat_struct *) batch->msgs, batch->n)
           : EVMLR_ERR_ARG;
        evmlr_voting_batch_clear(batch);
    }
    evmlr_voting_pp_clear(pp);
    return st;
}

#ifdef MAIN
#include "evmlr_main.h"

#define K 3

/* An election run end to end through the byte API, each party seeing only
 * what it would receive: k servers, an authority, voters, a bulletin board
 * and an auditor. */
static void test(void) {
    evmlr_bytes pk[K], kp[K], pp = { 0 }, ballots[6], batch[K + 1], mixes[K], tr = { 0 };
    const uint8_t *pk_ptr[K], *ballot_ptr[6], *mix_ptr[K];
    size_t pk_len[K], ballot_len[6], mix_len[K], n_votes = 0;
    uint8_t id[EVMLR_ELECTION_ID_BYTES];
    const int64_t voters[6] = { 11, 12, 13, 12, 14, 15 };
    const uint16_t choice[6] = { 1, 2, 3, 3, 2, 1 };
    uint16_t votes[6 * DEGREE_N];
    int n_ballots = 0, n_batches = 0, n_mixes = 0;

    memset(pk, 0, sizeof(pk));
    memset(kp, 0, sizeof(kp));
    memset(batch, 0, sizeof(batch));
    memset(mixes, 0, sizeof(mixes));
    for (int i = 0; i < EVMLR_ELECTION_ID_BYTES; i++) {
        id[i] = (uint8_t) (0xa0 + i);
    }

    TEST_ONCE("api: servers generate their keys and the authority the parameters") {
        for (int j = 0; j < K; j++) {
            TEST_ASSERT(evmlr_keygen(&pk[j], &kp[j]) == EVMLR_OK, end);
            pk_ptr[j] = pk[j].data;
            pk_len[j] = pk[j].len;
        }
        TEST_ASSERT(evmlr_params(&pp, id, K, pk_ptr, pk_len) == EVMLR_OK, end);
        uint32_t k = 0;
        uint8_t id2[EVMLR_ELECTION_ID_BYTES];
        TEST_ASSERT(evmlr_params_info(pp.data, pp.len, &k, id2) == EVMLR_OK && k == K && memcmp(id, id2, sizeof(id)) == 0, end);
    } TEST_END;

    TEST_ONCE("api: voters cast, and the board checks every ballot") {
        for (int i = 0; i < 6; i++) {
            int64_t who = 0;
            TEST_ASSERT(evmlr_cast(&ballots[i], pp.data, pp.len, voters[i], &choice[i], 1) == EVMLR_OK, end);
            n_ballots = i + 1;
            TEST_ASSERT(evmlr_check_ballot(pp.data, pp.len, ballots[i].data, ballots[i].len, &who) == EVMLR_OK, end);
            TEST_ASSERT(who == voters[i], end);
            ballot_ptr[i] = ballots[i].data;
            ballot_len[i] = ballots[i].len;
        }
    } TEST_END;

    TEST_ONCE("api: the board refuses a ballot copied under another voter id") {
        uint8_t *copy = (uint8_t *) malloc(ballots[0].len);
        int st;
        memcpy(copy, ballots[0].data, ballots[0].len);
        copy[EVMLR_SER_HEADER_LEN] ^= 0x40; /* the voter id, a u64 right after the header */
        st = evmlr_check_ballot(pp.data, pp.len, copy, ballots[0].len, NULL);
        free(copy);
        TEST_ASSERT(st == EVMLR_ERR_INVALID, end);
    } TEST_END;

    TEST_ONCE("api: the servers mix in turn, each with its own key") {
        TEST_ASSERT(evmlr_first_batch(&batch[0], pp.data, pp.len, ballot_ptr, ballot_len, 6) == EVMLR_OK, end);
        n_batches = 1;
        /* Wrong key, and a missing previous record. */
        TEST_ASSERT(evmlr_mix(&batch[1], &mixes[0], pp.data, pp.len, kp[1].data, kp[1].len,
                              batch[0].data, batch[0].len, NULL, 0) == EVMLR_ERR_KEY, end);
        for (int j = 0; j < K; j++) {
            const uint8_t *prev = j ? mixes[j - 1].data : NULL;
            size_t prev_len = j ? mixes[j - 1].len : 0;
            TEST_ASSERT(evmlr_mix(&batch[j + 1], &mixes[j], pp.data, pp.len, kp[j].data, kp[j].len,
                                  batch[j].data, batch[j].len, prev, prev_len) == EVMLR_OK, end);
            n_batches = j + 2;
            n_mixes = j + 1;
            mix_ptr[j] = mixes[j].data;
            mix_len[j] = mixes[j].len;
            if (j == 1) {
                /* Server 2 refuses server 1's output passed with server 0's record. */
                evmlr_bytes x = { 0 }, y = { 0 };
                TEST_ASSERT(evmlr_mix(&x, &y, pp.data, pp.len, kp[2].data, kp[2].len,
                                      batch[2].data, batch[2].len, mixes[0].data, mixes[0].len) == EVMLR_ERR_INVALID, end);
            }
        }
        uint32_t server, n_in, n_out, n_rej;
        TEST_ASSERT(evmlr_mix_info(pp.data, pp.len, mixes[0].data, mixes[0].len, &server, &n_in, &n_out, &n_rej,
                                   NULL, NULL, 0) == EVMLR_OK, end);
        TEST_ASSERT(server == 0 && n_in == 5 && n_out == 5 && n_rej == 0, end);
    } TEST_END;

    TEST_ONCE("api: the last batch holds the votes, and the auditor gets the same") {
        uint16_t check[6 * DEGREE_N];
        size_t n_check = 0;
        int seen[4] = { 0, 0, 0, 0 };
        TEST_ASSERT(evmlr_results(votes, 6, &n_votes, pp.data, pp.len, batch[K].data, batch[K].len) == EVMLR_OK, end);
        TEST_ASSERT(n_votes == 5, end);
        for (size_t i = 0; i < n_votes; i++) {
            seen[votes[i * DEGREE_N] & 3]++;
        }
        TEST_ASSERT(seen[1] == 2 && seen[2] == 1 && seen[3] == 2, end); /* voter 12 changed 2 -> 3 */
        TEST_ASSERT(evmlr_transcript(&tr, pp.data, pp.len, mix_ptr, mix_len, K) == EVMLR_OK, end);
        TEST_ASSERT(evmlr_verify(check, 6, &n_check, pp.data, pp.len, ballot_ptr, ballot_len, 6, tr.data, tr.len) == EVMLR_OK, end);
        TEST_ASSERT(n_check == n_votes && memcmp(check, votes, n_votes * DEGREE_N * sizeof(uint16_t)) == 0, end);
        TEST_ASSERT(evmlr_verify(check, 1, &n_check, pp.data, pp.len, ballot_ptr, ballot_len, 6, tr.data, tr.len) == EVMLR_ERR_SPACE
                    && n_check == 5, end);
        /* Not the ballots that were counted. */
        TEST_ASSERT(evmlr_verify(check, 6, &n_check, pp.data, pp.len, ballot_ptr, ballot_len, 5, tr.data, tr.len) == EVMLR_ERR_INVALID, end);
    } TEST_END;

end:
    for (int j = 0; j < K; j++) {
        evmlr_bytes_free(&pk[j]);
        evmlr_bytes_free(&kp[j]);
    }
    for (int i = 0; i < n_ballots; i++) evmlr_bytes_free(&ballots[i]);
    for (int i = 0; i < n_batches; i++) evmlr_bytes_free(&batch[i]);
    for (int i = 0; i < n_mixes; i++) evmlr_bytes_free(&mixes[i]);
    evmlr_bytes_free(&pp);
    evmlr_bytes_free(&tr);
}

int main(int argc, char *argv[]) {
    evmlr_mode_t mode = evmlr_mode(argc, argv);

    if (evmlr_runs_tests(mode)) test();
    return test_status();
}
#endif
