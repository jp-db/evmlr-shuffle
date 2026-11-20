/**
 * @file
 *
 * JavaScript-facing entry points of the WebAssembly build.
 *
 * The FLINT types cannot cross into JavaScript, so everything here works on
 * opaque handles (pointers into the WASM heap, numbers on the JS side), plain
 * integer arrays and byte strings in the format of evmlr_serialize.h.
 * wasm/evmlr.mjs wraps these in a small JS API; see wasm/README.md.
 *
 * An election handle holds the public parameters and, optionally, the
 * servers' secret keys:
 *
 *  - authority:  setup(k), then export the public parameters (and the keys,
 *                for the servers);
 *  - voter page: import the public parameters, cast, export the ballot;
 *  - counting:   import parameters and keys, import the ballots, count.
 *
 * Exported byte strings are malloc'ed in the WASM heap: the caller copies
 * them out and frees them with free().
 */

#include <stdint.h>
#include <stdlib.h>
#include <sys/random.h>

#include <emscripten.h>

#include "flint/flint.h"

#include "../evmlr_serialize.h"
#include "../evmlr_voting.h"

typedef struct {
    flint_rand_t state;
    evmlr_voting_pp_t pp;
    evmlr_voting_sk_t sk;
    int has_sk;
} evmlr_wasm_election;

typedef struct {
    evmlr_voting_submission_t sub;
} evmlr_wasm_ballot;

static evmlr_wasm_election *election_new(void) {
    ulong seed[2];
    evmlr_wasm_election *e = (evmlr_wasm_election *) calloc(1, sizeof(*e));

    if (e == NULL) {
        return NULL;
    }
    flint_rand_init(e->state);
    if (getentropy(seed, sizeof(seed)) != 0) {
        abort();
    }
    flint_rand_set_seed(e->state, seed[0], seed[1]);
    return e;
}

/* Hands a finished buffer to the caller: its data, with the length in *len. */
static uint8_t *buf_release(evmlr_buf_t *b, int status, uint32_t *len) {
    uint8_t *data;

    if (status != EVMLR_SER_OK || b->failed || b->len > UINT32_MAX) {
        evmlr_buf_clear(b);
        *len = 0;
        return NULL;
    }
    data = b->data;
    *len = (uint32_t) b->len;
    return data; /* b is not cleared: the caller owns data now */
}

/* Every ballot in the demo is a 1 x 1 matrix: one polynomial of R_q. */
static void vote_from_coeffs(nmod_poly_mat_t vote, const uint32_t *coeffs, int len) {
    nmod_poly_mat_init(vote, 1, 1, MOD_Q);
    for (int i = 0; i < len && i < DEGREE_N; i++) {
        nmod_poly_set_coeff_ui(nmod_poly_mat_entry(vote, 0, 0), i, coeffs[i] % MOD_Q);
    }
}

/** Constants of the build, so the JS side can check what it talks to. */
EMSCRIPTEN_KEEPALIVE int evmlr_wasm_param(int which) {
    switch (which) {
        case 0: return MOD_Q;
        case 1: return DEGREE_N;
        case 2: return K_LWE;
        case 3: return K_LWR;
        default: return -1;
    }
}

EMSCRIPTEN_KEEPALIVE const char *evmlr_wasm_strerror(int status) {
    return evmlr_ser_strerror(status);
}

/**
 * Generates the parameters and the keys of @p k mixing servers.
 * Returns an election handle, or 0 if @p k < 1.
 */
EMSCRIPTEN_KEEPALIVE evmlr_wasm_election *evmlr_wasm_setup(int k) {
    evmlr_wasm_election *e;

    if (k < 1 || (e = election_new()) == NULL) {
        return NULL;
    }
    evmlr_voting_setup(e->pp, e->sk, k, e->state);
    e->has_sk = 1;
    return e;
}

EMSCRIPTEN_KEEPALIVE void evmlr_wasm_election_free(evmlr_wasm_election *e) {
    if (e == NULL) {
        return;
    }
    evmlr_voting_pp_clear(e->pp);
    if (e->has_sk) {
        evmlr_voting_sk_clear(e->sk);
    }
    flint_rand_clear(e->state);
    free(e);
}

EMSCRIPTEN_KEEPALIVE int evmlr_wasm_has_sk(const evmlr_wasm_election *e) {
    return e != NULL && e->has_sk;
}

/** Encodes the public parameters; returns the bytes and their length in *len. */
EMSCRIPTEN_KEEPALIVE uint8_t *evmlr_wasm_export_pp(const evmlr_wasm_election *e, uint32_t *len) {
    evmlr_buf_t b;

    evmlr_buf_init(&b);
    return buf_release(&b, e ? evmlr_serialize_pp(&b, e->pp) : EVMLR_SER_ERR_ARG, len);
}

/**
 * Creates an election handle from encoded public parameters, without keys.
 * Returns 0 on failure, with the evmlr_ser_status_t in *status.
 */
EMSCRIPTEN_KEEPALIVE evmlr_wasm_election *evmlr_wasm_import_pp(const uint8_t *in, uint32_t len, int *status) {
    evmlr_wasm_election *e = election_new();

    if (e == NULL) {
        *status = EVMLR_SER_ERR_NOMEM;
        return NULL;
    }
    *status = evmlr_deserialize_pp(e->pp, in, len);
    if (*status != EVMLR_SER_OK) {
        flint_rand_clear(e->state);
        free(e);
        return NULL;
    }
    return e;
}

EMSCRIPTEN_KEEPALIVE uint8_t *evmlr_wasm_export_sk(const evmlr_wasm_election *e, uint32_t *len) {
    evmlr_buf_t b;

    evmlr_buf_init(&b);
    return buf_release(&b, (e && e->has_sk) ? evmlr_serialize_sk(&b, e->sk) : EVMLR_SER_ERR_ARG, len);
}

/** Attaches encoded secret keys to an election; they must be for its k servers. */
EMSCRIPTEN_KEEPALIVE int evmlr_wasm_import_sk(evmlr_wasm_election *e, const uint8_t *in, uint32_t len) {
    evmlr_voting_sk_t sk;
    int st;

    if (e == NULL) {
        return EVMLR_SER_ERR_ARG;
    }
    st = evmlr_deserialize_sk(sk, in, len);
    if (st != EVMLR_SER_OK) {
        return st;
    }
    if (sk->k != e->pp->k) {
        evmlr_voting_sk_clear(sk);
        return EVMLR_SER_ERR_SHAPE;
    }
    if (e->has_sk) {
        evmlr_voting_sk_clear(e->sk);
    }
    e->sk[0] = sk[0];
    e->has_sk = 1;
    return EVMLR_SER_OK;
}

/**
 * Casts a ballot: encrypts the vote in one layer per server and proves each
 * layer's encryption correct. @p coeffs holds @p len coefficients of the vote
 * polynomial (reduced mod q; at most DEGREE_N are used). Needs only the
 * public parameters. Returns a ballot handle.
 */
EMSCRIPTEN_KEEPALIVE evmlr_wasm_ballot *evmlr_wasm_cast(evmlr_wasm_election *e, int voter_id,
                                                        const uint32_t *coeffs, int len) {
    evmlr_wasm_ballot *b;
    nmod_poly_mat_t vote;

    if (e == NULL || (b = (evmlr_wasm_ballot *) malloc(sizeof(*b))) == NULL) {
        return NULL;
    }
    vote_from_coeffs(vote, coeffs, len);
    evmlr_voting_submission_init(b->sub, e->pp);
    evmlr_voting_casting(b->sub, voter_id, vote, e->pp, e->state);
    nmod_poly_mat_clear(vote);
    return b;
}

EMSCRIPTEN_KEEPALIVE void evmlr_wasm_ballot_free(evmlr_wasm_election *e, evmlr_wasm_ballot *b) {
    if (e == NULL || b == NULL) {
        return;
    }
    evmlr_voting_submission_clear(b->sub, e->pp);
    free(b);
}

EMSCRIPTEN_KEEPALIVE int evmlr_wasm_ballot_voter(const evmlr_wasm_ballot *b) {
    return b ? (int) b->sub->voter_id : 0;
}

EMSCRIPTEN_KEEPALIVE uint8_t *evmlr_wasm_export_ballot(const evmlr_wasm_election *e, const evmlr_wasm_ballot *b,
                                                      uint32_t *len) {
    evmlr_buf_t buf;

    evmlr_buf_init(&buf);
    return buf_release(&buf, (e && b) ? evmlr_serialize_submission(&buf, b->sub, e->pp) : EVMLR_SER_ERR_ARG, len);
}

/** Decodes a ballot for this election. Returns 0 on failure, with *status set. */
EMSCRIPTEN_KEEPALIVE evmlr_wasm_ballot *evmlr_wasm_import_ballot(const evmlr_wasm_election *e, const uint8_t *in,
                                                                uint32_t len, int *status) {
    evmlr_wasm_ballot *b;

    if (e == NULL) {
        *status = EVMLR_SER_ERR_ARG;
        return NULL;
    }
    if ((b = (evmlr_wasm_ballot *) malloc(sizeof(*b))) == NULL) {
        *status = EVMLR_SER_ERR_NOMEM;
        return NULL;
    }
    *status = evmlr_deserialize_submission(b->sub, e->pp, in, len);
    if (*status != EVMLR_SER_OK) {
        free(b);
        return NULL;
    }
    return b;
}

/* counting takes an array of submissions by value; the structs are copied
 * shallowly, so the ballots keep owning their memory. */
static evmlr_voting_submission_struct *ballot_array(evmlr_wasm_ballot **ballots, int n) {
    evmlr_voting_submission_struct *subs = (evmlr_voting_submission_struct *) malloc(n * sizeof(*subs));

    for (int i = 0; subs != NULL && i < n; i++) {
        subs[i] = ballots[i]->sub[0];
    }
    return subs;
}

static void write_votes(uint32_t *out, nmod_poly_mat_t *results, slong counted) {
    for (slong i = 0; i < counted; i++) {
        out[i] = (uint32_t) nmod_poly_get_coeff_ui(nmod_poly_mat_entry(results[i], 0, 0), 0);
        nmod_poly_mat_clear(results[i]);
    }
}

/**
 * Runs the mixnet over @p n ballots (the last ballot of each voter counts):
 * every server checks the encryption proofs, decrypts its layer, shuffles and
 * proves the shuffle. Writes the constant coefficient of each decrypted vote
 * to @p out (room for @p n entries; the order is the shuffled one) and
 * returns how many votes were counted, or -1 without the secret keys.
 *
 * If @p transcript is non-NULL, it receives the encoded counting transcript
 * (every server's shuffle statement and proof, for evmlr_wasm_verify_counting),
 * malloc'ed, with its length in *transcript_len.
 */
EMSCRIPTEN_KEEPALIVE int evmlr_wasm_count(evmlr_wasm_election *e, evmlr_wasm_ballot **ballots,
                                          int n, uint32_t *out,
                                          uint8_t **transcript, uint32_t *transcript_len) {
    evmlr_voting_submission_struct *subs;
    nmod_poly_mat_t *results;
    slong counted;

    if (transcript != NULL) {
        *transcript = NULL;
        *transcript_len = 0;
    }
    if (e == NULL || !e->has_sk || n < 0) {
        return -1;
    }
    if (n == 0) {
        return 0;
    }
    subs = ballot_array(ballots, n);
    results = (nmod_poly_mat_t *) malloc(n * sizeof(*results));
    if (subs == NULL || results == NULL) {
        free(subs);
        free(results);
        return -1;
    }

    if (transcript != NULL) {
        evmlr_voting_transcript_t tr;
        evmlr_buf_t b;

        counted = evmlr_voting_counting_transcript(results, tr, (const evmlr_voting_submission_t *) subs,
                                                   n, e->pp, e->sk, e->state);
        evmlr_buf_init(&b);
        *transcript = buf_release(&b, evmlr_serialize_transcript(&b, tr, e->pp), transcript_len);
        evmlr_voting_transcript_clear(tr);
    } else {
        evmlr_shuffle_proof_t *proofs = (evmlr_shuffle_proof_t *) malloc(e->pp->k * sizeof(*proofs));

        counted = evmlr_voting_counting(results, proofs, (const evmlr_voting_submission_t *) subs,
                                        n, e->pp, e->sk, e->state);
        if (counted > 0) {
            for (slong j = 0; j < e->pp->k; j++) {
                evmlr_proof_clear(proofs[j], e->pp->shuf_ctx);
            }
        }
        free(proofs);
    }
    write_votes(out, results, counted);
    free(results);
    free(subs);
    return (int) counted;
}

/**
 * Checks a counting transcript against the ballots, with the public
 * parameters only (see evmlr_voting_verify_counting). On success writes the
 * votes to @p out (room for @p n entries) and returns their number. Returns
 * -1 if the count does not verify, with *status 0, or if the transcript does
 * not decode, with its evmlr_ser_status_t in *status.
 */
EMSCRIPTEN_KEEPALIVE int evmlr_wasm_verify_counting(const evmlr_wasm_election *e, evmlr_wasm_ballot **ballots,
                                                    int n, const uint8_t *transcript, uint32_t transcript_len,
                                                    uint32_t *out, int *status) {
    evmlr_voting_submission_struct *subs;
    evmlr_voting_transcript_t tr;
    nmod_poly_mat_t *results;
    slong verified;

    *status = EVMLR_SER_OK;
    if (e == NULL || n < 0) {
        *status = EVMLR_SER_ERR_ARG;
        return -1;
    }
    *status = evmlr_deserialize_transcript(tr, e->pp, transcript, transcript_len);
    if (*status != EVMLR_SER_OK) {
        return -1;
    }
    subs = ballot_array(ballots, n > 0 ? n : 0);
    results = (nmod_poly_mat_t *) malloc((n > 0 ? n : 1) * sizeof(*results));
    if ((n > 0 && subs == NULL) || results == NULL) {
        free(subs);
        free(results);
        evmlr_voting_transcript_clear(tr);
        *status = EVMLR_SER_ERR_NOMEM;
        return -1;
    }
    verified = evmlr_voting_verify_counting(results, tr, (const evmlr_voting_submission_t *) subs, n, e->pp);
    if (verified > 0) {
        write_votes(out, results, verified);
    }
    free(results);
    free(subs);
    evmlr_voting_transcript_clear(tr);
    return (int) verified;
}
