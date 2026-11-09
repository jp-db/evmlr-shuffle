/**
 * @file
 *
 * JavaScript-facing entry points of the WebAssembly build.
 *
 * The FLINT types cannot cross into JavaScript, so everything here works on
 * opaque handles (pointers into the WASM heap, numbers on the JS side) and
 * plain integer arrays. wasm/evmlr.mjs wraps these in a small JS class; see
 * wasm/README.md.
 *
 * An election handle owns the public parameters and, because this build has
 * no serialisation of keys yet, also the servers' secret keys: setup, casting
 * and counting all run in the same module. In a deployment the voter's page
 * would only hold the public parameters and call cast.
 */

#include <stdint.h>
#include <stdlib.h>
#include <sys/random.h>

#include <emscripten.h>

#include "flint/flint.h"

#include "../evmlr_voting.h"

typedef struct {
    flint_rand_t state;
    evmlr_voting_pp_t pp;
    evmlr_voting_sk_t sk;
} evmlr_wasm_election;

typedef struct {
    evmlr_voting_submission_t sub;
} evmlr_wasm_ballot;

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

/**
 * Generates the parameters and the keys of @p k mixing servers.
 * Returns an election handle, or 0 if @p k < 1.
 */
EMSCRIPTEN_KEEPALIVE evmlr_wasm_election *evmlr_wasm_setup(int k) {
    ulong seed[2];
    evmlr_wasm_election *e;

    if (k < 1) {
        return NULL;
    }
    e = (evmlr_wasm_election *) malloc(sizeof(*e));
    if (e == NULL) {
        return NULL;
    }
    flint_rand_init(e->state);
    if (getentropy(seed, sizeof(seed)) != 0) {
        abort();
    }
    flint_rand_set_seed(e->state, seed[0], seed[1]);
    evmlr_voting_setup(e->pp, e->sk, k, e->state);
    return e;
}

EMSCRIPTEN_KEEPALIVE void evmlr_wasm_election_free(evmlr_wasm_election *e) {
    if (e == NULL) {
        return;
    }
    evmlr_voting_pp_clear(e->pp);
    evmlr_voting_sk_clear(e->sk);
    flint_rand_clear(e->state);
    free(e);
}

/**
 * Casts a ballot: encrypts the vote in one layer per server and proves each
 * layer's encryption correct. @p coeffs holds @p len coefficients of the vote
 * polynomial (reduced mod q; at most DEGREE_N are used).
 * Returns a ballot handle.
 */
EMSCRIPTEN_KEEPALIVE evmlr_wasm_ballot *evmlr_wasm_cast(evmlr_wasm_election *e, int voter_id,
                                                        const uint32_t *coeffs, int len) {
    evmlr_wasm_ballot *b;
    nmod_poly_mat_t vote;

    if (e == NULL) {
        return NULL;
    }
    b = (evmlr_wasm_ballot *) malloc(sizeof(*b));
    if (b == NULL) {
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

/**
 * Runs the mixnet over @p n ballots (the last ballot of each voter counts):
 * every server checks the encryption proofs, decrypts its layer, shuffles and
 * proves the shuffle. Writes the constant coefficient of each decrypted vote
 * to @p out (room for @p n entries; the order is the shuffled one) and
 * returns how many votes were counted.
 */
EMSCRIPTEN_KEEPALIVE int evmlr_wasm_count(evmlr_wasm_election *e, evmlr_wasm_ballot **ballots,
                                          int n, uint32_t *out) {
    evmlr_voting_submission_struct *subs;
    nmod_poly_mat_t *results;
    evmlr_shuffle_proof_t *proofs;
    slong counted;

    if (e == NULL || n < 1) {
        return 0;
    }
    /* counting takes an array of submissions by value; the structs are
     * copied shallowly, so the ballots keep owning their memory. */
    subs = (evmlr_voting_submission_struct *) malloc(n * sizeof(*subs));
    results = (nmod_poly_mat_t *) malloc(n * sizeof(*results));
    proofs = (evmlr_shuffle_proof_t *) malloc(e->pp->k * sizeof(*proofs));
    for (int i = 0; i < n; i++) {
        subs[i] = ballots[i]->sub[0];
    }

    counted = evmlr_voting_counting(results, proofs, (const evmlr_voting_submission_t *) subs,
                                    n, e->pp, e->sk, e->state);

    for (slong i = 0; i < counted; i++) {
        out[i] = (uint32_t) nmod_poly_get_coeff_ui(nmod_poly_mat_entry(results[i], 0, 0), 0);
        nmod_poly_mat_clear(results[i]);
    }
    if (counted > 0) {
        for (slong j = 0; j < e->pp->k; j++) {
            evmlr_proof_clear(proofs[j], e->pp->shuf_ctx);
        }
    }
    free(proofs);
    free(results);
    free(subs);
    return (int) counted;
}
