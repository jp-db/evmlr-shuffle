#ifndef EVMLR_VOTING_H
#define EVMLR_VOTING_H

#include <stdint.h>

#include "evmlr_shuffle.h"
#include "evmlr_enc_proof.h"
#include "sha.h"

/*
 * The voting protocol, as separate parties run it:
 *
 *   1. every mixing server j generates its own key pair (evmlr_voting_keygen)
 *      and publishes the public key;
 *   2. anyone builds the public parameters from the election id and the k
 *      public keys (evmlr_voting_pp_init); the rest of the parameters is
 *      derived from a hash, so nobody chooses it;
 *   3. voters cast ballots (evmlr_voting_casting); the bulletin board accepts
 *      a ballot only if evmlr_voting_check_ballot holds;
 *   4. the ballots become the first batch (evmlr_voting_batch_from_submissions),
 *      and server j = 0, ..., k - 1 in turn runs evmlr_voting_mix_step on the
 *      batch server j - 1 output, with its own secret key only. Each step
 *      publishes a mix record (what it rejected and why, the shuffle
 *      statement and its proof) and passes a batch to the next server;
 *   5. anyone checks the k mix records against the ballots
 *      (evmlr_voting_verify_counting); the last server's outputs are the votes.
 *
 * evmlr_voting_setup and evmlr_voting_counting_transcript run all of this in
 * one process with every secret key. They are for tests and benchmarks: in
 * an election, whoever holds all k keys can decrypt every ballot.
 */

#define EVMLR_ELECTION_ID_LEN 32

// Public parameters for the voting protocol
typedef struct {
    slong k;                    // Number of mixnet servers (k >= 1)
    slong L_max;                // Maximum vector length for the outermost layer
    uint8_t election_id[EVMLR_ELECTION_ID_LEN];
    evmlr_hpke_ctx_t hpke_ctx;  // HPKE context; its OTSE matrices are derived from (election_id, k)
    evmlr_mlpke_pk_t* pks;      // Public keys for servers 0 to k - 1
} evmlr_voting_pp_struct;
typedef evmlr_voting_pp_struct evmlr_voting_pp_t[1];

// Secret keys of all the servers, as evmlr_voting_setup deals them (tests only)
typedef struct {
    slong k;
    evmlr_mlpke_sk_t* sks;      // Secret keys for servers 0 to k - 1
} evmlr_voting_sk_struct;
typedef evmlr_voting_sk_struct evmlr_voting_sk_t[1];

// A ciphertext layer (holds a single HPKE cipher representation)
typedef struct {
    evmlr_mlpke_cipher_t enc_cipher[K_LWR];
    evmlr_otse_ciphertext_t otse_cipher;
} evmlr_voting_layer_struct;
typedef evmlr_voting_layer_struct evmlr_voting_layer_t[1];

// Ballot ciphertext containing the outermost encrypted layer
typedef struct {
    evmlr_voting_layer_t layer; // The outermost ciphertext layer
} evmlr_voting_ciphertext_struct;
typedef evmlr_voting_ciphertext_struct evmlr_voting_ciphertext_t[1];

// ZK Proof of correct encryption for a layer
typedef struct {
    evmlr_enc_proof_struct proof;
} evmlr_voting_proof_enc_struct;
typedef evmlr_voting_proof_enc_struct evmlr_voting_proof_enc_t[1];

// Submission from a voter containing voter ID, ballot ciphertext, and encryption proofs
typedef struct {
    int64_t voter_id;                 // the same width in every build (slong is 32 bits on wasm32)
    evmlr_voting_ciphertext_t ballot;
    evmlr_voting_proof_enc_t* proofs; // proofs[j] proves correct encryption of layers[j]
} evmlr_voting_submission_struct;
typedef evmlr_voting_submission_struct evmlr_voting_submission_t[1];

/* ---------------------------------------------------------------------------
 * Keys and parameters
 * ------------------------------------------------------------------------- */

// A mixing server's key pair. It depends on nothing but the scheme
// parameters, so each server generates its own; kp is initialised here and
// cleared with evmlr_mlpke_keypair_clear.
void evmlr_voting_keygen(evmlr_mlpke_keypair_t kp, flint_rand_t state);

// The public parameters of an election with k servers whose public keys are
// pks[0..k-1] (copied), in mixing order. The OTSE matrices are expanded from
// a hash of the election id and k. Returns 0, or -1 if k < 1.
int evmlr_voting_pp_init(evmlr_voting_pp_t pp, const uint8_t election_id[EVMLR_ELECTION_ID_LEN],
                         slong k, const evmlr_mlpke_pk_t* pks);

// Generates all k key pairs and the parameters of an election with a random
// id in one call. Tests and benchmarks only: the caller learns every key.
void evmlr_voting_setup(evmlr_voting_pp_t pp, evmlr_voting_sk_t sk, slong k, flint_rand_t state);

// The shuffle context of a mix of N ciphertexts. Its commitment key is
// expanded from a hash of pp and N, so every party derives the same one.
void evmlr_voting_shuffle_ctx_init(evmlr_shuffle_ctx_t ctx, const evmlr_voting_pp_t pp, slong N);

/* ---------------------------------------------------------------------------
 * Casting
 * ------------------------------------------------------------------------- */

// Encrypts vote in k layers and proves every layer. Every proof is bound to
// the election id and its layer, and the outermost one also to voter_id and
// to the proofs of the inner layers: a ballot copied under another voter id
// or into another election, or changed anywhere, fails
// evmlr_voting_check_ballot. Checks every layer's proof against its
// plaintext before returning, and returns 1 if all hold (0 means a bug, or
// a fault, and the ballot must not be cast).
int evmlr_voting_casting(evmlr_voting_submission_t submission, int64_t voter_id, const nmod_poly_mat_t vote, const evmlr_voting_pp_t pp, flint_rand_t state);

// What can be checked of a ballot without a secret key: the ML-PKE half of
// the outermost proof, under the ballot's voter id. It is bound to the rest
// of the ballot, so it fails if any part of the ballot was changed. It does
// not show that the inner layers are well formed, which only the servers can
// check as they decrypt; a voter who makes them badly loses their own vote.
// Returns 1 if it holds. A bulletin board should refuse ballots that fail it.
int evmlr_voting_check_ballot(const evmlr_voting_pp_t pp, const evmlr_voting_submission_t sub);

// The whole proof of encryption of a ballot's outermost layer, against its
// decryption (what server 0 checks). Returns 1 if it holds.
int evmlr_voting_verify_ballot(const evmlr_voting_pp_t pp, const evmlr_voting_submission_t sub,
                               const nmod_poly_mat_t plaintext);

/* ---------------------------------------------------------------------------
 * Mixing, one server at a time
 * ------------------------------------------------------------------------- */

// What a server takes in and hands on: n ciphertexts of layer `layer`, each
// flattened to a ((k - layer) K_LWR (K_LWE + 1) + 1) x 1 message, with the
// proofs of encryption of its remaining layers. Server j consumes the batch
// of layer j and outputs the batch of layer j + 1; the batch of layer k holds
// the votes (1 x 1 messages) and no proofs.
typedef struct {
    slong k;
    slong layer;
    slong n;
    nmod_poly_mat_t* msgs;            // n messages
    evmlr_enc_proof_struct* proofs;   // n * (k - layer): item i, layer + l at [i * (k - layer) + l]
    int64_t* voter_ids;               // layer 0 only (the outer proofs bind the voter); else NULL
} evmlr_voting_batch_struct;
typedef evmlr_voting_batch_struct evmlr_voting_batch_t[1];

// Why a server dropped a ciphertext.
typedef enum {
    EVMLR_VOTING_REJECT_HPKE = 1, // the ML-PKE half of its proof of encryption fails
    EVMLR_VOTING_REJECT_OTSE = 2  // the OTSE half fails against the decryption
} evmlr_voting_reject_t;

// What server j publishes: of the n_in ciphertexts it received, which it
// rejected (indices into its input, increasing) and why; the statement of its
// shuffle of the N = n_in - n_rejected others, and the proof. statement and
// proof are set only when N > 0, and have evmlr_voting_shuffle_len(N)
// entries: when one ciphertext is left, it is shuffled with a copy of itself.
typedef struct {
    slong k;
    slong j;
    slong n_in;
    slong n_rejected;
    slong* rejected;
    uint8_t* reasons;                 // evmlr_voting_reject_t, one per rejected index
    slong N;
    evmlr_shuffle_pp_t statement;     // D, c_hat (outputs), c_star (inputs)
    evmlr_shuffle_proof_t proof;
} evmlr_voting_mix_struct;
typedef evmlr_voting_mix_struct evmlr_voting_mix_t[1];

typedef enum {
    EVMLR_VOTING_OK = 0,
    EVMLR_VOTING_ERR_ARG = -1,        // the batch is not for this election, or of the wrong layer
    EVMLR_VOTING_ERR_KEY = -2,        // the key pair is not server j's in pp
    EVMLR_VOTING_ERR_NOMEM = -3
} evmlr_voting_status_t;

// The first batch: the ballots that count under the revote policy (the last
// one of each voter), in the order they were cast. Their proofs are copied.
int evmlr_voting_batch_from_submissions(evmlr_voting_batch_t batch, const evmlr_voting_submission_t* submissions,
                                        slong num_submissions, const evmlr_voting_pp_t pp);

// Server in->layer's step, with its key pair kp only: decrypts every
// ciphertext, checks its proof of encryption, drops those that fail, shuffles
// the rest and proves the shuffle. Initialises out (the batch of the next
// layer) and mix (to publish) on EVMLR_VOTING_OK; on an error, neither.
int evmlr_voting_mix_step(evmlr_voting_batch_t out, evmlr_voting_mix_t mix, const evmlr_voting_batch_t in,
                          const evmlr_mlpke_keypair_t kp, const evmlr_voting_pp_t pp, flint_rand_t state);

// Whether batch is what mix says its server output: same layer and count,
// and every message equal to the published output. The next server should
// check this before mixing; the proofs a batch carries are checked when they
// are used.
int evmlr_voting_batch_matches_mix(const evmlr_voting_batch_t batch, const evmlr_voting_mix_t mix,
                                   const evmlr_voting_pp_t pp);

// Allocates a batch of n zero messages of layer `layer`, with proofs for
// layers layer..k-1 and (layer 0) voter ids, for the caller to fill in.
int evmlr_voting_batch_init(evmlr_voting_batch_t batch, slong k, slong layer, slong n);
void evmlr_voting_batch_clear(evmlr_voting_batch_t batch);
void evmlr_voting_mix_clear(evmlr_voting_mix_t mix);

/* ---------------------------------------------------------------------------
 * The count and its verification
 * ------------------------------------------------------------------------- */

// The k mix records of a count, in server order.
typedef struct {
    slong k;
    evmlr_voting_mix_struct* mixes;   // k records
} evmlr_voting_transcript_struct;
typedef evmlr_voting_transcript_struct evmlr_voting_transcript_t[1];

// All k mix steps in one process, with every secret key (tests only).
// Initialises results[0..N-1] with the votes and tr with the mix records;
// returns N, the number of votes, or -1 on an error.
slong evmlr_voting_counting_transcript(nmod_poly_mat_t* results, evmlr_voting_transcript_t tr,
                                       const evmlr_voting_submission_t* submissions, slong num_submissions,
                                       const evmlr_voting_pp_t pp, const evmlr_voting_sk_t sk, flint_rand_t state);

// Checks a count without any secret: for every server j, that its inputs are
// the ballots that count (j = 0) or what server j - 1 output (j > 0); that
// its rejections are well formed, and at server 0 that every ballot rejected
// for its ML-PKE proof does fail that public check; and that the proof of
// its shuffle of the others verifies under the derived commitment key.
// Returns the number of votes and, if results is non-NULL, initialises
// results[0..N-1] with them (in the published order); returns -1, leaving
// results alone, if anything fails. Neither correct decryption nor rejections
// for the OTSE proof (which need the decryption) are proven yet.
slong evmlr_voting_verify_counting(nmod_poly_mat_t* results, const evmlr_voting_transcript_t tr,
                                   const evmlr_voting_submission_t* submissions, slong num_submissions,
                                   const evmlr_voting_pp_t pp);

// The checks of evmlr_voting_verify_counting on the first m (0..k) records
// only: mixes[j] is server j's, for j < m. Server m runs this before it
// mixes, to refuse a bad chain instead of extending it, and the board can
// run it after every step. Returns the number of ciphertexts record m - 1
// says it output (for m = 0: the number of ballots that count), or -1 if
// anything fails.
slong evmlr_voting_verify_mixes(const evmlr_voting_mix_struct* mixes, slong m,
                                const evmlr_voting_submission_t* submissions, slong num_submissions,
                                const evmlr_voting_pp_t pp);
void evmlr_voting_transcript_clear(evmlr_voting_transcript_t tr);

// Memory Management
void evmlr_voting_pp_clear(evmlr_voting_pp_t pp);
void evmlr_voting_sk_clear(evmlr_voting_sk_t sk);
void evmlr_voting_submission_init(evmlr_voting_submission_t sub, const evmlr_voting_pp_t pp);
void evmlr_voting_submission_clear(evmlr_voting_submission_t sub, const evmlr_voting_pp_t pp);

void evmlr_voting_layer_init(evmlr_voting_layer_t layer, slong otse_len);
void evmlr_voting_layer_clear(evmlr_voting_layer_t layer);

// Lengths fixed by the protocol: the OTSE part of a layer-`layer` ciphertext,
// and the flattened message server `layer` takes in (1 for layer k: a vote).
slong evmlr_voting_otse_len(slong k, slong layer);
slong evmlr_voting_msg_len(slong k, slong layer);
// The number of ciphertexts a shuffle of N proves (N, or 2 when N = 1).
slong evmlr_voting_shuffle_len(slong N);

#endif // EVMLR_VOTING_H
