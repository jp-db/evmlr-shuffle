/**
 * @file
 *
 * The voting protocol for an external application: a C ABI on byte strings.
 *
 * Nothing here exposes a FLINT type, so the header needs only the C standard
 * library, and any language with a C FFI (or JavaScript, through the
 * WebAssembly build) can use it. Every object crosses the boundary in the
 * encoding of evmlr_serialize.h: the same bytes on every platform, decoded
 * strictly, so a party never has to trust another's in-memory state.
 *
 * Who calls what:
 *
 *     each mixing server j   evmlr_keygen            -> public key (publish), key pair (keep)
 *     election authority     evmlr_params            -> public parameters, from k public keys
 *     voter                  evmlr_cast              -> ballot
 *     bulletin board         evmlr_check_ballot      -> accept or refuse a ballot
 *     server 0               evmlr_first_batch       -> its input, from the ballots on the board
 *     server j, in turn      evmlr_verify_mixes      -> whether records 0..j-1 hold, before mixing
 *                            evmlr_mix               -> batch for server j + 1, mix record (publish)
 *     anyone                 evmlr_mix_info          -> what a server rejected, and why
 *     anyone                 evmlr_transcript        -> the k mix records as one transcript
 *     anyone                 evmlr_verify            -> the votes, if the whole count verifies
 *     server k - 1           evmlr_results           -> the votes in its output batch
 *
 * No party ever holds more than one secret key. The election authority only
 * assembles public keys; the hash-derived parts of the parameters cannot be
 * chosen by anyone.
 *
 * Voter authentication is the application's job. The library binds the
 * voter id it is given into the ballot and counts the last valid ballot of
 * each id, but cannot tell who may post under an id: the bulletin board must
 * authenticate the voter before it accepts a ballot.
 *
 * Conventions:
 *  - Functions return an evmlr_status (0 on success).
 *  - Outputs go into evmlr_bytes the library allocates; release them with
 *    evmlr_bytes_free. They are left empty on failure.
 *  - A vote is EVMLR_VOTE_LEN coefficients mod EVMLR_MODULUS (a polynomial
 *    of R_q); shorter inputs are zero-padded. Most applications use only the
 *    first coefficient, as a candidate number.
 *  - Lists of byte strings are passed as an array of pointers and an array
 *    of lengths.
 *  - Calls are not thread-safe: the library keeps one global random stream.
 */

#ifndef EVMLR_API_H
#define EVMLR_API_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#define EVMLR_ELECTION_ID_BYTES 32

typedef enum {
    EVMLR_OK = 0,
    EVMLR_ERR_TRUNCATED = 1,  /**< input ends early */
    EVMLR_ERR_MAGIC = 2,      /**< not an evmlr encoding */
    EVMLR_ERR_VERSION = 3,    /**< unknown format version */
    EVMLR_ERR_TYPE = 4,       /**< an encoding of another object type */
    EVMLR_ERR_PARAMS = 5,     /**< made by a build with other scheme parameters */
    EVMLR_ERR_SHAPE = 6,      /**< a dimension differs from the one expected */
    EVMLR_ERR_VALUE = 7,      /**< unreduced or non-canonical value */
    EVMLR_ERR_TRAILING = 8,   /**< bytes left after the object */
    EVMLR_ERR_NOMEM = 9,      /**< allocation failed */
    EVMLR_ERR_ARG = 10,       /**< invalid argument, or an object for another election or step */
    EVMLR_ERR_KEY = 11,       /**< a key pair that is not server j's, or whose halves differ */
    EVMLR_ERR_INVALID = 12,   /**< decodes, but a proof or a count does not verify */
    /* 13 is no longer returned (it was "one ciphertext left to mix") */
    EVMLR_ERR_SPACE = 14      /**< an output array is too small; the count needed is returned */
} evmlr_status;

/** A byte string the library allocated. */
typedef struct {
    uint8_t *data;
    size_t len;
} evmlr_bytes;

void evmlr_bytes_free(evmlr_bytes *b);

const char *evmlr_strerror(int status);

/** Scheme constants: EVMLR_MODULUS (q) and EVMLR_VOTE_LEN (n, coefficients per vote). */
uint32_t evmlr_modulus(void);
uint32_t evmlr_vote_len(void);

/** A mixing server's key pair: the public key to publish, the pair to keep secret. */
int evmlr_keygen(evmlr_bytes *pk, evmlr_bytes *keypair);

/**
 * The public parameters of an election with @p k servers, whose public keys
 * are pks[0..k-1] in mixing order. @p election_id identifies the election
 * (e.g. a hash of its description); ballots of one election do not verify
 * in another.
 */
int evmlr_params(evmlr_bytes *pp, const uint8_t election_id[EVMLR_ELECTION_ID_BYTES], uint32_t k,
                 const uint8_t *const *pks, const size_t *pk_lens);

/** Reads k and the election id out of public parameters (either may be NULL). */
int evmlr_params_info(const uint8_t *pp, size_t pp_len, uint32_t *k, uint8_t election_id[EVMLR_ELECTION_ID_BYTES]);

/**
 * Encrypts @p vote (vote_len <= EVMLR_VOTE_LEN coefficients, each < q) for
 * every server and proves every layer, bound to the election, to
 * @p voter_id (any int64_t, in every build) and, through the outermost
 * proof, to the whole ballot. Checks every layer before returning
 * (EVMLR_ERR_INVALID if one fails, which means a bug or a fault: do not
 * cast). Slow (seconds); in a browser, call it from a Web Worker.
 */
int evmlr_cast(evmlr_bytes *ballot, const uint8_t *pp, size_t pp_len, int64_t voter_id,
               const uint16_t *vote, size_t vote_len);

/**
 * What a bulletin board checks before accepting a ballot: that it decodes
 * for this election and that its public proof holds for its voter id
 * (EVMLR_ERR_INVALID otherwise). The proof is bound to every byte of the
 * ballot, so a ballot copied under another id, or changed anywhere, fails.
 * It cannot show that the inner layers were made correctly (only the
 * servers can, as they decrypt): a voter whose client makes them wrong
 * loses their own vote, and nobody else's. Writes the voter id to
 * @p voter_id if non-NULL.
 */
int evmlr_check_ballot(const uint8_t *pp, size_t pp_len, const uint8_t *ballot, size_t ballot_len,
                       int64_t *voter_id);

/**
 * Server 0's input: of the ballots that pass evmlr_check_ballot, the last of
 * each voter, in the order given (the board's order). Ballots that do not
 * decode make the call fail; ballots that fail the check are skipped, as
 * evmlr_verify will skip them.
 */
int evmlr_first_batch(evmlr_bytes *batch, const uint8_t *pp, size_t pp_len,
                      const uint8_t *const *ballots, const size_t *ballot_lens, size_t n);

/**
 * One server's step. The batch says which server it is for; @p keypair must
 * be that server's. For server j > 0, @p prev_mix is server j - 1's mix
 * record, and the batch must be what it says was output (EVMLR_ERR_INVALID
 * otherwise); for server 0 it must be NULL. That checks the hand-over only:
 * call evmlr_verify_mixes first to check the proofs of servers 0..j-1. A
 * single ciphertext left is mixed too (shuffled with a copy of itself),
 * though nothing can hide whose it is. Writes the batch for the next
 * server (for the last server: the votes, see evmlr_results) and the mix
 * record to publish.
 */
int evmlr_mix(evmlr_bytes *out_batch, evmlr_bytes *mix, const uint8_t *pp, size_t pp_len,
              const uint8_t *keypair, size_t keypair_len, const uint8_t *in_batch, size_t in_len,
              const uint8_t *prev_mix, size_t prev_len);

/**
 * What a mix record says: the server, how many ciphertexts it received and
 * passed on, and which it rejected (indices into its input) with the
 * reasons (1: the ML-PKE proof fails, 2: the OTSE proof fails against the
 * decryption). @p rejected and @p reasons have room for @p cap entries
 * (both may be NULL with cap 0); EVMLR_ERR_SPACE if there are more, with
 * the number in *n_rejected.
 */
int evmlr_mix_info(const uint8_t *pp, size_t pp_len, const uint8_t *mix, size_t mix_len,
                   uint32_t *server, uint32_t *n_in, uint32_t *n_out, uint32_t *n_rejected,
                   uint32_t *rejected, uint8_t *reasons, size_t cap);

/**
 * Checks the first @p m (0..k) mix records, mixes[j] being server j's,
 * against the ballots on the board (in the board's order): everything
 * evmlr_verify checks of those servers. Server m calls it before mixing,
 * so a bad record is refused instead of extended, and the board can call
 * it after every step. Writes the number of ciphertexts record m - 1 says
 * were output (for m = 0: the ballots that count) to @p n_out if non-NULL.
 * EVMLR_ERR_INVALID if anything fails; EVMLR_ERR_ARG if the records are
 * out of order or more than k.
 */
int evmlr_verify_mixes(uint32_t *n_out, const uint8_t *pp, size_t pp_len,
                       const uint8_t *const *ballots, const size_t *ballot_lens, size_t n,
                       const uint8_t *const *mixes, const size_t *mix_lens, size_t m);

/** The k mix records (in server order) as one counting transcript. */
int evmlr_transcript(evmlr_bytes *transcript, const uint8_t *pp, size_t pp_len,
                     const uint8_t *const *mixes, const size_t *mix_lens, size_t k);

/**
 * Checks a count from public data: the parameters, every ballot on the board
 * (in the board's order) and the transcript. On success writes the votes,
 * EVMLR_VOTE_LEN coefficients each, to @p votes (room for @p cap votes) and
 * their number to *n_votes. EVMLR_ERR_INVALID if the count does not verify;
 * EVMLR_ERR_SPACE if cap is too small, with the number in *n_votes.
 */
int evmlr_verify(uint16_t *votes, size_t cap, size_t *n_votes, const uint8_t *pp, size_t pp_len,
                 const uint8_t *const *ballots, const size_t *ballot_lens, size_t n,
                 const uint8_t *transcript, size_t transcript_len);

/** The votes in the last server's output batch, as evmlr_verify writes them. */
int evmlr_results(uint16_t *votes, size_t cap, size_t *n_votes, const uint8_t *pp, size_t pp_len,
                  const uint8_t *batch, size_t batch_len);

#ifdef __cplusplus
}
#endif

#endif /* EVMLR_API_H */
