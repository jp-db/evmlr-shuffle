/**
 * @file
 *
 * Byte encodings of the voting protocol's keys, parameters and proofs, so
 * they can move between processes and platforms (a native server and a
 * WebAssembly voter page produce and accept the same bytes).
 *
 * Every encoding starts with a 20-byte header:
 *
 *     "EVML"  magic
 *     u8      format version (EVMLR_SER_VERSION)
 *     u8      object type (evmlr_ser_type_t)
 *     14 B    parameter fingerprint: MOD_Q, DEGREE_N, ALPHA, NONZERO (u16 each),
 *             K_LWE, K_LWR, K_SIS, ETA, ZETA, LOG_Q_CEIL (u8 each)
 *
 * so bytes from a build with other parameters are refused instead of being
 * misread. All integers are little-endian and fixed width, whatever the
 * platform's long is (64 bits on x86-64, 32 on wasm32).
 *
 * A matrix is written as u32 rows, u32 cols and then its entries row by row.
 * A polynomial is a u16 length followed by its coefficients packed at
 * ceil(log2(modulus)) bits each (12 bits mod q = 3109, 6 bits mod 2^ZETA),
 * zero-padded to a whole byte.
 *
 * Decoding is strict, because ballots and proofs come from untrusted parties.
 * It checks every dimension against the one the protocol fixes, requires
 * every coefficient to be reduced, and rejects polynomials with a zero
 * leading coefficient or nonzero padding bits, and trailing bytes. So each
 * object has exactly one encoding, and an object read back hashes into the
 * Fiat-Shamir transcript exactly as the original did. Sizes are checked
 * against the input length before anything is allocated.
 *
 * Every deserialize function initialises its output only on success
 * (EVMLR_SER_OK); the caller then owns it and clears it with the usual clear
 * function. On failure nothing is left to clear.
 */

#ifndef EVMLR_SERIALIZE_H
#define EVMLR_SERIALIZE_H

#include <stddef.h>
#include <stdint.h>

#include "evmlr_voting.h"

#define EVMLR_SER_VERSION 1
#define EVMLR_SER_HEADER_LEN 20

typedef enum {
    EVMLR_SER_PP              = 1, /**< voting public parameters (incl. the servers' public keys) */
    EVMLR_SER_SK              = 2, /**< the servers' secret keys */
    EVMLR_SER_SUBMISSION      = 3, /**< a ballot with its proofs of correct encryption */
    EVMLR_SER_SHUFFLE_PROOF   = 4, /**< a proof of shuffle */
    EVMLR_SER_SHUFFLE_STMT    = 5  /**< what a proof of shuffle is about: D, c_hat, c_star */
} evmlr_ser_type_t;

typedef enum {
    EVMLR_SER_OK = 0,
    EVMLR_SER_ERR_TRUNCATED,  /**< input ends early */
    EVMLR_SER_ERR_MAGIC,      /**< not an evmlr encoding */
    EVMLR_SER_ERR_VERSION,    /**< unknown format version */
    EVMLR_SER_ERR_TYPE,       /**< an encoding of another object type */
    EVMLR_SER_ERR_PARAMS,     /**< made by a build with other scheme parameters */
    EVMLR_SER_ERR_SHAPE,      /**< a dimension differs from the one expected */
    EVMLR_SER_ERR_VALUE,      /**< unreduced or non-canonical value */
    EVMLR_SER_ERR_TRAILING,   /**< bytes left after the object */
    EVMLR_SER_ERR_NOMEM,      /**< allocation failed */
    EVMLR_SER_ERR_ARG         /**< invalid argument (e.g. an unencodable object) */
} evmlr_ser_status_t;

/** A growable output buffer. Zero-initialise it or call evmlr_buf_init. */
typedef struct {
    uint8_t *data;
    size_t len;
    size_t cap;
    int failed; /**< set if an allocation failed; the content is then unusable */
} evmlr_buf_t;

void evmlr_buf_init(evmlr_buf_t *buf);
void evmlr_buf_clear(evmlr_buf_t *buf);

const char *evmlr_ser_strerror(int status);

/** Reads the type out of an encoding's header, or returns a negative status. */
int evmlr_ser_peek_type(const uint8_t *in, size_t len);

/* Public parameters: k, the HPKE/OTSE matrices, the servers' public keys and
 * the current shuffle (commitment) context. */
int evmlr_serialize_pp(evmlr_buf_t *out, const evmlr_voting_pp_t pp);
int evmlr_deserialize_pp(evmlr_voting_pp_t pp, const uint8_t *in, size_t len);

/* Secret keys of all k servers. */
int evmlr_serialize_sk(evmlr_buf_t *out, const evmlr_voting_sk_t sk);
int evmlr_deserialize_sk(evmlr_voting_sk_t sk, const uint8_t *in, size_t len);

/* A submission: voter id, the outermost ciphertext layer and the k proofs of
 * correct encryption. Its dimensions follow from pp. */
int evmlr_serialize_submission(evmlr_buf_t *out, const evmlr_voting_submission_t sub,
                               const evmlr_voting_pp_t pp);
int evmlr_deserialize_submission(evmlr_voting_submission_t sub, const evmlr_voting_pp_t pp,
                                 const uint8_t *in, size_t len);

/* A proof of shuffle, and the statement it proves, for shuffle context ctx
 * (in the voting protocol, pp->shuf_ctx as counting left it). */
int evmlr_serialize_shuffle_proof(evmlr_buf_t *out, const evmlr_shuffle_proof_t proof,
                                  const evmlr_shuffle_ctx_t ctx);
int evmlr_deserialize_shuffle_proof(evmlr_shuffle_proof_t proof, const evmlr_shuffle_ctx_t ctx,
                                    const uint8_t *in, size_t len);
int evmlr_serialize_shuffle_statement(evmlr_buf_t *out, const evmlr_shuffle_pp_t stmt,
                                      const evmlr_shuffle_ctx_t ctx);
int evmlr_deserialize_shuffle_statement(evmlr_shuffle_pp_t stmt, const evmlr_shuffle_ctx_t ctx,
                                        const uint8_t *in, size_t len);

#endif /* EVMLR_SERIALIZE_H */
