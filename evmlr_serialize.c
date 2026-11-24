#include "evmlr_serialize.h"
#include "evmlr_utils.h"

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#ifdef MAIN
#include "test.h"
#include "bench.h"
#endif

/* Polynomial lengths are written as u16; elements of R_q never exceed n. */
_Static_assert(DEGREE_N <= 65535, "polynomial lengths are encoded in 16 bits");

/* Upper bounds on the counts read from untrusted input, so that the size
 * arithmetic below cannot overflow; the input length is the real bound. */
#define MAX_SERVERS  (1 << 16)
#define MAX_MESSAGES (1 << 24)

static const uint8_t MAGIC[4] = { 'E', 'V', 'M', 'L' };

/* ---------------------------------------------------------------------------
 * Output buffer
 * ------------------------------------------------------------------------- */

void evmlr_buf_init(evmlr_buf_t *buf) {
    memset(buf, 0, sizeof(*buf));
}

void evmlr_buf_clear(evmlr_buf_t *buf) {
    free(buf->data);
    memset(buf, 0, sizeof(*buf));
}

typedef struct {
    evmlr_buf_t *buf;
    int err;
} writer_t;

static void put_bytes(writer_t *w, const void *p, size_t n) {
    evmlr_buf_t *b = w->buf;

    if (w->err || n == 0) {
        return;
    }
    if (b->failed || n > SIZE_MAX - b->len) {
        w->err = EVMLR_SER_ERR_NOMEM;
        b->failed = 1;
        return;
    }
    if (b->len + n > b->cap) {
        size_t cap = b->cap ? b->cap : 4096;
        uint8_t *d;

        while (cap < b->len + n) {
            cap = (cap > SIZE_MAX / 2) ? b->len + n : 2 * cap;
        }
        d = (uint8_t *) realloc(b->data, cap);
        if (d == NULL) {
            w->err = EVMLR_SER_ERR_NOMEM;
            b->failed = 1;
            return;
        }
        b->data = d;
        b->cap = cap;
    }
    memcpy(b->data + b->len, p, n);
    b->len += n;
}

static void put_uint(writer_t *w, uint64_t v, int nbytes) {
    uint8_t b[8];

    for (int i = 0; i < nbytes; i++) {
        b[i] = (uint8_t) (v >> (8 * i));
    }
    put_bytes(w, b, (size_t) nbytes);
}

static void put_u16(writer_t *w, uint64_t v) { put_uint(w, v, 2); }
static void put_u32(writer_t *w, uint64_t v) { put_uint(w, v, 4); }
static void put_u64(writer_t *w, uint64_t v) { put_uint(w, v, 8); }

/* Bits per coefficient for a modulus: enough for mod - 1. */
static unsigned coeff_bits(ulong mod) {
    unsigned bits = 0;

    for (ulong m = mod - 1; m != 0; m >>= 1) {
        bits++;
    }
    return bits ? bits : 1;
}

static void put_poly(writer_t *w, const nmod_poly_struct *poly) {
    const ulong mod = poly->mod.n;
    const unsigned bits = coeff_bits(mod);
    const slong len = poly->length;
    uint8_t tmp[256];
    size_t t = 0;
    uint64_t acc = 0;
    unsigned nacc = 0;

    if (w->err) {
        return;
    }
    if (len > DEGREE_N) {
        w->err = EVMLR_SER_ERR_ARG;
        return;
    }
    put_u16(w, (uint64_t) len);
    for (slong i = 0; i < len; i++) {
        acc |= (uint64_t) poly->coeffs[i] << nacc;
        nacc += bits;
        while (nacc >= 8) {
            tmp[t++] = (uint8_t) acc;
            acc >>= 8;
            nacc -= 8;
            if (t == sizeof(tmp)) {
                put_bytes(w, tmp, t);
                t = 0;
            }
        }
    }
    if (nacc > 0) {
        tmp[t++] = (uint8_t) acc;
    }
    put_bytes(w, tmp, t);
}

/* rows x cols entries of mat starting at (r0, c0), with their dimensions. */
static void put_mat_block(writer_t *w, const nmod_poly_mat_t mat, slong r0, slong c0, slong rows, slong cols) {
    put_u32(w, (uint64_t) rows);
    put_u32(w, (uint64_t) cols);
    for (slong i = 0; i < rows && !w->err; i++) {
        for (slong j = 0; j < cols && !w->err; j++) {
            put_poly(w, nmod_poly_mat_entry(mat, r0 + i, c0 + j));
        }
    }
}

static void put_mat(writer_t *w, const nmod_poly_mat_t mat) {
    put_mat_block(w, mat, 0, 0, nmod_poly_mat_nrows(mat), nmod_poly_mat_ncols(mat));
}

static void header_bytes(uint8_t h[EVMLR_SER_HEADER_LEN], uint8_t type) {
    uint8_t *p = h;

    memcpy(p, MAGIC, 4);
    p += 4;
    *p++ = EVMLR_SER_VERSION;
    *p++ = type;
    *p++ = (uint8_t) MOD_Q;
    *p++ = (uint8_t) (MOD_Q >> 8);
    *p++ = (uint8_t) DEGREE_N;
    *p++ = (uint8_t) (DEGREE_N >> 8);
    *p++ = (uint8_t) ALPHA;
    *p++ = (uint8_t) (ALPHA >> 8);
    *p++ = (uint8_t) NONZERO;
    *p++ = (uint8_t) (NONZERO >> 8);
    *p++ = K_LWE;
    *p++ = K_LWR;
    *p++ = K_SIS;
    *p++ = ETA;
    *p++ = ZETA;
    *p++ = LOG_Q_CEIL;
}

static void put_header(writer_t *w, uint8_t type) {
    uint8_t h[EVMLR_SER_HEADER_LEN];

    header_bytes(h, type);
    put_bytes(w, h, sizeof(h));
}

/* ---------------------------------------------------------------------------
 * Input
 * ------------------------------------------------------------------------- */

typedef struct {
    const uint8_t *p;
    size_t left;
    int err;
} reader_t;

static int need(reader_t *r, size_t n) {
    if (r->err) {
        return 0;
    }
    if (r->left < n) {
        r->err = EVMLR_SER_ERR_TRUNCATED;
        return 0;
    }
    return 1;
}

static uint64_t get_uint(reader_t *r, int nbytes) {
    uint64_t v = 0;

    if (!need(r, (size_t) nbytes)) {
        return 0;
    }
    for (int i = 0; i < nbytes; i++) {
        v |= (uint64_t) r->p[i] << (8 * i);
    }
    r->p += nbytes;
    r->left -= (size_t) nbytes;
    return v;
}

static uint64_t get_u16(reader_t *r) { return get_uint(r, 2); }
static uint64_t get_u32(reader_t *r) { return get_uint(r, 4); }
static uint64_t get_u64(reader_t *r) { return get_uint(r, 8); }

/* Reads a polynomial into poly, which carries the expected modulus. */
static void get_poly(reader_t *r, nmod_poly_struct *poly) {
    const ulong mod = poly->mod.n;
    const unsigned bits = coeff_bits(mod);
    const uint64_t mask = (bits >= 64) ? UINT64_MAX : ((uint64_t) 1 << bits) - 1;
    uint64_t len, acc = 0;
    size_t nbytes, k = 0;
    unsigned nacc = 0;

    len = get_u16(r);
    if (r->err) {
        return;
    }
    if (len > DEGREE_N) {
        r->err = EVMLR_SER_ERR_SHAPE;
        return;
    }
    nbytes = (size_t) ((len * bits + 7) / 8);
    if (!need(r, nbytes)) {
        return;
    }

    nmod_poly_fit_length(poly, (slong) len);
    for (uint64_t i = 0; i < len; i++) {
        uint64_t c;

        while (nacc < bits) {
            acc |= (uint64_t) r->p[k++] << nacc;
            nacc += 8;
        }
        c = acc & mask;
        acc >>= bits;
        nacc -= bits;
        if (c >= mod) {
            r->err = EVMLR_SER_ERR_VALUE;
            break;
        }
        poly->coeffs[i] = (ulong) c;
    }
    /* Canonical: zero padding bits, nonzero leading coefficient. */
    if (!r->err && (acc != 0 || (len > 0 && poly->coeffs[len - 1] == 0))) {
        r->err = EVMLR_SER_ERR_VALUE;
    }
    if (r->err) {
        poly->length = 0;
        return;
    }
    poly->length = (slong) len;
    r->p += nbytes;
    r->left -= nbytes;
}

/* Reads rows x cols entries into mat at (r0, c0); the dimensions written in
 * the input must be exactly rows x cols. */
static void get_mat_block(reader_t *r, nmod_poly_mat_t mat, slong r0, slong c0, slong rows, slong cols) {
    uint64_t rr = get_u32(r), cc = get_u32(r);

    if (r->err) {
        return;
    }
    if (rr != (uint64_t) rows || cc != (uint64_t) cols) {
        r->err = EVMLR_SER_ERR_SHAPE;
        return;
    }
    for (slong i = 0; i < rows && !r->err; i++) {
        for (slong j = 0; j < cols && !r->err; j++) {
            get_poly(r, nmod_poly_mat_entry(mat, r0 + i, c0 + j));
        }
    }
}

static void get_mat(reader_t *r, nmod_poly_mat_t mat) {
    get_mat_block(r, mat, 0, 0, nmod_poly_mat_nrows(mat), nmod_poly_mat_ncols(mat));
}

static int check_header(const uint8_t *in, size_t len, int type) {
    uint8_t h[EVMLR_SER_HEADER_LEN];

    if (in == NULL || len < EVMLR_SER_HEADER_LEN) {
        return (in != NULL && len >= 4 && memcmp(in, MAGIC, 4) != 0)
            ? EVMLR_SER_ERR_MAGIC : EVMLR_SER_ERR_TRUNCATED;
    }
    header_bytes(h, in[5]);
    if (memcmp(in, h, 4) != 0) return EVMLR_SER_ERR_MAGIC;
    if (in[4] != h[4]) return EVMLR_SER_ERR_VERSION;
    if (memcmp(in + 6, h + 6, EVMLR_SER_HEADER_LEN - 6) != 0) return EVMLR_SER_ERR_PARAMS;
    if (type >= 0 && in[5] != type) return EVMLR_SER_ERR_TYPE;
    return EVMLR_SER_OK;
}

static reader_t open_reader(const uint8_t *in, size_t len, int type) {
    reader_t r = { in, len, check_header(in, len, type) };

    if (!r.err) {
        r.p += EVMLR_SER_HEADER_LEN;
        r.left -= EVMLR_SER_HEADER_LEN;
    }
    return r;
}

static int close_reader(reader_t *r) {
    if (!r->err && r->left != 0) {
        r->err = EVMLR_SER_ERR_TRAILING;
    }
    return r->err;
}

int evmlr_ser_peek_type(const uint8_t *in, size_t len) {
    int st = check_header(in, len, -1);

    return st ? -st : in[5];
}

const char *evmlr_ser_strerror(int status) {
    switch (status) {
        case EVMLR_SER_OK:             return "ok";
        case EVMLR_SER_ERR_TRUNCATED:  return "input is truncated";
        case EVMLR_SER_ERR_MAGIC:      return "not an evmlr encoding";
        case EVMLR_SER_ERR_VERSION:    return "unsupported format version";
        case EVMLR_SER_ERR_TYPE:       return "encoding of a different object type";
        case EVMLR_SER_ERR_PARAMS:     return "encoded with different scheme parameters";
        case EVMLR_SER_ERR_SHAPE:      return "unexpected dimension";
        case EVMLR_SER_ERR_VALUE:      return "unreduced or non-canonical value";
        case EVMLR_SER_ERR_TRAILING:   return "trailing bytes after the object";
        case EVMLR_SER_ERR_NOMEM:      return "out of memory";
        case EVMLR_SER_ERR_ARG:        return "invalid argument";
        case EVMLR_SER_ERR_KEY:        return "key does not match";
        case EVMLR_SER_ERR_INVALID:    return "does not verify";
        default:                       return "unknown error";
    }
}

/* ---------------------------------------------------------------------------
 * Proof components
 * ------------------------------------------------------------------------- */

static void put_lin_proof(writer_t *w, const evmlr_lin_proof_struct *p) {
    put_mat(w, p->w);
    put_mat(w, p->z_1);
    put_mat(w, p->z_2);
}

static void get_lin_proof(reader_t *r, evmlr_lin_proof_struct *p) {
    get_mat(r, p->w);
    get_mat(r, p->z_1);
    get_mat(r, p->z_2);
}

static void put_bin_proof(writer_t *w, const evmlr_bin_proof_struct *p) {
    put_mat(w, p->c_g->c);
    put_mat(w, p->c_g1->c);
    put_mat(w, p->w);
    put_mat(w, p->w_g);
    put_mat(w, p->w_g1);
    put_poly(w, p->h);
    put_poly(w, p->v);
    put_mat(w, p->z_m);
    put_mat(w, p->z_r);
    put_poly(w, p->z_g);
    put_mat(w, p->z_rg);
    put_poly(w, p->z_g1);
    put_mat(w, p->z_rg1);
}

static void get_bin_proof(reader_t *r, evmlr_bin_proof_struct *p) {
    get_mat(r, p->c_g->c);
    get_mat(r, p->c_g1->c);
    get_mat(r, p->w);
    get_mat(r, p->w_g);
    get_mat(r, p->w_g1);
    get_poly(r, p->h);
    get_poly(r, p->v);
    get_mat(r, p->z_m);
    get_mat(r, p->z_r);
    get_poly(r, p->z_g);
    get_mat(r, p->z_rg);
    get_poly(r, p->z_g1);
    get_mat(r, p->z_rg1);
}

/* ---------------------------------------------------------------------------
 * Sizes fixed by the protocol
 * ------------------------------------------------------------------------- */

static uint64_t voting_L_max(uint64_t k) {
    return (k - 1) * K_LWR * (K_LWE + 1) + 1;
}

/* ---------------------------------------------------------------------------
 * Public parameters
 * ------------------------------------------------------------------------- */

int evmlr_serialize_pp(evmlr_buf_t *out, const evmlr_voting_pp_t pp) {
    writer_t w = { out, 0 };

    if (pp->k < 1 || pp->L_max != (slong) voting_L_max((uint64_t) pp->k)) {
        return EVMLR_SER_ERR_ARG;
    }
    put_header(&w, EVMLR_SER_PP);
    put_bytes(&w, pp->election_id, EVMLR_ELECTION_ID_LEN);
    put_u32(&w, (uint64_t) pp->k);
    for (slong j = 0; j < pp->k; j++) {
        put_mat(&w, pp->pks[j]->A);
        put_mat(&w, pp->pks[j]->t);
    }
    /* The OTSE matrices are derived from the election id and k (see
     * evmlr_voting_pp_init), so they are not encoded. */
    return w.err;
}

static void pk_init(evmlr_mlpke_pk_struct *pk) {
    nmod_poly_mat_init(pk->A, K_LWE, K_LWE, MOD_Q);
    nmod_poly_mat_init(pk->t, K_LWE, 1, MOD_Q);
}

static void pk_clear(evmlr_mlpke_pk_struct *pk) {
    nmod_poly_mat_clear(pk->A);
    nmod_poly_mat_clear(pk->t);
}

int evmlr_deserialize_pp(evmlr_voting_pp_t pp, const uint8_t *in, size_t len) {
    reader_t r = open_reader(in, len, EVMLR_SER_PP);
    uint8_t id[EVMLR_ELECTION_ID_LEN];
    evmlr_mlpke_pk_t *pks;
    uint64_t k;

    if (need(&r, sizeof(id))) {
        memcpy(id, r.p, sizeof(id));
        r.p += sizeof(id);
        r.left -= sizeof(id);
    }
    k = get_u32(&r);
    if (r.err) {
        return r.err;
    }
    if (k < 1 || k > MAX_SERVERS) {
        return EVMLR_SER_ERR_SHAPE;
    }
    /* Every polynomial takes at least its 2-byte length. */
    if (k * (K_LWE * K_LWE + K_LWE) > r.left / 2) {
        return EVMLR_SER_ERR_TRUNCATED;
    }
    pks = (evmlr_mlpke_pk_t *) malloc((size_t) k * sizeof(evmlr_mlpke_pk_t));
    if (pks == NULL) {
        return EVMLR_SER_ERR_NOMEM;
    }
    for (uint64_t j = 0; j < k; j++) {
        pk_init(pks[j]);
    }
    for (uint64_t j = 0; j < k && !r.err; j++) {
        get_mat(&r, pks[j]->A);
        get_mat(&r, pks[j]->t);
    }
    if (!close_reader(&r)) {
        evmlr_voting_pp_init(pp, id, (slong) k, (const evmlr_mlpke_pk_t *) pks);
    }
    for (uint64_t j = 0; j < k; j++) {
        pk_clear(pks[j]);
    }
    free(pks);
    return r.err;
}

/* ---------------------------------------------------------------------------
 * A server's keys
 * ------------------------------------------------------------------------- */

int evmlr_serialize_pk(evmlr_buf_t *out, const evmlr_mlpke_pk_t pk) {
    writer_t w = { out, 0 };

    put_header(&w, EVMLR_SER_PK);
    put_mat(&w, pk->A);
    put_mat(&w, pk->t);
    return w.err;
}

int evmlr_deserialize_pk(evmlr_mlpke_pk_t pk, const uint8_t *in, size_t len) {
    reader_t r = open_reader(in, len, EVMLR_SER_PK);

    if (r.err) {
        return r.err;
    }
    pk_init(pk);
    get_mat(&r, pk->A);
    get_mat(&r, pk->t);
    if (close_reader(&r)) {
        pk_clear(pk);
    }
    return r.err;
}

int evmlr_serialize_keypair(evmlr_buf_t *out, const evmlr_mlpke_keypair_t kp) {
    writer_t w = { out, 0 };

    put_header(&w, EVMLR_SER_KEYPAIR);
    put_mat(&w, kp->sk->s);
    put_mat(&w, kp->pk->A);
    put_mat(&w, kp->pk->t);
    return w.err;
}

/* Whether s is the secret of (A, t): e = t - A s, and s itself, must be
 * small, as keygen samples them from B_ETA. */
static int keypair_consistent(const evmlr_mlpke_keypair_struct *kp) {
    evmlr_mlpke_ctx_t ctx;
    nmod_poly_mat_t e;
    int ok = 1;

    evmlr_mlpke_ctx_init(ctx);
    nmod_poly_mat_init(e, K_LWE, 1, MOD_Q);
    nmod_poly_mat_mulmod(e, kp->pk->A, kp->sk->s, ctx->cyclo_poly);
    nmod_poly_mat_sub(e, kp->pk->t, e);
    for (slong i = 0; i < K_LWE && ok; i++) {
        const nmod_poly_struct *ei = nmod_poly_mat_entry(e, i, 0), *si = nmod_poly_mat_entry(kp->sk->s, i, 0);
        for (slong c = 0; c < DEGREE_N && ok; c++) {
            ulong x = nmod_poly_get_coeff_ui(ei, c), y = nmod_poly_get_coeff_ui(si, c);
            ok = (x <= ETA || x >= MOD_Q - ETA) && (y <= ETA || y >= MOD_Q - ETA);
        }
    }
    nmod_poly_mat_clear(e);
    evmlr_mlpke_ctx_clear(ctx);
    return ok;
}

int evmlr_deserialize_keypair(evmlr_mlpke_keypair_t kp, const uint8_t *in, size_t len) {
    reader_t r = open_reader(in, len, EVMLR_SER_KEYPAIR);

    if (r.err) {
        return r.err;
    }
    nmod_poly_mat_init(kp->sk->s, K_LWE, 1, MOD_Q);
    pk_init(kp->pk);
    get_mat(&r, kp->sk->s);
    get_mat(&r, kp->pk->A);
    get_mat(&r, kp->pk->t);
    if (!close_reader(&r) && !keypair_consistent(kp)) {
        r.err = EVMLR_SER_ERR_KEY;
    }
    if (r.err) {
        evmlr_mlpke_keypair_clear(kp);
    }
    return r.err;
}

/* ---------------------------------------------------------------------------
 * Secret keys
 * ------------------------------------------------------------------------- */

int evmlr_serialize_sk(evmlr_buf_t *out, const evmlr_voting_sk_t sk) {
    writer_t w = { out, 0 };

    if (sk->k < 1) {
        return EVMLR_SER_ERR_ARG;
    }
    put_header(&w, EVMLR_SER_SK);
    put_u32(&w, (uint64_t) sk->k);
    for (slong j = 0; j < sk->k; j++) {
        put_mat(&w, sk->sks[j]->s);
    }
    return w.err;
}

int evmlr_deserialize_sk(evmlr_voting_sk_t sk, const uint8_t *in, size_t len) {
    reader_t r = open_reader(in, len, EVMLR_SER_SK);
    uint64_t k = get_u32(&r);

    if (r.err) {
        return r.err;
    }
    if (k < 1 || k > MAX_SERVERS) {
        return EVMLR_SER_ERR_SHAPE;
    }
    if (k * K_LWE > r.left / 2) {
        return EVMLR_SER_ERR_TRUNCATED;
    }
    sk->sks = (evmlr_mlpke_sk_t *) malloc((size_t) k * sizeof(evmlr_mlpke_sk_t));
    if (sk->sks == NULL) {
        return EVMLR_SER_ERR_NOMEM;
    }
    sk->k = (slong) k;
    for (uint64_t j = 0; j < k; j++) {
        nmod_poly_mat_init(sk->sks[j]->s, K_LWE, 1, MOD_Q);
    }
    for (uint64_t j = 0; j < k && !r.err; j++) {
        get_mat(&r, sk->sks[j]->s);
    }
    if (close_reader(&r)) {
        evmlr_voting_sk_clear(sk);
        return r.err;
    }
    return EVMLR_SER_OK;
}

/* ---------------------------------------------------------------------------
 * Submissions
 * ------------------------------------------------------------------------- */

int evmlr_serialize_submission(evmlr_buf_t *out, const evmlr_voting_submission_t sub,
                               const evmlr_voting_pp_t pp) {
    writer_t w = { out, 0 };
    const evmlr_voting_layer_struct *layer = sub->ballot->layer;

    if (layer->otse_cipher->c->r != pp->L_max) {
        return EVMLR_SER_ERR_ARG;
    }
    put_header(&w, EVMLR_SER_SUBMISSION);
    put_u64(&w, (uint64_t) sub->voter_id);
    for (int i = 0; i < K_LWR; i++) {
        put_mat(&w, layer->enc_cipher[i]->uT);
        put_poly(&w, layer->enc_cipher[i]->v);
    }
    put_mat(&w, layer->otse_cipher->c);
    for (slong j = 0; j < pp->k; j++) {
        put_lin_proof(&w, &sub->proofs[j]->proof.proof_hpke);
        put_lin_proof(&w, &sub->proofs[j]->proof.proof_otse);
    }
    return w.err;
}

int evmlr_deserialize_submission(evmlr_voting_submission_t sub, const evmlr_voting_pp_t pp,
                                 const uint8_t *in, size_t len) {
    reader_t r = open_reader(in, len, EVMLR_SER_SUBMISSION);
    evmlr_voting_layer_struct *layer;
    int64_t voter_id;

    voter_id = (int64_t) get_u64(&r);
    if (r.err) {
        return r.err;
    }
    evmlr_voting_submission_init(sub, pp);
    sub->voter_id = voter_id;
    layer = sub->ballot->layer;
    for (int i = 0; i < K_LWR; i++) {
        get_mat(&r, layer->enc_cipher[i]->uT);
        get_poly(&r, layer->enc_cipher[i]->v);
    }
    get_mat(&r, layer->otse_cipher->c);
    for (slong j = 0; j < pp->k && !r.err; j++) {
        get_lin_proof(&r, &sub->proofs[j]->proof.proof_hpke);
        get_lin_proof(&r, &sub->proofs[j]->proof.proof_otse);
    }
    if (close_reader(&r)) {
        evmlr_voting_submission_clear(sub, pp);
        return r.err;
    }
    return EVMLR_SER_OK;
}

/* ---------------------------------------------------------------------------
 * Shuffle proofs, statements and counting transcripts
 * ------------------------------------------------------------------------- */

static void put_statement(writer_t *w, const evmlr_shuffle_pp_struct *stmt, slong N) {
    put_mat(w, stmt->D->c);
    for (slong i = 0; i < N; i++) {
        put_mat(w, stmt->c_hat[i]);
        put_mat(w, stmt->c_star[i]->c);
    }
}

/* Allocates stmt for ctx's dimensions and reads it; frees it again on
 * failure. Returns r->err. */
static int read_statement(reader_t *r, evmlr_shuffle_pp_struct *stmt, const evmlr_shuffle_ctx_struct *ctx) {
    const slong n = ctx->N;

    if (r->err) {
        return r->err;
    }
    stmt->c_hat = (nmod_poly_mat_t *) malloc((size_t) n * sizeof(nmod_poly_mat_t));
    stmt->c_star = (evmlr_otse_ciphertext_t *) malloc((size_t) n * sizeof(evmlr_otse_ciphertext_t));
    if (stmt->c_hat == NULL || stmt->c_star == NULL) {
        free(stmt->c_hat);
        free(stmt->c_star);
        return r->err = EVMLR_SER_ERR_NOMEM;
    }
    evmlr_commit_init(stmt->D);
    for (slong i = 0; i < n; i++) {
        nmod_poly_mat_init(stmt->c_hat[i], ctx->L, 1, MOD_Q);
        evmlr_otse_ciphertext_init(stmt->c_star[i], ctx->L);
    }
    get_mat(r, stmt->D->c);
    for (slong i = 0; i < n && !r->err; i++) {
        get_mat(r, stmt->c_hat[i]);
        get_mat(r, stmt->c_star[i]->c);
    }
    if (r->err) {
        evmlr_shuffle_pp_clear(stmt, ctx);
    }
    return r->err;
}

static void put_shuffle_proof(writer_t *w, const evmlr_shuffle_proof_struct *proof) {
    put_mat(w, proof->u);
    put_mat(w, proof->P->c);
    put_mat(w, proof->W->c);
    put_bin_proof(w, proof->proof_D);
    put_bin_proof(w, proof->proof_P);
    put_bin_proof(w, proof->proof_W);
    put_lin_proof(w, proof->proof_u);
}

static int read_shuffle_proof(reader_t *r, evmlr_shuffle_proof_struct *proof, const evmlr_shuffle_ctx_struct *ctx) {
    if (r->err) {
        return r->err;
    }
    /* evmlr_proof_init leaves P and W to the prover; prove initialises them. */
    evmlr_proof_init(proof, ctx);
    evmlr_commit_init(proof->P);
    evmlr_commit_init(proof->W);
    get_mat(r, proof->u);
    get_mat(r, proof->P->c);
    get_mat(r, proof->W->c);
    get_bin_proof(r, proof->proof_D);
    get_bin_proof(r, proof->proof_P);
    get_bin_proof(r, proof->proof_W);
    get_lin_proof(r, proof->proof_u);
    if (r->err) {
        evmlr_proof_clear(proof, ctx);
    }
    return r->err;
}

/* The dimensions of a shuffle of N messages of length L, without the
 * matrices: what proof_init, statements and the clear functions read. */
static void shape_ctx(evmlr_shuffle_ctx_struct *shape, slong N, slong L) {
    memset(shape, 0, sizeof(*shape));
    shape->N = N;
    shape->L = L;
    shape->com_ctx->N = evmlr_shuffle_com_len(N, L);
}

/* How many polynomials a statement and a proof for N messages of length L
 * hold; each takes at least two bytes of input. */
static uint64_t shuffle_polys(uint64_t N, uint64_t L) {
    const uint64_t com_N = (uint64_t) evmlr_shuffle_com_len((slong) N, (slong) L);
    const uint64_t bin = 5 * K_SIS + 4 + com_N + 3 * 2 * K_SIS;
    const uint64_t lin = 2 * N + N * (LOG_Q_CEIL + 1 + 2 * ETA * (K_LWE + L));

    return K_SIS + 2 * N * L                     /* statement */
         + (N - 1) + 2 * K_SIS + 3 * bin + lin;  /* proof */
}

int evmlr_serialize_shuffle_proof(evmlr_buf_t *out, const evmlr_shuffle_proof_t proof,
                                  const evmlr_shuffle_ctx_t ctx) {
    writer_t w = { out, 0 };

    put_header(&w, EVMLR_SER_SHUFFLE_PROOF);
    put_u32(&w, (uint64_t) ctx->N);
    put_shuffle_proof(&w, proof);
    return w.err;
}

int evmlr_deserialize_shuffle_proof(evmlr_shuffle_proof_t proof, const evmlr_shuffle_ctx_t ctx,
                                    const uint8_t *in, size_t len) {
    reader_t r = open_reader(in, len, EVMLR_SER_SHUFFLE_PROOF);
    uint64_t N = get_u32(&r);

    if (r.err) {
        return r.err;
    }
    if (N != (uint64_t) ctx->N) {
        return EVMLR_SER_ERR_SHAPE;
    }
    if (read_shuffle_proof(&r, proof, ctx) || close_reader(&r)) {
        if (r.err == EVMLR_SER_ERR_TRAILING) {
            evmlr_proof_clear(proof, ctx);
        }
        return r.err;
    }
    return EVMLR_SER_OK;
}

int evmlr_serialize_shuffle_statement(evmlr_buf_t *out, const evmlr_shuffle_pp_t stmt,
                                      const evmlr_shuffle_ctx_t ctx) {
    writer_t w = { out, 0 };

    put_header(&w, EVMLR_SER_SHUFFLE_STMT);
    put_u32(&w, (uint64_t) ctx->N);
    put_statement(&w, stmt, ctx->N);
    return w.err;
}

int evmlr_deserialize_shuffle_statement(evmlr_shuffle_pp_t stmt, const evmlr_shuffle_ctx_t ctx,
                                        const uint8_t *in, size_t len) {
    reader_t r = open_reader(in, len, EVMLR_SER_SHUFFLE_STMT);
    uint64_t N = get_u32(&r);

    if (r.err) {
        return r.err;
    }
    if (N != (uint64_t) ctx->N) {
        return EVMLR_SER_ERR_SHAPE;
    }
    if (read_statement(&r, stmt, ctx) || close_reader(&r)) {
        if (r.err == EVMLR_SER_ERR_TRAILING) {
            evmlr_shuffle_pp_clear(stmt, ctx);
        }
        return r.err;
    }
    return EVMLR_SER_OK;
}

/* ---------------------------------------------------------------------------
 * Batches, mix records and transcripts
 * ------------------------------------------------------------------------- */

/* Polynomials in a proof of encryption of a layer with OTSE length L. */
static uint64_t enc_proof_polys(uint64_t L) {
    const uint64_t hdr = K_LWR * (K_LWE + 1);

    return 3 * hdr + 2 * L + 2 * ETA * (K_LWE + L);
}

int evmlr_serialize_batch(evmlr_buf_t *out, const evmlr_voting_batch_t batch, const evmlr_voting_pp_t pp) {
    writer_t w = { out, 0 };
    const slong k = pp->k, nprf = k - batch->layer;

    if (batch->k != k || batch->layer < 0 || batch->layer > k || batch->n < 0
        || (batch->layer == 0 && batch->n > 0 && batch->voter_ids == NULL)) {
        return EVMLR_SER_ERR_ARG;
    }
    put_header(&w, EVMLR_SER_BATCH);
    put_u32(&w, (uint64_t) k);
    put_u32(&w, (uint64_t) batch->layer);
    put_u32(&w, (uint64_t) batch->n);
    for (slong i = 0; i < batch->n && !w.err; i++) {
        if (batch->layer == 0) {
            put_u64(&w, (uint64_t) batch->voter_ids[i]);
        }
        put_mat(&w, batch->msgs[i]);
        for (slong l = 0; l < nprf; l++) {
            put_lin_proof(&w, &batch->proofs[i * nprf + l].proof_hpke);
            put_lin_proof(&w, &batch->proofs[i * nprf + l].proof_otse);
        }
    }
    return w.err;
}

int evmlr_deserialize_batch(evmlr_voting_batch_t batch, const evmlr_voting_pp_t pp, const uint8_t *in, size_t len) {
    reader_t r = open_reader(in, len, EVMLR_SER_BATCH);
    uint64_t k, layer, n, per_item;
    slong nprf;

    k = get_u32(&r);
    layer = get_u32(&r);
    n = get_u32(&r);
    if (r.err) {
        return r.err;
    }
    if (k != (uint64_t) pp->k || layer > k || n > MAX_MESSAGES) {
        return EVMLR_SER_ERR_SHAPE;
    }
    nprf = (slong) (k - layer);
    per_item = (uint64_t) evmlr_voting_msg_len((slong) k, (slong) layer);
    for (slong l = 0; l < nprf; l++) {
        per_item += enc_proof_polys((uint64_t) evmlr_voting_otse_len((slong) k, (slong) layer + l));
    }
    if (n * per_item > r.left / 2) {
        return EVMLR_SER_ERR_TRUNCATED;
    }
    if (evmlr_voting_batch_init(batch, (slong) k, (slong) layer, (slong) n) != EVMLR_VOTING_OK) {
        return EVMLR_SER_ERR_NOMEM;
    }
    for (slong i = 0; i < (slong) n && !r.err; i++) {
        if (layer == 0) {
            batch->voter_ids[i] = (int64_t) get_u64(&r);
        }
        get_mat(&r, batch->msgs[i]);
        for (slong l = 0; l < nprf; l++) {
            get_lin_proof(&r, &batch->proofs[i * nprf + l].proof_hpke);
            get_lin_proof(&r, &batch->proofs[i * nprf + l].proof_otse);
        }
    }
    if (close_reader(&r)) {
        evmlr_voting_batch_clear(batch);
    }
    return r.err;
}

static void put_mix_body(writer_t *w, const evmlr_voting_mix_struct *mix) {
    put_u32(w, (uint64_t) mix->k);
    put_u32(w, (uint64_t) mix->j);
    put_u32(w, (uint64_t) mix->n_in);
    put_u32(w, (uint64_t) mix->n_rejected);
    for (slong r = 0; r < mix->n_rejected; r++) {
        put_u32(w, (uint64_t) mix->rejected[r]);
        put_uint(w, mix->reasons[r], 1);
    }
    if (mix->N > 0) {
        put_statement(w, mix->statement, evmlr_voting_shuffle_len(mix->N));
        put_shuffle_proof(w, mix->proof);
    }
}

static int check_mix(const evmlr_voting_mix_struct *mix, const evmlr_voting_pp_t pp) {
    if (mix->k != pp->k || mix->j < 0 || mix->j >= pp->k || mix->n_in < 0 || mix->n_rejected < 0
        || mix->N != mix->n_in - mix->n_rejected
        || (mix->n_rejected > 0 && (mix->rejected == NULL || mix->reasons == NULL))) {
        return EVMLR_SER_ERR_ARG;
    }
    return EVMLR_SER_OK;
}

/* Reads a mix record into mix (all fields set, even on failure, so that
 * evmlr_voting_mix_clear can free what was read). */
static int read_mix_body(reader_t *r, evmlr_voting_mix_struct *mix, const evmlr_voting_pp_t pp) {
    uint64_t k, j, n_in, n_rej, prev = 0;
    evmlr_shuffle_ctx_struct shape;

    memset(mix, 0, sizeof(*mix));
    k = get_u32(r);
    j = get_u32(r);
    n_in = get_u32(r);
    n_rej = get_u32(r);
    if (r->err) {
        return r->err;
    }
    if (k != (uint64_t) pp->k || j >= k || n_in > MAX_MESSAGES || n_rej > n_in) {
        return r->err = EVMLR_SER_ERR_SHAPE;
    }
    if (n_rej > r->left / 5) {
        return r->err = EVMLR_SER_ERR_TRUNCATED;
    }
    mix->k = (slong) k;
    mix->j = (slong) j;
    mix->n_in = (slong) n_in;
    if (n_rej > 0) {
        mix->rejected = (slong *) malloc((size_t) n_rej * sizeof(slong));
        mix->reasons = (uint8_t *) malloc((size_t) n_rej);
        if (mix->rejected == NULL || mix->reasons == NULL) {
            return r->err = EVMLR_SER_ERR_NOMEM;
        }
    }
    for (uint64_t i = 0; i < n_rej && !r->err; i++) {
        uint64_t idx = get_u32(r), why = get_uint(r, 1);

        if (r->err) {
            break;
        }
        /* Strictly increasing, in range, with a known reason: one encoding. */
        if (idx >= n_in || (i > 0 && idx <= prev)
            || (why != EVMLR_VOTING_REJECT_HPKE && why != EVMLR_VOTING_REJECT_OTSE)) {
            r->err = EVMLR_SER_ERR_VALUE;
            break;
        }
        prev = idx;
        mix->rejected[i] = (slong) idx;
        mix->reasons[i] = (uint8_t) why;
        mix->n_rejected = (slong) i + 1;
    }
    if (r->err) {
        return r->err;
    }
    mix->N = (slong) (n_in - n_rej);
    if (mix->N == 0) {
        return EVMLR_SER_OK;
    }
    if (shuffle_polys((uint64_t) evmlr_voting_shuffle_len(mix->N), (uint64_t) pp->L_max) > r->left / 2) {
        mix->N = 0;
        return r->err = EVMLR_SER_ERR_TRUNCATED;
    }
    shape_ctx(&shape, evmlr_voting_shuffle_len(mix->N), pp->L_max);
    if (read_statement(r, mix->statement, &shape)) {
        mix->N = 0;
        return r->err;
    }
    if (read_shuffle_proof(r, mix->proof, &shape)) {
        evmlr_shuffle_pp_clear(mix->statement, &shape);
        mix->N = 0;
        return r->err;
    }
    return EVMLR_SER_OK;
}

int evmlr_serialize_mix(evmlr_buf_t *out, const evmlr_voting_mix_t mix, const evmlr_voting_pp_t pp) {
    writer_t w = { out, 0 };
    int st = check_mix(mix, pp);

    if (st) {
        return st;
    }
    put_header(&w, EVMLR_SER_MIX);
    put_mix_body(&w, mix);
    return w.err;
}

int evmlr_deserialize_mix(evmlr_voting_mix_t mix, const evmlr_voting_pp_t pp, const uint8_t *in, size_t len) {
    reader_t r = open_reader(in, len, EVMLR_SER_MIX);

    if (r.err) {
        return r.err;
    }
    if (read_mix_body(&r, mix, pp) || close_reader(&r)) {
        evmlr_voting_mix_clear(mix);
    }
    return r.err;
}

int evmlr_serialize_transcript(evmlr_buf_t *out, const evmlr_voting_transcript_t tr,
                               const evmlr_voting_pp_t pp) {
    writer_t w = { out, 0 };

    if (tr->k != pp->k || tr->mixes == NULL) {
        return EVMLR_SER_ERR_ARG;
    }
    for (slong j = 0; j < tr->k; j++) {
        int st = check_mix(&tr->mixes[j], pp);
        if (st) {
            return st;
        }
    }
    put_header(&w, EVMLR_SER_TRANSCRIPT);
    put_u32(&w, (uint64_t) tr->k);
    for (slong j = 0; j < tr->k; j++) {
        put_mix_body(&w, &tr->mixes[j]);
    }
    return w.err;
}

int evmlr_deserialize_transcript(evmlr_voting_transcript_t tr, const evmlr_voting_pp_t pp,
                                 const uint8_t *in, size_t len) {
    reader_t r = open_reader(in, len, EVMLR_SER_TRANSCRIPT);
    uint64_t k = get_u32(&r);

    if (r.err) {
        return r.err;
    }
    if (k != (uint64_t) pp->k) {
        return EVMLR_SER_ERR_SHAPE;
    }
    tr->k = (slong) k;
    tr->mixes = (evmlr_voting_mix_struct *) calloc((size_t) k, sizeof(evmlr_voting_mix_struct));
    if (tr->mixes == NULL) {
        return EVMLR_SER_ERR_NOMEM;
    }
    for (slong j = 0; j < tr->k && !r.err; j++) {
        read_mix_body(&r, &tr->mixes[j], pp);
    }
    if (close_reader(&r)) {
        evmlr_voting_transcript_clear(tr);
    }
    return r.err;
}

#ifdef MAIN
#include "evmlr_main.h"
#include "evmlr_utils.h"

/* Defined in evmlr_shuffle.c, which keeps it out of its header. */
void evmlr_shuffle_clear_sp(evmlr_shuffle_sp_t sp, const evmlr_shuffle_ctx_t ctx);

/* Serialises with f and checks that it succeeded. */
#define ENCODE(buf, call) (evmlr_buf_init(&(buf)), (call) == EVMLR_SER_OK && !(buf).failed)

static int same_bytes(const evmlr_buf_t *a, const evmlr_buf_t *b) {
    return a->len == b->len && memcmp(a->data, b->data, a->len) == 0;
}

/* Decrypts the outermost layer of a submission and checks its proof of
 * correct encryption, as the first server does during counting. */
static int verify_outer_layer(const evmlr_voting_submission_t sub, const evmlr_voting_pp_t pp,
                              const evmlr_voting_sk_t sk) {
    evmlr_hpke_cipher_t cipher;
    nmod_poly_mat_t msg;
    int ok;

    evmlr_hpke_cipher_init(cipher, pp->L_max);
    for (int i = 0; i < K_LWR; i++) {
        nmod_poly_mat_set(cipher->enc_cipher[i]->uT, sub->ballot->layer->enc_cipher[i]->uT);
        nmod_poly_set(cipher->enc_cipher[i]->v, sub->ballot->layer->enc_cipher[i]->v);
    }
    nmod_poly_mat_set(cipher->otse_cipher->c, sub->ballot->layer->otse_cipher->c);
    nmod_poly_mat_init(msg, pp->L_max, 1, MOD_Q);
    evmlr_hpke_decrypt(msg, cipher, sk->sks[0], pp->hpke_ctx);
    ok = evmlr_voting_verify_ballot(pp, sub, msg);
    nmod_poly_mat_clear(msg);
    evmlr_hpke_cipher_clear(cipher);
    return ok;
}

/* A statement and a witness for a shuffle of ctx->N random messages, as in
 * evmlr_shuffle.c's test. */
static void shuffle_instance(evmlr_shuffle_sp_t sp, evmlr_shuffle_pp_t st, const evmlr_shuffle_ctx_t ctx,
                             flint_rand_t state) {
    slong N = ctx->N, L = ctx->L, row_len = (2 * ETA) * (K_LWE + L);
    evmlr_hpke_keypair_t keypair;
    nmod_poly_mat_t *m = (nmod_poly_mat_t *) malloc(N * sizeof(nmod_poly_mat_t));
    nmod_poly_mat_t d_flat;

    evmlr_hpke_keypair_gen(keypair, state, ctx->hpke_ctx);
    for (slong i = 0; i < N; i++) {
        nmod_poly_mat_init(m[i], L, 1, MOD_Q);
        nmod_poly_mat_randtest(m[i], state, DEGREE_N);
    }
    sp->pi = (size_t *) malloc(N * sizeof(size_t));
    evmlr_utils_new_perm(sp->pi, N);
    st->c_hat = (nmod_poly_mat_t *) malloc(N * sizeof(nmod_poly_mat_t));
    st->c_star = (evmlr_otse_ciphertext_t *) malloc(N * sizeof(evmlr_otse_ciphertext_t));
    sp->d_dagger = (nmod_poly_mat_t *) malloc(N * sizeof(nmod_poly_mat_t));
    for (slong i = 0; i < N; i++) {
        evmlr_hpke_cipher_t cipher;

        nmod_poly_mat_init(st->c_hat[i], L, 1, MOD_Q);
        nmod_poly_mat_set(st->c_hat[i], m[sp->pi[i]]);
        evmlr_hpke_cipher_init(cipher, L);
        nmod_poly_mat_init(sp->d_dagger[i], row_len, 1, MOD_Q);
        evmlr_hpke_encrypt(cipher, sp->d_dagger[i], m[i], keypair->enc_keypair->pk, ctx->hpke_ctx, state);
        evmlr_otse_ciphertext_init(st->c_star[i], L);
        nmod_poly_mat_set(st->c_star[i]->c, cipher->otse_cipher->c);
        evmlr_hpke_cipher_clear(cipher);
    }
    nmod_poly_mat_init(sp->r_D, 2 * K_SIS, 1, MOD_Q);
    evmlr_commit_sample_r(sp->r_D);
    nmod_poly_mat_init(d_flat, ctx->com_ctx->N, 1, MOD_Q);
    for (slong i = 0; i < N; i++) {
        for (slong j = 0; j < row_len; j++) {
            nmod_poly_set(nmod_poly_mat_entry(d_flat, i * row_len + j, 0),
                          nmod_poly_mat_entry(sp->d_dagger[i], j, 0));
        }
    }
    evmlr_commit_init(st->D);
    evmlr_commit(st->D, d_flat, sp->r_D, ctx->com_ctx);

    nmod_poly_mat_clear(d_flat);
    evmlr_hpke_keypair_clear(keypair);
    for (slong i = 0; i < N; i++) {
        nmod_poly_mat_clear(m[i]);
    }
    free(m);
}

static void test(flint_rand_t state) {
    const slong k = 2, n_voters = 3;
    const ulong votes[3] = { 7, 11, 13 };
    evmlr_voting_pp_t pp, pp2;
    evmlr_voting_sk_t sk, sk2;
    evmlr_voting_submission_t subs[3], subs2[3];
    evmlr_buf_t b_pp, b_pp2, b_sk, b_sk2, b_sub[3], b_sub2;
    int have_pp2 = 0, have_sk2 = 0, have_subs2 = 0;

    evmlr_voting_setup(pp, sk, k, state);
    for (slong i = 0; i < n_voters; i++) {
        nmod_poly_mat_t vote;

        nmod_poly_mat_init(vote, 1, 1, MOD_Q);
        nmod_poly_set_coeff_ui(nmod_poly_mat_entry(vote, 0, 0), 0, votes[i]);
        evmlr_voting_submission_init(subs[i], pp);
        evmlr_voting_casting(subs[i], i + 1, vote, pp, state);
        nmod_poly_mat_clear(vote);
    }
    evmlr_buf_init(&b_pp2);
    evmlr_buf_init(&b_sk2);
    evmlr_buf_init(&b_sub2);

    TEST_ONCE("parameters, keys and ballots serialize") {
        TEST_ASSERT(ENCODE(b_pp, evmlr_serialize_pp(&b_pp, pp)), end);
        TEST_ASSERT(ENCODE(b_sk, evmlr_serialize_sk(&b_sk, sk)), end);
        for (slong i = 0; i < n_voters; i++) {
            TEST_ASSERT(ENCODE(b_sub[i], evmlr_serialize_submission(&b_sub[i], subs[i], pp)), end);
        }
        TEST_ASSERT(evmlr_ser_peek_type(b_sub[0].data, b_sub[0].len) == EVMLR_SER_SUBMISSION, end);
    } TEST_END;

    TEST_ONCE("they deserialize and re-encode to the same bytes") {
        TEST_ASSERT(evmlr_deserialize_pp(pp2, b_pp.data, b_pp.len) == EVMLR_SER_OK, end);
        have_pp2 = 1;
        TEST_ASSERT(evmlr_deserialize_sk(sk2, b_sk.data, b_sk.len) == EVMLR_SER_OK, end);
        have_sk2 = 1;
        for (slong i = 0; i < n_voters; i++) {
            TEST_ASSERT(evmlr_deserialize_submission(subs2[i], pp2, b_sub[i].data, b_sub[i].len) == EVMLR_SER_OK, end);
            have_subs2 = (int) i + 1;
        }
        TEST_ASSERT(ENCODE(b_pp2, evmlr_serialize_pp(&b_pp2, pp2)) && same_bytes(&b_pp, &b_pp2), end);
        TEST_ASSERT(ENCODE(b_sk2, evmlr_serialize_sk(&b_sk2, sk2)) && same_bytes(&b_sk, &b_sk2), end);
        for (slong i = 0; i < n_voters; i++) {
            evmlr_buf_clear(&b_sub2);
            TEST_ASSERT(ENCODE(b_sub2, evmlr_serialize_submission(&b_sub2, subs2[i], pp2)), end);
            TEST_ASSERT(same_bytes(&b_sub[i], &b_sub2), end);
        }
    } TEST_END;

    TEST_ONCE("imported ballots verify under imported parameters") {
        for (slong i = 0; i < n_voters; i++) {
            TEST_ASSERT(verify_outer_layer(subs2[i], pp2, sk2) == 1, end);
        }
    } TEST_END;

    TEST_ONCE("a ballot with a changed coefficient fails verification") {
        evmlr_voting_submission_t bad;
        /* The v polynomial of the first ML-PKE ciphertext: after the header,
         * the voter id and the uT matrix. Change its constant coefficient
         * to another reduced value, so the encoding stays valid. */
        nmod_poly_struct *v;
        ulong c;
        int ok;

        TEST_ASSERT(evmlr_deserialize_submission(bad, pp2, b_sub[0].data, b_sub[0].len) == EVMLR_SER_OK, end);
        v = bad->ballot->layer->enc_cipher[0]->v;
        c = nmod_poly_get_coeff_ui(v, 0);
        nmod_poly_set_coeff_ui(v, 0, (c + 1) % MOD_Q);
        ok = verify_outer_layer(bad, pp2, sk2);
        evmlr_voting_submission_clear(bad, pp2);
        TEST_ASSERT(ok == 0, end);
    } TEST_END;

    TEST_ONCE("malformed encodings are rejected") {
        const evmlr_buf_t *b = &b_sub[0];
        uint8_t *mut = (uint8_t *) malloc(b->len + 1);
        evmlr_voting_submission_t tmp;
        evmlr_voting_pp_t pp_tmp;
        int st, all = 1;
        size_t off;

        /* Truncations: every header length, then a spread of prefixes. */
        for (size_t cut = 0; cut <= EVMLR_SER_HEADER_LEN + 16 && all; cut++) {
            all &= evmlr_deserialize_submission(tmp, pp2, b->data, cut) != EVMLR_SER_OK;
        }
        for (size_t cut = EVMLR_SER_HEADER_LEN; cut < b->len && all; cut += 1 + b->len / 97) {
            all &= evmlr_deserialize_submission(tmp, pp2, b->data, cut) == EVMLR_SER_ERR_TRUNCATED;
        }
        all &= evmlr_deserialize_submission(tmp, pp2, b->data, b->len - 1) == EVMLR_SER_ERR_TRUNCATED;
        /* A trailing byte. */
        memcpy(mut, b->data, b->len);
        mut[b->len] = 0;
        all &= evmlr_deserialize_submission(tmp, pp2, mut, b->len + 1) == EVMLR_SER_ERR_TRAILING;
        /* Header fields. */
        mut[0] ^= 1;
        all &= evmlr_deserialize_submission(tmp, pp2, mut, b->len) == EVMLR_SER_ERR_MAGIC;
        mut[0] ^= 1;
        mut[4] ^= 1;
        all &= evmlr_deserialize_submission(tmp, pp2, mut, b->len) == EVMLR_SER_ERR_VERSION;
        mut[4] ^= 1;
        all &= evmlr_deserialize_pp(pp_tmp, mut, b->len) == EVMLR_SER_ERR_TYPE;
        mut[6] ^= 1;
        all &= evmlr_deserialize_submission(tmp, pp2, mut, b->len) == EVMLR_SER_ERR_PARAMS;
        mut[6] ^= 1;
        /* The first coefficient of uT[0][0] set to 4095 >= q. Its bytes
         * follow the header, the voter id, the 8-byte matrix dimensions and
         * the 2-byte polynomial length. */
        off = EVMLR_SER_HEADER_LEN + 8 + 8 + 2;
        mut[off] = 0xff;
        mut[off + 1] |= 0x0f;
        st = evmlr_deserialize_submission(tmp, pp2, mut, b->len);
        all &= st == EVMLR_SER_ERR_VALUE;
        memcpy(mut, b->data, b->len);
        /* A matrix dimension: the row count of uT. */
        mut[EVMLR_SER_HEADER_LEN + 8] ^= 1;
        all &= evmlr_deserialize_submission(tmp, pp2, mut, b->len) == EVMLR_SER_ERR_SHAPE;
        free(mut);
        TEST_ASSERT(all, end);
    } TEST_END;

    TEST_ONCE("shuffle proofs and statements round-trip and verify") {
        evmlr_shuffle_sp_t sp;
        evmlr_shuffle_pp_t st, st2;
        evmlr_shuffle_proof_t proof, proof2;
        evmlr_buf_t b_proof, b_st, b_again;
        evmlr_shuffle_ctx_t ctx, ctx2;
        int ok_rt, ok_ver, ok_bad;

        evmlr_voting_shuffle_ctx_init(ctx, pp, SHUFFLE_N_MSGS);
        evmlr_voting_shuffle_ctx_init(ctx2, pp2, SHUFFLE_N_MSGS);

        /* Proved under the original shuffle context, verified under the one
         * imported with the public parameters. */
        shuffle_instance(sp, st, ctx, state);
        evmlr_proof_init(proof, ctx);
        evmlr_shuffle_prove(proof, sp, st, ctx, state);

        ok_rt = ENCODE(b_proof, evmlr_serialize_shuffle_proof(&b_proof, proof, ctx))
             && ENCODE(b_st, evmlr_serialize_shuffle_statement(&b_st, st, ctx))
             && evmlr_deserialize_shuffle_proof(proof2, ctx2, b_proof.data, b_proof.len) == EVMLR_SER_OK
             && evmlr_deserialize_shuffle_statement(st2, ctx2, b_st.data, b_st.len) == EVMLR_SER_OK;
        ok_ver = ok_bad = 0;
        if (ok_rt) {
            ok_rt = ENCODE(b_again, evmlr_serialize_shuffle_proof(&b_again, proof2, ctx2))
                 && same_bytes(&b_proof, &b_again);
            evmlr_buf_clear(&b_again);
            ok_ver = evmlr_shuffle_verify(proof2, ctx2, st2);
            /* The proof must not verify for another statement. */
            nmod_poly_set_coeff_ui(nmod_poly_mat_entry(st2->c_hat[0], 0, 0), 0,
                                   (nmod_poly_get_coeff_ui(nmod_poly_mat_entry(st2->c_hat[0], 0, 0), 0) + 1) % MOD_Q);
            ok_bad = !evmlr_shuffle_verify(proof2, ctx2, st2);
            evmlr_proof_clear(proof2, ctx2);
            evmlr_shuffle_pp_clear(st2, ctx2);
        }
        evmlr_buf_clear(&b_proof);
        evmlr_buf_clear(&b_st);
        evmlr_proof_clear(proof, ctx);
        evmlr_shuffle_clear_sp(sp, ctx);
        evmlr_shuffle_pp_clear(st, ctx);
        evmlr_shuffle_ctx_clear(ctx);
        evmlr_shuffle_ctx_clear(ctx2);
        TEST_ASSERT(ok_rt && ok_ver && ok_bad, end);
    } TEST_END;

    TEST_ONCE("imported ballots are counted, and the transcript verifies") {
        nmod_poly_mat_t results[3], checked[3];
        evmlr_voting_transcript_t tr, tr2;
        evmlr_voting_submission_struct some[3];
        evmlr_buf_t b_tr, b_tr2;
        slong counted, verified = -1, n_checked = 0;
        int seen[3] = { 0, 0, 0 }, ok, same = 1, ok_rt, bad_stmt, bad_missing, bad_order, bad_trunc;

        counted = evmlr_voting_counting_transcript(results, tr, (const evmlr_voting_submission_t *) subs2,
                                                   n_voters, pp2, sk2, state);
        ok = (counted == n_voters);
        for (slong i = 0; i < counted; i++) {
            ulong v = nmod_poly_get_coeff_ui(nmod_poly_mat_entry(results[i], 0, 0), 0);
            for (int j = 0; j < 3; j++) {
                if (v == votes[j]) seen[j]++;
            }
        }
        ok = ok && seen[0] == 1 && seen[1] == 1 && seen[2] == 1;

        /* A third party: the public parameters, the ballots and the
         * transcript, all from bytes. */
        ok_rt = ENCODE(b_tr, evmlr_serialize_transcript(&b_tr, tr, pp2))
             && evmlr_deserialize_transcript(tr2, pp, b_tr.data, b_tr.len) == EVMLR_SER_OK
             && ENCODE(b_tr2, evmlr_serialize_transcript(&b_tr2, tr2, pp))
             && same_bytes(&b_tr, &b_tr2);
        if (ok_rt) {
            verified = evmlr_voting_verify_counting(checked, tr2, (const evmlr_voting_submission_t *) subs2, n_voters, pp);
            n_checked = verified > 0 ? verified : 0;
            for (slong i = 0; i < n_checked && i < counted; i++) {
                same &= nmod_poly_mat_equal(checked[i], results[i]);
            }

            /* Ballots missing, or in another order than the one shuffled. */
            some[0] = subs2[0][0];
            some[1] = subs2[1][0];
            bad_missing = evmlr_voting_verify_counting(NULL, tr2, (const evmlr_voting_submission_t *) some, 2, pp) == -1;
            some[0] = subs2[1][0];
            some[1] = subs2[0][0];
            some[2] = subs2[2][0];
            bad_order = evmlr_voting_verify_counting(NULL, tr2, (const evmlr_voting_submission_t *) some, 3, pp) == -1;

            /* One output of the first server changed. */
            nmod_poly_struct *c = nmod_poly_mat_entry(tr2->mixes[0].statement->c_hat[0], 0, 0);
            nmod_poly_set_coeff_ui(c, 0, (nmod_poly_get_coeff_ui(c, 0) + 1) % MOD_Q);
            bad_stmt = evmlr_voting_verify_counting(NULL, tr2, (const evmlr_voting_submission_t *) subs2, n_voters, pp) == -1;
            evmlr_voting_transcript_clear(tr2);

            bad_trunc = evmlr_deserialize_transcript(tr2, pp, b_tr.data, b_tr.len - 1) == EVMLR_SER_ERR_TRUNCATED;
        } else {
            bad_stmt = bad_missing = bad_order = bad_trunc = 0;
        }
        evmlr_buf_clear(&b_tr);
        evmlr_buf_clear(&b_tr2);
        for (slong i = 0; i < counted; i++) {
            nmod_poly_mat_clear(results[i]);
        }
        for (slong i = 0; i < n_checked; i++) {
            nmod_poly_mat_clear(checked[i]);
        }
        evmlr_voting_transcript_clear(tr);
        TEST_ASSERT(ok && ok_rt && verified == n_voters && same, end);
        TEST_ASSERT(bad_stmt && bad_missing && bad_order && bad_trunc, end);
    } TEST_END;

    TEST_ONCE("a server's keys round-trip, and a mismatched key pair is refused") {
        evmlr_mlpke_keypair_t kp, kp2, other;
        evmlr_mlpke_pk_t pk2;
        evmlr_buf_t b_pk, b_kp, b_again;
        int ok;

        evmlr_voting_keygen(kp, state);
        evmlr_voting_keygen(other, state);
        ok = ENCODE(b_pk, evmlr_serialize_pk(&b_pk, kp->pk))
          && ENCODE(b_kp, evmlr_serialize_keypair(&b_kp, kp))
          && evmlr_deserialize_pk(pk2, b_pk.data, b_pk.len) == EVMLR_SER_OK;
        if (ok) {
            ok = nmod_poly_mat_equal(pk2->A, kp->pk->A) && nmod_poly_mat_equal(pk2->t, kp->pk->t);
            nmod_poly_mat_clear(pk2->A);
            nmod_poly_mat_clear(pk2->t);
        }
        if (ok && evmlr_deserialize_keypair(kp2, b_kp.data, b_kp.len) == EVMLR_SER_OK) {
            ok = ENCODE(b_again, evmlr_serialize_keypair(&b_again, kp2)) && same_bytes(&b_kp, &b_again);
            evmlr_buf_clear(&b_again);
            evmlr_mlpke_keypair_clear(kp2);
        } else {
            ok = 0;
        }
        /* Another server's secret with this public key. */
        nmod_poly_mat_swap(kp->sk->s, other->sk->s);
        evmlr_buf_clear(&b_again);
        ok = ok && ENCODE(b_again, evmlr_serialize_keypair(&b_again, kp))
                && evmlr_deserialize_keypair(kp2, b_again.data, b_again.len) == EVMLR_SER_ERR_KEY;
        evmlr_buf_clear(&b_again);
        evmlr_buf_clear(&b_pk);
        evmlr_buf_clear(&b_kp);
        evmlr_mlpke_keypair_clear(kp);
        evmlr_mlpke_keypair_clear(other);
        TEST_ASSERT(ok, end);
    } TEST_END;

    TEST_ONCE("the parameters depend only on the election id and the public keys") {
        evmlr_mlpke_pk_t *pks = (evmlr_mlpke_pk_t *) malloc(k * sizeof(evmlr_mlpke_pk_t));
        evmlr_voting_pp_t pp3;
        evmlr_buf_t b_pp3;
        int ok;

        for (slong j = 0; j < k; j++) {
            pks[j][0] = pp->pks[j][0];
        }
        evmlr_voting_pp_init(pp3, pp->election_id, k, (const evmlr_mlpke_pk_t *) pks);
        ok = ENCODE(b_pp3, evmlr_serialize_pp(&b_pp3, pp3)) && same_bytes(&b_pp, &b_pp3)
          && nmod_poly_mat_equal(pp3->hpke_ctx->otse_ctx->H, pp->hpke_ctx->otse_ctx->H)
          && nmod_poly_mat_equal(pp3->hpke_ctx->otse_ctx->H_prime, pp->hpke_ctx->otse_ctx->H_prime);
        evmlr_buf_clear(&b_pp3);
        evmlr_voting_pp_clear(pp3);
        free(pks);
        TEST_ASSERT(ok, end);
    } TEST_END;

    TEST_ONCE("batches and mix records round-trip between servers") {
        evmlr_voting_batch_t in, in2, out, out2, tmp_b;
        evmlr_voting_mix_t mix, mix2, tmp_m;
        evmlr_mlpke_keypair_struct kp;
        evmlr_buf_t b_in, b_out, b_mix, b_again;
        int have_in = 0, have_in2 = 0, have_out = 0, have_out2 = 0, have_mix2 = 0, ok = 0;

        evmlr_buf_init(&b_in);
        evmlr_buf_init(&b_out);
        evmlr_buf_init(&b_mix);
        evmlr_buf_init(&b_again);
        kp.sk[0] = sk->sks[0][0];
        kp.pk[0] = pp->pks[0][0];
        do {
            if (evmlr_voting_batch_from_submissions(in, (const evmlr_voting_submission_t *) subs, n_voters, pp)) break;
            have_in = 1;
            if (evmlr_serialize_batch(&b_in, in, pp) || evmlr_deserialize_batch(in2, pp2, b_in.data, b_in.len)) break;
            have_in2 = 1;
            /* Server 0 mixes the batch it read from bytes, under imported parameters. */
            if (evmlr_voting_mix_step(out, mix, in2, &kp, pp2, state)) break;
            have_out = 1;
            if (evmlr_serialize_batch(&b_out, out, pp2) || evmlr_serialize_mix(&b_mix, mix, pp2)) break;
            if (evmlr_deserialize_batch(out2, pp, b_out.data, b_out.len)) break;
            have_out2 = 1;
            if (evmlr_deserialize_mix(mix2, pp, b_mix.data, b_mix.len)) break;
            have_mix2 = 1;
            if (!evmlr_voting_batch_matches_mix(out2, mix2, pp)) break;
            if (evmlr_serialize_mix(&b_again, mix2, pp) || !same_bytes(&b_mix, &b_again)) break;
            if (evmlr_deserialize_mix(tmp_m, pp, b_mix.data, b_mix.len - 1) != EVMLR_SER_ERR_TRUNCATED) break;
            if (evmlr_deserialize_batch(tmp_b, pp, b_out.data, b_out.len - 1) != EVMLR_SER_ERR_TRUNCATED) break;
            ok = 1;
        } while (0);
        if (have_in) evmlr_voting_batch_clear(in);
        if (have_in2) evmlr_voting_batch_clear(in2);
        if (have_out) {
            evmlr_voting_batch_clear(out);
            evmlr_voting_mix_clear(mix);
        }
        if (have_out2) evmlr_voting_batch_clear(out2);
        if (have_mix2) evmlr_voting_mix_clear(mix2);
        evmlr_buf_clear(&b_in);
        evmlr_buf_clear(&b_out);
        evmlr_buf_clear(&b_mix);
        evmlr_buf_clear(&b_again);
        TEST_ASSERT(ok, end);
    } TEST_END;

end:
    evmlr_buf_clear(&b_pp);
    evmlr_buf_clear(&b_pp2);
    evmlr_buf_clear(&b_sk);
    evmlr_buf_clear(&b_sk2);
    evmlr_buf_clear(&b_sub2);
    for (slong i = 0; i < n_voters; i++) {
        evmlr_buf_clear(&b_sub[i]);
        evmlr_voting_submission_clear(subs[i], pp);
    }
    for (int i = 0; i < have_subs2; i++) {
        evmlr_voting_submission_clear(subs2[i], pp2);
    }
    if (have_sk2) evmlr_voting_sk_clear(sk2);
    if (have_pp2) evmlr_voting_pp_clear(pp2);
    evmlr_voting_sk_clear(sk);
    evmlr_voting_pp_clear(pp);
}

/* BENCH_ADD repeats its argument, so each timed call cleans up after itself. */
static void bench_import_pp(const evmlr_buf_t *b) {
    evmlr_voting_pp_t pp;
    if (evmlr_deserialize_pp(pp, b->data, b->len) == EVMLR_SER_OK) evmlr_voting_pp_clear(pp);
}

static void bench_export_sub(const evmlr_voting_submission_t sub, const evmlr_voting_pp_t pp) {
    evmlr_buf_t b;
    evmlr_buf_init(&b);
    evmlr_serialize_submission(&b, sub, pp);
    evmlr_buf_clear(&b);
}

static void bench_import_sub(const evmlr_buf_t *b, const evmlr_voting_pp_t pp) {
    evmlr_voting_submission_t sub;
    if (evmlr_deserialize_submission(sub, pp, b->data, b->len) == EVMLR_SER_OK) evmlr_voting_submission_clear(sub, pp);
}

static void bench(flint_rand_t state) {
    evmlr_voting_pp_t pp;
    evmlr_voting_sk_t sk;
    evmlr_voting_submission_t sub;
    evmlr_buf_t b_pp, b_sub;
    nmod_poly_mat_t vote;

    evmlr_voting_setup(pp, sk, 2, state);
    nmod_poly_mat_init(vote, 1, 1, MOD_Q);
    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(vote, 0, 0), 0, 5);
    evmlr_voting_submission_init(sub, pp);
    evmlr_voting_casting(sub, 1, vote, pp, state);

    evmlr_buf_init(&b_pp);
    evmlr_buf_init(&b_sub);
    evmlr_serialize_pp(&b_pp, pp);
    evmlr_serialize_submission(&b_sub, sub, pp);
    printf("encoded public parameters (k = 2, N = %d): %zu bytes\n", SHUFFLE_N_MSGS, b_pp.len);
    printf("encoded submission (k = 2): %zu bytes\n", b_sub.len);

    BENCH_BEGIN("evmlr_deserialize_pp") {
        BENCH_ADD(bench_import_pp(&b_pp));
    } BENCH_END;
    BENCH_BEGIN("evmlr_serialize_submission") {
        BENCH_ADD(bench_export_sub(sub, pp));
    } BENCH_END;
    BENCH_BEGIN("evmlr_deserialize_submission") {
        BENCH_ADD(bench_import_sub(&b_sub, pp));
    } BENCH_END;

    evmlr_buf_clear(&b_pp);
    evmlr_buf_clear(&b_sub);
    nmod_poly_mat_clear(vote);
    evmlr_voting_submission_clear(sub, pp);
    evmlr_voting_sk_clear(sk);
    evmlr_voting_pp_clear(pp);
}

int main(int argc, char *argv[]) {
    evmlr_mode_t mode = evmlr_mode(argc, argv);
    flint_rand_t state;

    evmlr_rand_init(state);
    if (evmlr_runs_tests(mode))   test(state);
    if (evmlr_runs_benches(mode)) bench(state);
    flint_rand_clear(state);
    return test_status();
}
#endif
