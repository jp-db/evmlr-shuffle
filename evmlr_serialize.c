#include "evmlr_serialize.h"

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
 * Contexts without sampling
 *
 * The *_ctx_init functions sample their matrices; on import the matrices come
 * from the input instead. These mirror those functions' allocations and
 * fixed parts (keep them in sync with evmlr_otse_ctx_init,
 * evmlr_commit_ctx_init and evmlr_shuffle_ctx_init).
 * ------------------------------------------------------------------------- */

static void set_cyclo(nmod_poly_t cyclo) {
    nmod_poly_init(cyclo, MOD_Q);
    nmod_poly_set_coeff_ui(cyclo, DEGREE_N, 1);
    nmod_poly_set_coeff_ui(cyclo, 0, 1);
}

static void hpke_ctx_init_empty(evmlr_hpke_ctx_t ctx, slong L) {
    evmlr_otse_ctx_struct *otse = ctx->otse_ctx;

    evmlr_mlpke_ctx_init(ctx->enc_ctx);
    otse->L = L;
    nmod_poly_mat_init(otse->H, K_LWE + L, K_LWR, 1UL << ZETA);
    /* H' = (H'' | I_L); only H'' is encoded. */
    nmod_poly_mat_init(otse->H_prime, L, K_LWE + L, MOD_Q);
    for (slong i = 0; i < L; i++) {
        nmod_poly_one(nmod_poly_mat_entry(otse->H_prime, i, K_LWE + i));
    }
    set_cyclo(otse->cyclo_poly);
}

static void commit_ctx_init_empty(evmlr_commit_ctx_t ctx, slong N) {
    ctx->b_sqr = compute_b(N);
    ctx->N = N;
    nmod_poly_mat_init(ctx->A_1, K_SIS, N, MOD_Q);
    nmod_poly_mat_init(ctx->A_2, K_SIS, 2 * K_SIS, MOD_Q);
    set_cyclo(ctx->cyclo_poly);
}

static uint64_t voting_L_max(uint64_t k) {
    return (k - 1) * K_LWR * (K_LWE + 1) + 1;
}

static uint64_t shuffle_com_N(uint64_t N, uint64_t L) {
    uint64_t row = (uint64_t) (2 * ETA) * (K_LWE + L);

    return N * (LOG_Q_CEIL > row ? LOG_Q_CEIL : row);
}

/* ---------------------------------------------------------------------------
 * Public parameters
 * ------------------------------------------------------------------------- */

int evmlr_serialize_pp(evmlr_buf_t *out, const evmlr_voting_pp_t pp) {
    writer_t w = { out, 0 };
    const evmlr_otse_ctx_struct *otse = pp->hpke_ctx->otse_ctx;
    const evmlr_shuffle_ctx_struct *shuf = pp->shuf_ctx;

    if (pp->k < 1 || pp->L_max != (slong) voting_L_max((uint64_t) pp->k) || otse->L != pp->L_max
        || shuf->L != pp->L_max
        || shuf->com_ctx->N != (slong) shuffle_com_N((uint64_t) shuf->N, (uint64_t) shuf->L)) {
        return EVMLR_SER_ERR_ARG;
    }
    put_header(&w, EVMLR_SER_PP);
    put_u32(&w, (uint64_t) pp->k);
    put_u32(&w, (uint64_t) shuf->N);
    put_mat(&w, otse->H);
    put_mat_block(&w, otse->H_prime, 0, 0, otse->L, K_LWE);
    for (slong j = 0; j < pp->k; j++) {
        put_mat(&w, pp->pks[j]->A);
        put_mat(&w, pp->pks[j]->t);
    }
    put_mat(&w, shuf->com_ctx->A_1);
    put_mat(&w, shuf->com_ctx->A_2);
    return w.err;
}

int evmlr_deserialize_pp(evmlr_voting_pp_t pp, const uint8_t *in, size_t len) {
    reader_t r = open_reader(in, len, EVMLR_SER_PP);
    uint64_t k, N, L, com_N, polys;
    evmlr_mlpke_pk_t *pks;

    k = get_u32(&r);
    N = get_u32(&r);
    if (r.err) {
        return r.err;
    }
    if (k < 1 || k > MAX_SERVERS || N < 1 || N > MAX_MESSAGES) {
        return EVMLR_SER_ERR_SHAPE;
    }
    L = voting_L_max(k);
    com_N = shuffle_com_N(N, L);
    /* Every polynomial takes at least its 2-byte length. */
    polys = (K_LWE + L) * K_LWR + L * K_LWE + k * (K_LWE * K_LWE + K_LWE)
            + K_SIS * com_N + K_SIS * 2 * K_SIS;
    if (polys > r.left / 2) {
        return EVMLR_SER_ERR_TRUNCATED;
    }

    pks = (evmlr_mlpke_pk_t *) malloc((size_t) k * sizeof(evmlr_mlpke_pk_t));
    if (pks == NULL) {
        return EVMLR_SER_ERR_NOMEM;
    }
    pp->k = (slong) k;
    pp->L_max = (slong) L;
    pp->pks = pks;
    hpke_ctx_init_empty(pp->hpke_ctx, (slong) L);
    commit_ctx_init_empty(pp->shuf_ctx->com_ctx, (slong) com_N);
    hpke_ctx_init_empty(pp->shuf_ctx->hpke_ctx, (slong) L);
    pp->shuf_ctx->N = (slong) N;
    pp->shuf_ctx->L = (slong) L;
    for (uint64_t j = 0; j < k; j++) {
        nmod_poly_mat_init(pks[j]->A, K_LWE, K_LWE, MOD_Q);
        nmod_poly_mat_init(pks[j]->t, K_LWE, 1, MOD_Q);
    }

    get_mat(&r, pp->hpke_ctx->otse_ctx->H);
    get_mat_block(&r, pp->hpke_ctx->otse_ctx->H_prime, 0, 0, (slong) L, K_LWE);
    for (uint64_t j = 0; j < k && !r.err; j++) {
        get_mat(&r, pks[j]->A);
        get_mat(&r, pks[j]->t);
    }
    get_mat(&r, pp->shuf_ctx->com_ctx->A_1);
    get_mat(&r, pp->shuf_ctx->com_ctx->A_2);

    if (close_reader(&r)) {
        evmlr_voting_pp_clear(pp);
        return r.err;
    }
    /* The shuffle runs over the same OTSE matrices (see sync_hpke_ctx). */
    nmod_poly_mat_set(pp->shuf_ctx->hpke_ctx->otse_ctx->H, pp->hpke_ctx->otse_ctx->H);
    nmod_poly_mat_set(pp->shuf_ctx->hpke_ctx->otse_ctx->H_prime, pp->hpke_ctx->otse_ctx->H_prime);
    return EVMLR_SER_OK;
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
    put_u64(&w, (uint64_t) (int64_t) sub->voter_id);
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
    if ((int64_t) (slong) voter_id != voter_id) { /* slong is 32 bits on wasm32 */
        return EVMLR_SER_ERR_VALUE;
    }

    evmlr_voting_submission_init(sub, pp);
    sub->voter_id = (slong) voter_id;
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
 * Shuffle proofs and statements
 * ------------------------------------------------------------------------- */

int evmlr_serialize_shuffle_proof(evmlr_buf_t *out, const evmlr_shuffle_proof_t proof,
                                  const evmlr_shuffle_ctx_t ctx) {
    writer_t w = { out, 0 };

    put_header(&w, EVMLR_SER_SHUFFLE_PROOF);
    put_u32(&w, (uint64_t) ctx->N);
    put_mat(&w, proof->u);
    put_mat(&w, proof->P->c);
    put_mat(&w, proof->W->c);
    put_bin_proof(&w, proof->proof_D);
    put_bin_proof(&w, proof->proof_P);
    put_bin_proof(&w, proof->proof_W);
    put_lin_proof(&w, proof->proof_u);
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
    /* evmlr_proof_init leaves P and W to the prover; prove initialises them. */
    evmlr_proof_init(proof, ctx);
    evmlr_commit_init(proof->P);
    evmlr_commit_init(proof->W);

    get_mat(&r, proof->u);
    get_mat(&r, proof->P->c);
    get_mat(&r, proof->W->c);
    get_bin_proof(&r, proof->proof_D);
    get_bin_proof(&r, proof->proof_P);
    get_bin_proof(&r, proof->proof_W);
    get_lin_proof(&r, proof->proof_u);
    if (close_reader(&r)) {
        evmlr_proof_clear(proof, ctx);
        return r.err;
    }
    return EVMLR_SER_OK;
}

int evmlr_serialize_shuffle_statement(evmlr_buf_t *out, const evmlr_shuffle_pp_t stmt,
                                      const evmlr_shuffle_ctx_t ctx) {
    writer_t w = { out, 0 };

    put_header(&w, EVMLR_SER_SHUFFLE_STMT);
    put_u32(&w, (uint64_t) ctx->N);
    put_mat(&w, stmt->D->c);
    for (slong i = 0; i < ctx->N; i++) {
        put_mat(&w, stmt->c_hat[i]);
        put_mat(&w, stmt->c_star[i]->c);
    }
    return w.err;
}

int evmlr_deserialize_shuffle_statement(evmlr_shuffle_pp_t stmt, const evmlr_shuffle_ctx_t ctx,
                                        const uint8_t *in, size_t len) {
    reader_t r = open_reader(in, len, EVMLR_SER_SHUFFLE_STMT);
    uint64_t N = get_u32(&r);
    slong n;

    if (r.err) {
        return r.err;
    }
    if (N != (uint64_t) ctx->N) {
        return EVMLR_SER_ERR_SHAPE;
    }
    n = ctx->N;
    stmt->c_hat = (nmod_poly_mat_t *) malloc((size_t) n * sizeof(nmod_poly_mat_t));
    stmt->c_star = (evmlr_otse_ciphertext_t *) malloc((size_t) n * sizeof(evmlr_otse_ciphertext_t));
    if (stmt->c_hat == NULL || stmt->c_star == NULL) {
        free(stmt->c_hat);
        free(stmt->c_star);
        return EVMLR_SER_ERR_NOMEM;
    }
    evmlr_commit_init(stmt->D);
    for (slong i = 0; i < n; i++) {
        nmod_poly_mat_init(stmt->c_hat[i], ctx->L, 1, MOD_Q);
        evmlr_otse_ciphertext_init(stmt->c_star[i], ctx->L);
    }

    get_mat(&r, stmt->D->c);
    for (slong i = 0; i < n && !r.err; i++) {
        get_mat(&r, stmt->c_hat[i]);
        get_mat(&r, stmt->c_star[i]->c);
    }
    if (close_reader(&r)) {
        evmlr_shuffle_pp_clear(stmt, ctx);
        return r.err;
    }
    return EVMLR_SER_OK;
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
    /* GCC >= 11 at -O2 reports a 128-byte read from a 40-byte region here
     * (-Wstringop-overread); the argument is the 128-byte layer, and the
     * test is clean under AddressSanitizer. A false positive of its
     * array-parameter analysis. */
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstringop-overread"
#endif
    ok = evmlr_voting_verify_enc(pp, sub->proofs[0], sub->ballot->layer, msg, 0);
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic pop
#endif
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
        int ok_rt, ok_ver, ok_bad;

        /* Proved under the original shuffle context, verified under the one
         * imported with the public parameters. */
        shuffle_instance(sp, st, pp->shuf_ctx, state);
        evmlr_proof_init(proof, pp->shuf_ctx);
        evmlr_shuffle_prove(proof, sp, st, pp->shuf_ctx, state);

        ok_rt = ENCODE(b_proof, evmlr_serialize_shuffle_proof(&b_proof, proof, pp->shuf_ctx))
             && ENCODE(b_st, evmlr_serialize_shuffle_statement(&b_st, st, pp->shuf_ctx))
             && evmlr_deserialize_shuffle_proof(proof2, pp2->shuf_ctx, b_proof.data, b_proof.len) == EVMLR_SER_OK
             && evmlr_deserialize_shuffle_statement(st2, pp2->shuf_ctx, b_st.data, b_st.len) == EVMLR_SER_OK;
        ok_ver = ok_bad = 0;
        if (ok_rt) {
            ok_rt = ENCODE(b_again, evmlr_serialize_shuffle_proof(&b_again, proof2, pp2->shuf_ctx))
                 && same_bytes(&b_proof, &b_again);
            evmlr_buf_clear(&b_again);
            ok_ver = evmlr_shuffle_verify(proof2, pp2->shuf_ctx, st2);
            /* The proof must not verify for another statement. */
            nmod_poly_set_coeff_ui(nmod_poly_mat_entry(st2->c_hat[0], 0, 0), 0,
                                   (nmod_poly_get_coeff_ui(nmod_poly_mat_entry(st2->c_hat[0], 0, 0), 0) + 1) % MOD_Q);
            ok_bad = !evmlr_shuffle_verify(proof2, pp2->shuf_ctx, st2);
            evmlr_proof_clear(proof2, pp2->shuf_ctx);
            evmlr_shuffle_pp_clear(st2, pp2->shuf_ctx);
        }
        evmlr_buf_clear(&b_proof);
        evmlr_buf_clear(&b_st);
        evmlr_proof_clear(proof, pp->shuf_ctx);
        evmlr_shuffle_clear_sp(sp, pp->shuf_ctx);
        evmlr_shuffle_pp_clear(st, pp->shuf_ctx);
        TEST_ASSERT(ok_rt && ok_ver && ok_bad, end);
    } TEST_END;

    TEST_ONCE("imported ballots are counted under imported keys") {
        nmod_poly_mat_t results[3];
        evmlr_shuffle_proof_t shuffle_proofs[2];
        slong counted = evmlr_voting_counting(results, shuffle_proofs, (const evmlr_voting_submission_t *) subs2,
                                              n_voters, pp2, sk2, state);
        int seen[3] = { 0, 0, 0 }, ok = (counted == n_voters);

        for (slong i = 0; i < counted; i++) {
            ulong v = nmod_poly_get_coeff_ui(nmod_poly_mat_entry(results[i], 0, 0), 0);
            for (int j = 0; j < 3; j++) {
                if (v == votes[j]) seen[j]++;
            }
            nmod_poly_mat_clear(results[i]);
        }
        for (slong j = 0; j < k && counted > 0; j++) {
            evmlr_proof_clear(shuffle_proofs[j], pp2->shuf_ctx);
        }
        TEST_ASSERT(ok && seen[0] == 1 && seen[1] == 1 && seen[2] == 1, end);
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
