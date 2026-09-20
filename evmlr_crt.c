#include "evmlr_crt.h"
#include "flint/nmod_poly.h"
#include "flint/ulong_extras.h"
#include "flint/nmod.h"

#define HALF (DEGREE_N / 2)

static int  crt_ready = 0;   // constants computed?
static int  crt_usable = 0;  // does the decomposition apply to these parameters?
static nmod_t crt_mod;
static ulong crt_r;      // r with r^2 = -1 mod q
static ulong crt_inv2;   // 2^-1 mod q
static ulong crt_inv2r;  // (2r)^-1 mod q

static void crt_setup(void) {
    crt_ready = 1;
    crt_usable = 0;
    if (DEGREE_N % 2 != 0) return;

    nmod_init(&crt_mod, MOD_Q);
    for (ulong x = 2; x < MOD_Q; x++) {
        if (nmod_mul(x, x, crt_mod) == (ulong) MOD_Q - 1) { crt_r = x; break; }
        if (x == (ulong) MOD_Q - 1) return; // no square root of -1
    }
    if (crt_r == 0) return;

    crt_inv2  = n_invmod(2 % MOD_Q, MOD_Q);
    crt_inv2r = n_invmod(nmod_mul(2, crt_r, crt_mod), MOD_Q);
    if (crt_inv2 == 0 || crt_inv2r == 0) return;
    crt_usable = 1;
}

int evmlr_crt_available(void) {
    if (!crt_ready) crt_setup();
    return crt_usable;
}

// Split a into its two residues: a mod (x^HALF - r) and a mod (x^HALF + r).
// Writing a = lo + x^HALF * hi, those are lo + r*hi and lo - r*hi.
static void crt_split(ulong *res1, ulong *res2, const nmod_poly_t a) {
    const slong len = a->length;
    for (slong i = 0; i < HALF; i++) {
        ulong lo = (i < len) ? a->coeffs[i] : 0;
        ulong hi = (i + HALF < len) ? a->coeffs[i + HALF] : 0;
        ulong rh = nmod_mul(hi, crt_r, crt_mod);
        res1[i] = nmod_add(lo, rh, crt_mod);
        res2[i] = nmod_sub(lo, rh, crt_mod);
    }
}

// c = a * b mod (x^HALF - s), for operands already reduced to degree < HALF.
static void crt_mul_half(ulong *c, const ulong *a, const ulong *b, ulong s, ulong *scratch) {
    _nmod_poly_mul(scratch, a, HALF, b, HALF, crt_mod);
    const slong plen = 2 * HALF - 1;
    for (slong i = 0; i < HALF; i++) {
        ulong hi = (i + HALF < plen) ? scratch[i + HALF] : 0;
        c[i] = nmod_add(scratch[i], nmod_mul(hi, s, crt_mod), crt_mod);
    }
}

void evmlr_crt_mulmod(nmod_poly_t res, const nmod_poly_t a, const nmod_poly_t b,
                      const nmod_poly_t cyclo) {
    if (!crt_ready) crt_setup();

    if (!crt_usable
        || a->mod.n != (ulong) MOD_Q || b->mod.n != (ulong) MOD_Q
        || cyclo->mod.n != (ulong) MOD_Q || res->mod.n != (ulong) MOD_Q
        || nmod_poly_degree(cyclo) != DEGREE_N
        || a->length > DEGREE_N || b->length > DEGREE_N) {
        nmod_poly_mulmod(res, a, b, cyclo);
        return;
    }

    ulong a1[HALF], a2[HALF], b1[HALF], b2[HALF];
    ulong c1[HALF], c2[HALF], scratch[2 * HALF];

    // a and b are fully read here, so res may alias either of them.
    crt_split(a1, a2, a);
    crt_split(b1, b2, b);

    crt_mul_half(c1, a1, b1, crt_r, scratch);
    crt_mul_half(c2, a2, b2, MOD_Q - crt_r, scratch);

    // Recombine: lo = (c1 + c2)/2, hi = (c1 - c2)/(2r).
    nmod_poly_fit_length(res, DEGREE_N);
    for (slong i = 0; i < HALF; i++) {
        ulong x = c1[i], y = c2[i];
        res->coeffs[i]        = nmod_mul(nmod_add(x, y, crt_mod), crt_inv2,  crt_mod);
        res->coeffs[i + HALF] = nmod_mul(nmod_sub(x, y, crt_mod), crt_inv2r, crt_mod);
    }
    res->length = DEGREE_N;
    _nmod_poly_normalise(res);
}

// Is cyclo exactly x^n + 1?
static int is_negacyclic(const nmod_poly_t cyclo, slong *n_out) {
    const slong n = nmod_poly_degree(cyclo);
    if (n <= 0) return 0;
    if (nmod_poly_get_coeff_ui(cyclo, n) != 1) return 0;
    if (nmod_poly_get_coeff_ui(cyclo, 0) != 1) return 0;
    for (slong i = 1; i < n; i++)
        if (nmod_poly_get_coeff_ui(cyclo, i) != 0) return 0;
    *n_out = n;
    return 1;
}

void evmlr_crt_negacyclic_rem(nmod_poly_t res, const nmod_poly_t a, const nmod_poly_t cyclo) {
    slong n;
    if (!is_negacyclic(cyclo, &n)) {
        nmod_poly_rem(res, a, cyclo);
        return;
    }

    const slong alen = a->length;
    if (alen <= n) { // already reduced
        if (res != a) nmod_poly_set(res, a);
        return;
    }

    nmod_t mod = a->mod;
    ulong acc[DEGREE_N];
    if (n > DEGREE_N) { // no room in the scratch buffer; let FLINT handle it
        nmod_poly_rem(res, a, cyclo);
        return;
    }

    for (slong i = 0; i < n; i++) acc[i] = (i < alen) ? a->coeffs[i] : 0;
    // x^(i + k*n) = (-1)^k x^i
    for (slong k = 1; k * n < alen; k++) {
        const slong base = k * n;
        for (slong i = 0; i < n && base + i < alen; i++) {
            ulong v = a->coeffs[base + i];
            acc[i] = (k % 2) ? nmod_sub(acc[i], v, mod) : nmod_add(acc[i], v, mod);
        }
    }

    nmod_poly_fit_length(res, n);
    for (slong i = 0; i < n; i++) res->coeffs[i] = acc[i];
    res->length = n;
    _nmod_poly_normalise(res);
}
