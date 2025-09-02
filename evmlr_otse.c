#include "evmlr_otse.h"
#include "evmlr_utils.h"
#include "flint/fmpz.h"
#include "flint/fmpz_poly.h"
#ifdef MAIN
#include "test.h"
#include "bench.h"
#endif

void calc_a(nmod_poly_t a[], const evmlr_otse_key_t key, const evmlr_otse_ctx_t ctx);

void evmlr_otse_ctx_init(evmlr_otse_ctx_t ctx, size_t L, flint_rand_t state) {
    uint alloc; // prevent compiler warning
    ctx->L = L;
    // H \sample R_{2^ZETA}^{(K_LWE + L) x K_LWR}
    alloc = (K_LWE + L) * sizeof(nmod_poly_t*);
    ctx->H = (nmod_poly_t**) malloc(alloc);
    for (int i = 0; i < K_LWE + L; i++) {
        ctx->H[i] = (nmod_poly_t*) malloc(K_LWR * sizeof(nmod_poly_t));
        for (int j = 0; j < K_LWR; j++) {
            nmod_poly_init(ctx->H[i][j], 2*ZETA);
            nmod_poly_randtest(ctx->H[i][j], state, DEGREE_N);
        }
    }

    // H' = (H'' | I_{L}) \in R_q^{L x (K_LWE + L)}
    alloc = L * sizeof(nmod_poly_t*);
    ctx->H_prime = (nmod_poly_t**) malloc(alloc);
    for (int i = 0; i < L; i++) {
        ctx->H_prime[i] = (nmod_poly_t*) malloc((K_LWE + L) * sizeof(nmod_poly_t));
        for (int j = 0; j < K_LWE + L; j++) {
            nmod_poly_init(ctx->H_prime[i][j], MOD_Q);
            if (j < K_LWE) { // H'' \sample R_q^{L x K_LWE}
                nmod_poly_randtest(ctx->H_prime[i][j], state, DEGREE_N);
            } else if (j == K_LWE + i) { // Identity matrix where row=col
                nmod_poly_one(ctx->H_prime[i][j]);
            } else { // Identity matrix where row!=col
                nmod_poly_zero(ctx->H_prime[i][j]);
            }
        }
    }

    nmod_poly_init(ctx->cyclo_poly, MOD_Q);
    nmod_poly_set_coeff_ui(ctx->cyclo_poly, DEGREE_N, 1);
    nmod_poly_set_coeff_ui(ctx->cyclo_poly, 0, 1);
}

void evmlr_otse_ctx_clear(evmlr_otse_ctx_t ctx) {
    for (int i = 0; i < K_LWE + ctx->L; i++) {
        for (int j = 0; j < K_LWR; j++) {
            nmod_poly_clear(ctx->H[i][j]);
        }
        free(ctx->H[i]);
    }
    free(ctx->H);
    for (int i = 0; i < ctx->L; i++) {
        for (int j = 0; j < K_LWE + ctx->L; j++) {
            nmod_poly_clear(ctx->H_prime[i][j]);
        }
        free(ctx->H_prime[i]);
    }
    free(ctx->H_prime);
    nmod_poly_clear(ctx->cyclo_poly);
}

void evmlr_otse_keygen(evmlr_otse_key_t key, flint_rand_t state) {
    for (int i = 0; i < K_LWR; i++) {
        nmod_poly_init(key->s[i], MOD_Q);
        evmlr_utils_sample_binary_poly(key->s[i], DEGREE_N, state);
    }
}

void evmlr_otse_keyclear(evmlr_otse_key_t key) {
    for (int i = 0; i < K_LWR; i++) {
        nmod_poly_clear(key->s[i]);
    }
}

void evmlr_otse_encrypt(evmlr_otse_ciphertext_t ct, const nmod_poly_t m[], const evmlr_otse_key_t key, const evmlr_otse_ctx_t ctx) {
    nmod_poly_t a[ctx->L], tmp;
    nmod_poly_init(tmp, MOD_Q);
    calc_a(a, key, ctx);
    // c = m + a
    ct->c = (nmod_poly_t*) malloc(ctx->L * sizeof(nmod_poly_t));
    for (int i = 0; i < ctx->L; i++) {
        nmod_poly_init(ct->c[i], MOD_Q);
        nmod_poly_add(ct->c[i], m[i], a[i]);
    }

    nmod_poly_clear(tmp);
    for (int i = 0; i < ctx->L; i++) {
        nmod_poly_clear(a[i]);
    }
}

void evmlr_otse_decrypt(nmod_poly_t m[], const evmlr_otse_ciphertext_t ct, const evmlr_otse_key_t key, const evmlr_otse_ctx_t ctx) {
    nmod_poly_t a[ctx->L], tmp;
    nmod_poly_init(tmp, MOD_Q);
    calc_a(a, key, ctx);
    // m = c - a
    for (int i = 0; i < ctx->L; i++) {
        nmod_poly_init(m[i], MOD_Q);
        nmod_poly_sub(m[i], ct->c[i], a[i]);
    }
    nmod_poly_clear(tmp);
    for (int i = 0; i < ctx->L; ++i) {
        nmod_poly_clear(a[i]);
    }
}

void evmlr_otse_ciphertext_clear(evmlr_otse_ciphertext_t ct, const evmlr_otse_ctx_t ctx) {
    for (int i = 0; i < ctx->L; i++) {
        nmod_poly_clear(ct->c[i]);
    }
    free(ct->c);
}

void round_poly(nmod_poly_t poly, const nmod_poly_t hs, const fmpz_t q1, const fmpz_t q2) {
    fmpz_t coeff, half_q1;
    fmpz_init(coeff);
    // coeff for rounding equivalent to adding 0.5 before flooring so we round to nearest integer with ties rounding up
    fmpz_init(half_q1);
    fmpz_fdiv_q_2exp(half_q1, q1, 1); // half_q1 = q1 / 2

    slong deg = nmod_poly_degree(hs);
    nmod_poly_fit_length(poly, deg + 1);
    for (slong i = 0; i <= deg; i++) {
        ulong coeff_hs = nmod_poly_get_coeff_ui(hs, i);
        fmpz_set_ui(coeff, coeff_hs);
        // round(coeff * q2 / q1)
        fmpz_mul(coeff, coeff, q2);
        fmpz_add(coeff, coeff, half_q1);
        fmpz_fdiv_q(coeff, coeff, q1);

        fmpz_mod(coeff, coeff, q2); // coeff mod q2
        nmod_poly_set_coeff_ui(poly, i, fmpz_get_ui(coeff));
    }

    fmpz_clear(coeff);
    fmpz_clear(half_q1);
}

void calc_d(ulong L, nmod_poly_t d[K_LWE + L], nmod_poly_t** H, const nmod_poly_t s[K_LWR], const nmod_poly_t cyclo_poly) {
    // ⌈Hs⌋_{2^ζ →2^{2η}} := ⌈(Hs * 2^{2η})/2^ζ⌋ mod 2^{2η} := Hs * q2/q1 mod q2
    fmpz_t q1, q2;
    fmpz_init(q1);
    fmpz_init(q2);
    fmpz_set_ui(q1, 1UL << ZETA);      // q1 = 2^ζ
    fmpz_set_ui(q2, 1UL << (2 * ETA)); // q2 = 2^{2η}

    // Hs
    nmod_poly_t Hs[K_LWE + L], tmp;
    nmod_poly_init(tmp, MOD_Q);
    for (int i = 0; i < K_LWE + L; i++) {
        nmod_poly_init(Hs[i], MOD_Q);
        nmod_poly_zero(Hs[i]);
        for (int j = 0; j < K_LWR; j++) {
            nmod_poly_mulmod(tmp, H[i][j], s[j], cyclo_poly);
            nmod_poly_add(Hs[i], Hs[i], tmp);
        }
    }
    // d = Hs * q2/q1 mod q2
    for (int i = 0; i < K_LWE + L; i++) {
        nmod_poly_init(d[i], MOD_Q);
        round_poly(d[i], Hs[i], q1, q2);
    }

    nmod_poly_clear(tmp);
    fmpz_clear(q1);
    fmpz_clear(q2);
    for (int i = 0; i < K_LWE + L; i++) {
        nmod_poly_clear(Hs[i]);
    }
}

void multiply_by_B_eta(ulong L, nmod_poly_t w[K_LWE + L], const nmod_poly_t d[(2*ETA)* (K_LWE + L)], ulong mod) {
    // Accumulate the first eta vectors
    for (int i = 0; i < K_LWE + L; i++) {
        nmod_poly_init(w[i], mod);
        nmod_poly_zero(w[i]);
        for (int b = 0; b < ETA; b++) {
            nmod_poly_add(w[i], w[i], d[b * (K_LWE + L) + i]);
        }
    }

    // Subtract the last eta vectors
    for (int b = ETA; b < 2 * ETA; b++) {
        for (int i = 0; i < K_LWE + L; i++) {
            nmod_poly_sub(w[i], w[i], d[b * (K_LWE + L) + i]);
        }
    }
}

void calc_a(nmod_poly_t a[], const evmlr_otse_key_t key, const evmlr_otse_ctx_t ctx) {
    nmod_poly_t d[K_LWE + ctx->L], d_bin[2 * ETA][K_LWE + ctx->L];
    fmpz_poly_t d_fmpz[K_LWE + ctx->L];

    calc_d(ctx->L, d, ctx->H, key->s, ctx->cyclo_poly);
    for (int i = 0; i < K_LWE + ctx->L; i++) {
        fmpz_poly_init(d_fmpz[i]);
        fmpz_poly_set_nmod_poly(d_fmpz[i], d[i]);
    }

    slong len = K_LWE + ctx->L;
    int bits = 2 * ETA;
    evmlr_utils_ring_to_bin(len, bits, d_bin, d_fmpz, MOD_Q);

    // a = H' B^n stack(d)
    evmlr_calc_a(a, ctx, (nmod_poly_t *) d_bin);

    for (int i = 0; i < K_LWE + ctx->L; i++) {
        nmod_poly_clear(d[i]);
        fmpz_poly_clear(d_fmpz[i]);
        for (int j = 0; j < 2*ETA; ++j) {
            nmod_poly_clear(d_bin[j][i]);
        }
    }
}

void evmlr_calc_a(nmod_poly_t a[], const evmlr_otse_ctx_t ctx, const nmod_poly_t d_stack[2 * ETA * (K_LWE + ctx->L)]) {
    nmod_poly_t w[K_LWE + ctx->L], tmp;
    nmod_poly_init(tmp, MOD_Q);
    // w = B_eta stack(d)
    multiply_by_B_eta(ctx->L, w, d_stack, MOD_Q);

    // a = H' w
    for (int i = 0; i < ctx->L; i++) {
        nmod_poly_init(a[i], MOD_Q);
        nmod_poly_zero(a[i]);
        for (int j = 0; j < K_LWE + ctx->L; j++) {
            nmod_poly_init(tmp, MOD_Q);
            nmod_poly_mulmod(tmp, ctx->H_prime[i][j], w[j], ctx->cyclo_poly);
            nmod_poly_add(a[i], a[i], tmp);
        }
    }

    nmod_poly_clear(tmp);
    for (int i = 0; i < K_LWE + ctx->L; i++) {
        nmod_poly_clear(w[i]);
    }
}

#ifdef MAIN
void test(flint_rand_t rand, evmlr_otse_ctx_t ctx) {
    evmlr_otse_key_t key;
    evmlr_otse_keygen(key, rand);

    nmod_poly_t msg[ctx->L], decrypted_msg[ctx->L];
    for (int i = 0; i < ctx->L; i++) {
        nmod_poly_init(msg[i], MOD_Q);
        nmod_poly_randtest(msg[i], rand, DEGREE_N);
    }

    evmlr_otse_ciphertext_t ct;
    TEST_BEGIN("encryption and decryption are consistent") {
        evmlr_otse_encrypt(ct, msg, key, ctx);
        evmlr_otse_decrypt(decrypted_msg, ct, key, ctx);
        for (int i = 0; i < ctx->L; i++) {
            TEST_ASSERT(nmod_poly_equal(msg[i], decrypted_msg[i]) == 1, end)
        }
    } TEST_END;

    end:
        evmlr_otse_ciphertext_clear(ct, ctx);
        evmlr_otse_keyclear(key);
        for (int i = 0; i < ctx->L; i++) {
            nmod_poly_clear(msg[i]);
            nmod_poly_clear(decrypted_msg[i]);
        }
}

void bench(flint_rand_t rand, evmlr_otse_ctx_t ctx) {
    evmlr_otse_key_t key;
    BENCH_BEGIN("evmlr_otse_keygen") {
        BENCH_ADD(evmlr_otse_keygen(key, rand))
    } BENCH_END;

    nmod_poly_t msg[ctx->L];
    for (int i = 0; i < ctx->L; i++) {
        nmod_poly_init(msg[i], MOD_Q);
        nmod_poly_randtest(msg[i], rand, DEGREE_N);
    }

    evmlr_otse_ciphertext_t ct;
    BENCH_BEGIN("evmlr_otse_encrypt") {
        BENCH_ADD(evmlr_otse_encrypt(ct, msg, key, ctx))
    } BENCH_END;

    nmod_poly_t decrypted_msg[ctx->L];
    BENCH_BEGIN("evmlr_otse_decrypt") {
        BENCH_ADD(evmlr_otse_decrypt(decrypted_msg, ct, key, ctx))
    } BENCH_END;

    evmlr_otse_keyclear(key);
    for (int i = 0; i < ctx->L; i++) {
        nmod_poly_clear(msg[i]);
        nmod_poly_clear(decrypted_msg[i]);
    }
}


int main() {
    flint_rand_t state;
    flint_rand_init(state);
    ulong seed[2];
    getrandom(seed, sizeof(ulong)*2, 0);
    flint_rand_set_seed(state, seed[0], seed[1]);

    evmlr_otse_ctx_t ctx;
    evmlr_otse_ctx_init(ctx, M_LEN, state);

    test(state, ctx);
    bench(state, ctx);

    evmlr_otse_ctx_clear(ctx);
    flint_rand_clear(state);
    return 0;
}
#endif


