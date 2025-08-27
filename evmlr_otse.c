#include "evmlr_otse.h"
#include "evmlr_utils.h"
#include "flint/fmpz.h"
#include "flint/fmpz_poly.h"
#ifdef MAIN
#include "test.h"
#include "bench.h"
#endif

void calc_a(nmod_poly_t a[M_LEN], const evmlr_otse_key_t key, const evmlr_otse_ctx_t ctx);

void evmlr_otse_ctx_init(evmlr_otse_ctx_t ctx, flint_rand_t state) {
    // H \sample R_{2^ZETA}^{(K_LWE + M_LEN) x K_LWR}
    for (int i = 0; i < K_LWE + M_LEN; i++) {
        for (int j = 0; j < K_LWR; j++) {
            nmod_poly_init(ctx->H[i][j], 2*ZETA);
            nmod_poly_randtest(ctx->H[i][j], state, DEGREE_N);
        }
    }

    // H'' \sample R_q^{M_LEN x K_LWE}
    nmod_poly_t H_dprime[M_LEN][K_LWE];
    for (int i = 0; i < M_LEN; i++) {
        for (int j = 0; j < K_LWE; j++) {
            nmod_poly_init(H_dprime[i][j], MOD_Q);
            nmod_poly_randtest(H_dprime[i][j], state, DEGREE_N);
        }
    }

    // I_{M_LEN}
    nmod_poly_t I[M_LEN][M_LEN];
    for (int i = 0; i < M_LEN; i++) {
        for (int j = 0; j < M_LEN; j++) {
            nmod_poly_init(I[i][j], MOD_Q);
            if (i == j) nmod_poly_one(I[i][j]);
            else        nmod_poly_zero(I[i][j]);
        }
    }

    // H' = (H'' | I_{M_LEN}) \in R_q^{M_LEN x (K_LWE + M_LEN)}
    for (int i = 0; i < M_LEN; i++) {
        for (int j = 0; j < K_LWE + M_LEN; j++) {
            nmod_poly_init(ctx->H_prime[i][j], MOD_Q);
            if (j < K_LWE) {
                nmod_poly_set(ctx->H_prime[i][j], H_dprime[i][j]);
            } else {
                nmod_poly_set(ctx->H_prime[i][j], I[i][j - K_LWE]);
            }
        }
    }

    nmod_poly_init(ctx->cyclo_poly, MOD_Q);
    nmod_poly_set_coeff_ui(ctx->cyclo_poly, DEGREE_N, 1);
    nmod_poly_set_coeff_ui(ctx->cyclo_poly, 0, 1);

    // Clear temporary variables
    for (int i = 0; i < M_LEN; i++) {
        for (int j = 0; j < K_LWE; j++) {
            nmod_poly_clear(H_dprime[i][j]);
        }
        for (int j = 0; j < M_LEN; ++j) {
            nmod_poly_clear(I[i][j]);
        }
    }
}

void evmlr_otse_ctx_clear(evmlr_otse_ctx_t ctx) {
    for (int i = 0; i < K_LWE + M_LEN; i++) {
        for (int j = 0; j < K_LWR; j++) {
            nmod_poly_clear(ctx->H[i][j]);
        }
    }
    for (int i = 0; i < M_LEN; i++) {
        for (int j = 0; j < K_LWE + M_LEN; j++) {
            nmod_poly_clear(ctx->H_prime[i][j]);
        }
    }
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

void evmlr_otse_encrypt(evmlr_otse_ciphertext_t ct, const nmod_poly_t m[M_LEN], const evmlr_otse_key_t key, const evmlr_otse_ctx_t ctx) {
    nmod_poly_t a[M_LEN], tmp;
    nmod_poly_init(tmp, MOD_Q);
    calc_a(a, key, ctx);
    // c = m + a
    for (int i = 0; i < M_LEN; i++) {
        nmod_poly_init(ct->c[i], MOD_Q);
        nmod_poly_add(ct->c[i], m[i], a[i]);
    }

    nmod_poly_clear(tmp);
    for (int i = 0; i < M_LEN; i++) {
        nmod_poly_clear(a[i]);
    }
}

void evmlr_otse_decrypt(nmod_poly_t m[M_LEN], const evmlr_otse_ciphertext_t ct, const evmlr_otse_key_t key, const evmlr_otse_ctx_t ctx) {
    nmod_poly_t a[M_LEN], tmp;
    nmod_poly_init(tmp, MOD_Q);
    calc_a(a, key, ctx);
    // m = c - a
    for (int i = 0; i < M_LEN; i++) {
        nmod_poly_init(m[i], MOD_Q);
        nmod_poly_sub(m[i], ct->c[i], a[i]);
    }
    nmod_poly_clear(tmp);
    for (int i = 0; i < M_LEN; ++i) {
        nmod_poly_clear(a[i]);
    }
}

void evmlr_otse_ciphertext_clear(evmlr_otse_ciphertext_t ct) {
    for (int i = 0; i < M_LEN; i++) {
        nmod_poly_clear(ct->c[i]);
    }
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

void calc_d(nmod_poly_t d[K_LWE + M_LEN], const nmod_poly_t H[K_LWE + M_LEN][K_LWR], const nmod_poly_t s[K_LWR], const nmod_poly_t cyclo_poly) {
    // ⌈Hs⌋_{2^ζ →2^{2η}} := ⌈(Hs * 2^{2η})/2^ζ⌋ mod 2^{2η} := Hs * q2/q1 mod q2
    fmpz_t q1, q2;
    fmpz_init(q1);
    fmpz_init(q2);
    fmpz_set_ui(q1, 1UL << ZETA);      // q1 = 2^ζ
    fmpz_set_ui(q2, 1UL << (2 * ETA)); // q2 = 2^{2η}

    // Hs
    nmod_poly_t Hs[K_LWE + M_LEN], tmp;
    nmod_poly_init(tmp, MOD_Q);
    for (int i = 0; i < K_LWE + M_LEN; i++) {
        nmod_poly_init(Hs[i], MOD_Q);
        nmod_poly_zero(Hs[i]);
        for (int j = 0; j < K_LWR; j++) {
            nmod_poly_mulmod(tmp, H[i][j], s[j], cyclo_poly);
            nmod_poly_add(Hs[i], Hs[i], tmp);
        }
    }
    // d = Hs * q2/q1 mod q2
    for (int i = 0; i < K_LWE + M_LEN; i++) {
        nmod_poly_init(d[i], MOD_Q);
        round_poly(d[i], Hs[i], q1, q2);
    }

    nmod_poly_clear(tmp);
    fmpz_clear(q1);
    fmpz_clear(q2);
    for (int i = 0; i < K_LWE + M_LEN; i++) {
        nmod_poly_clear(Hs[i]);
    }
}

void multiply_by_B_eta(nmod_poly_t w[K_LWE + M_LEN], const nmod_poly_t d[2*ETA][K_LWE + M_LEN], ulong mod) {
    for (int i = 0; i < K_LWE + M_LEN; i++) {
        nmod_poly_init(w[i], mod);
        nmod_poly_zero(w[i]);
    }

    // Accumulate the first eta vectors
    for (int b = 0; b < ETA; b++) {
        for (int i = 0; i < K_LWE + M_LEN; i++) {
            nmod_poly_add(w[i], w[i], d[b][i]);
        }
    }

    // Subtract the last eta vectors
    for (int b = ETA; b < 2 * ETA; b++) {
        for (int i = 0; i < K_LWE + M_LEN; i++) {
            nmod_poly_sub(w[i], w[i], d[b][i]);
        }
    }
}

void calc_a(nmod_poly_t a[M_LEN], const evmlr_otse_key_t key, const evmlr_otse_ctx_t ctx) {
    nmod_poly_t d[K_LWE + M_LEN], d_bin[2 * ETA][K_LWE + M_LEN];
    fmpz_poly_t d_fmpz[K_LWE + M_LEN];

    calc_d(d, ctx->H, key->s, ctx->cyclo_poly);
    for (int i = 0; i < K_LWE + M_LEN; i++) {
        fmpz_poly_init(d_fmpz[i]);
        fmpz_poly_set_nmod_poly(d_fmpz[i], d[i]);
    }

    slong len = K_LWE + M_LEN;
    int bits = 2 * ETA;
    evmlr_utils_ring_to_bin(len, bits, d_bin, d_fmpz, MOD_Q);

    // a = H' B^n stack(d)
    nmod_poly_t w[K_LWE + M_LEN], tmp;
    nmod_poly_init(tmp, MOD_Q);
    // w = B_eta stack(d)
    multiply_by_B_eta(w, d_bin, MOD_Q);

    // a = H' w
    for (int i = 0; i < M_LEN; i++) {
        nmod_poly_init(a[i], MOD_Q);
        nmod_poly_zero(a[i]);
        for (int j = 0; j < K_LWE + M_LEN; j++) {
            nmod_poly_init(tmp, MOD_Q);
            nmod_poly_mulmod(tmp, ctx->H_prime[i][j], w[j], ctx->cyclo_poly);
            nmod_poly_add(a[i], a[i], tmp);
        }
    }

    nmod_poly_clear(tmp);
    for (int i = 0; i < K_LWE + M_LEN; i++) {
        nmod_poly_clear(d[i]);
        fmpz_poly_clear(d_fmpz[i]);
        for (int j = 0; j < 2*ETA; ++j) {
            nmod_poly_clear(d_bin[j][i]);
        }
        nmod_poly_clear(w[i]);
    }
}

#ifdef MAIN
void test(flint_rand_t rand, evmlr_otse_ctx_t ctx) {
    evmlr_otse_key_t key;
    evmlr_otse_keygen(key, rand);

    nmod_poly_t msg[M_LEN], decrypted_msg[M_LEN];
    for (int i = 0; i < M_LEN; i++) {
        nmod_poly_init(msg[i], MOD_Q);
        nmod_poly_randtest(msg[i], rand, DEGREE_N);
    }

    evmlr_otse_ciphertext_t ct;
    TEST_BEGIN("encryption and decryption are consistent") {
        evmlr_otse_encrypt(ct, msg, key, ctx);
        evmlr_otse_decrypt(decrypted_msg, ct, key, ctx);
        for (int i = 0; i < M_LEN; i++) {
            TEST_ASSERT(nmod_poly_equal(msg[i], decrypted_msg[i]) == 1, end)
        }
    } TEST_END;

    end:
        evmlr_otse_ciphertext_clear(ct);
        evmlr_otse_keyclear(key);
        for (int i = 0; i < M_LEN; i++) {
            nmod_poly_clear(msg[i]);
            nmod_poly_clear(decrypted_msg[i]);
        }
}

void bench(flint_rand_t rand, evmlr_otse_ctx_t ctx) {
    evmlr_otse_key_t key;
    BENCH_BEGIN("evmlr_otse_keygen") {
        BENCH_ADD(evmlr_otse_keygen(key, rand))
    } BENCH_END;

    nmod_poly_t msg[M_LEN];
    for (int i = 0; i < M_LEN; i++) {
        nmod_poly_init(msg[i], MOD_Q);
        nmod_poly_randtest(msg[i], rand, DEGREE_N);
    }

    evmlr_otse_ciphertext_t ct;
    BENCH_BEGIN("evmlr_otse_encrypt") {
        BENCH_ADD(evmlr_otse_encrypt(ct, msg, key, ctx))
    } BENCH_END;

    nmod_poly_t decrypted_msg[M_LEN];
    BENCH_BEGIN("evmlr_otse_decrypt") {
        BENCH_ADD(evmlr_otse_decrypt(decrypted_msg, ct, key, ctx))
    } BENCH_END;

    evmlr_otse_keyclear(key);
    for (int i = 0; i < M_LEN; i++) {
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
    evmlr_otse_ctx_init(ctx, state);

    test(state, ctx);
    bench(state, ctx);

    evmlr_otse_ctx_clear(ctx);
    flint_rand_clear(state);
    return 0;
}
#endif


