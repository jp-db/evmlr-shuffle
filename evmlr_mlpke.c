#include "evmlr_mlpke.h"
#include "evmlr_utils.h"
#ifdef MAIN
#include "test.h"
#include "bench.h"
#endif


void evmlr_mlpke_ctx_init(evmlr_mlpke_ctx_t ctx) {
    nmod_poly_init(ctx->cyclo_poly, MOD_Q);
    // Initialize cyclotomic polynomial x^N + 1
    nmod_poly_set_coeff_ui(ctx->cyclo_poly, DEGREE_N, 1);
    nmod_poly_set_coeff_ui(ctx->cyclo_poly, 0, 1);
}

void evmlr_mlpke_ctx_clear(evmlr_mlpke_ctx_t ctx) {
    nmod_poly_clear(ctx->cyclo_poly);
}

void evmlr_mlpke_keypair_gen(evmlr_mlpke_keypair_t keypair, flint_rand_t state, const evmlr_mlpke_ctx_t ctx) {
    // Generate secret key s and error e with coefficients from B_eta
    nmod_poly_t e[K_LWE];
    for (size_t i = 0; i < K_LWE; i++) {
        nmod_poly_init(keypair->sk->s[i], MOD_Q);
        nmod_poly_init(e[i], MOD_Q);
        evmlr_utils_binom_sample_ring(keypair->sk->s[i], ETA);
        evmlr_utils_binom_sample_ring(e[i], ETA);
    }

    // Generate public matrix A with uniform random coefficients
    for (size_t i = 0; i < K_LWE; i++) {
        for (size_t j = 0; j < K_LWE; j++) {
            nmod_poly_init(keypair->pk->A[i][j], MOD_Q);
            nmod_poly_randtest(keypair->pk->A[i][j], state, DEGREE_N);
        }
    }

    // Compute t = A * s + e
    nmod_poly_t tmp;
    nmod_poly_init(tmp, MOD_Q);
    for (size_t i = 0; i < K_LWE; i++) {
        nmod_poly_init(keypair->pk->t[i], MOD_Q);
        nmod_poly_zero(keypair->pk->t[i]);
        for (size_t j = 0; j < K_LWE; j++) {
            nmod_poly_mulmod(tmp, keypair->pk->A[i][j], keypair->sk->s[j], ctx->cyclo_poly);
            nmod_poly_add(keypair->pk->t[i], keypair->pk->t[i], tmp);
        }
        nmod_poly_add(keypair->pk->t[i], keypair->pk->t[i], e[i]);
    }

    nmod_poly_clear(tmp);
    for (size_t i = 0; i < K_LWE; i++) {
        nmod_poly_clear(e[i]);
    }
}

void evmlr_mlpke_keypair_clear(evmlr_mlpke_keypair_t keypair) {
    // Clear secret key
    for (size_t i = 0; i < K_LWE; i++) {
        nmod_poly_clear(keypair->sk->s[i]);
    }
    // Clear public key
    for (size_t i = 0; i < K_LWE; i++) {
        for (size_t j = 0; j < K_LWE; j++) {
            nmod_poly_clear(keypair->pk->A[i][j]);
        }
        nmod_poly_clear(keypair->pk->t[i]);
    }
}

void evmlr_mlpke_enc(evmlr_mlpke_cipher_t cipher, const nmod_poly_t msg, const evmlr_mlpke_pk_t pk, const evmlr_mlpke_ctx_t ctx) {
    // assert(evmlr_utils_is_poly_bin(msg));
    nmod_poly_t r[K_LWE], e2[K_LWE], e3, tmp, msg_scaled;
    // Sample r and e2 from B^{k_lwe}_{eta}
    for (size_t i = 0; i < K_LWE; i++) {
        nmod_poly_init(r[i], MOD_Q);
        nmod_poly_init(e2[i], MOD_Q);
        evmlr_utils_binom_sample_ring(r[i], ETA);
        evmlr_utils_binom_sample_ring(e2[i], ETA);
    }

    // Sample e3 from B_{eta}
    nmod_poly_init(e3, MOD_Q);
    evmlr_utils_binom_sample_ring(e3, ETA);

    nmod_poly_init(tmp, MOD_Q);
    // Compute u^T = r^T A + e_2^T
    for (size_t i = 0; i < K_LWE; i++) {
        nmod_poly_init(cipher->uT[i], MOD_Q);
        nmod_poly_zero(cipher->uT[i]);
        for (size_t j = 0; j < K_LWE; j++) {
            nmod_poly_mulmod(tmp, r[j], pk->A[j][i], ctx->cyclo_poly);
            nmod_poly_add(cipher->uT[i], cipher->uT[i], tmp);
        }
        nmod_poly_add(cipher->uT[i], cipher->uT[i], e2[i]);
    }
    // c = q/2 getting the closest integer with ties rounded up
    ulong c = (MOD_Q + 1) / 2;
    // Compute v = r^T t + e_3 + c * msg
    nmod_poly_init(cipher->v, MOD_Q);
    nmod_poly_zero(cipher->v);
    for (size_t i = 0; i < K_LWE; i++) {
        nmod_poly_mulmod(tmp, r[i], pk->t[i], ctx->cyclo_poly);
        nmod_poly_add(cipher->v, cipher->v, tmp);
    }
    nmod_poly_add(cipher->v, cipher->v, e3);

    nmod_poly_init(msg_scaled, MOD_Q);
    nmod_poly_scalar_mul_nmod(msg_scaled, msg, c);

    nmod_poly_add(cipher->v, cipher->v, msg_scaled);

    nmod_poly_clear(tmp);
    nmod_poly_clear(msg_scaled);
    for (size_t i = 0; i < K_LWE; i++) {
        nmod_poly_clear(r[i]);
        nmod_poly_clear(e2[i]);
    }
    nmod_poly_clear(e3);
}

void evmlr_mlpke_dec(nmod_poly_t msg, const evmlr_mlpke_cipher_t cipher, const evmlr_mlpke_sk_t sk, const evmlr_mlpke_ctx_t ctx) {
    // \tilde{m} = v - u^T s
    nmod_poly_t m_tilde, tmp;
    nmod_poly_init(m_tilde, MOD_Q);
    nmod_poly_init(tmp, MOD_Q);
    nmod_poly_set(m_tilde, cipher->v);
    for (size_t i = 0; i < K_LWE; i++) {
        nmod_poly_mulmod(tmp, cipher->uT[i], sk->s[i], ctx->cyclo_poly);
        nmod_poly_sub(m_tilde, m_tilde, tmp);
    }

    nmod_poly_zero(msg);
    // msg = 0 if \tilde{m} is closer to 0 than to q/2, and 1 otherwise
    ulong threshold = (MOD_Q) / 4;
    for (slong i = 0; i < DEGREE_N; i++) {
        ulong coeff = nmod_poly_get_coeff_ui(m_tilde, i);
        if (coeff > threshold && coeff < (MOD_Q - threshold)) {
            nmod_poly_set_coeff_ui(msg, i, 1); // closer to q/2
        } else {
            nmod_poly_set_coeff_ui(msg, i, 0); // closer to 0
        }
    }

    nmod_poly_clear(m_tilde);
    nmod_poly_clear(tmp);
}

void evmlr_mlpke_cipher_clear(evmlr_mlpke_cipher_t cipher) {
    for (size_t i = 0; i < K_LWE; i++) {
        nmod_poly_clear(cipher->uT[i]);
    }
    nmod_poly_clear(cipher->v);
}

#ifdef MAIN
static void test(flint_rand_t rand, evmlr_mlpke_ctx_t ctx) {
    evmlr_mlpke_keypair_t keypair;
    evmlr_mlpke_cipher_t cipher;
    evmlr_mlpke_keypair_gen(keypair, rand, ctx);

    nmod_poly_t msg, decrypted_msg;
    nmod_poly_init(msg, MOD_Q);
    nmod_poly_init(decrypted_msg, MOD_Q);

    ulong msg_val[1];
    getrandom(msg_val, sizeof(ulong), 0);
    evmlr_utils_int_to_bin(msg, msg_val[0]);

    TEST_BEGIN("encryption and decryption are consistent") {
        evmlr_mlpke_enc(cipher, msg, keypair->pk, ctx);
        evmlr_mlpke_dec(decrypted_msg, cipher, keypair->sk, ctx);
        TEST_ASSERT(nmod_poly_equal(msg, decrypted_msg) == 1, end)
    } TEST_END;

    // clear memory
    end:
        evmlr_mlpke_cipher_clear(cipher);
        evmlr_mlpke_keypair_clear(keypair);
        nmod_poly_clear(msg);
        nmod_poly_clear(decrypted_msg);
}

static void bench(flint_rand_t rand, evmlr_mlpke_ctx_t ctx) {
    evmlr_mlpke_keypair_t keypair;
    evmlr_mlpke_keypair_gen(keypair, rand, ctx);

    nmod_poly_t msg, decrypted_msg;
    nmod_poly_init(msg, MOD_Q);
    nmod_poly_init(decrypted_msg, MOD_Q);

    evmlr_mlpke_cipher_t cipher;

    BENCH_BEGIN("evmlr_mlpke_keypair_gen") {
        BENCH_ADD(evmlr_mlpke_keypair_gen(keypair, rand, ctx))
    } BENCH_END;

    ulong msg_val[1];
    getrandom(msg_val, sizeof(ulong), 0);
    evmlr_utils_int_to_bin(msg, msg_val[0]);

    BENCH_BEGIN("evmlr_mlpke_enc") {
        BENCH_ADD(evmlr_mlpke_enc(cipher, msg, keypair->pk, ctx))
    } BENCH_END;

    BENCH_BEGIN("evmlr_mlpke_dec") {
        BENCH_ADD(evmlr_mlpke_dec(decrypted_msg, cipher, keypair->sk, ctx))
    } BENCH_END;

    // clear memory
    evmlr_mlpke_cipher_clear(cipher);
    evmlr_mlpke_keypair_clear(keypair);
    nmod_poly_clear(msg);
    nmod_poly_clear(decrypted_msg);
}


int main() {
    flint_rand_t state;
    flint_rand_init(state);
    ulong seed[2];
    getrandom(seed, sizeof(ulong)*2, 0);
    flint_rand_set_seed(state, seed[0], seed[1]);

    evmlr_mlpke_ctx_t ctx;
    evmlr_mlpke_ctx_init(ctx);

    test(state, ctx);
    bench(state, ctx);

    evmlr_mlpke_ctx_clear(ctx);
    flint_rand_clear(state);

    return 0;
}
#endif