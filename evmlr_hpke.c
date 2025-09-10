#include "evmlr_hpke.h"
#ifdef MAIN
#include "test.h"
#include "bench.h"
#endif

void evmlr_hpke_ctx_init(evmlr_hpke_ctx_t ctx, slong L, flint_rand_t state) {
    evmlr_mlpke_ctx_init(ctx->enc_ctx);
    evmlr_otse_ctx_init(ctx->otse_ctx, L, state);
}

void evmlr_hpke_ctx_clear(evmlr_hpke_ctx_t ctx) {
    evmlr_mlpke_ctx_clear(ctx->enc_ctx);
    evmlr_otse_ctx_clear(ctx->otse_ctx);
}

void evmlr_hpke_keypair_gen(evmlr_hpke_keypair_t keypair, flint_rand_t state, const evmlr_hpke_ctx_t ctx) {
    evmlr_mlpke_keypair_gen(keypair->enc_keypair, state, ctx->enc_ctx);
}

void evmlr_hpke_keypair_clear(evmlr_hpke_keypair_t keypair) {
    evmlr_mlpke_keypair_clear(keypair->enc_keypair);
}

void evmlr_hpke_encrypt(evmlr_hpke_cipher_t cipher, nmod_poly_mat_t d_dagger, const nmod_poly_mat_t msg,
                        const evmlr_mlpke_pk_t pk, const evmlr_hpke_ctx_t ctx, flint_rand_t state) {
    evmlr_otse_key_t key;
    evmlr_otse_keygen(key, state);
    for (int i = 0; i < K_LWR; i++) {
        evmlr_mlpke_enc(cipher->enc_cipher[i], nmod_poly_mat_entry(key->s, i, 0), pk, ctx->enc_ctx);
    }
    evmlr_otse_encrypt(cipher->otse_cipher, d_dagger, msg, key, ctx->otse_ctx);
    evmlr_otse_keyclear(key);
}

void evmlr_hpke_decrypt(nmod_poly_mat_t msg, const evmlr_hpke_cipher_t cipher, const evmlr_mlpke_sk_t sk, const evmlr_hpke_ctx_t ctx) {
    evmlr_otse_key_t key;
    nmod_poly_t key_poly;
    nmod_poly_init(key_poly, MOD_Q);

    nmod_poly_mat_init(key->s, K_LWR, 1, MOD_Q);
    for (int i = 0; i < K_LWR; i++) {
        evmlr_mlpke_dec(key_poly, cipher->enc_cipher[i], sk, ctx->enc_ctx);
        nmod_poly_set(nmod_poly_mat_entry(key->s, i, 0), key_poly);
    }
    evmlr_otse_decrypt(msg, cipher->otse_cipher, key, ctx->otse_ctx);

    evmlr_otse_keyclear(key);
    nmod_poly_clear(key_poly);
}

void evmlr_hpke_cipher_clear(evmlr_hpke_cipher_t cipher) {
    for (int i = 0; i < K_LWR; i++) {
        evmlr_mlpke_cipher_clear(cipher->enc_cipher[i]);
    }
    evmlr_otse_ciphertext_clear(cipher->otse_cipher);
}

#ifdef MAIN
static void test(flint_rand_t rand, evmlr_hpke_ctx_t ctx) {
    evmlr_hpke_keypair_t keypair;
    evmlr_hpke_keypair_gen(keypair, rand, ctx);
    nmod_poly_mat_t msg, decrypted_msg;
    nmod_poly_mat_init(msg, ctx->otse_ctx->L, 1, MOD_Q);
    nmod_poly_mat_init(decrypted_msg, ctx->otse_ctx->L, 1, MOD_Q);
    nmod_poly_mat_randtest(msg, rand, DEGREE_N);

    evmlr_hpke_cipher_t cipher;
    TEST_BEGIN("encryption and decryption are consistent") {
        evmlr_hpke_encrypt(cipher, nullptr, msg, keypair->enc_keypair->pk, ctx, rand);
        evmlr_hpke_decrypt(decrypted_msg, cipher, keypair->enc_keypair->sk, ctx);
        TEST_ASSERT(nmod_poly_mat_equal(msg, decrypted_msg) == 1, end)
    } TEST_END;

    // clear memory
    end:
        evmlr_hpke_cipher_clear(cipher);
        evmlr_hpke_keypair_clear(keypair);
        nmod_poly_mat_clear(msg);
        nmod_poly_mat_clear(decrypted_msg);
}

static void bench(flint_rand_t rand, evmlr_hpke_ctx_t ctx) {
    evmlr_hpke_keypair_t keypair;
    evmlr_hpke_keypair_gen(keypair, rand, ctx);

    nmod_poly_mat_t msg, decrypted_msg;
    nmod_poly_mat_init(msg, ctx->otse_ctx->L, 1, MOD_Q);
    nmod_poly_mat_init(decrypted_msg, ctx->otse_ctx->L, 1, MOD_Q);
    nmod_poly_mat_randtest(msg, rand, DEGREE_N);

    evmlr_hpke_cipher_t cipher;

    BENCH_BEGIN("evmlr_hpke_keypair_gen") {
        BENCH_ADD(evmlr_hpke_keypair_gen(keypair, rand, ctx))
    } BENCH_END;

    BENCH_BEGIN("evmlr_hpke_encrypt") {
        BENCH_ADD(evmlr_hpke_encrypt(cipher, nullptr, msg, keypair->enc_keypair->pk, ctx, rand))
    } BENCH_END;

    BENCH_BEGIN("evmlr_hpke_decrypt") {
        BENCH_ADD(evmlr_hpke_decrypt(decrypted_msg, cipher, keypair->enc_keypair->sk, ctx))
    } BENCH_END;

    // clear memory
    evmlr_hpke_cipher_clear(cipher);
    evmlr_hpke_keypair_clear(keypair);
    nmod_poly_mat_clear(msg);
    nmod_poly_mat_clear(decrypted_msg);
}

int main() {
    flint_rand_t state;
    flint_rand_init(state);
    ulong seed[2];
    getrandom(seed, sizeof(ulong)*2, 0);
    flint_rand_set_seed(state, seed[0], seed[1]);

    evmlr_hpke_ctx_t ctx;
    evmlr_hpke_ctx_init(ctx, M_LEN, state);

    test(state, ctx);
    bench(state, ctx);

    evmlr_hpke_ctx_clear(ctx);
    flint_rand_clear(state);

    return 0;
}
#endif
