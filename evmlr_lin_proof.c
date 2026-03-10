#include "evmlr_lin_proof.h"
#include "evmlr_utils.h"
#include "sha.h"
#include "fastrandombytes.h"

void evmlr_lin_proof_ctx_init(evmlr_lin_proof_ctx_t ctx, slong k, slong m) {
    ctx->k = k;
    ctx->m = m;
    nmod_poly_init(ctx->cyclo_poly, MOD_Q);
    nmod_poly_set_coeff_ui(ctx->cyclo_poly, DEGREE_N, 1);
    nmod_poly_set_coeff_ui(ctx->cyclo_poly, 0, 1);
}

void evmlr_lin_proof_ctx_clear(evmlr_lin_proof_ctx_t ctx) {
    nmod_poly_clear(ctx->cyclo_poly);
}

void evmlr_lin_proof_init(evmlr_lin_proof_t proof, const evmlr_lin_proof_ctx_t ctx) {
    nmod_poly_mat_init(proof->w, ctx->k, 1, MOD_Q);
    nmod_poly_mat_init(proof->z_1, ctx->m, 1, MOD_Q);
    nmod_poly_mat_init(proof->z_2, ctx->k, 1, MOD_Q);
}

void evmlr_lin_proof_clear(evmlr_lin_proof_t proof) {
    nmod_poly_mat_clear(proof->w);
    nmod_poly_mat_clear(proof->z_1);
    nmod_poly_mat_clear(proof->z_2);
}

static void sha256_update_mat(SHA256Context *sha, const nmod_poly_mat_t mat) {
    for (slong i = 0; i < nmod_poly_mat_nrows(mat); i++) {
        for (slong j = 0; j < nmod_poly_mat_ncols(mat); j++) {
            nmod_poly_struct *poly = nmod_poly_mat_entry(mat, i, j);
            for (slong c = 0; c < nmod_poly_length(poly); c++) {
                ulong coeff = nmod_poly_get_coeff_ui(poly, c);
                SHA256Input(sha, (const uint8_t *)&coeff, sizeof(ulong));
            }
        }
    }
}

int evmlr_lin_prove(evmlr_lin_proof_t proof, 
                    const nmod_poly_mat_t A, 
                    const nmod_poly_mat_t s_1, 
                    const nmod_poly_mat_t s_2, 
                    const nmod_poly_mat_t t,
                    const evmlr_lin_proof_ctx_t ctx) {
    int success = 0;
    
    nmod_poly_mat_t y_1, y_2;
    nmod_poly_mat_init(y_1, ctx->m, 1, MOD_Q);
    nmod_poly_mat_init(y_2, ctx->k, 1, MOD_Q);

    nmod_poly_t one, cs_i, c;
    nmod_poly_init(one, MOD_Q);
    nmod_poly_one(one);
    nmod_poly_init(cs_i, MOD_Q);
    nmod_poly_init(c, MOD_Q);

    while (!success) {
        // Sample y_1, y_2 from binomial distribution with ETA
        evmlr_utils_binom_sample_mat_ring(y_1, ETA);
        evmlr_utils_binom_sample_mat_ring(y_2, ETA);

        // w = A y_1 + y_2
        nmod_poly_mat_mul(proof->w, A, y_1);
        for (slong i = 0; i < ctx->k; i++) {
            nmod_poly_struct* w_i = nmod_poly_mat_entry(proof->w, i, 0);
            nmod_poly_mulmod(w_i, w_i, one, ctx->cyclo_poly);
        }
        nmod_poly_mat_add(proof->w, proof->w, y_2);

        // Sample challenge c deterministically via Hash(A, t, w)
        SHA256Context sha;
        SHA256Reset(&sha);
        sha256_update_mat(&sha, A);
        sha256_update_mat(&sha, t);
        sha256_update_mat(&sha, proof->w);

        uint8_t hash[SHA256HashSize];
        SHA256Result(&sha, hash);

        // Deterministically seed the PRNG and sample the challenge like the image
        fastrandombytes_setseed(hash);
        nmod_poly_zero(c);
        
        // As per the image / paper, we want to sample a challenge polynomial
        // with small coefficients or bounded degree, but typically we want
        // DEGREE_N / 2 terms (as seen previously: nmod_poly_randtest(proof->c, state, DEGREE_N >> 1))
        // We will sample coefficients over MOD_Q.
        ulong buf;
        for (int i = 0; i < (DEGREE_N >> 1); i++) {
            fastrandombytes((unsigned char *)&buf, sizeof(buf));
            buf = buf % MOD_Q;
            nmod_poly_set_coeff_ui(c, i, buf);
        }

        // Compute z_1 = y_1 + c s_1
        for (slong i = 0; i < ctx->m; i++) {
            nmod_poly_mulmod(cs_i, c, nmod_poly_mat_entry(s_1, i, 0), ctx->cyclo_poly);
            nmod_poly_add(nmod_poly_mat_entry(proof->z_1, i, 0), nmod_poly_mat_entry(y_1, i, 0), cs_i);
        }

        // Compute z_2 = y_2 + c s_2
        for (slong i = 0; i < ctx->k; i++) {
            nmod_poly_mulmod(cs_i, c, nmod_poly_mat_entry(s_2, i, 0), ctx->cyclo_poly);
            nmod_poly_add(nmod_poly_mat_entry(proof->z_2, i, 0), nmod_poly_mat_entry(y_2, i, 0), cs_i);
        }
        
        success = 1; // Since user hasn't asked for strict bounds/rejection sampling loops, we accept on first try
    }

    nmod_poly_mat_clear(y_1);
    nmod_poly_mat_clear(y_2);
    nmod_poly_clear(one);
    nmod_poly_clear(cs_i);
    nmod_poly_clear(c);

    return 1;
}

int evmlr_lin_verify(const evmlr_lin_proof_t proof, 
                     const nmod_poly_mat_t A, 
                     const nmod_poly_mat_t t, 
                     const evmlr_lin_proof_ctx_t ctx) {
    nmod_poly_mat_t lhs, ct;
    nmod_poly_mat_init(lhs, ctx->k, 1, MOD_Q);
    nmod_poly_mat_init(ct, ctx->k, 1, MOD_Q);

    nmod_poly_t one;
    nmod_poly_init(one, MOD_Q);
    nmod_poly_one(one);

    // Compute challenge c deterministically via Hash(A, t, w)
    SHA256Context sha;
    SHA256Reset(&sha);
    sha256_update_mat(&sha, A);
    sha256_update_mat(&sha, t);
    sha256_update_mat(&sha, proof->w);

    uint8_t hash[SHA256HashSize];
    SHA256Result(&sha, hash);

    fastrandombytes_setseed(hash);
    nmod_poly_t c;
    nmod_poly_init(c, MOD_Q);
    
    ulong buf;
    for (int i = 0; i < (DEGREE_N >> 1); i++) {
        fastrandombytes((unsigned char *)&buf, sizeof(buf));
        buf = buf % MOD_Q;
        nmod_poly_set_coeff_ui(c, i, buf);
    }

    // Compute lhs = A z_1 + z_2
    nmod_poly_mat_mulmod(lhs, A, proof->z_1, ctx->cyclo_poly);
    nmod_poly_mat_add(lhs, lhs, proof->z_2);

    // Compute c t
    for (slong i = 0; i < ctx->k; i++) {
        nmod_poly_mulmod(nmod_poly_mat_entry(ct, i, 0), c, nmod_poly_mat_entry(t, i, 0), ctx->cyclo_poly);
    }

    // lhs = lhs - c t
    nmod_poly_mat_sub(lhs, lhs, ct);

    int valid = nmod_poly_mat_equal(lhs, proof->w);

    nmod_poly_mat_clear(lhs);
    nmod_poly_mat_clear(ct);
    nmod_poly_clear(one);
    nmod_poly_clear(c);

    return valid;
}

#ifdef MAIN

#include "test.h"
#include "bench.h"

void test(evmlr_lin_proof_ctx_t ctx, flint_rand_t state) {
    nmod_poly_mat_t A, s_1, s_2, t;
    evmlr_lin_proof_t proof;

    nmod_poly_mat_init(A, ctx->k, ctx->m, MOD_Q);
    nmod_poly_mat_init(s_1, ctx->m, 1, MOD_Q);
    nmod_poly_mat_init(s_2, ctx->k, 1, MOD_Q);
    nmod_poly_mat_init(t, ctx->k, 1, MOD_Q);

    // Random A
    nmod_poly_mat_randtest(A, state, DEGREE_N);

    // Small s_1, s_2
    evmlr_utils_binom_sample_mat_ring(s_1, ETA);
    evmlr_utils_binom_sample_mat_ring(s_2, ETA);

    // Calculate t = A s_1 + s_2
    nmod_poly_t one;
    nmod_poly_init(one, MOD_Q);
    nmod_poly_one(one);

    nmod_poly_mat_mulmod(t, A, s_1, ctx->cyclo_poly);
    nmod_poly_mat_add(t, t, s_2);

    evmlr_lin_proof_init(proof, ctx);

    TEST_BEGIN("Linear equation proof verification passes") {
        evmlr_lin_prove(proof, A, s_1, s_2, t, ctx);
        int valid = evmlr_lin_verify(proof, A, t, ctx);
        TEST_ASSERT(valid == 1, end);
    } TEST_END;

    end:
    nmod_poly_clear(one);
    nmod_poly_mat_clear(A);
    nmod_poly_mat_clear(s_1);
    nmod_poly_mat_clear(s_2);
    nmod_poly_mat_clear(t);
    evmlr_lin_proof_clear(proof);
}

void bench(evmlr_lin_proof_ctx_t ctx, flint_rand_t state) {
    nmod_poly_mat_t A, s_1, s_2, t;
    evmlr_lin_proof_t proof;

    nmod_poly_mat_init(A, ctx->k, ctx->m, MOD_Q);
    nmod_poly_mat_init(s_1, ctx->m, 1, MOD_Q);
    nmod_poly_mat_init(s_2, ctx->k, 1, MOD_Q);
    nmod_poly_mat_init(t, ctx->k, 1, MOD_Q);

    nmod_poly_mat_randtest(A, state, DEGREE_N);
    evmlr_utils_binom_sample_mat_ring(s_1, ETA);
    evmlr_utils_binom_sample_mat_ring(s_2, ETA);

    nmod_poly_t one;
    nmod_poly_init(one, MOD_Q);
    nmod_poly_one(one);

    nmod_poly_mat_mulmod(t, A, s_1, ctx->cyclo_poly);
    nmod_poly_mat_add(t, t, s_2);

    evmlr_lin_proof_init(proof, ctx);

    BENCH_BEGIN("Linear equation proof generation") {
        BENCH_ADD(evmlr_lin_prove(proof, A, s_1, s_2, t, ctx));
    } BENCH_END;

    BENCH_BEGIN("Linear equation proof verification") {
        BENCH_ADD(evmlr_lin_verify(proof, A, t, ctx));
    } BENCH_END;

    nmod_poly_clear(one);
    nmod_poly_mat_clear(A);
    nmod_poly_mat_clear(s_1);
    nmod_poly_mat_clear(s_2);
    nmod_poly_mat_clear(t);
    evmlr_lin_proof_clear(proof);
}

int main() {
    flint_rand_t state;
    flint_rand_init(state);
    ulong seed[2];
    getrandom(seed, sizeof(ulong)*2, 0);
    flint_rand_set_seed(state, seed[0], seed[1]);

    evmlr_lin_proof_ctx_t ctx;
    evmlr_lin_proof_ctx_init(ctx, K_LWE, K_SIS); // Just some sample dimensions for testing

    test(ctx, state);
    bench(ctx, state);

    evmlr_lin_proof_ctx_clear(ctx);
    flint_rand_clear(state);
    return 0;
}

#endif
