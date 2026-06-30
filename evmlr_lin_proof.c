#include "evmlr_lin_proof.h"
#include "evmlr_utils.h"
#include "evmlr_challenge.h"

static void evmlr_lin_sample_mask_binomial(nmod_poly_mat_t y, const struct evmlr_lin_proof_ctx_struct* ctx) {
    evmlr_utils_binom_sample_mat_ring(y, ctx->eta);
}

static int evmlr_lin_rejection_check_binomial(const nmod_poly_mat_t z_1, const nmod_poly_mat_t z_2, const struct evmlr_lin_proof_ctx_struct* ctx) {
    return 1;
}

static int evmlr_lin_norm_check_binomial(const nmod_poly_mat_t z_1, const nmod_poly_mat_t z_2, const struct evmlr_lin_proof_ctx_struct* ctx) {
    return 1;
}

static void evmlr_lin_sample_mask_linear(nmod_poly_mat_t y, const struct evmlr_lin_proof_ctx_struct* ctx) {
    evmlr_utils_linear_sample_mat_ring(y, ctx->beta);
}

static int evmlr_lin_rejection_check_linear(const nmod_poly_mat_t z_1, const nmod_poly_mat_t z_2, const struct evmlr_lin_proof_ctx_struct* ctx) {
    if (!evmlr_utils_is_bounded(z_1, ctx->beta)) return 0;
    if (z_2 != NULL && !evmlr_utils_is_bounded(z_2, ctx->beta)) return 0;
    return 1;
}

static int evmlr_lin_norm_check_linear(const nmod_poly_mat_t z_1, const nmod_poly_mat_t z_2, const struct evmlr_lin_proof_ctx_struct* ctx) {
    return evmlr_lin_rejection_check_linear(z_1, z_2, ctx);
}

void evmlr_lin_proof_ctx_init(evmlr_lin_proof_ctx_t ctx, slong k, slong m) {
    ctx->k = k;
    ctx->m = m;
    nmod_poly_init(ctx->cyclo_poly, MOD_Q);
    nmod_poly_set_coeff_ui(ctx->cyclo_poly, DEGREE_N, 1);
    nmod_poly_set_coeff_ui(ctx->cyclo_poly, 0, 1);
    
    ctx->eta = ETA;
    ctx->beta = 0;
    ctx->sample_mask = evmlr_lin_sample_mask_binomial;
    ctx->rejection_check = evmlr_lin_rejection_check_binomial;
    ctx->norm_check = evmlr_lin_norm_check_binomial;
}

void evmlr_lin_proof_ctx_set_linear(evmlr_lin_proof_ctx_t ctx, slong beta) {
    ctx->beta = beta;
    ctx->sample_mask = evmlr_lin_sample_mask_linear;
    ctx->rejection_check = evmlr_lin_rejection_check_linear;
    ctx->norm_check = evmlr_lin_norm_check_linear;
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
        // Sample y_1, y_2 using context's sample_mask
        ctx->sample_mask(y_1, ctx);
        ctx->sample_mask(y_2, ctx);
 
        // w = A y_1 + y_2
        nmod_poly_mat_mul(proof->w, A, y_1);
        for (slong i = 0; i < ctx->k; i++) {
            nmod_poly_struct* w_i = nmod_poly_mat_entry(proof->w, i, 0);
            nmod_poly_mulmod(w_i, w_i, one, ctx->cyclo_poly);
        }
        nmod_poly_mat_add(proof->w, proof->w, y_2);
 
        evmlr_challenge_t chal;
        evmlr_challenge_init(chal);
        evmlr_challenge_add_matrix(chal, A);
        evmlr_challenge_add_matrix(chal, t);
        evmlr_challenge_add_matrix(chal, proof->w);
        evmlr_challenge_get_poly_half(c, chal);

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
        
        success = ctx->rejection_check(proof->z_1, proof->z_2, ctx);
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
 
    nmod_poly_t c;
    nmod_poly_init(c, MOD_Q);

    evmlr_challenge_t chal;
    evmlr_challenge_init(chal);
    evmlr_challenge_add_matrix(chal, A);
    evmlr_challenge_add_matrix(chal, t);
    evmlr_challenge_add_matrix(chal, proof->w);
    evmlr_challenge_get_poly_half(c, chal);
 
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
    if (valid) {
        valid = ctx->norm_check(proof->z_1, proof->z_2, ctx);
    }
 
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

    // Test with linear sampling
    evmlr_lin_proof_ctx_set_linear(ctx, 2000);
    
    // Resample t = A s_1 + s_2
    nmod_poly_mat_mulmod(t, A, s_1, ctx->cyclo_poly);
    nmod_poly_mat_add(t, t, s_2);

    TEST_BEGIN("Linear equation proof verification (linear/uniform sampling) passes") {
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
