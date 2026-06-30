#include "evmlr_lin_proof.h"
#include "evmlr_utils.h"
#include "sha.h"
#include "fastrandombytes.h"

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
    nmod_poly_mat_init(proof->z_2, ctx->k, 1, MOD_Q); // z_2 stores the Hint matrix "h"
}

void evmlr_lin_proof_clear(evmlr_lin_proof_t proof) {
    nmod_poly_mat_clear(proof->w);
    nmod_poly_mat_clear(proof->z_1);
    nmod_poly_mat_clear(proof->z_2);
}

static void sha256_update_mat(SHA256Context *sha, const nmod_poly_mat_t mat) {
    for (slong i = 0; i < mat->r; i++) {
        for (slong j = 0; j < mat->c; j++) {
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
    
    nmod_poly_mat_t y, Ay;
    nmod_poly_mat_init(y, ctx->m, 1, MOD_Q);
    nmod_poly_mat_init(Ay, ctx->k, 1, MOD_Q);

    nmod_poly_t one, cs_i, c;
    nmod_poly_init(one, MOD_Q);
    nmod_poly_one(one);
    nmod_poly_init(cs_i, MOD_Q);
    nmod_poly_init(c, MOD_Q);
    
    nmod_poly_mat_t t_1, t_0, Az, ct_1, ct_0, Az_minus_ct1;
    nmod_poly_mat_init(t_1, ctx->k, 1, MOD_Q);
    nmod_poly_mat_init(t_0, ctx->k, 1, MOD_Q);
    nmod_poly_mat_init(Az, ctx->k, 1, MOD_Q);
    nmod_poly_mat_init(ct_1, ctx->k, 1, MOD_Q);
    nmod_poly_mat_init(ct_0, ctx->k, 1, MOD_Q);
    nmod_poly_mat_init(Az_minus_ct1, ctx->k, 1, MOD_Q);

    // Compute t_1 = HIGHS(t) and t_0 = LOWS(t)
    evmlr_utils_highs_mat(t_1, t);
    evmlr_utils_lows_mat(t_0, t);

    while (!success) {
        ctx->sample_mask(y, ctx);

        // w = HIGHS(A y)
        nmod_poly_mat_mul(Ay, A, y);
        for (slong i = 0; i < ctx->k; i++) {
            nmod_poly_mulmod(nmod_poly_mat_entry(Ay, i, 0), nmod_poly_mat_entry(Ay, i, 0), one, ctx->cyclo_poly);
        }
        evmlr_utils_highs_mat(proof->w, Ay);

        // Sample challenge c deterministically via Hash(A, t_1, w)
        SHA256Context sha;
        SHA256Reset(&sha);
        sha256_update_mat(&sha, A);
        sha256_update_mat(&sha, t_1);
        sha256_update_mat(&sha, proof->w);

        uint8_t hash[SHA256HashSize];
        SHA256Result(&sha, hash);

        fastrandombytes_setseed(hash);
        nmod_poly_zero(c);
        
        int count = 0;
        while (count < 1) {
            ulong idx;
            fastrandombytes((unsigned char *)&idx, sizeof(ulong));
            idx %= DEGREE_N;
            if (nmod_poly_get_coeff_ui(c, idx) == 0) {
                ulong sign;
                fastrandombytes((unsigned char *)&sign, 1);
                nmod_poly_set_coeff_ui(c, idx, (sign & 1) ? 1 : MOD_Q - 1);
                count++;
            }
        }

        // Compute z_1 = y + c s_1
        for (slong i = 0; i < ctx->m; i++) {
            nmod_poly_mulmod(cs_i, c, nmod_poly_mat_entry(s_1, i, 0), ctx->cyclo_poly);
            nmod_poly_add(nmod_poly_mat_entry(proof->z_1, i, 0), nmod_poly_mat_entry(y, i, 0), cs_i);
        }

        // Rejection Sampling directly testing bounded limits
        nmod_poly_mat_t Ay_minus_cs2, Ay_highs, Ay_sub_highs;
        nmod_poly_mat_init(Ay_minus_cs2, ctx->k, 1, MOD_Q);
        nmod_poly_mat_init(Ay_highs, ctx->k, 1, MOD_Q);
        nmod_poly_mat_init(Ay_sub_highs, ctx->k, 1, MOD_Q);
        
        for (slong i = 0; i < ctx->k; i++) {
            nmod_poly_t cs2_i;
            nmod_poly_init(cs2_i, MOD_Q);
            nmod_poly_mulmod(cs2_i, c, nmod_poly_mat_entry(s_2, i, 0), ctx->cyclo_poly);
            nmod_poly_sub(nmod_poly_mat_entry(Ay_minus_cs2, i, 0), nmod_poly_mat_entry(Ay, i, 0), cs2_i);
            nmod_poly_clear(cs2_i);
        }
        
        evmlr_utils_highs_mat(Ay_highs, Ay);
        evmlr_utils_highs_mat(Ay_sub_highs, Ay_minus_cs2);
        
        int reject = !nmod_poly_mat_equal(Ay_highs, Ay_sub_highs);

        nmod_poly_mat_clear(Ay_minus_cs2);
        nmod_poly_mat_clear(Ay_highs);
        nmod_poly_mat_clear(Ay_sub_highs);
        
        if (reject) {
            continue;
        }

        if (!ctx->rejection_check(proof->z_1, NULL, ctx)) {
            continue;
        }

        // Compute z_2 = HINT(Az_1 - ct_1, ct_0)
        nmod_poly_mat_mul(Az, A, proof->z_1);
        for (slong i = 0; i < ctx->k; i++) {
            nmod_poly_mulmod(nmod_poly_mat_entry(Az, i, 0), nmod_poly_mat_entry(Az, i, 0), one, ctx->cyclo_poly);
            
            nmod_poly_mulmod(nmod_poly_mat_entry(ct_1, i, 0), c, nmod_poly_mat_entry(t_1, i, 0), ctx->cyclo_poly);
            nmod_poly_mulmod(nmod_poly_mat_entry(ct_0, i, 0), c, nmod_poly_mat_entry(t_0, i, 0), ctx->cyclo_poly);
            
            nmod_poly_sub(nmod_poly_mat_entry(Az_minus_ct1, i, 0), nmod_poly_mat_entry(Az, i, 0), nmod_poly_mat_entry(ct_1, i, 0));
        }

        evmlr_utils_make_hint_mat(proof->z_2, Az_minus_ct1, ct_0);
        
        success = 1;
    }

    nmod_poly_mat_clear(y);
    nmod_poly_mat_clear(Ay);
    nmod_poly_clear(one);
    nmod_poly_clear(cs_i);
    nmod_poly_clear(c);
    nmod_poly_mat_clear(t_1);
    nmod_poly_mat_clear(t_0);
    nmod_poly_mat_clear(Az);
    nmod_poly_mat_clear(ct_1);
    nmod_poly_mat_clear(ct_0);
    nmod_poly_mat_clear(Az_minus_ct1);

    return 1;
}

int evmlr_lin_verify(const evmlr_lin_proof_t proof, 
                          const nmod_poly_mat_t A, 
                          const nmod_poly_mat_t t, 
                          const evmlr_lin_proof_ctx_t ctx) {
    nmod_poly_mat_t t_1;
    nmod_poly_mat_init(t_1, ctx->k, 1, MOD_Q);
    evmlr_utils_highs_mat(t_1, t);

    SHA256Context sha;
    SHA256Reset(&sha);
    sha256_update_mat(&sha, A);
    sha256_update_mat(&sha, t_1);
    sha256_update_mat(&sha, proof->w);

    uint8_t hash[SHA256HashSize];
    SHA256Result(&sha, hash);

    fastrandombytes_setseed(hash);
    
    nmod_poly_t c, one;
    nmod_poly_init(c, MOD_Q);
    nmod_poly_init(one, MOD_Q);
    nmod_poly_one(one);
    
    int count = 0;
    while (count < 1) {
        ulong idx;
        fastrandombytes((unsigned char *)&idx, sizeof(ulong));
        idx %= DEGREE_N;
        if (nmod_poly_get_coeff_ui(c, idx) == 0) {
            ulong sign;
            fastrandombytes((unsigned char *)&sign, 1);
            nmod_poly_set_coeff_ui(c, idx, (sign & 1) ? 1 : MOD_Q - 1);
            count++;
        }
    }

    nmod_poly_mat_t Az, ct_1, Az_minus_ct1, w_prime;
    nmod_poly_mat_init(Az, ctx->k, 1, MOD_Q);
    nmod_poly_mat_init(ct_1, ctx->k, 1, MOD_Q);
    nmod_poly_mat_init(Az_minus_ct1, ctx->k, 1, MOD_Q);
    nmod_poly_mat_init(w_prime, ctx->k, 1, MOD_Q);

    nmod_poly_mat_mul(Az, A, proof->z_1);
    
    for (slong i = 0; i < ctx->k; i++) {
        nmod_poly_mulmod(nmod_poly_mat_entry(Az, i, 0), nmod_poly_mat_entry(Az, i, 0), one, ctx->cyclo_poly);
        nmod_poly_mulmod(nmod_poly_mat_entry(ct_1, i, 0), c, nmod_poly_mat_entry(t_1, i, 0), ctx->cyclo_poly);
        nmod_poly_sub(nmod_poly_mat_entry(Az_minus_ct1, i, 0), nmod_poly_mat_entry(Az, i, 0), nmod_poly_mat_entry(ct_1, i, 0));
    }

    evmlr_utils_use_hint_mat(w_prime, proof->z_2, Az_minus_ct1);

    int valid = nmod_poly_mat_equal(w_prime, proof->w);
    if (valid) {
        valid = ctx->norm_check(proof->z_1, NULL, ctx);
    }

    nmod_poly_clear(c);
    nmod_poly_clear(one);
    nmod_poly_mat_clear(t_1);
    nmod_poly_mat_clear(Az);
    nmod_poly_mat_clear(ct_1);
    nmod_poly_mat_clear(Az_minus_ct1);
    nmod_poly_mat_clear(w_prime);

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

    nmod_poly_mat_randtest(A, state, DEGREE_N);

    evmlr_utils_binom_sample_mat_ring(s_1, ETA);
    evmlr_utils_binom_sample_mat_ring(s_2, ETA);

    nmod_poly_t one;
    nmod_poly_init(one, MOD_Q);
    nmod_poly_one(one);

    nmod_poly_mat_mulmod(t, A, s_1, ctx->cyclo_poly);
    nmod_poly_mat_add(t, t, s_2);
    
    evmlr_lin_proof_init(proof, ctx);

    TEST_BEGIN("Compressed Hint linear proof verification passes") {
        evmlr_lin_prove(proof, A, s_1, s_2, t, ctx);
        int valid = evmlr_lin_verify(proof, A, t, ctx);
        TEST_ASSERT(valid == 1, end_comp);
    } TEST_END;

    // Test with linear sampling
    evmlr_lin_proof_ctx_set_linear(ctx, 2000);

    // Resample t = A s_1 + s_2
    nmod_poly_mat_mulmod(t, A, s_1, ctx->cyclo_poly);
    nmod_poly_mat_add(t, t, s_2);

    TEST_BEGIN("Compressed Hint linear proof verification (linear/uniform sampling) passes") {
        evmlr_lin_prove(proof, A, s_1, s_2, t, ctx);
        int valid = evmlr_lin_verify(proof, A, t, ctx);
        TEST_ASSERT(valid == 1, end_comp);
    } TEST_END;

    end_comp:
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

    BENCH_BEGIN("Compressed Hint linear proof generation") {
        BENCH_ADD(evmlr_lin_prove(proof, A, s_1, s_2, t, ctx));
    } BENCH_END;

    BENCH_BEGIN("Compressed Hint linear proof verification") {
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
    evmlr_lin_proof_ctx_init(ctx, K_LWE, K_SIS);

    test(ctx, state);
    bench(ctx, state);

    evmlr_lin_proof_ctx_clear(ctx);
    flint_rand_clear(state);
    return 0;
}

#endif
