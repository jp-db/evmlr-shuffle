#include "evmlr_commit.h"
#include "evmlr_utils.h"
#include "evmlr_crt.h"
#ifdef MAIN
#include "test.h"
#include "bench.h"
#endif

// B^2 := nL + n 2 K_SIS ETA^2
ulong compute_b(ulong L) {
    return (DEGREE_N * L) + (DEGREE_N * 2 * K_SIS * ETA * ETA);
}

void evmlr_commit_ctx_init(evmlr_commit_ctx_t ctx, slong N, flint_rand_t state) {
    ctx->b_sqr = compute_b(N);
    ctx->N = N;

    // A_1 \sample R_q^{K_SIS x N}
    nmod_poly_mat_init(ctx->A_1, K_SIS, N, MOD_Q);
    nmod_poly_mat_randtest(ctx->A_1, state, DEGREE_N);

    // A_2 \sample R_q^{K_SIS x 2K_SIS}
    nmod_poly_mat_init(ctx->A_2, K_SIS, 2 * K_SIS, MOD_Q);
    nmod_poly_mat_randtest(ctx->A_2, state, DEGREE_N);

    nmod_poly_init(ctx->cyclo_poly, MOD_Q);
    nmod_poly_set_coeff_ui(ctx->cyclo_poly, DEGREE_N, 1);
    nmod_poly_set_coeff_ui(ctx->cyclo_poly, 0, 1);
    
    evmlr_crt_setup(ctx->crt_ctx, ctx->cyclo_poly);
}

void evmlr_commit_ctx_clear(evmlr_commit_ctx_t ctx) {
    evmlr_crt_clear(ctx->crt_ctx);
    nmod_poly_mat_clear(ctx->A_1);
    nmod_poly_mat_clear(ctx->A_2);
}

void evmlr_commit_sample_r(nmod_poly_mat_t r) {
    nmod_poly_mat_init(r, 2 * K_SIS, 1, MOD_Q);
    evmlr_utils_binom_sample_mat_ring(r, ETA);
}

void evmlr_commit(evmlr_commit_t com, const nmod_poly_mat_t msg, const nmod_poly_mat_t r, const evmlr_commit_ctx_t ctx) {
    // Initialize commitment polynomials
    nmod_poly_mat_init(com->c, K_SIS, 1, MOD_Q);

    evmlr_poly_mat_crt_t crt_A1, crt_msg, crt_A2, crt_r, crt_c, crt_tmp;
    evmlr_poly_mat_crt_init(crt_A1, K_SIS, ctx->N, MOD_Q);
    evmlr_poly_mat_crt_init(crt_msg, ctx->N, 1, MOD_Q);
    evmlr_poly_mat_crt_init(crt_A2, K_SIS, 2 * K_SIS, MOD_Q);
    evmlr_poly_mat_crt_init(crt_r, 2 * K_SIS, 1, MOD_Q);
    evmlr_poly_mat_crt_init(crt_c, K_SIS, 1, MOD_Q);
    evmlr_poly_mat_crt_init(crt_tmp, K_SIS, 1, MOD_Q);

    // Compute c = A_1 * msg + A_2 * r
    evmlr_poly_mat_to_crt(crt_A1, ctx->A_1, ctx->crt_ctx);
    evmlr_poly_mat_to_crt(crt_msg, msg, ctx->crt_ctx);
    evmlr_poly_mat_crt_mul(crt_c, crt_A1, crt_msg, ctx->crt_ctx);
    evmlr_poly_mat_from_crt(com->c, crt_c, ctx->crt_ctx);

    evmlr_poly_mat_to_crt(crt_A2, ctx->A_2, ctx->crt_ctx);
    evmlr_poly_mat_to_crt(crt_r, r, ctx->crt_ctx);
    evmlr_poly_mat_crt_mul(crt_tmp, crt_A2, crt_r, ctx->crt_ctx);

    nmod_poly_mat_t tmp;
    nmod_poly_mat_init(tmp, K_SIS, 1, MOD_Q);
    evmlr_poly_mat_from_crt(tmp, crt_tmp, ctx->crt_ctx);

    nmod_poly_mat_add(com->c, com->c, tmp);

    evmlr_poly_mat_crt_clear(crt_A1);
    evmlr_poly_mat_crt_clear(crt_msg);
    evmlr_poly_mat_crt_clear(crt_A2);
    evmlr_poly_mat_crt_clear(crt_r);
    evmlr_poly_mat_crt_clear(crt_c);
    evmlr_poly_mat_crt_clear(crt_tmp);
    nmod_poly_mat_clear(tmp);
}

void evmlr_commit_clear(evmlr_commit_t com) {
    nmod_poly_mat_clear(com->c);
}

slong nmod_poly_norm_sqr(const nmod_poly_t poly) {
   slong norm = 0;
   for (slong i = 0; i < nmod_poly_length(poly); i++) {
         slong coeff = nmod_poly_get_coeff_ui(poly, i);
         if (coeff > MOD_Q / 2) {
            coeff -= MOD_Q;
         }
         norm += coeff * coeff;
   }
   return norm;
}

int verify_norm(const ulong L, const nmod_poly_mat_t msg, const nmod_poly_mat_t r, ulong b_sqr) {
    ulong total_norm = 0;
    for (slong i = 0; i < L; i++) {
        total_norm += nmod_poly_hamming_weight(nmod_poly_mat_entry(msg, i, 0));
    }
    for (slong i = 0; i < 2 * K_SIS; i++) {
        total_norm += nmod_poly_norm_sqr(nmod_poly_mat_entry(r, i, 0));
    }
    return total_norm <= b_sqr;
}

int evmlr_commit_verify(const evmlr_commit_t com, const nmod_poly_mat_t msg, const nmod_poly_mat_t r, const evmlr_commit_ctx_t ctx) {
    // ||\binom{m}{r}|| <= B
    int valid = verify_norm(ctx->N, msg, r, ctx->b_sqr);
    if (!valid) {
        return 0;
    }
    // c = A_1 * msg + A_2 * r
    evmlr_commit_t recomputed_com;
    evmlr_commit(recomputed_com, msg, r, ctx);
    valid = nmod_poly_mat_equal(com->c, recomputed_com->c);

    evmlr_commit_clear(recomputed_com);
    return valid;
}

#ifdef MAIN

void test(evmlr_commit_ctx_t ctx) {
    ulong msg_value[ctx->N];
    nmod_poly_mat_t msg, r;
    evmlr_commit_t com;
    // Sample message
    getrandom(msg_value, sizeof(ulong) * ctx->N, 0);
    nmod_poly_mat_init(msg, ctx->N, 1, MOD_Q);
    for (slong i = 0; i < ctx->N; i++) {
        evmlr_utils_int_to_bin(nmod_poly_mat_entry(msg, i, 0), msg_value[i]);
    }
    // Sample r
    evmlr_commit_sample_r(r);


    TEST_BEGIN("commitments can be created and verified") {
        evmlr_commit(com, msg, r, ctx);
        int valid = evmlr_commit_verify(com, msg, r, ctx);
        TEST_ASSERT(valid == 1, end)
    } TEST_END;

    end:
        nmod_poly_mat_clear(msg);
        nmod_poly_mat_clear(r);
        evmlr_commit_clear(com);
}

void bench(evmlr_commit_ctx_t ctx) {
    ulong msg_value[ctx->N];
    nmod_poly_mat_t msg, r;
    evmlr_commit_t com;
    // Sample message
    getrandom(msg_value, sizeof(ulong) * ctx->N, 0);
    nmod_poly_mat_init(msg, ctx->N, 1, MOD_Q);
    for (slong i = 0; i < ctx->N; i++) {
        evmlr_utils_int_to_bin(nmod_poly_mat_entry(msg, i, 0), msg_value[i]);
    }
    // Sample r
    evmlr_commit_sample_r(r);

    BENCH_BEGIN("evmlr_commit") {
        BENCH_ADD(evmlr_commit(com, msg, r, ctx))
        evmlr_commit_clear(com);
    } BENCH_END;

    BENCH_BEGIN("evmlr_commit_verify") {
        evmlr_commit(com, msg, r, ctx);
        BENCH_ADD(evmlr_commit_verify(com,msg, r, ctx))
    } BENCH_END;

    nmod_poly_mat_clear(msg);
    nmod_poly_mat_clear(r);
    evmlr_commit_clear(com);
}

int main() {
    flint_rand_t state;
    flint_rand_init(state);
    ulong seed[2];
    getrandom(seed, sizeof(ulong)*2, 0);
    flint_rand_set_seed(state, seed[0], seed[1]);

    evmlr_commit_ctx_t ctx;
    evmlr_commit_ctx_init(ctx, M_LEN, state);

    test(ctx);
    bench(ctx);

    evmlr_commit_ctx_clear(ctx);
    flint_rand_clear(state);
    return 0;
}
#endif