#include "evmlr_commit.h"
#include "evmlr_utils.h"
#ifdef MAIN
#include "test.h"
#include "bench.h"
#endif

// B^2 := nL + n 2 K_SIS ETA^2
ulong compute_b(ulong L) {
    return (DEGREE_N * L) + (DEGREE_N * 2 * K_SIS * ETA * ETA);
}

void evmlr_commit_ctx_init(evmlr_commit_ctx_t ctx, ulong L, flint_rand_t state) {
    ctx->b_sqr = compute_b(L);
    ctx->L = L;

    // A_1 \sample R_q^{K_SIS x L}
    for (ulong i = 0; i < K_SIS; i++) {
        ctx->A_1[i] = (nmod_poly_t*) malloc(sizeof(nmod_poly_t) * L);
        for (ulong j = 0; j < L; j++) {
            nmod_poly_init(ctx->A_1[i][j], MOD_Q);
            nmod_poly_randtest(ctx->A_1[i][j], state, DEGREE_N);
        }
    }
    // A_2 \sample R_q^{K_SIS x 2K_SIS}
    for (ulong i = 0; i < K_SIS; i++) {
        for (ulong j = 0; j < 2 * K_SIS; j++) {
            nmod_poly_init(ctx->A_2[i][j], MOD_Q);
            nmod_poly_randtest(ctx->A_2[i][j], state, DEGREE_N);
        }
    }
    nmod_poly_init(ctx->cyclo_poly, MOD_Q);
    nmod_poly_set_coeff_ui(ctx->cyclo_poly, DEGREE_N, 1);
    nmod_poly_set_coeff_ui(ctx->cyclo_poly, 0, 1);
}

void evmlr_commit_ctx_clear(evmlr_commit_ctx_t ctx) {
    // Clear A_1
    for (ulong i = 0; i < K_SIS; i++) {
        for (ulong j = 0; j < ctx->L; j++) {
            nmod_poly_clear(ctx->A_1[i][j]);
        }
        free(ctx->A_1[i]);
    }
    // Clear A_2
    for (ulong i = 0; i < K_SIS; i++) {
        for (ulong j = 0; j < 2 * K_SIS; j++) {
            nmod_poly_clear(ctx->A_2[i][j]);
        }
    }
}

void evmlr_sample_r(nmod_poly_t r[2*K_SIS]) {
    for (ulong i = 0; i < 2 * K_SIS; i++) {
        nmod_poly_init(r[i], MOD_Q);
        evmlr_utils_binom_sample_ring(r[i], ETA);
    }
}

void evmlr_commit(evmlr_commit_t com, const evmlr_commit_ctx_t ctx, const nmod_poly_t msg[ctx->L], const nmod_poly_t r[2 * K_SIS]) {
    // Initialize commitment polynomials
    for (ulong i = 0; i < K_SIS; i++) {
        nmod_poly_init(com->c[i], MOD_Q);
        nmod_poly_zero(com->c[i]);
    }

    nmod_poly_t tmp;
    nmod_poly_init(tmp, MOD_Q);

    // Compute c = A_1 * msg + A_2 * r
    for (ulong i = 0; i < K_SIS; i++) {
        // A_1 * msg
        for (ulong j = 0; j < ctx->L; j++) {
            nmod_poly_mul(tmp, ctx->A_1[i][j], msg[j]);
            nmod_poly_add(com->c[i], com->c[i], tmp);
        }
        // A_2 * r
        for (ulong j = 0; j < 2 * K_SIS; j++) {
            nmod_poly_mul(tmp, ctx->A_2[i][j], r[j]);
            nmod_poly_add(com->c[i], com->c[i], tmp);
        }
    }
    nmod_poly_clear(tmp);
}

void evmlr_commit_clear(evmlr_commit_t com) {
    for (ulong i = 0; i < K_SIS; i++) {
        nmod_poly_clear(com->c[i]);
    }
}

ulong nmod_poly_norm_sqr(const nmod_poly_t poly) {
   ulong norm = 0;
   for (slong i = 0; i < nmod_poly_length(poly); i++) {
         ulong coeff = nmod_poly_get_coeff_ui(poly, i);
         if (coeff > MOD_Q / 2) {
            coeff -= MOD_Q;
         }
         norm += coeff * coeff;
   }
   return norm;
}

int verify_norm(const ulong L, const nmod_poly_t msg[L], const nmod_poly_t r[2 * K_SIS], ulong b_sqr) {
    ulong total_norm = 0;
    for (ulong i = 0; i < L; i++) {
        total_norm += nmod_poly_hamming_weight(msg[i]);
    }
    for (ulong i = 0; i < 2 * K_SIS; i++) {
        total_norm += nmod_poly_norm_sqr(r[i]);
    }
    return total_norm <= b_sqr;
}

int evmlr_commit_verify(const evmlr_commit_t com, const evmlr_commit_ctx_t ctx, const nmod_poly_t msg[ctx->L], const nmod_poly_t r[2 * K_SIS]) {
    // ||\binom{m}{r}|| <= B
    int valid = verify_norm(ctx->L, msg, r, ctx->b_sqr);
    if (!valid) {
        return 0;
    }
    // c = A_1 * msg + A_2 * r
    evmlr_commit_t recomputed_com;
    evmlr_commit(recomputed_com, ctx, msg, r);

    for (int i = 0; i < K_SIS && valid; i++) {
        valid &= nmod_poly_equal(com->c[i], recomputed_com->c[i]);
    }

    evmlr_commit_clear(recomputed_com);
    return valid;
}

#ifdef MAIN

void test(evmlr_commit_ctx_t ctx) {
    ulong msg_value[ctx->L];
    nmod_poly_t msg[ctx->L], r[2 * K_SIS];
    evmlr_commit_t com;
    // Sample message
    getrandom(msg_value, sizeof(ulong) * ctx->L, 0);
    for (ulong i = 0; i < ctx->L; i++) {
        nmod_poly_init(msg[i], MOD_Q);
        evmlr_utils_int_to_bin(msg[i], msg_value[i]);
    }
    // Sample r
    evmlr_sample_r(r);


    TEST_BEGIN("commitments can be created and verified") {
        evmlr_commit(com, ctx, msg, r);
        int valid = evmlr_commit_verify(com, ctx, msg, r);
        TEST_ASSERT(valid == 1, end)
    } TEST_END;

    end:
        for (ulong i = 0; i < ctx->L; i++) {
            nmod_poly_clear(msg[i]);
        }
        for (ulong i = 0; i < 2 * K_SIS; i++) {
            nmod_poly_clear(r[i]);
        }
        evmlr_commit_clear(com);
}

void bench(evmlr_commit_ctx_t ctx) {
    ulong msg_value[ctx->L];
    nmod_poly_t msg[ctx->L], r[2 * K_SIS];
    evmlr_commit_t com;
    // Sample message
    getrandom(msg_value, sizeof(ulong) * ctx->L, 0);
    for (ulong i = 0; i < ctx->L; i++) {
        nmod_poly_init(msg[i], MOD_Q);
        evmlr_utils_int_to_bin(msg[i], msg_value[i]);
    }
    // Sample r
    evmlr_sample_r(r);

    BENCH_BEGIN("evmlr_commit") {
        BENCH_ADD(evmlr_commit(com, ctx, msg, r))
        evmlr_commit_clear(com);
    } BENCH_END;

    BENCH_BEGIN("evmlr_commit_verify") {
        evmlr_commit(com, ctx, msg, r);
        BENCH_ADD(evmlr_commit_verify(com, ctx, msg, r))
    } BENCH_END;

    for (ulong i = 0; i < ctx->L; i++) {
        nmod_poly_clear(msg[i]);
    }
    for (ulong i = 0; i < 2 * K_SIS; i++) {
        nmod_poly_clear(r[i]);
    }
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