#include "evmlr_commit.h"
#include "evmlr_utils.h"
#ifdef MAIN
#include "test.h"
#include "bench.h"
#endif

// B^2 := nL + n 2 K_SIS ETA^2
ulong compute_b() {
    return (DEGREE_N * M_LEN) + (DEGREE_N * 2 * K_SIS * ETA * ETA);
}

void evmlr_commit_scheme_init(evmlr_commit_ctx_t ctx, flint_rand_t state) {
    ctx->b_sqr = compute_b();

    // A_1 \sample R_q^{K_SIS x M_LEN}
    for (size_t i = 0; i < K_SIS; i++) {
        for (size_t j = 0; j < M_LEN; j++) {
            nmod_poly_init(ctx->A_1[i][j], MOD_Q);
            nmod_poly_randtest(ctx->A_1[i][j], state, DEGREE_N);
        }
    }
    // A_2 \sample R_q^{K_SIS x 2K_SIS}
    for (size_t i = 0; i < K_SIS; i++) {
        for (size_t j = 0; j < 2 * K_SIS; j++) {
            nmod_poly_init(ctx->A_2[i][j], MOD_Q);
            nmod_poly_randtest(ctx->A_2[i][j], state, DEGREE_N);
        }
    }
    nmod_poly_init(ctx->cyclo_poly, MOD_Q);
    nmod_poly_set_coeff_ui(ctx->cyclo_poly, DEGREE_N, 1);
    nmod_poly_set_coeff_ui(ctx->cyclo_poly, 0, 1);
}

void evmlr_commit_scheme_clear(evmlr_commit_ctx_t ctx) {
    // Clear A_1
    for (size_t i = 0; i < K_SIS; i++) {
        for (size_t j = 0; j < M_LEN; j++) {
            nmod_poly_clear(ctx->A_1[i][j]);
        }
    }
    // Clear A_2
    for (size_t i = 0; i < K_SIS; i++) {
        for (size_t j = 0; j < 2 * K_SIS; j++) {
            nmod_poly_clear(ctx->A_2[i][j]);
        }
    }
}

void evmlr_sample_r(nmod_poly_t r[2*K_SIS]) {
    for (size_t i = 0; i < 2 * K_SIS; i++) {
        nmod_poly_init(r[i], MOD_Q);
        evmlr_utils_binom_sample_ring(r[i], ETA);
    }
}

void evmlr_commit(evmlr_commit_t com, const nmod_poly_t msg[M_LEN], const nmod_poly_t r[2 * K_SIS], const evmlr_commit_ctx_t ctx) {
    // Initialize commitment polynomials
    for (size_t i = 0; i < K_SIS; i++) {
        nmod_poly_init(com->c[i], MOD_Q);
        nmod_poly_zero(com->c[i]);
    }

    nmod_poly_t tmp;
    nmod_poly_init(tmp, MOD_Q);

    // Compute c = A_1 * msg + A_2 * r
    for (size_t i = 0; i < K_SIS; i++) {
        // A_1 * msg
        for (size_t j = 0; j < M_LEN; j++) {
            nmod_poly_mul(tmp, ctx->A_1[i][j], msg[j]);
            nmod_poly_add(com->c[i], com->c[i], tmp);
        }
        // A_2 * r
        for (size_t j = 0; j < 2 * K_SIS; j++) {
            nmod_poly_mul(tmp, ctx->A_2[i][j], r[j]);
            nmod_poly_add(com->c[i], com->c[i], tmp);
        }
    }
    nmod_poly_clear(tmp);
}

void evmlr_commit_clear(evmlr_commit_t com) {
    for (size_t i = 0; i < K_SIS; i++) {
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

int verify_norm(const nmod_poly_t msg[M_LEN], const nmod_poly_t r[2 * K_SIS], ulong b_sqr) {
    ulong total_norm = 0;
    for (size_t i = 0; i < M_LEN; i++) {
        total_norm += nmod_poly_hamming_weight(msg[i]);
    }
    for (size_t i = 0; i < 2 * K_SIS; i++) {
        total_norm += nmod_poly_norm_sqr(r[i]);
    }
    return total_norm <= b_sqr;
}

int evmlr_commit_verify(const evmlr_commit_t com, const nmod_poly_t msg[M_LEN], const nmod_poly_t r[2 * K_SIS], const evmlr_commit_ctx_t ctx) {
    // ||\binom{m}{r}|| <= B
    int valid = verify_norm(msg, r, ctx->b_sqr);
    if (!valid) {
        return 0;
    }
    // c = A_1 * msg + A_2 * r
    evmlr_commit_t recomputed_com;
    evmlr_commit(recomputed_com, msg, r, ctx);

    for (int i = 0; i < K_SIS && valid; i++) {
        valid &= nmod_poly_equal(com->c[i], recomputed_com->c[i]);
    }

    evmlr_commit_clear(recomputed_com);
    return valid;
}

#ifdef MAIN

void test(evmlr_commit_ctx_t ctx) {
    ulong msg_value[M_LEN];
    nmod_poly_t msg[M_LEN], r[2 * K_SIS];
    evmlr_commit_t com;
    // Sample message
    getrandom(msg_value, sizeof(ulong) * M_LEN, 0);
    for (size_t i = 0; i < M_LEN; i++) {
        nmod_poly_init(msg[i], MOD_Q);
        evmlr_utils_int_to_bin(msg[i], msg_value[i]);
    }
    // Sample r
    evmlr_sample_r(r);


    TEST_BEGIN("commitments can be created and verified") {
        evmlr_commit(com, msg, r, ctx);
        int valid = evmlr_commit_verify(com, msg, r, ctx);
        TEST_ASSERT(valid == 1, end)
    } TEST_END;

    end:
        for (size_t i = 0; i < M_LEN; i++) {
            nmod_poly_clear(msg[i]);
        }
        for (size_t i = 0; i < 2 * K_SIS; i++) {
            nmod_poly_clear(r[i]);
        }
        evmlr_commit_clear(com);
}

void bench(evmlr_commit_ctx_t ctx) {
    ulong msg_value[M_LEN];
    nmod_poly_t msg[M_LEN], r[2 * K_SIS];
    evmlr_commit_t com;
    // Sample message
    getrandom(msg_value, sizeof(ulong) * M_LEN, 0);
    for (size_t i = 0; i < M_LEN; i++) {
        nmod_poly_init(msg[i], MOD_Q);
        evmlr_utils_int_to_bin(msg[i], msg_value[i]);
    }
    // Sample r
    evmlr_sample_r(r);

    BENCH_BEGIN("evmlr_commit") {
        BENCH_ADD(evmlr_commit(com, msg, r, ctx))
        evmlr_commit_clear(com);
    } BENCH_END;

    BENCH_BEGIN("evmlr_commit_verify") {
        evmlr_commit(com, msg, r, ctx);
        BENCH_ADD(evmlr_commit_verify(com, msg, r, ctx))
    } BENCH_END;

    for (size_t i = 0; i < M_LEN; i++) {
        nmod_poly_clear(msg[i]);
    }
    for (size_t i = 0; i < 2 * K_SIS; i++) {
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
    evmlr_commit_scheme_init(ctx, state);

    test(ctx);
    bench(ctx);

    evmlr_commit_scheme_clear(ctx);
    flint_rand_clear(state);
    return 0;
}
#endif