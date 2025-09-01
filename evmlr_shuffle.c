#include "evmlr_shuffle.h"
#include "evmlr_utils.h"

void evmlr_shuffle_ctx_init(evmlr_shuffle_ctx_t ctx, ulong L, flint_rand_t state) {
    evmlr_commit_ctx_init(ctx->com_ctx, L, state);
    evmlr_hpke_ctx_init(ctx->hpke_ctx, state);
}

void evmlr_shuffle_ctx_clear(evmlr_shuffle_ctx_t ctx) {
    evmlr_commit_ctx_clear(ctx->com_ctx);
    evmlr_hpke_ctx_clear(ctx->hpke_ctx);
}

void evmlr_shuffle_sample_challenge(nmod_poly_t chal, flint_rand_t state) {
    nmod_poly_init(chal, MOD_Q);
    nmod_poly_randtest(chal, state, DEGREE_N >> 1);
}

void evmlr_shuffle_phase_1(nmod_poly_t sigma[], evmlr_commit_t com, nmod_poly_t r[2*K_SIS], const evmlr_shuffle_ctx_t ctx, const evmlr_shuffle_sp_t sp) {
    // Sample r from B^{2K_SIS}
    evmlr_sample_r(r);
    for (ulong i = 0; i < sp->N; i++) {
        evmlr_utils_int_to_bin(sigma[i], sp->pi[i]);
    }
    evmlr_commit(com, ctx->com_ctx, sigma, r);
}

void evmlr_shuffle_phase_2(evmlr_commit_t com, nmod_poly_t r[2*K_SIS], const evmlr_shuffle_ctx_t ctx, const evmlr_shuffle_sp_t sp) {

}
