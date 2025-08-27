#include "evmlr_shuffle.h"

void evmlr_shuffle_ctx_init(evmlr_shuffle_ctx_t ctx, flint_rand_t state) {
    evmlr_commit_ctx_init(ctx->com_ctx, state);
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
