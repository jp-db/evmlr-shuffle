#include "evmlr_crt.h"
#include "evmlr_utils.h"
#include "flint/ulong_extras.h"

void evmlr_crt_setup(evmlr_crt_ctx_t ctx, const nmod_poly_t cyclo_poly) {
    ctx->mod = cyclo_poly->mod.n;

    nmod_poly_factor_init(ctx->factors);
    // Factor X^512 + 1 mod Q
    nmod_poly_factor(ctx->factors, cyclo_poly);

    ctx->C_0 = nmod_poly_get_coeff_ui(ctx->factors->p + 0, 0);
    ctx->C_1 = nmod_poly_get_coeff_ui(ctx->factors->p + 1, 0);

    ulong diff = n_submod(ctx->C_0, ctx->C_1, ctx->mod);
    ctx->inv_C0_sub_C1 = n_invmod(diff, ctx->mod);

    // C0_inv_diff = C_0 * inv_C0_sub_C1 mod Q
    ctx->C0_inv_diff = (ctx->C_0 * ctx->inv_C0_sub_C1) % ctx->mod;
    // C1_inv_diff = C_1 * inv_C0_sub_C1 mod Q
    ctx->C1_inv_diff = (ctx->C_1 * ctx->inv_C0_sub_C1) % ctx->mod;
}

void evmlr_crt_clear(evmlr_crt_ctx_t ctx) {
    nmod_poly_factor_clear(ctx->factors);
}

void evmlr_poly_crt_init(evmlr_poly_crt_t crt, ulong n) {
    nmod_poly_init(crt->poly + 0, n);
    nmod_poly_init(crt->poly + 1, n);
}

void evmlr_poly_crt_clear(evmlr_poly_crt_t crt) {
    nmod_poly_clear(crt->poly + 0);
    nmod_poly_clear(crt->poly + 1);
}

void evmlr_poly_mat_crt_init(evmlr_poly_mat_crt_t crt, slong rows, slong cols, ulong n) {
    nmod_poly_mat_init(crt->mat + 0, rows, cols, n);
    nmod_poly_mat_init(crt->mat + 1, rows, cols, n);
}

void evmlr_poly_mat_crt_clear(evmlr_poly_mat_crt_t crt) {
    nmod_poly_mat_clear(crt->mat + 0);
    nmod_poly_mat_clear(crt->mat + 1);
}

static inline void _evmlr_poly_to_crt(nmod_poly_struct* dest0, nmod_poly_struct* dest1, const nmod_poly_struct* src, const evmlr_crt_ctx_t ctx) {
    ulong Q = ctx->mod;
    nmod_poly_zero(dest0);
    nmod_poly_zero(dest1);
    nmod_poly_fit_length(dest0, 256);
    nmod_poly_fit_length(dest1, 256);
    for (int i = 0; i < 256; i++) {
        ulong a_i = nmod_poly_get_coeff_ui(src, i);
        ulong a_i_256 = nmod_poly_get_coeff_ui(src, i + 256);
        
        ulong C0_a = (ctx->C_0 * a_i_256) % Q;
        ulong r0 = (a_i >= C0_a) ? (a_i - C0_a) : (a_i + Q - C0_a);
        
        ulong C1_a = (ctx->C_1 * a_i_256) % Q;
        ulong r1 = (a_i >= C1_a) ? (a_i - C1_a) : (a_i + Q - C1_a);

        nmod_poly_set_coeff_ui(dest0, i, r0);
        nmod_poly_set_coeff_ui(dest1, i, r1);
    }
}

void evmlr_poly_to_crt(evmlr_poly_crt_t dest, const nmod_poly_t src, const evmlr_crt_ctx_t ctx) {
    _evmlr_poly_to_crt(dest->poly + 0, dest->poly + 1, src, ctx);
}

static inline void _evmlr_poly_from_crt(nmod_poly_struct* dest, const nmod_poly_struct* src0, const nmod_poly_struct* src1, const evmlr_crt_ctx_t ctx) {
    ulong Q = ctx->mod;
    nmod_poly_zero(dest);
    nmod_poly_fit_length(dest, 512);
    for (int i = 0; i < 256; i++) {
        ulong R0 = nmod_poly_get_coeff_ui(src0, i);
        ulong R1 = nmod_poly_get_coeff_ui(src1, i);
        
        ulong diff = (R1 >= R0) ? (R1 - R0) : (R1 + Q - R0);
        ulong a_i_256 = (diff * ctx->inv_C0_sub_C1) % Q;
        
        ulong term1 = (R1 * ctx->C0_inv_diff) % Q;
        ulong term2 = (R0 * ctx->C1_inv_diff) % Q;
        ulong a_i = (term1 >= term2) ? (term1 - term2) : (term1 + Q - term2);
        
        nmod_poly_set_coeff_ui(dest, i, a_i);
        nmod_poly_set_coeff_ui(dest, i + 256, a_i_256);
    }
}

void evmlr_poly_from_crt(nmod_poly_t dest, const evmlr_poly_crt_t src, const evmlr_crt_ctx_t ctx) {
    _evmlr_poly_from_crt(dest, src->poly + 0, src->poly + 1, ctx);
}

void evmlr_poly_mat_to_crt(evmlr_poly_mat_crt_t dest, const nmod_poly_mat_t src, const evmlr_crt_ctx_t ctx) {
    for (slong i = 0; i < src->r; i++) {
        for (slong j = 0; j < src->c; j++) {
            nmod_poly_struct* src_ij = nmod_poly_mat_entry(src, i, j);
            nmod_poly_struct* dest_0_ij = nmod_poly_mat_entry(dest->mat + 0, i, j);
            nmod_poly_struct* dest_1_ij = nmod_poly_mat_entry(dest->mat + 1, i, j);

            _evmlr_poly_to_crt(dest_0_ij, dest_1_ij, src_ij, ctx);
        }
    }
}

void evmlr_poly_mat_from_crt(nmod_poly_mat_t dest, const evmlr_poly_mat_crt_t src, const evmlr_crt_ctx_t ctx) {
    for (slong i = 0; i < src->mat[0].r; i++) {
        for (slong j = 0; j < src->mat[0].c; j++) {
            nmod_poly_struct* src_0_ij = nmod_poly_mat_entry(src->mat + 0, i, j);
            nmod_poly_struct* src_1_ij = nmod_poly_mat_entry(src->mat + 1, i, j);
            nmod_poly_struct* dest_ij = nmod_poly_mat_entry(dest, i, j);

            _evmlr_poly_from_crt(dest_ij, src_0_ij, src_1_ij, ctx);
        }
    }
}

void evmlr_poly_crt_mul(evmlr_poly_crt_t res, const evmlr_poly_crt_t op1, const evmlr_poly_crt_t op2, const evmlr_crt_ctx_t ctx) {
    nmod_poly_mulmod(res->poly + 0, op1->poly + 0, op2->poly + 0, ctx->factors->p + 0);
    nmod_poly_mulmod(res->poly + 1, op1->poly + 1, op2->poly + 1, ctx->factors->p + 1);
}

void evmlr_poly_mat_crt_mul(evmlr_poly_mat_crt_t res, const evmlr_poly_mat_crt_t op1, const evmlr_poly_mat_crt_t op2, const evmlr_crt_ctx_t ctx) {
    nmod_poly_mat_mulmod(res->mat + 0, op1->mat + 0, op2->mat + 0, ctx->factors->p + 0);
    nmod_poly_mat_mulmod(res->mat + 1, op1->mat + 1, op2->mat + 1, ctx->factors->p + 1);
}
