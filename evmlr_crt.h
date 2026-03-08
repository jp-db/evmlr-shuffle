#ifndef EVMLR_SHUFFLE_EVMLR_CRT_H
#define EVMLR_SHUFFLE_EVMLR_CRT_H

#include "flint/nmod_poly.h"
#include "flint/nmod_poly_mat.h"
#include "flint/nmod_poly_factor.h"

// Struct for a CRT-transformed polynomial
typedef struct {
    nmod_poly_struct poly[2]; // Two degree-256 remainders
} evmlr_poly_crt_struct;

typedef evmlr_poly_crt_struct evmlr_poly_crt_t[1];

// Struct for a CRT-transformed matrix of polynomials
typedef struct {
    nmod_poly_mat_struct mat[2];
} evmlr_poly_mat_crt_struct;

typedef evmlr_poly_mat_crt_struct evmlr_poly_mat_crt_t[1];

// CRT Setup and Precomputation Context
typedef struct {
    nmod_poly_factor_t factors;
    ulong C_0;
    ulong C_1;
    ulong inv_C0_sub_C1; // (C0 - C1)^{-1} mod Q
    ulong C0_inv_diff;   // C0 * inv_C0_sub_C1 mod Q
    ulong C1_inv_diff;   // C1 * inv_C0_sub_C1 mod Q
    ulong mod;
} evmlr_crt_ctx_struct;

typedef evmlr_crt_ctx_struct evmlr_crt_ctx_t[1];

// Global Initialization
void evmlr_crt_setup(evmlr_crt_ctx_t ctx, const nmod_poly_t cyclo_poly);
void evmlr_crt_clear(evmlr_crt_ctx_t ctx);

// Memory management
void evmlr_poly_crt_init(evmlr_poly_crt_t crt, ulong n);
void evmlr_poly_crt_clear(evmlr_poly_crt_t crt);

void evmlr_poly_mat_crt_init(evmlr_poly_mat_crt_t crt, slong rows, slong cols, ulong n);
void evmlr_poly_mat_crt_clear(evmlr_poly_mat_crt_t crt);

// Transform operations
void evmlr_poly_to_crt(evmlr_poly_crt_t dest, const nmod_poly_t src, const evmlr_crt_ctx_t ctx);
void evmlr_poly_from_crt(nmod_poly_t dest, const evmlr_poly_crt_t src, const evmlr_crt_ctx_t ctx);

void evmlr_poly_mat_to_crt(evmlr_poly_mat_crt_t dest, const nmod_poly_mat_t src, const evmlr_crt_ctx_t ctx);
void evmlr_poly_mat_from_crt(nmod_poly_mat_t dest, const evmlr_poly_mat_crt_t src, const evmlr_crt_ctx_t ctx);

// Arithmetic in CRT domain
// point-wise multiplication
void evmlr_poly_crt_mul(evmlr_poly_crt_t res, const evmlr_poly_crt_t op1, const evmlr_poly_crt_t op2, const evmlr_crt_ctx_t ctx);

// matrix multiplication in CRT domain
void evmlr_poly_mat_crt_mul(evmlr_poly_mat_crt_t res, const evmlr_poly_mat_crt_t op1, const evmlr_poly_mat_crt_t op2, const evmlr_crt_ctx_t ctx);

#endif // EVMLR_SHUFFLE_EVMLR_CRT_H
