#include "evmlr_shuffle.h"
#include "evmlr_utils.h"

#ifdef MAIN
#include "test.h"
#include "bench.h"
#endif

void evmlr_shuffle_ctx_init(evmlr_shuffle_ctx_t ctx, slong N, slong L, flint_rand_t state) {
    slong com_N = N * (LOG_Q_CEIL > (2*ETA) * (K_LWE + L)
            ? LOG_Q_CEIL
            : (2*ETA) * (K_LWE + L)); // initialize with the biggest vector len used
    evmlr_commit_ctx_init(ctx->com_ctx, com_N, state);
    evmlr_hpke_ctx_init(ctx->hpke_ctx, L, state);
    ctx->N = N;
    ctx->L = L;
}

void evmlr_shuffle_ctx_clear(evmlr_shuffle_ctx_t ctx) {
    evmlr_commit_ctx_clear(ctx->com_ctx);
    evmlr_hpke_ctx_clear(ctx->hpke_ctx);
}

void evmlr_proof_init(evmlr_shuffle_proof_t proof, const evmlr_shuffle_ctx_t ctx) {
    nmod_poly_init(proof->alpha, MOD_Q);
    nmod_poly_init(proof->beta, MOD_Q);
    nmod_poly_init(proof->lambda, MOD_Q);
    nmod_poly_init(proof->gamma, MOD_Q);
    nmod_poly_mat_init(proof->u, ctx->N - 1, 1, MOD_Q);
}

void evmlr_proof_clear(evmlr_shuffle_proof_t proof, const evmlr_shuffle_ctx_t ctx) {
    nmod_poly_clear(proof->alpha);
    nmod_poly_clear(proof->beta);
    nmod_poly_clear(proof->lambda);
    nmod_poly_clear(proof->gamma);
    nmod_poly_mat_clear(proof->u);
    evmlr_commit_clear(proof->P);
    evmlr_commit_clear(proof->W);
}

void evmlr_shuffle_sp_clear(evmlr_shuffle_sp_t sp, const evmlr_shuffle_ctx_t ctx) {
    free(sp->pi);
    nmod_poly_mat_clear(sp->r_D);
    for (size_t i = 0; i < ctx->N; i++) {
        nmod_poly_mat_clear(sp->d_dagger[i]);
    }
    free(sp->d_dagger);
}

void evmlr_shuffle_pp_clear(evmlr_shuffle_pp_t pp, const evmlr_shuffle_ctx_t ctx) {
    evmlr_commit_clear(pp->D);
    for (size_t i = 0; i < ctx->N; i++) {
        nmod_poly_mat_clear(pp->c_hat[i]);
        evmlr_otse_ciphertext_clear(pp->c_star[i]);
    }
    free(pp->c_hat);
    free(pp->c_star);
}

void evmlr_shuffle_sample_challenge(nmod_poly_t chal, flint_rand_t state) {
    nmod_poly_init(chal, MOD_Q);
    nmod_poly_randtest(chal, state, DEGREE_N >> 1);
}

void evmlr_shuffle_phase_1(nmod_poly_mat_t sigma, evmlr_commit_t com, nmod_poly_mat_t r, const evmlr_shuffle_ctx_t ctx, const evmlr_shuffle_sp_t sp) {
    // Sample r from B^{2K_SIS}
    evmlr_commit_sample_r(r);
    nmod_poly_mat_init(sigma, ctx->com_ctx->N, 1, MOD_Q);
    nmod_poly_mat_zero(sigma);
    for (slong i = 0; i < ctx->N; i++) {
        evmlr_utils_int_to_bin(nmod_poly_mat_entry(sigma, i, 0), sp->pi[i]);
    }
    evmlr_commit(com, sigma, r, ctx->com_ctx);
}

void lambda_setup(nmod_poly_mat_t Lambda, const nmod_poly_t lambda, const nmod_poly_t cyclo_poly, slong L);

void calc_z(nmod_poly_t z, const evmlr_otse_ciphertext_t cipher, const nmod_poly_mat_t d_dagger,
            const nmod_poly_mat_t Lambda, const evmlr_shuffle_ctx_t ctx);

void calc_z_hat(nmod_poly_t z_hat, const nmod_poly_mat_t c_hat, const nmod_poly_mat_t Lambda, const evmlr_shuffle_ctx_t ctx);

void calc_h(nmod_poly_t h, const nmod_poly_t z, const nmod_poly_t bin_poly, const nmod_poly_t alpha, const nmod_poly_t beta, const nmod_poly_t cyclo_poly);

void evmlr_shuffle_phase_2(nmod_poly_mat_t h, nmod_poly_mat_t h_hat, evmlr_shuffle_proof_t proof, nmod_poly_mat_t r_W,
                           nmod_poly_mat_t Theta, nmod_poly_mat_t w_dagger[],
                           const nmod_poly_mat_t sigma, const evmlr_shuffle_ctx_t ctx, const evmlr_shuffle_sp_t sp,
                           const evmlr_shuffle_pp_t pp, flint_rand_t state) {
    slong L = ctx->L;
    slong N = ctx->N;
    nmod_poly_struct *alpha = proof->alpha, *beta = proof->beta, *lambda = proof->lambda;

    nmod_poly_mat_init(h, N, 1, MOD_Q);
    nmod_poly_mat_init(h_hat, N, 1, MOD_Q);

    nmod_poly_mat_t Lambda;
    lambda_setup(Lambda, lambda, ctx->com_ctx->cyclo_poly, L);

    nmod_poly_mat_t z, z_hat;
    nmod_poly_mat_init(z, N, 1, MOD_Q);
    nmod_poly_mat_init(z_hat, N, 1, MOD_Q);
    nmod_poly_mat_zero(z);
    nmod_poly_mat_zero(z_hat);

    nmod_poly_t tmp;
    nmod_poly_init(tmp, MOD_Q);
    for (slong i = 0; i < N; i++) {
        calc_z(nmod_poly_mat_entry(z, i, 0), pp->c_star[i], sp->d_dagger[i], Lambda, ctx); // z_i = ⟨c^⋆_i − H'B(η)d^dagger_i , Λ⟩
        calc_z_hat(nmod_poly_mat_entry(z_hat, i, 0), pp->c_hat[i], Lambda, ctx); // ẑ_i = ⟨c^hat_i, Λ⟩
        evmlr_utils_int_to_bin(tmp, i);
        calc_h(nmod_poly_mat_entry(h, i, 0), nmod_poly_mat_entry(z, i, 0), tmp, alpha, beta, ctx->com_ctx->cyclo_poly); // h_i = z_i + int_to_bin(i) * beta - alpha
        calc_h(nmod_poly_mat_entry(h_hat, i, 0), nmod_poly_mat_entry(z_hat, i, 0), nmod_poly_mat_entry(sigma, i, 0), alpha, beta, ctx->com_ctx->cyclo_poly); // ĥ_i = ẑ_i + sigma_i * beta - alpha
    }

    nmod_poly_mat_init(Theta, N - 1, 1, MOD_Q);
    nmod_poly_mat_randtest(Theta, state, DEGREE_N);

    nmod_poly_mat_t w;
    nmod_poly_mat_init(w, N, 1, MOD_Q);

    nmod_poly_mulmod(nmod_poly_vec_entry(w, 0), nmod_poly_vec_entry(Theta, 0), nmod_poly_vec_entry(h, 0), ctx->com_ctx->cyclo_poly);
    for (slong i = 1; i < N - 1; i++) {
        // w_i = Θ_{i-1} * h^hat_i + Θ_i * h_i
        nmod_poly_mulmod(nmod_poly_vec_entry(w, i), nmod_poly_vec_entry(Theta, i-1), nmod_poly_vec_entry(h_hat, i), ctx->com_ctx->cyclo_poly);
        nmod_poly_mulmod(tmp, nmod_poly_vec_entry(Theta, i), nmod_poly_vec_entry(h, i), ctx->com_ctx->cyclo_poly);
        nmod_poly_add(nmod_poly_vec_entry(w, i), nmod_poly_vec_entry(w, i), tmp);
    }
    nmod_poly_mulmod(nmod_poly_vec_entry(w, N-1), nmod_poly_vec_entry(Theta, N-2),
                     nmod_poly_vec_entry(h_hat, N-1), ctx->com_ctx->cyclo_poly);

    // TODO see coefficient range
    // w_i = G^{LOG_Q_CEIL} w_dagger_i
    nmod_poly_mat_t tmp_mat;
    nmod_poly_mat_init(tmp_mat, 1, 1, MOD_Q);
    for (int i = 0; i < N; i++) {
        nmod_poly_set(nmod_poly_vec_entry(tmp_mat, 0), nmod_poly_vec_entry(w, i));
        evmlr_utils_ring_to_bin(w_dagger[i], tmp_mat, LOG_Q_CEIL);
    }
    nmod_poly_mat_clear(tmp_mat);

    nmod_poly_mat_t w_flat;
    nmod_poly_mat_init(w_flat, ctx->com_ctx->N, 1, MOD_Q);
    nmod_poly_mat_zero(w_flat);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < LOG_Q_CEIL; j++) {
            nmod_poly_struct* w_poly = nmod_poly_vec_entry(w_dagger[i], j);
            nmod_poly_set(nmod_poly_vec_entry(w_flat, i * LOG_Q_CEIL + j), w_poly);
        }
    }

    evmlr_commit_sample_r(r_W);
    evmlr_commit(proof->W, w_flat, r_W, ctx->com_ctx);

    // cleanup
    nmod_poly_clear(tmp);
    nmod_poly_mat_clear(Lambda);
    nmod_poly_mat_clear(z);
    nmod_poly_mat_clear(z_hat);
    nmod_poly_mat_clear(w);
    nmod_poly_mat_clear(w_flat);
}

void evmlr_shuffle_phase_3(evmlr_shuffle_proof_t proof, const nmod_poly_mat_t h, const nmod_poly_mat_t h_hat,
                           const nmod_poly_mat_t Theta, const evmlr_shuffle_ctx_t ctx) {
    slong N = ctx->N;

    // Temporary polynomials for calculations
    nmod_poly_t running_prod, tmp;
    nmod_poly_init(running_prod, MOD_Q);
    nmod_poly_init(tmp, MOD_Q);

    // Initialize the running product with gamma
    nmod_poly_set(running_prod, proof->gamma);

    // Loop to calculate u_1 through u_{N-1} (indexed 0 to N-2 in code)
    for (slong i = 0; i < N - 1; i++) {
        // h_hat[i] / h[i]
        nmod_poly_invmod(tmp, nmod_poly_vec_entry(h, i), ctx->com_ctx->cyclo_poly);
        nmod_poly_mulmod(tmp, nmod_poly_vec_entry(h_hat, i), tmp, ctx->com_ctx->cyclo_poly);

        // update gamma * \prod_{j=1}^{i} (h_hat[j] / h[j])
        nmod_poly_mulmod(running_prod, running_prod, tmp, ctx->com_ctx->cyclo_poly);

        // 3. (-1)^(i+1). If i is even, negate the running product
        if (i % 2 == 0) nmod_poly_neg(tmp, running_prod);
        else            nmod_poly_set(tmp, running_prod);

        // +- gamma * prod + Θ_i
        nmod_poly_add(nmod_poly_vec_entry(proof->u, i), tmp, nmod_poly_vec_entry(Theta, i));
    }

    // Cleanup
    nmod_poly_clear(running_prod);
    nmod_poly_clear(tmp);
}

// check if everything is correct. TEST code, not zero-knowledge
int evmlr_shuffle_proof_ver(evmlr_shuffle_proof_t proof, const nmod_poly_mat_t h, const nmod_poly_mat_t h_hat,
                            const nmod_poly_mat_t r_W, nmod_poly_mat_t w_dagger[], const nmod_poly_mat_t r_P,
                            const nmod_poly_mat_t sigma, const evmlr_shuffle_ctx_t ctx,
                            const evmlr_shuffle_sp_t sp, const evmlr_shuffle_pp_t pp) {
    nmod_poly_struct *gamma = proof->gamma;
    nmod_poly_mat_struct *u = proof->u;
    evmlr_commit_struct *P = proof->P, *W = proof->W;

    nmod_poly_mat_t tmp1, tmp2;
    nmod_poly_mat_init(tmp1, 1, 1, MOD_Q);
    nmod_poly_mat_init(tmp2, 1, 1, MOD_Q);

    nmod_poly_mat_t G;
    evmlr_utils_gadget_matrix(G, 1, LOG_Q_CEIL, MOD_Q);

    nmod_poly_mat_t mat_flattened;
    nmod_poly_mat_init(mat_flattened, ctx->com_ctx->N, 1, MOD_Q);

    int valid = 1;
    // \gamma h_hat_1 + u_1 h_1 = G^{LOG_Q_CEIL} w_dagger_1
    // gamma * h_hat_1
    nmod_poly_mulmod(nmod_poly_vec_entry(tmp1, 0), gamma, nmod_poly_vec_entry(h_hat, 0), ctx->com_ctx->cyclo_poly);
    // u_1 * h_1
    nmod_poly_mulmod(nmod_poly_vec_entry(tmp2, 0), nmod_poly_vec_entry(u, 0), nmod_poly_vec_entry(h, 0), ctx->com_ctx->cyclo_poly);
    // gamma * h_hat_1 + u_1 * h_1
    nmod_poly_mat_add(tmp1, tmp1, tmp2);
    // G^{LOG_Q_CEIL} w_dagger_1
    nmod_poly_mat_mulmod(tmp2, G, w_dagger[0], ctx->com_ctx->cyclo_poly);
    valid &= nmod_poly_mat_equal(tmp1, tmp2);
    if (!valid) goto end;

    // u_{i−1} h_hat_i + u_i * h_i = G^{LOG_Q_CEIL} w_dagger_i for i = 2,...,N-1
    for (slong i = 1; i < ctx->N - 1; i++) {
        // u_{i-1} * h_hat_i
        nmod_poly_mulmod(nmod_poly_vec_entry(tmp1, 0), nmod_poly_vec_entry(u,i - 1), nmod_poly_vec_entry(h_hat, i), ctx->com_ctx->cyclo_poly);
        // u_i * h_i
        nmod_poly_mulmod(nmod_poly_vec_entry(tmp2, 0), nmod_poly_vec_entry(u, i), nmod_poly_vec_entry(h, i), ctx->com_ctx->cyclo_poly);
        // u_{i- 1} * h_hat_i + u_i * h_i
        nmod_poly_mat_add(tmp1, tmp1, tmp2);
        // G^{LOG_Q_CEIL} w_dagger_i
        nmod_poly_mat_mulmod(tmp2, G, w_dagger[i], ctx->com_ctx->cyclo_poly);
        valid &= nmod_poly_mat_equal(tmp1, tmp2);
        if (!valid) goto end;
    }

    // u_{N−1} h_hat_N + (-1)^N \gamma h_n = G^{LOG_Q_CEIL} w_dagger_N
    // u_{N-1} * h_hat_N
    nmod_poly_mulmod(nmod_poly_vec_entry(tmp1, 0), nmod_poly_vec_entry(u, ctx->N - 2), nmod_poly_vec_entry(h_hat, ctx->N - 1), ctx->com_ctx->cyclo_poly);
    // gamma * h_N
    nmod_poly_mulmod(nmod_poly_vec_entry(tmp2, 0), gamma, nmod_poly_vec_entry(h, ctx->N - 1), ctx->com_ctx->cyclo_poly);
    // (-1)^N * gamma * h_N
    if (ctx->N % 2 != 0) nmod_poly_neg(nmod_poly_vec_entry(tmp2, 0), nmod_poly_vec_entry(tmp2, 0));
    // u_{N-1} * h_hat_N + (-1)^N * gamma * h_N
    nmod_poly_mat_add(tmp1, tmp1, tmp2);

    // G^{LOG_Q_CEIL} w_dagger_N
    nmod_poly_mat_mulmod(tmp2, G, w_dagger[ctx->N - 1], ctx->com_ctx->cyclo_poly);
    valid &= nmod_poly_mat_equal(tmp1, tmp2);
    if (!valid) goto end;

    valid &= evmlr_utils_is_mat_bin(sigma);
    for (ulong i = 0; i < ctx->N && valid; i++) {
        valid &= evmlr_utils_is_mat_bin(sp->d_dagger[i]) && evmlr_utils_is_mat_bin(w_dagger[i]);
    }
    if (!valid) goto end;

    nmod_poly_mat_zero(mat_flattened);
    for (slong i = 0; i < ctx->N; i++) {
        for (slong j = 0; j < LOG_Q_CEIL; j++) {
            nmod_poly_struct* w_poly = nmod_poly_vec_entry(w_dagger[i], j);
            nmod_poly_set(nmod_poly_vec_entry(mat_flattened, i * LOG_Q_CEIL + j), w_poly);
        }
    }

    valid &= evmlr_commit_verify(W, mat_flattened, r_W, ctx->com_ctx);
    if (!valid) goto end;

    nmod_poly_mat_zero(mat_flattened);
    for (slong i = 0; i < ctx->N; i++) {
        nmod_poly_struct* sigma_i = nmod_poly_mat_entry(sigma, i, 0);
        nmod_poly_set(nmod_poly_vec_entry(mat_flattened, i), sigma_i);
    }

    valid &= evmlr_commit_verify(P, mat_flattened, r_P, ctx->com_ctx);
    if (!valid) goto end;

    nmod_poly_mat_zero(mat_flattened);
    for (slong i = 0; i < ctx->N; i++) {
        for (slong j = 0; j < (2*ETA) * (K_LWE + ctx->L); j++) {
            nmod_poly_struct* d_poly = nmod_poly_mat_entry(sp->d_dagger[i], j, 0);
            nmod_poly_set(nmod_poly_vec_entry(mat_flattened, i * (2*ETA) * (K_LWE + ctx->L) + j), d_poly);
        }
    }
    valid = evmlr_commit_verify(pp->D, mat_flattened, sp->r_D, ctx->com_ctx);

    end:
        nmod_poly_mat_clear(tmp1);
        nmod_poly_mat_clear(tmp2);
        nmod_poly_mat_clear(G);
        nmod_poly_mat_clear(mat_flattened);
        return valid;
}

void evmlr_shuffle_clear_sp(evmlr_shuffle_sp_t sp, const evmlr_shuffle_ctx_t ctx) {
    free(sp->pi);
    for (size_t i = 0; i < ctx->N; i++) {
        nmod_poly_mat_clear(sp->d_dagger[i]);
    }
    free(sp->d_dagger);
    nmod_poly_mat_clear(sp->r_D);
}

void evmlr_shuffle_clear_pp(evmlr_shuffle_pp_t pp, const evmlr_shuffle_ctx_t ctx) {
    evmlr_commit_clear(pp->D);
    for (size_t i = 0; i < ctx->N; i++) {
        nmod_poly_mat_clear(pp->c_hat[i]);
        evmlr_otse_ciphertext_clear(pp->c_star[i]);
    }
    free(pp->c_hat);
    free(pp->c_star);
}

void lambda_setup(nmod_poly_mat_t Lambda, const nmod_poly_t lambda, const nmod_poly_t cyclo_poly, slong L) {
    nmod_poly_mat_init(Lambda, L, 1, MOD_Q);
    for (int i = 0; i < L; i++) {
        nmod_poly_struct* Lambda_i = nmod_poly_mat_entry(Lambda, i, 0);
        switch (i) {
            case 0: nmod_poly_one(Lambda_i);           break;
            case 1: nmod_poly_set(Lambda_i, lambda); break;
            default:
                nmod_poly_mulmod(Lambda_i, nmod_poly_mat_entry(Lambda, i-1, 0), lambda, cyclo_poly);
                break;
        }
    }
}

void calc_z(nmod_poly_t z, const evmlr_otse_ciphertext_t cipher, const nmod_poly_mat_t d_dagger,
            const nmod_poly_mat_t Lambda, const evmlr_shuffle_ctx_t ctx) {
    slong L = ctx->L;
    nmod_poly_mat_t a;
    nmod_poly_mat_init(a, L, 1, MOD_Q);
    evmlr_calc_a(a, d_dagger, ctx->hpke_ctx->otse_ctx);

    for (int j = 0; j < L; j++) {
        nmod_poly_struct* a_j = nmod_poly_mat_entry(a, j, 0);
        nmod_poly_sub(a_j, nmod_poly_mat_entry(cipher->c, j, 0), a_j);
        nmod_poly_mulmod(a_j, a_j, nmod_poly_mat_entry(Lambda, j, 0), ctx->com_ctx->cyclo_poly);
        nmod_poly_add(z, z, a_j);
    }

    nmod_poly_mat_clear(a);
}

void calc_z_hat(nmod_poly_t z_hat, const nmod_poly_mat_t c_hat, const nmod_poly_mat_t Lambda, const evmlr_shuffle_ctx_t ctx) {
    slong L = ctx->L;
    nmod_poly_t tmp;
    nmod_poly_init(tmp, MOD_Q);

    for (int j = 0; j < L; j++) {
        nmod_poly_struct* c_hat_j = nmod_poly_mat_entry(c_hat, j, 0);
        nmod_poly_struct* Lambda_j = nmod_poly_mat_entry(Lambda, j, 0);
        nmod_poly_mulmod(tmp, c_hat_j, Lambda_j, ctx->com_ctx->cyclo_poly);
        nmod_poly_add(z_hat, z_hat, tmp);
    }
    nmod_poly_clear(tmp);
}

void calc_h(nmod_poly_t h, const nmod_poly_t z, const nmod_poly_t bin_poly, const nmod_poly_t alpha,
            const nmod_poly_t beta, const nmod_poly_t cyclo_poly) {
    nmod_poly_set(h, bin_poly);
    nmod_poly_mulmod(h, h, beta, cyclo_poly);
    nmod_poly_add(h, z, h);
    nmod_poly_sub(h, h, alpha);
}

int evmlr_shuffle_run(evmlr_shuffle_sp_t sp, evmlr_shuffle_pp_t pp, const evmlr_shuffle_ctx_t ctx, flint_rand_t state) {
    evmlr_shuffle_proof_t proof;
    evmlr_proof_init(proof, ctx);

    size_t N = ctx->N;
    nmod_poly_mat_t sigma;
    nmod_poly_mat_t r_P;
    // Prover computes P, sigma, r_P
    evmlr_shuffle_phase_1(sigma, proof->P, r_P, ctx, sp);
    // Verifier samples challenges
    evmlr_shuffle_sample_challenge(proof->alpha, state);
    evmlr_shuffle_sample_challenge(proof->beta, state);
    evmlr_shuffle_sample_challenge(proof->lambda, state);

    nmod_poly_mat_t h, h_hat;
    nmod_poly_mat_t r_W;
    nmod_poly_mat_t Theta;
    nmod_poly_mat_t w_dagger[N];

    // Prover computes W, h, h_hat, Theta, w_dagger, r_W
    evmlr_shuffle_phase_2(h, h_hat, proof, r_W, Theta, w_dagger, sigma, ctx, sp, pp, state);
    // Verifier samples gamma
    evmlr_shuffle_sample_challenge(proof->gamma, state);

    // Prover computes u
    evmlr_shuffle_phase_3(proof, h, h_hat, Theta, ctx);

    int valid = evmlr_shuffle_proof_ver(proof, h, h_hat, r_W, w_dagger, r_P, sigma, ctx, sp, pp);

    // cleanup
    evmlr_proof_clear(proof, ctx);
    nmod_poly_mat_clear(sigma);
    nmod_poly_mat_clear(r_P);

    nmod_poly_mat_clear(h);
    nmod_poly_mat_clear(h_hat);
    nmod_poly_mat_clear(r_W);
    nmod_poly_mat_clear(Theta);
    for (size_t i = 0; i < N; i++) {
        nmod_poly_mat_clear(w_dagger[i]);
    }
    return valid;

}

#ifdef MAIN
void setup(evmlr_shuffle_sp_t sp, evmlr_shuffle_pp_t pp, const evmlr_shuffle_ctx_t ctx, flint_rand_t state) {
    slong N = ctx->N;
    slong L = ctx->L;
    evmlr_hpke_keypair_t keypair;
    evmlr_hpke_keypair_gen(keypair, state, ctx->hpke_ctx);

    nmod_poly_mat_t m[N];
    for (size_t i = 0; i < N; i++) {
        nmod_poly_mat_init(m[i], L, 1, MOD_Q);
        nmod_poly_mat_randtest(m[i], state, DEGREE_N);
    }

    sp->pi = (size_t *) malloc(N * sizeof(size_t));
    evmlr_utils_new_perm(sp->pi, N);

    // shuffle messages
    pp->c_hat = (nmod_poly_mat_t *) malloc(N * sizeof(nmod_poly_mat_t));
    for (size_t i = 0; i < N; i++) {
        nmod_poly_mat_init(pp->c_hat[i], L, 1, MOD_Q);
        nmod_poly_mat_set(pp->c_hat[i], m[sp->pi[i]]);
    }

    pp->c_star = (evmlr_otse_ciphertext_t *) malloc(N * sizeof(evmlr_otse_ciphertext_t));
    sp->d_dagger = (nmod_poly_mat_t *) malloc(N * sizeof(nmod_poly_mat_t));
    for (slong i = 0; i < N; i++) {
        evmlr_hpke_cipher_t cipher;
        nmod_poly_mat_struct* d_dagger_i = sp->d_dagger[i];
        evmlr_hpke_encrypt(cipher, d_dagger_i, m[i], keypair->enc_keypair->pk, ctx->hpke_ctx, state);
        nmod_poly_mat_init(pp->c_star[i]->c, L, 1, MOD_Q);
        nmod_poly_mat_set(pp->c_star[i]->c, cipher->otse_cipher->c);
        evmlr_hpke_cipher_clear(cipher);
    }

    evmlr_commit_sample_r(sp->r_D);

    nmod_poly_mat_t d_flat;
    nmod_poly_mat_init(d_flat, ctx->com_ctx->N, 1, MOD_Q);
    nmod_poly_mat_zero(d_flat);
    slong row_len = (2*ETA) * (K_LWE + L);
    for (slong i = 0; i < N; i++) {
        for (slong j = 0; j < row_len; j++) {
            nmod_poly_struct* d_dag_ij = nmod_poly_mat_entry(sp->d_dagger[i], j, 0);
            nmod_poly_struct* poly = nmod_poly_mat_entry(d_flat, i * row_len + j, 0);
            nmod_poly_set(poly, d_dag_ij);
        }
    }
    evmlr_commit(pp->D, d_flat, sp->r_D, ctx->com_ctx);

    nmod_poly_mat_clear(d_flat);
    evmlr_hpke_keypair_clear(keypair);
    for (int i = 0; i < N; ++i) {
        nmod_poly_mat_clear(m[i]);
    }
}

void test(evmlr_shuffle_ctx_t ctx, flint_rand_t state) {
    evmlr_shuffle_sp_t sp;
    evmlr_shuffle_pp_t pp;
    setup(sp, pp, ctx, state);

    TEST_BEGIN("shuffle proof works") {
        TEST_ASSERT(evmlr_shuffle_run(sp, pp, ctx, state) == 1, end)
    } TEST_END;

    end:
        evmlr_shuffle_clear_sp(sp, ctx);
        evmlr_shuffle_clear_pp(pp, ctx);
}

void bench(evmlr_shuffle_ctx_t ctx, flint_rand_t state) {
    evmlr_shuffle_sp_t sp;
    evmlr_shuffle_pp_t pp;
    setup(sp, pp, ctx, state);

    BENCH_BEGIN("evmlr_shuffle_run") {
        BENCH_ADD(evmlr_shuffle_run(sp, pp, ctx, state))
    } BENCH_END;

    evmlr_shuffle_clear_sp(sp, ctx);
    evmlr_shuffle_clear_pp(pp, ctx);
}

int main() {
    slong n_messages = SHUFFLE_N_MSGS;

    flint_rand_t state;
    flint_rand_init(state);
    ulong seed[2];
    getrandom(seed, sizeof(ulong)*2, 0);
    flint_rand_set_seed(state, seed[0], seed[1]);

    evmlr_shuffle_ctx_t ctx;
    evmlr_shuffle_ctx_init(ctx, n_messages, M_LEN, state);

    test(ctx, state);
    bench(ctx, state);

    evmlr_shuffle_ctx_clear(ctx);
    flint_rand_clear(state);
    return 0;
}
#endif