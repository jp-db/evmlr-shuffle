#include "evmlr_shuffle.h"
#include "evmlr_utils.h"

#include "test.h"
#include "bench.h"

int checkpoint=0;

void evmlr_shuffle_ctx_init(evmlr_shuffle_ctx_t ctx, size_t N, size_t L, flint_rand_t state) {
    size_t com_N = N * (LOG_Q_CEIL > (2*ETA) * (K_LWE + L)
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

void evmlr_shuffle_sp_init(evmlr_shuffle_sp_t sp, size_t pi[], nmod_poly_t** d_dagger,
                           nmod_poly_t r_D[2 * K_SIS]) {
    sp->pi = pi;
    sp->d_dagger = d_dagger;
    for (ulong i = 0; i < 2 * K_SIS; i++) {
        nmod_poly_init(sp->r_D[i], MOD_Q);
        nmod_poly_set(sp->r_D[i], r_D[i]);
    }
}

void evmlr_shuffle_sp_clear(evmlr_shuffle_sp_t sp) {
    for (ulong i = 0; i < 2 * K_SIS; i++) {
        nmod_poly_clear(sp->r_D[i]);
    }
}

void evmlr_shuffle_pp_init(evmlr_shuffle_pp_t pp, evmlr_otse_ciphertext_t c_star[], nmod_poly_t** c_hat,
                           evmlr_commit_t D) {
    // set c_star and c_hat as the memory address passed
    pp->c_star = c_star;
    pp->c_hat = c_hat;
    for (size_t i = 0; i < K_SIS; i++) {
        nmod_poly_init(pp->D->c[i], MOD_Q);
        nmod_poly_set(pp->D->c[i], D->c[i]);
    }
}

void evmlr_shuffle_pp_clear(evmlr_shuffle_pp_t pp) {
    evmlr_commit_clear(pp->D);
}

void evmlr_shuffle_sample_challenge(nmod_poly_t chal, flint_rand_t state) {
    nmod_poly_init(chal, MOD_Q);
    nmod_poly_randtest(chal, state, DEGREE_N >> 1);
}

void evmlr_shuffle_phase_1(nmod_poly_t sigma[], evmlr_commit_t com, nmod_poly_t r[2*K_SIS], const evmlr_shuffle_ctx_t ctx, const evmlr_shuffle_sp_t sp) {
    // Sample r from B^{2K_SIS}
    evmlr_commit_sample_r(r);
    for (size_t i = 0; i < ctx->N; i++) {
        nmod_poly_init(sigma[i], MOD_Q);
        evmlr_utils_int_to_bin(sigma[i], sp->pi[i]);
    }
    for (size_t i = ctx->N; i < ctx->com_ctx->N; i++) {
        nmod_poly_init(sigma[i], MOD_Q);
        nmod_poly_zero(sigma[i]);
    }
    evmlr_commit(com, ctx->com_ctx, sigma, r);
}

void lambda_setup(nmod_poly_t Lambda[], const nmod_poly_t lambda, const nmod_poly_t cyclo_poly, ulong L);

void calc_z(nmod_poly_t z, const evmlr_otse_ciphertext_t cipher, const nmod_poly_t* d_dagger,
            const nmod_poly_t Lambda[], const evmlr_shuffle_ctx_t ctx);

void calc_z_hat(nmod_poly_t z_hat, const nmod_poly_t* c_hat, const nmod_poly_t Lambda[], const evmlr_shuffle_ctx_t ctx);

void calc_h(nmod_poly_t h, const nmod_poly_t z, const nmod_poly_t bin_poly, const nmod_poly_t alpha, const nmod_poly_t beta, const nmod_poly_t cyclo_poly);

void evmlr_shuffle_phase_2(nmod_poly_t h[], nmod_poly_t h_hat[], evmlr_commit_t W, nmod_poly_t r_W[2*K_SIS],
                           nmod_poly_t Theta[], nmod_poly_t w_dagger[][LOG_Q_CEIL][1],
                           const nmod_poly_t sigma[], const evmlr_shuffle_ctx_t ctx, const evmlr_shuffle_sp_t sp,
                           const evmlr_shuffle_pp_t pp, const nmod_poly_t alpha, const nmod_poly_t beta,
                           const nmod_poly_t lambda, flint_rand_t state) {
    ulong L = ctx->L;
    ulong N = ctx->N;
    nmod_poly_t Lambda[L];
    lambda_setup(Lambda, lambda, ctx->com_ctx->cyclo_poly, L);
    printf("Checkpoint %d: Lambda setup done.\n", ++checkpoint);

    nmod_poly_t z[N], z_hat[N], tmp;
    nmod_poly_init(tmp, MOD_Q);
    for (int i = 0; i < N; i++) {
        calc_z(z[i], pp->c_star[i], sp->d_dagger[i], Lambda, ctx); // z_i = ⟨c^⋆_i − H'B(η)d^dagger_i , Λ⟩
        calc_z_hat(z_hat[i], pp->c_hat[i], Lambda, ctx); // ẑ_i = ⟨c^hat_i, Λ⟩
        evmlr_utils_int_to_bin(tmp, i);
        calc_h(h[i], z[i], tmp, alpha, beta, ctx->com_ctx->cyclo_poly); // h_i = z_i + int_to_bin(i) * beta - alpha
        calc_h(h_hat[i], z_hat[i], sigma[i], beta, alpha, ctx->com_ctx->cyclo_poly); // ĥ_i = ẑ_i + sigma_i * beta - alpha
    }
    printf("Checkpoint %d: h and h_hat calculation done.\n", ++checkpoint);

    for (int i = 0; i < N - 1; i++) {
        nmod_poly_init(Theta[i], MOD_Q);
        nmod_poly_randtest(Theta[i], state, DEGREE_N);
    }

    nmod_poly_t w[N];
    for (int i = 0; i < N; i++) {
        nmod_poly_init(w[i], MOD_Q);
    }
    nmod_poly_mulmod(w[0], Theta[0], h[0], ctx->com_ctx->cyclo_poly);

    for (ulong i = 1; i < N - 1; i++) {
        // w_i = Θ_{i-1} * h^hat_i + Θ_i * ĥ_i
        nmod_poly_mulmod(w[i], Theta[i - 1], h_hat[i], ctx->com_ctx->cyclo_poly);
        nmod_poly_mulmod(tmp, Theta[i], h[i], ctx->com_ctx->cyclo_poly);
        nmod_poly_add(w[i], w[i], tmp);
    }

    nmod_poly_mulmod(w[N-1], Theta[N-2], h_hat[N-1], ctx->com_ctx->cyclo_poly);

    // TODO see coefficient range
    fmpz_poly_t w_fmpz[N][1];
    for (int i = 0; i < N; ++i) {
        fmpz_poly_init(w_fmpz[i][0]);
        fmpz_poly_set_nmod_poly(w_fmpz[i][0], w[i]);
    }

    // w_i = G^{LOG_Q_CEIL} w_dagger_i
    for (int i = 0; i < N; i++) {
        evmlr_utils_ring_to_bin(1, LOG_Q_CEIL, w_dagger[i], w_fmpz[i], MOD_Q);
    }

    evmlr_commit_sample_r(r_W);
    evmlr_commit(W, ctx->com_ctx, (nmod_poly_t *) w_dagger, r_W);

    // cleanup
    nmod_poly_clear(tmp);
    for (int i = 0; i < L; i++) {
        nmod_poly_clear(Lambda[i]);
    }
    for (int i = 0; i < N; i++) {
        nmod_poly_clear(z[i]);
        nmod_poly_clear(z_hat[i]);
        nmod_poly_clear(w[i]);
        fmpz_poly_clear(w_fmpz[i][0]);
        for (int j = 0; j < LOG_Q_CEIL; j++) {
            nmod_poly_clear(w_dagger[i][j][0]);
        }
    }
}

void evmlr_shuffle_phase_3(nmod_poly_t u[], const nmod_poly_t h[], const nmod_poly_t h_hat[], const nmod_poly_t Theta[],
                           const nmod_poly_t gamma, const evmlr_shuffle_ctx_t ctx) {
    slong N = ctx->N;

    // Temporary polynomials for calculations
    nmod_poly_t running_prod, tmp;
    nmod_poly_init(running_prod, MOD_Q);
    nmod_poly_init(tmp, MOD_Q);

    // Initialize the running product with gamma
    nmod_poly_set(running_prod, gamma);

    // Loop to calculate u_1 through u_{N-1} (indexed 0 to N-2 in code)
    for (slong i = 0; i < N - 1; i++) {
        nmod_poly_init(u[i], MOD_Q);

        // h_hat[i] / h[i]
        nmod_poly_invmod(tmp, h[i], ctx->com_ctx->cyclo_poly);
        nmod_poly_mulmod(tmp, h_hat[i], tmp, ctx->com_ctx->cyclo_poly);

        // update gamma * \prod_{j=1}^{i} (h_hat[j] / h[j])
        nmod_poly_mulmod(running_prod, running_prod, tmp, ctx->com_ctx->cyclo_poly);

        // 3. (-1)^(i+1). If i is even, negate the running product
        if (i % 2 == 0) nmod_poly_neg(tmp, running_prod);
        else            nmod_poly_set(tmp, running_prod);

        // +- gamma * prod + Θ_i
        nmod_poly_add(u[i], tmp, Theta[i]);
    }

    // Cleanup
    nmod_poly_clear(running_prod);
    nmod_poly_clear(tmp);
}

// check if everything is correct. TEST code, not zero-knowledge
int evmlr_shuffle_proof_ver(const nmod_poly_t u[], const nmod_poly_t h[], const nmod_poly_t h_hat[],
                            const evmlr_commit_t W, const nmod_poly_t r_W[2*K_SIS], nmod_poly_t*** w_dagger,
                            const evmlr_commit_t P, const nmod_poly_t r_P[2*K_SIS], const nmod_poly_t sigma[],
                            const nmod_poly_t gamma, const evmlr_shuffle_ctx_t ctx, const evmlr_shuffle_sp_t sp,
                            const evmlr_shuffle_pp_t pp) {
    nmod_poly_t tmp1, tmp2, tmp3;
    nmod_poly_init(tmp1, MOD_Q);
    nmod_poly_init(tmp2, MOD_Q);
    nmod_poly_init(tmp3, MOD_Q);

    nmod_poly_t G[LOG_Q_CEIL][1];
    for (int i = 0; i < LOG_Q_CEIL; i++) {
        nmod_poly_init(G[i][0], MOD_Q);
        int coeff = 1 << i; // 2^i
        if (i == LOG_Q_CEIL - 1) coeff = MOD_Q - coeff; // -2^(b-1) mod q
        nmod_poly_set_coeff_ui(G[i][0], 0, coeff);
    }

    int valid = 1;
    // \gamma h_hat_1 + u_1 h_1 = G^{LOG_Q_CEIL} w_dagger_1
    nmod_poly_mulmod(tmp1, gamma, h_hat[0], ctx->com_ctx->cyclo_poly); // gamma * h_hat_1
    nmod_poly_mulmod(tmp2, u[0], h[0], ctx->com_ctx->cyclo_poly); // u_1 * h_1
    nmod_poly_add(tmp1, tmp1, tmp2); // gamma * h_hat_1 + u_1 * h_1

    // G^{LOG_Q_CEIL} w_dagger_1
    nmod_poly_zero(tmp2);
    for (int i = 0; i < LOG_Q_CEIL; i++) {
        nmod_poly_mulmod(tmp3, G[i][0], w_dagger[0][i][0], ctx->com_ctx->cyclo_poly);
        nmod_poly_add(tmp2, tmp2, tmp3);
    }

    valid &= nmod_poly_equal(tmp1, tmp2);
    if (!valid) goto end;

    // u_{i−1} h_hat_i + u_i * h_i = G^{LOG_Q_CEIL} w_dagger_i for i = 2,...,N-1
    for (ulong i = 1; i < ctx->N - 1; i++) {
        nmod_poly_mulmod(tmp1, u[i - 1], h_hat[i], ctx->com_ctx->cyclo_poly); // u_{i-1} * h_hat_i
        nmod_poly_mulmod(tmp2, u[i], h[i], ctx->com_ctx->cyclo_poly); // u_i * h_i
        nmod_poly_add(tmp1, tmp1, tmp2); // u_{i- 1} * h_hat_i + u_i * h_i
        // G^{LOG_Q_CEIL} w_dagger_i
        nmod_poly_zero(tmp2);
        for (int j = 0; j < LOG_Q_CEIL; j++) {
            nmod_poly_mulmod(tmp3, G[j][0], w_dagger[i][j][0], ctx->com_ctx->cyclo_poly);
            nmod_poly_add(tmp2, tmp2, tmp3);
        }
        valid &= nmod_poly_equal(tmp1, tmp2);
        if (!valid) goto end;
    }

    // u_{N−1} h_hat_N + (-1)^N \gamma h_n = G^{LOG_Q_CEIL} w_dagger_N
    nmod_poly_mulmod(tmp1, u[ctx->N - 2], h_hat[ctx->N - 1], ctx->com_ctx->cyclo_poly); // u_{N-1} * h_hat_N
    nmod_poly_mulmod(tmp2, gamma, h[ctx->N - 1], ctx->com_ctx->cyclo_poly); // gamma * h_N
    if (ctx->N % 2 != 0) nmod_poly_neg(tmp2, tmp2); // (-1)^N * gamma * h_N
    nmod_poly_add(tmp1, tmp1, tmp2); // u_{N-1} * h_hat_N + (-1)^N

    // G^{LOG_Q_CEIL} w_dagger_N
    nmod_poly_zero(tmp2);
    for (int j = 0; j < LOG_Q_CEIL; j++) {
        nmod_poly_mulmod(tmp3, G[j][0], w_dagger[ctx->N - 1][j][0], ctx->com_ctx->cyclo_poly);
        nmod_poly_add(tmp2, tmp2, tmp3);
    }
    valid &= nmod_poly_equal(tmp1, tmp2);
    if (!valid) goto end;

    for (ulong i = 0; i < ctx->N; i++) {
        valid &= evmlr_utils_is_bin(sp->d_dagger[i], (2*ETA) * (K_LWE + M_LEN))
                && evmlr_utils_is_poly_bin(sigma[i])
                && evmlr_utils_is_bin((nmod_poly_t*) w_dagger[i], LOG_Q_CEIL);
        if (!valid) goto end;
    }

    valid &= evmlr_commit_verify(W, ctx->com_ctx, (nmod_poly_t *) w_dagger, r_W);
    if (!valid) goto end;
    valid &= evmlr_commit_verify(P, ctx->com_ctx, sigma, r_P);
    if (!valid) goto end;
    valid = evmlr_commit_verify(pp->D, ctx->com_ctx, (nmod_poly_t *) sp->d_dagger, sp->r_D);

    end:
        nmod_poly_clear(tmp1);
        nmod_poly_clear(tmp2);
        nmod_poly_clear(tmp3);
        for (int i = 0; i < LOG_Q_CEIL; i++)
            nmod_poly_clear(G[i][0]);

        return valid;
}

void evmlr_shuffle_clear_sp(evmlr_shuffle_sp_t sp, const evmlr_shuffle_ctx_t ctx) {
    for (size_t i = 0; i < ctx->N; i++) {
        for (ulong j = 0; j < (2*ETA) * (K_LWE + M_LEN); j++) {
            nmod_poly_clear(sp->d_dagger[i][j]);
            free(sp->d_dagger[i][j]);
        }
        free(sp->d_dagger[i]);
    }
    for (ulong i = 0; i < 2 * K_SIS; i++) {
        nmod_poly_clear(sp->r_D[i]);
    }
}

void evmlr_shuffle_clear_pp(evmlr_shuffle_pp_t pp, const evmlr_shuffle_ctx_t ctx) {
    evmlr_commit_clear(pp->D);
    for (size_t i = 0; i < ctx->N; i++) {
        for (ulong j = 0; j < LOG_Q_CEIL; j++) {
            nmod_poly_clear(pp->c_hat[i][j]);
        }
        free(pp->c_hat[i]);
        evmlr_otse_ciphertext_clear(pp->c_star[i], ctx->hpke_ctx->otse_ctx);
    }
    free(pp->c_hat);
    free(pp->c_star);
}

void lambda_setup(nmod_poly_t Lambda[], const nmod_poly_t lambda, const nmod_poly_t cyclo_poly, ulong L) {
    for (int i = 0; i < L; i++) {
        nmod_poly_init(Lambda[i], MOD_Q);
        switch (i) {
            case 0: nmod_poly_one(Lambda[i]);         break;
            case 1: nmod_poly_set(Lambda[i], lambda); break;
            default:
                nmod_poly_mulmod(Lambda[i], Lambda[i - 1], lambda, cyclo_poly);
                break;
        }
    }
}

void calc_z(nmod_poly_t z, const evmlr_otse_ciphertext_t cipher, const nmod_poly_t* d_dagger,
            const nmod_poly_t Lambda[], const evmlr_shuffle_ctx_t ctx) {
    ulong L = ctx->hpke_ctx->otse_ctx->L;
    nmod_poly_t a[L];
    printf("Checkpoint %d: Calculating z ...\n", ++checkpoint);
    evmlr_calc_a(a, ctx->hpke_ctx->otse_ctx, d_dagger);
    printf("Checkpoint %d: a calculation done.\n", ++checkpoint);

    nmod_poly_init(z, MOD_Q);
    nmod_poly_zero(z);

    for (int j = 0; j < L; j++) {
        nmod_poly_sub(a[j], cipher->c[j], a[j]);
        nmod_poly_mulmod(a[j], a[j], Lambda[j], ctx->com_ctx->cyclo_poly);
        nmod_poly_add(z, z, a[j]);
    }

    for (int j = 0; j < L; j++) {
        nmod_poly_clear(a[j]);
    }
}

void calc_z_hat(nmod_poly_t z_hat, const nmod_poly_t* c_hat, const nmod_poly_t Lambda[], const evmlr_shuffle_ctx_t ctx) {
    ulong L = ctx->hpke_ctx->otse_ctx->L;
    nmod_poly_t tmp;

    nmod_poly_init(z_hat, MOD_Q);
    nmod_poly_zero(z_hat);
    nmod_poly_init(tmp, MOD_Q);

    for (int j = 0; j < L; j++) {
        nmod_poly_mulmod(tmp, c_hat[j], Lambda[j], ctx->com_ctx->cyclo_poly);
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
    size_t N = ctx->N;
    evmlr_commit_t P;
    nmod_poly_t sigma[ctx->com_ctx->N];
    nmod_poly_t r_P[2 * K_SIS];
    evmlr_shuffle_phase_1(sigma, P, r_P, ctx, sp);
    printf("Checkpoint %d\n", checkpoint++);

    nmod_poly_t alpha, beta, lambda, gamma;
    nmod_poly_init(alpha, MOD_Q);
    nmod_poly_init(beta, MOD_Q);
    nmod_poly_init(lambda, MOD_Q);
    nmod_poly_init(gamma, MOD_Q);
    evmlr_shuffle_sample_challenge(alpha, state);
    evmlr_shuffle_sample_challenge(beta, state);
    evmlr_shuffle_sample_challenge(lambda, state);
    evmlr_shuffle_sample_challenge(gamma, state);

    printf("Checkpoint %d\n", checkpoint++);

    nmod_poly_t h[N], h_hat[N];
    evmlr_commit_t W;
    nmod_poly_t r_W[2 * K_SIS];
    nmod_poly_t Theta[N - 1];
    nmod_poly_t w_dagger[N][LOG_Q_CEIL][1];

    evmlr_shuffle_phase_2(h, h_hat, W, r_W, Theta, w_dagger, sigma, ctx, sp, pp, alpha, beta, lambda, state);

    printf("Checkpoint %d\n", checkpoint++);

    nmod_poly_t u[N - 1];
    evmlr_shuffle_phase_3(u, h, h_hat, Theta, gamma, ctx);

    int valid = evmlr_shuffle_proof_ver(u, h, h_hat, W, r_W, (nmod_poly_t (***)) w_dagger, P, r_P, sigma, gamma, ctx, sp, pp);

    printf("Checkpoint %d\n", checkpoint++);

    // cleanup
    evmlr_commit_clear(P);
    for (size_t i = 0; i < N; i++) {
        nmod_poly_clear(sigma[i]);
        nmod_poly_clear(h[i]);
        nmod_poly_clear(h_hat[i]);
    }

    for (size_t i = 0; i < N - 1; i++) {
        nmod_poly_clear(Theta[i]);
        nmod_poly_clear(u[i]);
    }

    for (size_t i = 0; i < N; i++) {
        for (ulong j = 0; j < LOG_Q_CEIL; j++) {
            nmod_poly_clear(w_dagger[i][j][0]);
            free(w_dagger[i][j]);
        }
    }

    nmod_poly_clear(alpha);
    nmod_poly_clear(beta);
    nmod_poly_clear(lambda);
    nmod_poly_clear(gamma);

    return valid;

}

void test(evmlr_shuffle_ctx_t ctx, flint_rand_t state) {
    size_t N = ctx->N;
    size_t L = ctx->L;

    evmlr_hpke_keypair_t keypair;
    evmlr_hpke_keypair_gen(keypair, state, ctx->hpke_ctx);

    nmod_poly_t m[N][L], m_shuffle[N][L];
    for (size_t i = 0; i < N; i++) {
        for (size_t j = 0; j < L; j++) {
            nmod_poly_init(m[i][j], MOD_Q);
            nmod_poly_randtest(m[i][j], state, DEGREE_N);
            nmod_poly_init(m_shuffle[i][j], MOD_Q);
        }
    }

    ulong perm[N];
    evmlr_utils_new_perm(perm, N);

    evmlr_hpke_cipher_t ciphers[N];
    evmlr_otse_ciphertext_t c_star[N];

    nmod_poly_t d_dagger[N][(2*ETA) * (K_LWE + ctx->L)];

    for (size_t i = 0; i < N; i++) {
        evmlr_hpke_encrypt(ciphers[i], d_dagger[i], m[i], keypair->enc_keypair->pk, ctx->hpke_ctx, state);
        c_star[i]->c = ciphers[i]->otse_cipher->c;
    }

    // shuffle
    for (size_t i = 0; i < N; i++) {
        for (size_t j = 0; j < L; j++) {
            nmod_poly_set(m_shuffle[i][j], m[perm[i]][j]);
        }
    }

    nmod_poly_t r_D[2 * K_SIS];
    evmlr_commit_sample_r(r_D);

    evmlr_commit_t D;
    evmlr_commit(D, ctx->com_ctx, (nmod_poly_t *) d_dagger, r_D);

    evmlr_shuffle_sp_t sp;
    evmlr_shuffle_sp_init(sp, perm, (nmod_poly_t **) d_dagger, r_D);

    evmlr_shuffle_pp_t pp;
    evmlr_shuffle_pp_init(pp, c_star, (nmod_poly_t **) m_shuffle, D);

    TEST_BEGIN("shuffle proof works") {
        TEST_ASSERT(evmlr_shuffle_run(sp, pp, ctx, state) == 1, end)
    } TEST_END;

    end:
        evmlr_hpke_keypair_clear(keypair);
        for (size_t i = 0; i < N; i++) {
            for (size_t j = 0; j < L; j++) {
                nmod_poly_clear(m[i][j]);
                nmod_poly_clear(m_shuffle[i][j]);
            }
            evmlr_hpke_cipher_clear(ciphers[i], ctx->hpke_ctx);
            for (int j = 0; j < (2*ETA) * (K_LWE + ctx->L); ++j) {
                nmod_poly_clear(d_dagger[i][j]);
            }
        }

        for (int i = 0; i < 2 * K_SIS; ++i) {
            nmod_poly_clear(r_D[i]);
        }
        evmlr_shuffle_clear_sp(sp, ctx);
}

int main() {
    size_t n_messages = 10;

    flint_rand_t state;
    flint_rand_init(state);
    ulong seed[2];
    getrandom(seed, sizeof(ulong)*2, 0);
    flint_rand_set_seed(state, seed[0], seed[1]);

    evmlr_shuffle_ctx_t ctx;
    evmlr_shuffle_ctx_init(ctx, n_messages, M_LEN, state);

    test(ctx, state);

    evmlr_shuffle_ctx_clear(ctx);
    flint_rand_clear(state);
    return 0;
}
