#include "evmlr_commit.h"
#include "math.h"
#include "evmlr_utils.h"

// B := \sqrt{nL + n 2 K_SIS ETA^2}
double compute_b() {
    return sqrt(MOD_N * L + MOD_N * 2 * K_SIS * ETA * ETA);
}

void evmlr_commit_scheme_init(evmlr_commit_scheme_t pp, flint_rand_t state) {
    pp->b = compute_b();

    // A_1 \sample R_q^{K_SIS x L}
    for (size_t i = 0; i < K_SIS; i++) {
        for (size_t j = 0; j < L; j++) {
            nmod_poly_init(pp->A_1[i][j], MOD_Q);
            nmod_poly_randtest(pp->A_1[i][j], state, MOD_N);
        }
    }
    // A_2 \sample R_q^{K_SIS x 2K_SIS}
    for (size_t i = 0; i < K_SIS; i++) {
        for (size_t j = 0; j < 2 * K_SIS; j++) {
            nmod_poly_init(pp->A_2[i][j], MOD_Q);
            nmod_poly_randtest(pp->A_2[i][j], state, MOD_N);
        }
    }
}

void evmlr_commit_scheme_clear(evmlr_commit_scheme_t pp) {
    // Clear A_1
    for (size_t i = 0; i < K_SIS; i++) {
        for (size_t j = 0; j < L; j++) {
            nmod_poly_clear(pp->A_1[i][j]);
        }
    }
    // Clear A_2
    for (size_t i = 0; i < K_SIS; i++) {
        for (size_t j = 0; j < 2 * K_SIS; j++) {
            nmod_poly_clear(pp->A_2[i][j]);
        }
    }
}

void evmlr_sample_r(nmod_poly_t r[2*K_SIS]) {
    for (size_t i = 0; i < 2 * K_SIS; i++) {
        nmod_poly_init(r[i], MOD_Q);
        evmlr_utils_binom_sample_ring(r[i], ETA);
    }
}

void evmlr_commit(evmlr_commit_t com, const nmod_poly_t msg[L], const nmod_poly_t r[2*K_SIS], const evmlr_commit_scheme_t pp) {
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
        for (size_t j = 0; j < L; j++) {
            nmod_poly_mul(tmp, pp->A_1[i][j], msg[j]);
            nmod_poly_add(com->c[i], com->c[i], tmp);
        }
        // A_2 * r
        for (size_t j = 0; j < 2 * K_SIS; j++) {
            nmod_poly_mul(tmp, pp->A_2[i][j], r[j]);
            nmod_poly_add(com->c[i], com->c[i], tmp);
        }
    }
    nmod_poly_clear(tmp);
}

void evmlr_commit_clear(evmlr_commit_t com, const evmlr_commit_scheme_t pp) {
    for (size_t i = 0; i < K_SIS; i++) {
        nmod_poly_clear(com->c[i]);
    }
}

int evmlr_commit_verify(const evmlr_commit_t com, const nmod_poly_t msg[L], const nmod_poly_t r[2*K_SIS], const evmlr_commit_scheme_t pp) {
    int valid = 1;
    // ||\binom{m}{r}|| <= B
    // TODO: Implement norm check if needed

    if (!valid) {
        return 0;
    }
    // c = A_1 * msg + A_2 * r
    evmlr_commit_t recomputed_com;
    evmlr_commit(recomputed_com, msg, r, pp);

    for (int i = 0; i < K_SIS && valid; i++) {
        valid &= nmod_poly_equal(com->c[i], recomputed_com->c[i]);
    }

    evmlr_commit_clear(recomputed_com, pp);
    return valid;
}