#include <assert.h>
#include "sys/random.h"
#include "evmlr_enc.h"
#include "evmlr_utils.h"

void evmlr_enc_keypair_gen(evmlr_enc_keypair_t keypair, flint_rand_t state) {
    // Generate secret key s and error e with coefficients from B_eta
    nmod_poly_t e[K_LWE];
    for (size_t i = 0; i < K_LWE; i++) {
        nmod_poly_init(keypair->sk->s[i], MOD_Q);
        nmod_poly_init(e[i], MOD_Q);
        evmlr_utils_binom_sample_ring(keypair->sk->s[i], ETA);
        evmlr_utils_binom_sample_ring(e[i], ETA);
    }

    // Generate public matrix A with uniform random coefficients
    for (size_t i = 0; i < K_LWE; i++) {
        for (size_t j = 0; j < K_LWE; j++) {
            nmod_poly_init(keypair->pk->A[i][j], MOD_Q);
            nmod_poly_randtest(keypair->pk->A[i][j], state, MOD_N);
        }
    }

    // Compute t = A * s + e
    nmod_poly_t tmp;
    nmod_poly_init(tmp, MOD_Q);
    for (size_t i = 0; i < K_LWE; i++) {
        nmod_poly_init(keypair->pk->t[i], MOD_Q);
        nmod_poly_zero(keypair->pk->t[i]);
        for (size_t j = 0; j < K_LWE; j++) {
            nmod_poly_mul(tmp, keypair->pk->A[i][j], keypair->sk->s[j]);
            nmod_poly_add(keypair->pk->t[i], keypair->pk->t[i], tmp);
        }
        nmod_poly_add(keypair->pk->t[i], keypair->pk->t[i], e[i]);
    }
}

void evmlr_enc_clear_keypair(evmlr_enc_keypair_t keypair) {
    // Clear secret key
    for (size_t i = 0; i < K_LWE; i++) {
        nmod_poly_clear(keypair->sk->s[i]);
    }
    // Clear public key
    for (size_t i = 0; i < K_LWE; i++) {
        for (size_t j = 0; j < K_LWE; j++) {
            nmod_poly_clear(keypair->pk->A[i][j]);
        }
        nmod_poly_clear(keypair->pk->t[i]);
    }
}

void evmlr_enc(evmlr_enc_cipher_t cipher, const nmod_poly_t msg, const evmlr_enc_pk_t pk) {
    assert(evmlr_utils_is_poly_bin(msg));
    nmod_poly_t r[K_LWE], e2[K_LWE], e3, tmp, msg_scaled;
    // Sample r and e2 from B_eta
    for (size_t i = 0; i < K_LWE; i++) {
        nmod_poly_init(r[i], MOD_Q);
        nmod_poly_init(e2[i], MOD_Q);
        evmlr_utils_binom_sample_ring(r[i], ETA);
        evmlr_utils_binom_sample_ring(e2[i], ETA);
    }
    nmod_poly_init(e3, MOD_Q);
    evmlr_utils_binom_sample_ring(e3, ETA);

    nmod_poly_init(tmp, MOD_Q);
    // Compute u^T = r^T A + e_2^T
    for (size_t i = 0; i < K_LWE; i++) {
        nmod_poly_init(cipher->uT[i], MOD_Q);
        nmod_poly_zero(cipher->uT[i]);
        for (size_t j = 0; j < K_LWE; j++) {
            nmod_poly_mul(tmp, r[j], pk->A[j][i]);
            nmod_poly_add(cipher->uT[i], cipher->uT[i], tmp);
        }
        nmod_poly_add(cipher->uT[i], cipher->uT[i], e2[i]);
    }
    // c = q/2 getting the closest integer with ties rounded up
    ulong c = (MOD_Q + 1) / 2;
    // Compute v = r^T t + e_3 + c * msg
    nmod_poly_init(cipher->v, MOD_Q);
    nmod_poly_zero(cipher->v);
    for (size_t i = 0; i < K_LWE; i++) {
        nmod_poly_mul(tmp, r[i], pk->t[i]);
        nmod_poly_add(cipher->v, cipher->v, tmp);
    }
    nmod_poly_add(cipher->v, cipher->v, e3);

    nmod_poly_init(msg_scaled, MOD_Q);
    nmod_poly_scalar_mul_nmod(msg_scaled, msg, c);

    nmod_poly_add(cipher->v, cipher->v, msg_scaled);

    nmod_poly_clear(tmp);
    nmod_poly_clear(msg_scaled);
    for (size_t i = 0; i < K_LWE; i++) {
        nmod_poly_clear(r[i]);
        nmod_poly_clear(e2[i]);
    }
    nmod_poly_clear(e3);
}

void evmlr_dec(nmod_poly_t msg, const evmlr_enc_cipher_t cipher, const evmlr_enc_sk_t sk) {
    // \tilde{m} = v - u^T s
    nmod_poly_t m_tilde, tmp;
    nmod_poly_init(m_tilde, MOD_Q);
    nmod_poly_init(tmp, MOD_Q);
    nmod_poly_set(m_tilde, cipher->v);
    for (size_t i = 0; i < K_LWE; i++) {
        nmod_poly_mul(tmp, cipher->uT[i], sk->s[i]);
        nmod_poly_sub(m_tilde, m_tilde, tmp);
    }

    nmod_poly_zero(msg);
    // msg = 0 if \tilde{m} is closer to 0 than to q/2, and 1 otherwise
    // TODO check if this is correct
    for (slong i = 0; i < MOD_N; i++) {
        ulong coeff = nmod_poly_get_coeff_ui(m_tilde, i);
        if (coeff < (MOD_Q + 1) / 2) {
            nmod_poly_set_coeff_ui(msg, i, 0);
        } else {
            nmod_poly_set_coeff_ui(msg, i, 1);
        }
    }

    nmod_poly_clear(m_tilde);
    nmod_poly_clear(tmp);
}

void evmlr_enc_clear_cipher(evmlr_enc_cipher_t cipher) {
    for (size_t i = 0; i < K_LWE; i++) {
        nmod_poly_clear(cipher->uT[i]);
    }
    nmod_poly_clear(cipher->v);
}

int main() {
    printf("EVMLR Encryption Scheme Test\n");
    flint_rand_t state;
    flint_rand_init(state);
    ulong seed[2];
    getrandom(seed, sizeof(ulong)*2, 0);
    flint_rand_set_seed(state, seed[0], seed[1]);

    // Keypair generation
    printf("Generating keypair...\n");
    evmlr_enc_keypair_t keypair;
    evmlr_enc_keypair_gen(keypair, state);

    // Message to be encrypted (binary polynomial)
    printf("Encrypting message...\n");
    nmod_poly_t msg;
    nmod_poly_init(msg, MOD_Q);
    evmlr_utils_int_to_bin(msg, 5); // Example message

    // Encryption
    evmlr_enc_cipher_t cipher;
    evmlr_enc(cipher, msg, keypair->pk);

    // Decryption
    printf("Decrypting message...\n");
    nmod_poly_t decrypted_msg;
    nmod_poly_init(decrypted_msg, MOD_Q);
    evmlr_dec(decrypted_msg, cipher, keypair->sk);

    // Check if decryption is correct
    assert(evmlr_utils_is_poly_bin(decrypted_msg));
    assert(nmod_poly_equal(msg, decrypted_msg));
    printf("Decryption successful!\n");

    // Clear memory
    evmlr_enc_clear_cipher(cipher);
    evmlr_enc_clear_keypair(keypair);
    nmod_poly_clear(msg);
    nmod_poly_clear(decrypted_msg);

    flint_rand_clear(state);
    return 0;
}