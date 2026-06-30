#ifndef EVMLR_SHUFFLE_EVMLR_CHALLENGE_H
#define EVMLR_SHUFFLE_EVMLR_CHALLENGE_H

#include "flint/nmod_poly.h"
#include "flint/nmod_poly_mat.h"
#include "sha.h"
#include "evmlr_params.h"

typedef struct {
    SHA256Context sha;
    uint8_t hash[SHA256HashSize];
    int finalized;
} evmlr_challenge_struct;
typedef evmlr_challenge_struct evmlr_challenge_t[1];


void evmlr_challenge_init(evmlr_challenge_t chal);
void evmlr_challenge_add_matrix(evmlr_challenge_t chal, const nmod_poly_mat_t mat);
void evmlr_challenge_add_poly(evmlr_challenge_t chal, const nmod_poly_t poly);
void evmlr_challenge_add_bytes(evmlr_challenge_t chal, const uint8_t* bytes, size_t len);

void evmlr_challenge_get_poly_half(nmod_poly_t out, evmlr_challenge_t chal);
void evmlr_challenge_get_poly_one_ternary(nmod_poly_t out, evmlr_challenge_t chal);
void evmlr_challenge_get_poly_symmetric(nmod_poly_t out, evmlr_challenge_t chal);
void evmlr_challenge_get_val(ulong* val, evmlr_challenge_t chal);
void evmlr_challenge_get_vector(ulong* vec, slong len, evmlr_challenge_t chal);
void evmlr_challenge_get_hash(uint8_t hash[SHA256HashSize], evmlr_challenge_t chal);

#endif // EVMLR_SHUFFLE_EVMLR_CHALLENGE_H
