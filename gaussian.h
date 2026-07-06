#ifndef GAUSSIAN_H
#define GAUSSIAN_H

#include <stdint.h>
#include <stddef.h>
#include <emmintrin.h> // for __m128d

#ifdef __cplusplus
extern "C" {
#endif

#define NORM_BATCH 8

struct discrete_gaussian_ctx_struct;

typedef struct discrete_gaussian_ctx_struct {
    double sigma;
    double discrete_normalisation;
    double sigma_inv;
    __m128d v_sigma_inv;
    __m128d v_discrete_normalisation;
    double norm[NORM_BATCH];
    uint64_t head;
} discrete_gaussian_ctx_struct;

typedef discrete_gaussian_ctx_struct discrete_gaussian_ctx_t[1];

void discrete_gaussian_ctx_init(discrete_gaussian_ctx_t ctx, double sigma);
void discrete_gaussian_ctx_clear(discrete_gaussian_ctx_t ctx);

int64_t discrete_gaussian(const double center, discrete_gaussian_ctx_t ctx);

void discrete_gaussian_vec(int64_t *samples, const double center, const size_t size, discrete_gaussian_ctx_t ctx);

#ifdef __cplusplus
}
#endif

#endif // GAUSSIAN_H
