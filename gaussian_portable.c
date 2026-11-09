/* Portable discrete Gaussian sampler, used instead of gaussian.cpp when
 * EVMLR_PORTABLE is set (every WebAssembly build).
 *
 * gaussian.cpp is Raymond K. Zhao's FACCT sampler (Zhao, Steinfeld and Sakzad,
 * "FACCT: FAst, Compact, and Constant-Time Discrete Gaussian Sampler over
 * Integers", IEEE TC 2020). It needs SSE4.1 (_mm_round_pd, _mm_floor_pd) and
 * the VCL vector library, which Emscripten cannot translate (see
 * https://emscripten.org/docs/porting/simd.html). This is the same algorithm,
 * step for step, in scalar C:
 *
 *  - every __m128d operation used only lane 0, so each becomes a double
 *    operation, and each bit-level cast becomes a memcpy;
 *  - _mm_round_pd(_MM_FROUND_TO_NEAREST_INT) is rint() (ties to even, under
 *    the default rounding mode, which WebAssembly cannot change);
 *  - the four-lane VCL Box-Muller becomes a loop over libm log/sqrt/sin/cos.
 *
 * Given the same random stream it returns the same samples as gaussian.cpp,
 * except where VCL's and libm's log/sin/cos round differently in the last bit
 * and that flips a rounding or rejection decision; the distribution is the
 * same. It consumes the same number of bytes from fastrandombytes() too.
 *
 * Build this file without floating-point contraction (-ffp-contract=off):
 * the SSE code does not fuse multiplies and adds. WebAssembly has no FMA, so
 * there it makes no difference. */

#include "gaussian.h"

#ifdef EVMLR_PORTABLE

#include <math.h>
#include <stdint.h>
#include <string.h>

#include "fastrandombytes.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef M_LN2
#define M_LN2 0.69314718055994530942
#endif
#ifndef M_SQRT1_2
#define M_SQRT1_2 0.70710678118654752440
#endif
#ifndef M_2_SQRTPI
#define M_2_SQRTPI 1.12837916709551257390
#endif

#define BOX_MULLER_BYTES (8 * 8)

/* Box-Muller sampler precision: 53 bits */
#define FP_PRECISION 53
static const double INT64_DOUBLE_FACTOR = 1.0 / (1LL << FP_PRECISION);
#define FP_MASK ((1LL << FP_PRECISION) - 1)

/* Parameters used by FACCT */
#define COMP_ENTRY_SIZE 9
#define EXP_MANTISSA_PRECISION 52
#define EXP_MANTISSA_MASK ((1LL << EXP_MANTISSA_PRECISION) - 1)
#define R_MANTISSA_PRECISION (EXP_MANTISSA_PRECISION + 1)
#define R_MANTISSA_MASK ((1LL << R_MANTISSA_PRECISION) - 1)
#define R_EXPONENT_L (8 * COMP_ENTRY_SIZE - R_MANTISSA_PRECISION)
#define DOUBLE_ONE (1023LL << 52)

#define DISCRETE_BYTES (1 + NORM_BATCH * COMP_ENTRY_SIZE)

static const double PORTABLE_M_1_LN2 = 1.0 / M_LN2;

/* coefficients of the FACCT exp(x) evaluation polynomial */
static const uint64_t EXP_COFF[] = {0x3e9e44a0d1ea10b4,
									0x3ec36ed25eaa5962,
									0x3efaa0f31065471a,
									0x3f29f143793315db,
									0x3f56c27c690d439b,
									0x3f811105d90e4b7e,
									0x3fa55555e692a339,
									0x3fc555555124dad1,
									0x3fe00000000f3594,
									0x3fefffffffffb3f3,
									0x3ff0000000000000};

static inline uint64_t dbl_bits(double x) {
	uint64_t u;
	memcpy(&u, &x, sizeof(u));
	return u;
}

static inline double bits_dbl(uint64_t u) {
	double x;
	memcpy(&x, &u, sizeof(x));
	return x;
}

/* _mm_cvtpd_epi32 on one lane: out-of-range and NaN inputs give the x86
 * "integer indefinite" value INT32_MIN rather than undefined behaviour. */
static inline int32_t cvt_to_i32(double x) {
	if (x >= -2147483648.0 && x < 2147483648.0) {
		return (int32_t) rint(x);
	}
	return INT32_MIN;
}

void discrete_gaussian_ctx_init(discrete_gaussian_ctx_t ctx, double sigma) {
	ctx->sigma = sigma;
	ctx->discrete_normalisation = 0.5 * M_SQRT1_2 * M_2_SQRTPI * (1.0 / sigma);
	ctx->sigma_inv = -1.0 / (2 * sigma * sigma);
	ctx->head = NORM_BATCH;
	memset(ctx->norm, 0, sizeof(ctx->norm));
}

void discrete_gaussian_ctx_clear(discrete_gaussian_ctx_t ctx) {
	memset(ctx, 0, sizeof(discrete_gaussian_ctx_struct));
}

/* Box-Muller sampler: eight continuous N(0, sigma^2) samples into ctx->norm. */
static inline void box_muller(discrete_gaussian_ctx_struct *ctx) {
	unsigned char r[BOX_MULLER_BYTES];
	uint64_t l1, l2;
	double r1, r2;

	fastrandombytes(r, BOX_MULLER_BYTES);

	for (int i = 0; i < 4; i++) {
		memcpy(&l1, r + 8 * i, sizeof(l1));
		memcpy(&l2, r + 32 + 8 * i, sizeof(l2));

		r1 = (double) (int64_t) ((l1 & FP_MASK) + 1) * INT64_DOUBLE_FACTOR;
		r2 = (double) (int64_t) ((l2 & FP_MASK) + 1) * INT64_DOUBLE_FACTOR;

		r1 = sqrt(-2.0 * log(r1)) * ctx->sigma;
		r2 = 2.0 * M_PI * r2;

		ctx->norm[i] = r1 * cos(r2);
		ctx->norm[i + 4] = r1 * sin(r2);
	}
}

/* FACCT exp(x * sigma_inv), returned as the bits of a double */
static inline uint64_t exp_constant(const double x, const discrete_gaussian_ctx_struct *ctx) {
	double vx, vx1, vx2, vsum;
	uint64_t vt;

	vx = x * ctx->sigma_inv;

	/* exp(x)=exp(floor(x/ln2)*ln2+a)=2^(floor(x/ln2))*exp(a), where a is in [0,ln2]
	 * we only evaluate exp(a) by using a polynomial */
	vx1 = floor(vx * PORTABLE_M_1_LN2);
	/* _mm_slli_epi64(_mm_cvtpd_epi32(.), 52): the 32-bit integer sits
	 * zero-extended in the low 64-bit lane before the shift. */
	vt = (uint64_t) (uint32_t) cvt_to_i32(vx1) << 52;
	vx2 = vx - vx1 * M_LN2;

	/* evaluate 2^a */
	vsum = bits_dbl(EXP_COFF[0]) * vx2 + bits_dbl(EXP_COFF[1]);
	for (int i = 2; i <= 10; i++) {
		vsum = vsum * vx2 + bits_dbl(EXP_COFF[i]);
	}

	/* combine to compute exp(x) */
	return vt + dbl_bits(vsum);
}

/* FACCT rejection check */
static inline int64_t comp(const unsigned char *r, const uint64_t res) {
	uint64_t res_mantissa, res_exponent;
	uint64_t r1, r2;
	uint64_t r_mantissa, r_exponent;

	res_mantissa = (res & EXP_MANTISSA_MASK) | (1LL << EXP_MANTISSA_PRECISION);
	res_exponent = R_EXPONENT_L - 1023 + 1 + (res >> EXP_MANTISSA_PRECISION);

	memcpy(&r1, r, sizeof(r1)); /* r is not guaranteed to be 8-byte aligned */
	r2 = (uint64_t)(r[8]);

	r_mantissa = r1 & R_MANTISSA_MASK;
	r_exponent = (r1 >> R_MANTISSA_PRECISION) | (r2 << (64 - R_MANTISSA_PRECISION));

	/* gaussian.cpp shifts by res_exponent unmasked, which is undefined in C
	 * once res_exponent leaves [0, 63] (see the README's TODO on the shift
	 * exponent); the x86 shl it compiles to masks the count to 6 bits. Mask
	 * it explicitly so that this file computes what the native build does,
	 * without the undefined behaviour. */
	return (((UINT64_C(1) << 63) ^ ((res - DOUBLE_ONE) | (DOUBLE_ONE - res))) | ((r_mantissa - res_mantissa) & (r_exponent - (UINT64_C(1) << (res_exponent & 63))))) >> 63;
}

/* See gaussian.cpp for the derivation of the merged branches. */
int64_t discrete_gaussian(const double center, discrete_gaussian_ctx_t ctx) {
	unsigned char r[DISCRETE_BYTES];

	uint64_t b, i, dcmp, dcmp1, drej, drc, dcmprc;
	double vc, vcr, vy, vyr, vyrc;

	vcr = rint(center);
	vc = vcr - center;

	drc = dbl_bits(bits_dbl(exp_constant(vc * vc, ctx)) * ctx->discrete_normalisation);

	fastrandombytes(r, COMP_ENTRY_SIZE);
	dcmprc = -comp(r, drc);

	while (1)
	{
		fastrandombytes(r, DISCRETE_BYTES);

		for (i = 0; i < 8; i++, ctx->head++)
		{
			if (ctx->head >= NORM_BATCH)
			{
				ctx->head = 0;
				box_muller(ctx);
			}

			vy = ctx->norm[ctx->head];
			b = (r[DISCRETE_BYTES - 1] >> i) & 0x01;

			b = -((-b | b) >> 63);
			dcmp = (vy <= 0.5);
			dcmp1 = dcmp ^ (vy >= -0.5);
			dcmp ^= b & dcmp1;

			/* round(y) - 1 when b == 0, round(y) + 1 when b == 1 */
			vyr = rint(vy) + bits_dbl(dbl_bits(-1.0) ^ (b & (dbl_bits(-1.0) ^ dbl_bits(1.0))));

			if (dcmp)
			{
				vyrc = vyr + vc;
				drej = exp_constant((vyrc + vy) * (vyrc - vy), ctx);

				if (comp(r + i * COMP_ENTRY_SIZE, drej))
				{
					ctx->head++;
					return (int64_t) rint(bits_dbl(dbl_bits(vyr) ^ (dcmprc & dbl_bits(vyr))) + vcr);
				}
			}
		}
	}
}

void discrete_gaussian_vec(int64_t *samples, const double center, const size_t size, discrete_gaussian_ctx_t ctx) {
	for (size_t i = 0; i < size; i++) {
		samples[i] = discrete_gaussian(center, ctx);
	}
}

#endif /* EVMLR_PORTABLE */
