/* Portable AES-256-CTR, used instead of the AES-NI code in fastrandombytes.c
 * when EVMLR_PORTABLE is set (every WebAssembly build).
 *
 * It produces the same bytes as the AES-NI version for the same seed: the
 * same key schedule, the counter block starts at zero and its low 64 bits
 * are incremented as a little-endian integer once per 16-byte block, and the
 * unused tail of a final partial block is discarded.
 *
 * The block cipher is the 32-bit bitsliced AES from BearSSL (aes_ct, by
 * Thomas Pornin, MIT licence), which uses no secret-dependent table lookups
 * or branches. A table-based AES would leak the key, which seeds the stream
 * the secret proof masks come from, through cache timing. Two blocks are
 * encrypted per call of the core. */

#include "fastrandombytes.h"

#ifdef EVMLR_PORTABLE

#include <stdint.h>
#include <string.h>

#define AES256_ROUNDS 14

static inline uint32_t load_le32(const unsigned char *p) {
	return (uint32_t) p[0] | ((uint32_t) p[1] << 8)
		| ((uint32_t) p[2] << 16) | ((uint32_t) p[3] << 24);
}

static inline void store_le32(unsigned char *p, uint32_t x) {
	p[0] = (unsigned char) x;
	p[1] = (unsigned char) (x >> 8);
	p[2] = (unsigned char) (x >> 16);
	p[3] = (unsigned char) (x >> 24);
}

/* The AES S-box on eight bitsliced words, as a Boyar-Peralta circuit. */
static void bitslice_sbox(uint32_t *q) {
	uint32_t x0, x1, x2, x3, x4, x5, x6, x7;
	uint32_t y1, y2, y3, y4, y5, y6, y7, y8, y9;
	uint32_t y10, y11, y12, y13, y14, y15, y16, y17, y18, y19;
	uint32_t y20, y21;
	uint32_t z0, z1, z2, z3, z4, z5, z6, z7, z8, z9;
	uint32_t z10, z11, z12, z13, z14, z15, z16, z17;
	uint32_t t0, t1, t2, t3, t4, t5, t6, t7, t8, t9;
	uint32_t t10, t11, t12, t13, t14, t15, t16, t17, t18, t19;
	uint32_t t20, t21, t22, t23, t24, t25, t26, t27, t28, t29;
	uint32_t t30, t31, t32, t33, t34, t35, t36, t37, t38, t39;
	uint32_t t40, t41, t42, t43, t44, t45, t46, t47, t48, t49;
	uint32_t t50, t51, t52, t53, t54, t55, t56, t57, t58, t59;
	uint32_t t60, t61, t62, t63, t64, t65, t66, t67;
	uint32_t s0, s1, s2, s3, s4, s5, s6, s7;

	x0 = q[7];
	x1 = q[6];
	x2 = q[5];
	x3 = q[4];
	x4 = q[3];
	x5 = q[2];
	x6 = q[1];
	x7 = q[0];

	/* Top linear transformation. */
	y14 = x3 ^ x5;
	y13 = x0 ^ x6;
	y9 = x0 ^ x3;
	y8 = x0 ^ x5;
	t0 = x1 ^ x2;
	y1 = t0 ^ x7;
	y4 = y1 ^ x3;
	y12 = y13 ^ y14;
	y2 = y1 ^ x0;
	y5 = y1 ^ x6;
	y3 = y5 ^ y8;
	t1 = x4 ^ y12;
	y15 = t1 ^ x5;
	y20 = t1 ^ x1;
	y6 = y15 ^ x7;
	y10 = y15 ^ t0;
	y11 = y20 ^ y9;
	y7 = x7 ^ y11;
	y17 = y10 ^ y11;
	y19 = y10 ^ y8;
	y16 = t0 ^ y11;
	y21 = y13 ^ y16;
	y18 = x0 ^ y16;

	/* Non-linear section. */
	t2 = y12 & y15;
	t3 = y3 & y6;
	t4 = t3 ^ t2;
	t5 = y4 & x7;
	t6 = t5 ^ t2;
	t7 = y13 & y16;
	t8 = y5 & y1;
	t9 = t8 ^ t7;
	t10 = y2 & y7;
	t11 = t10 ^ t7;
	t12 = y9 & y11;
	t13 = y14 & y17;
	t14 = t13 ^ t12;
	t15 = y8 & y10;
	t16 = t15 ^ t12;
	t17 = t4 ^ t14;
	t18 = t6 ^ t16;
	t19 = t9 ^ t14;
	t20 = t11 ^ t16;
	t21 = t17 ^ y20;
	t22 = t18 ^ y19;
	t23 = t19 ^ y21;
	t24 = t20 ^ y18;

	t25 = t21 ^ t22;
	t26 = t21 & t23;
	t27 = t24 ^ t26;
	t28 = t25 & t27;
	t29 = t28 ^ t22;
	t30 = t23 ^ t24;
	t31 = t22 ^ t26;
	t32 = t31 & t30;
	t33 = t32 ^ t24;
	t34 = t23 ^ t33;
	t35 = t27 ^ t33;
	t36 = t24 & t35;
	t37 = t36 ^ t34;
	t38 = t27 ^ t36;
	t39 = t29 & t38;
	t40 = t25 ^ t39;

	t41 = t40 ^ t37;
	t42 = t29 ^ t33;
	t43 = t29 ^ t40;
	t44 = t33 ^ t37;
	t45 = t42 ^ t41;
	z0 = t44 & y15;
	z1 = t37 & y6;
	z2 = t33 & x7;
	z3 = t43 & y16;
	z4 = t40 & y1;
	z5 = t29 & y7;
	z6 = t42 & y11;
	z7 = t45 & y17;
	z8 = t41 & y10;
	z9 = t44 & y12;
	z10 = t37 & y3;
	z11 = t33 & y4;
	z12 = t43 & y13;
	z13 = t40 & y5;
	z14 = t29 & y2;
	z15 = t42 & y9;
	z16 = t45 & y14;
	z17 = t41 & y8;

	/* Bottom linear transformation. */
	t46 = z15 ^ z16;
	t47 = z10 ^ z11;
	t48 = z5 ^ z13;
	t49 = z9 ^ z10;
	t50 = z2 ^ z12;
	t51 = z2 ^ z5;
	t52 = z7 ^ z8;
	t53 = z0 ^ z3;
	t54 = z6 ^ z7;
	t55 = z16 ^ z17;
	t56 = z12 ^ t48;
	t57 = t50 ^ t53;
	t58 = z4 ^ t46;
	t59 = z3 ^ t54;
	t60 = t46 ^ t57;
	t61 = z14 ^ t57;
	t62 = t52 ^ t58;
	t63 = t49 ^ t58;
	t64 = z4 ^ t59;
	t65 = t61 ^ t62;
	t66 = z1 ^ t63;
	s0 = t59 ^ t63;
	s6 = t56 ^ ~t62;
	s7 = t48 ^ ~t60;
	t67 = t64 ^ t65;
	s3 = t53 ^ t66;
	s4 = t51 ^ t66;
	s5 = t47 ^ t65;
	s1 = t64 ^ ~s3;
	s2 = t55 ^ ~t67;

	q[7] = s0;
	q[6] = s1;
	q[5] = s2;
	q[4] = s3;
	q[3] = s4;
	q[2] = s5;
	q[1] = s6;
	q[0] = s7;
}

/* Converts eight words between the natural and the bitsliced representation
 * (the transform is an involution). */
static void ortho(uint32_t *q) {
#define SWAPN(cl, ch, s, x, y) do { \
		uint32_t a_ = (x), b_ = (y); \
		(x) = (a_ & (uint32_t) (cl)) | ((b_ & (uint32_t) (cl)) << (s)); \
		(y) = ((a_ & (uint32_t) (ch)) >> (s)) | (b_ & (uint32_t) (ch)); \
	} while (0)
#define SWAP2(x, y) SWAPN(0x55555555, 0xAAAAAAAA, 1, x, y)
#define SWAP4(x, y) SWAPN(0x33333333, 0xCCCCCCCC, 2, x, y)
#define SWAP8(x, y) SWAPN(0x0F0F0F0F, 0xF0F0F0F0, 4, x, y)

	SWAP2(q[0], q[1]);
	SWAP2(q[2], q[3]);
	SWAP2(q[4], q[5]);
	SWAP2(q[6], q[7]);

	SWAP4(q[0], q[2]);
	SWAP4(q[1], q[3]);
	SWAP4(q[4], q[6]);
	SWAP4(q[5], q[7]);

	SWAP8(q[0], q[4]);
	SWAP8(q[1], q[5]);
	SWAP8(q[2], q[6]);
	SWAP8(q[3], q[7]);

#undef SWAP8
#undef SWAP4
#undef SWAP2
#undef SWAPN
}

/* Loads two 16-byte blocks into the bitsliced representation. */
static void load_blocks(uint32_t *q, const unsigned char *b0, const unsigned char *b1) {
	for (int i = 0; i < 4; i++) {
		q[2 * i] = load_le32(b0 + 4 * i);
		q[2 * i + 1] = load_le32(b1 + 4 * i);
	}
	ortho(q);
}

static void store_blocks(unsigned char *b0, unsigned char *b1, uint32_t *q) {
	ortho(q);
	for (int i = 0; i < 4; i++) {
		store_le32(b0 + 4 * i, q[2 * i]);
		store_le32(b1 + 4 * i, q[2 * i + 1]);
	}
}

/* SubWord of the key schedule, through the same constant-time S-box. */
static uint32_t sub_word(uint32_t x) {
	uint32_t q[8];

	memset(q, 0, sizeof(q));
	q[0] = x;
	ortho(q);
	bitslice_sbox(q);
	ortho(q);
	return q[0];
}

static inline void add_round_key(uint32_t *q, const uint32_t *sk) {
	for (int i = 0; i < 8; i++) {
		q[i] ^= sk[i];
	}
}

static inline void shift_rows(uint32_t *q) {
	for (int i = 0; i < 8; i++) {
		uint32_t x = q[i];

		q[i] = (x & 0x000000FF)
			| ((x & 0x0000FC00) >> 2) | ((x & 0x00000300) << 6)
			| ((x & 0x00F00000) >> 4) | ((x & 0x000F0000) << 4)
			| ((x & 0xC0000000) >> 6) | ((x & 0x3F000000) << 2);
	}
}

static inline uint32_t rotr16(uint32_t x) {
	return (x << 16) | (x >> 16);
}

static inline void mix_columns(uint32_t *q) {
	uint32_t q0, q1, q2, q3, q4, q5, q6, q7;
	uint32_t r0, r1, r2, r3, r4, r5, r6, r7;

	q0 = q[0];
	q1 = q[1];
	q2 = q[2];
	q3 = q[3];
	q4 = q[4];
	q5 = q[5];
	q6 = q[6];
	q7 = q[7];
	r0 = (q0 >> 8) | (q0 << 24);
	r1 = (q1 >> 8) | (q1 << 24);
	r2 = (q2 >> 8) | (q2 << 24);
	r3 = (q3 >> 8) | (q3 << 24);
	r4 = (q4 >> 8) | (q4 << 24);
	r5 = (q5 >> 8) | (q5 << 24);
	r6 = (q6 >> 8) | (q6 << 24);
	r7 = (q7 >> 8) | (q7 << 24);

	q[0] = q7 ^ r7 ^ r0 ^ rotr16(q0 ^ r0);
	q[1] = q0 ^ r0 ^ q7 ^ r7 ^ r1 ^ rotr16(q1 ^ r1);
	q[2] = q1 ^ r1 ^ r2 ^ rotr16(q2 ^ r2);
	q[3] = q2 ^ r2 ^ q7 ^ r7 ^ r3 ^ rotr16(q3 ^ r3);
	q[4] = q3 ^ r3 ^ q7 ^ r7 ^ r4 ^ rotr16(q4 ^ r4);
	q[5] = q4 ^ r4 ^ r5 ^ rotr16(q5 ^ r5);
	q[6] = q5 ^ r5 ^ r6 ^ rotr16(q6 ^ r6);
	q[7] = q6 ^ r6 ^ r7 ^ rotr16(q7 ^ r7);
}

/* Encrypts two blocks in place under the bitsliced round keys. */
static void aes256_encrypt2(const uint32_t *sk, unsigned char *b0, unsigned char *b1) {
	uint32_t q[8];

	load_blocks(q, b0, b1);
	add_round_key(q, sk);
	for (int u = 1; u < AES256_ROUNDS; u++) {
		bitslice_sbox(q);
		shift_rows(q);
		mix_columns(q);
		add_round_key(q, sk + 8 * u);
	}
	bitslice_sbox(q);
	shift_rows(q);
	add_round_key(q, sk + 8 * AES256_ROUNDS);
	store_blocks(b0, b1, q);
	memset(q, 0, sizeof(q));
}

/* round_key <-- aes256_key_expansion(seed), iv <-- 0 */
void prg_setseed(prg_t *prg, const unsigned char *seed) {
	static const uint32_t rcon[7] = { 0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40 };
	uint32_t w[4 * (AES256_ROUNDS + 1)];
	unsigned char rk[16];

	for (int i = 0; i < 8; i++) {
		w[i] = load_le32(seed + 4 * i);
	}
	for (int i = 8; i < 4 * (AES256_ROUNDS + 1); i++) {
		uint32_t t = w[i - 1];

		if (i % 8 == 0) {
			t = sub_word((t >> 8) | (t << 24)) ^ rcon[i / 8 - 1];
		} else if (i % 8 == 4) {
			t = sub_word(t);
		}
		w[i] = w[i - 8] ^ t;
	}

	/* Round key u, put through the same transform as the data. Both lanes
	 * hold the key, since both blocks are encrypted under it. */
	for (int u = 0; u <= AES256_ROUNDS; u++) {
		for (int i = 0; i < 4; i++) {
			store_le32(rk + 4 * i, w[4 * u + i]);
		}
		load_blocks(prg->round_key + 8 * u, rk, rk);
	}

	memset(prg->iv, 0, sizeof(prg->iv));
	memset(w, 0, sizeof(w));
	memset(rk, 0, sizeof(rk));
}

/* The low 64 bits of the counter block, as a little-endian integer, plus one;
 * the AES-NI version adds 1 to the low 64-bit lane in the same way. */
static inline void ctr_increment(uint8_t *iv) {
	for (int i = 0; i < 8; i++) {
		if (++iv[i] != 0) {
			break;
		}
	}
}

/* r <-- aes256_ctr(round_key, iv, rlen) */
void prg_bytes(prg_t *prg, unsigned char *r, unsigned long long rlen) {
	unsigned char b0[16], b1[16];

	while (rlen > 0) {
		memcpy(b0, prg->iv, 16);
		ctr_increment(prg->iv);
		memcpy(b1, prg->iv, 16);
		aes256_encrypt2(prg->round_key, b0, b1);

		/* Like the AES-NI version, a block whose output is not fully used
		 * still consumes its counter value, and its tail is dropped. */
		if (rlen <= 16) {
			memcpy(r, b0, rlen);
			break;
		}
		memcpy(r, b0, 16);
		r += 16;
		rlen -= 16;

		ctr_increment(prg->iv);
		if (rlen <= 16) {
			memcpy(r, b1, rlen);
			break;
		}
		memcpy(r, b1, 16);
		r += 16;
		rlen -= 16;
	}
	memset(b0, 0, sizeof(b0));
	memset(b1, 0, sizeof(b1));
}

#endif /* EVMLR_PORTABLE */
