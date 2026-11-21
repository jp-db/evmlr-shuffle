# WebAssembly build

The voting protocol (key generation, ballot encryption with its proofs of
correct encryption, and the mix steps with their proofs of shuffle) compiles
to WebAssembly with Emscripten. The native build (`make`) is unchanged; the WASM build is a
separate makefile that compiles the same sources with `emcc` and swaps in
portable versions of the code that needs x86.

## Building

```sh
source /path/to/emsdk/emsdk_env.sh   # emcc on PATH
wasm/build-deps.sh                   # once: GMP, MPFR, FLINT for wasm32 -> wasm/deps
make -f Makefile.wasm                # library + Node test binaries -> wasm/build
make -f Makefile.wasm test           # every test suite, under Node
make -f Makefile.wasm demo           # a 3-server election through the JS API
```

`build-deps.sh` downloads GMP 6.3.0, MPFR 4.2.2 and FLINT 3.6.0 (override with
`GMP_VERSION=` etc.) and builds them as static libraries, without assembly or
threads. FLINT takes a few minutes.

The outputs:

| File | What it is |
| --- | --- |
| `wasm/build/evmlr.mjs` + `evmlr.wasm` | ES module factory (`createEvmlrModule`) exporting the C API of `evmlr_api.h`, for browsers, workers and Node |
| `wasm/evmlr.mjs` | JS wrapper, see below |
| `wasm/demo.mjs` | servers, authority, voters, bulletin board and auditor exchanging only bytes (`make -f Makefile.wasm demo`) |
| `wasm/build/evmlr_*.js` | the test/benchmark binaries; `node wasm/build/evmlr_voting.js test` |

The module needs WASM SIMD (Chrome 91, Firefox 89, Safari 16.4, Node 16.4
or later), which the Gaussian sampler runs on (see below). It uses no
threads, so it needs no SharedArrayBuffer or cross-origin isolation. Casting
takes a while; in a page, run it in a Web Worker.

## JavaScript API

The module exports the functions of `evmlr_api.h` (the same API as the
native `libevmlr.so`); `wasm/evmlr.mjs` wraps them so that every value is a
`Uint8Array`:

```js
import { loadEvmlr } from './wasm/evmlr.mjs';
const evmlr = await loadEvmlr();

// Each mixing server, on its own
const { publicKey, keyPair } = evmlr.keygen();   // publish publicKey; keyPair stays with the server

// Election authority: public keys in, public parameters out
const id = await evmlr.electionId('Board election 2026');   // 32 bytes
const pp = evmlr.params(id, publicKeys);                     // in mixing order

// Voter page: the public parameters are all it needs
const ballot = evmlr.cast(pp, voterId, candidate);

// Bulletin board
evmlr.checkBallot(pp, ballot);       // the voter id, or throws EvmlrError

// Servers, one after the other, each with its own key pair
let batch = evmlr.firstBatch(pp, board);
const mixes = [];
for (let j = 0; j < k; j++) {
  const step = evmlr.mix(pp, keyPairOfServer[j], batch, j ? mixes[j - 1] : null);
  mixes.push(step.mix);              // publish
  batch = step.batch;                // send to server j + 1
}
evmlr.mixInfo(pp, mixes[0]);         // { server, received, passedOn, rejected: [{ index, reason }] }
const tally = evmlr.results(pp, batch).map((v) => v[0]);

// Anyone: check the count from public data (parameters, board, mix records)
const verified = evmlr.verify(pp, board, evmlr.transcript(pp, mixes));
```

A vote is a number or up to `evmlr.voteLen` (512) coefficients mod
`evmlr.modulus` (3109); votes come back as `Uint16Array`s of 512
coefficients, so `v[0]` is the candidate number. Errors throw `EvmlrError`,
whose `status` is an `evmlr_status`. No handles are kept: every call decodes
its inputs, and all WASM memory is freed before it returns.

`evmlr.mix` refuses a key pair that is not the batch's server's
(`EVMLR_ERR_KEY`), and for server j > 0 a batch that is not what server
j − 1's record says it output (`EVMLR_ERR_INVALID`).

## Serialization format

`evmlr_serialize.h` (also used natively) encodes the public parameters, a
server's public key and key pair, the secret keys of a one-process setup,
ballots with their proofs of correct encryption, batches, mix records,
proofs of shuffle with their statements, and counting transcripts. Every
encoding has:

- a 20-byte header: magic, format version (3), object type, and a
  fingerprint of the scheme parameters (`MOD_Q`, `DEGREE_N`, `K_LWE`, ...),
  so bytes from a build with different parameters are refused;
- fixed-width little-endian integers, the same on wasm32 and x86-64;
- polynomial coefficients bit-packed at 12 bits mod q (6 bits mod 2^zeta).

Decoding is strict: every dimension must be the one the protocol fixes,
coefficients must be reduced, and polynomials must be normalised. So every
object has exactly one encoding, and what is read back hashes into the
Fiat-Shamir transcripts like the original. A truncated, padded or tampered
ballot is refused before any proof is checked. A key pair is refused if its
secret does not match its public key.

The public parameters carry only the election id, k and the servers' public
keys. The OTSE matrices `H` and `H''` are expanded from a SHA-256 hash of
the election id and k, and the shuffle's commitment key from a hash of the
parameters and the number of ciphertexts mixed; both are rebuilt by every
party, and neither can be chosen by anyone.

| Object | Size, default parameters, 3 servers |
| --- | --- |
| public parameters | 14 KB |
| public key / key pair | 4.7 KB / 6.2 KB |
| ballot | 100 KB |
| mix record, 5 ciphertexts | 0.7 MB (grows linearly) |

### Mixing and verifiable counting

A batch is what server j hands to server j + 1: its outputs, each a
ciphertext of the next layer, with the proofs of encryption of the layers
left (server 0's batch, from `evmlr_voting_batch_from_submissions`, also
carries the voter ids). A mix record is what server j publishes:

- how many ciphertexts it received, and which it rejected, with the reason:
  the ML-PKE half of the proof of encryption fails, or the OTSE half fails
  against the decryption;
- the statement of its shuffle of the rest (the commitment `D`, the outputs
  `c_hat`, the inputs `c_star`) and the proof.

A transcript is the k mix records. `evmlr_voting_verify_counting` checks it
with public data only:

- of the ballots that pass the public check, the last of each voter makes
  server 0's input, in board order. A ballot that fails the check is skipped
  before the revote policy, so a copy posted under a voter's id cannot
  replace that voter's ballot. Who may post under an id is for the
  application to enforce (see the top-level README);
- server j's inputs, without the ones it rejected, are server j − 1's
  outputs (or those ballots' outer OTSE ciphertexts), in order;
- the rejections are well formed, and at server 0 every rejection for the
  ML-PKE proof really fails that public check;
- every proof of shuffle verifies.

It returns the tally the last server's outputs hold. It does not show that
the servers decrypted correctly, and so it cannot check rejections for the
OTSE proof either; that proof is not implemented yet (see the top-level
README). Neither can it stop the inner-layer proofs, which ballots publish,
from linking ballots through the shuffles. The paper avoids that with a
single proof per ballot over the whole layered encryption, checked by the
servers and then discarded (top-level README, TODO).

### Binding ballots to the election and the voter

The Fiat-Shamir challenges of a ballot's proofs hash a label: the election
id, k and the layer, and for the outermost layer the voter id. Without it a
ballot could be resubmitted under another voter id and counted twice (the
ballot-copying attack on Helios-style systems), or replayed in another
election. The outer proof's ML-PKE half needs no secret to check, so the
bulletin board checks it on arrival. Inner layers are not bound to the
voter: their ciphertexts stay hidden inside the outer layer until server 0
decrypts it, and binding them would make the servers carry voter ids past
the shuffles.

## What changes for WebAssembly

Native builds compile the same code as before: every WASM-only path is behind
`__EMSCRIPTEN__`, or behind `EVMLR_PORTABLE`, which `evmlr_portable.h`
defines only for Emscripten.

| Native | Problem in WASM | WASM replacement |
| --- | --- | --- |
| `getrandom(2)` (`evmlr_utils.c`, `evmlr_main.h`), `/dev/urandom` (`randombytes.c`) | no system calls | `getentropy()`, which Emscripten implements with `crypto.getRandomValues` (browser) or `crypto.randomFillSync` (Node), in chunks of 256 bytes |
| AES-256-CTR with AES-NI (`fastrandombytes.c`) | no AES-NI | `fastrandombytes_portable.c`: bitsliced, constant-time AES (BearSSL's `aes_ct` design) |
| FACCT Gaussian with SSE4.1 + VCL (`gaussian.cpp`) | no SSE natively; Emscripten has no `<x86intrin.h>` and no `_mm_setcsr` | none: `gaussian.cpp`, `gaussian.h` and `vcl/` are unchanged, compiled with `-msimd128 -msse4.1` on [Emscripten's SSE emulation](https://emscripten.org/docs/porting/simd.html), with `wasm/compat/x86intrin.h` filling the two gaps |
| `rdtsc` (`cpucycles.c`) | no cycle counter | nanoseconds from `performance.now()` |

### Proofs made in WASM verify natively

A ballot proven in a browser has to verify on a native server. The
Fiat-Shamir challenges are expanded from the transcript hash by the AES-CTR
stream, so both builds must produce the same stream and read it the same way:

- the portable AES produces the same bytes as the AES-NI version for the
  same seed, including the counter layout and how a partial final block is
  handled;
- `evmlr_challenge.c` reads challenge coins as `uint64_t`, not `ulong`.
  `ulong` is 64 bits on x86-64 but 32 bits on wasm32, so the old code would
  have expanded the same hash into different challenges. On x86-64 this
  changes nothing.

Also for 32-bit `ulong`: `evmlr_commit.c` sums the commitment-opening norm in
64 bits. With a 32-bit sum, a large `r` could wrap it and pass `verify_norm`.

### The Gaussian sampler on WASM SIMD

Emscripten implements the x86 SSE–AVX2 intrinsics on top of WASM SIMD.
`_mm_round_pd` and `_mm_floor_pd` map to native opcodes. Others the sampler
uses, such as `_mm_comile_sd`, `_mm_cvtpd_epi32` and `_mm_cvtsd_si64`, are
emulated one lane at a time. So `gaussian.cpp` compiles unchanged, given two
fixes in `wasm/compat/x86intrin.h`:

- Emscripten ships no `<x86intrin.h>`, and clang's own one pulls in x86-only
  builtins. The compat header maps it to Emscripten's `<immintrin.h>`.
- VCL references `_mm_setcsr`, which WebAssembly cannot implement (it has no
  rounding-mode or denormal control). With SSE4.1 enabled, the sampler never
  calls it, so the compat header defines it to trap.

Because the same code runs on both platforms, the WASM sampler returns the
same samples as the native one, bit for bit. That was checked on 1.4 million
samples, σ from 1 to 23456.7, from the same seed. It is not faster than a
scalar port (FACCT works on one lane at a time), but it leaves a single
sampler to maintain.

The AES stream cannot follow this route. Emscripten has no AES-NI intrinsics,
and WebAssembly has no AES instruction to map them to. A software
`_mm_aesenc_si128` would let `fastrandombytes.c` compile unchanged, but each
call encrypts a single block. Measured, that halved the stream's throughput
and made the encryption proof 65% slower than `fastrandombytes_portable.c`,
which encrypts two blocks per pass.
