# WebAssembly build

The voting protocol (setup, ballot encryption with its proofs of correct
encryption, and counting with the proofs of shuffle) compiles to WebAssembly
with Emscripten. The native build (`make`) is unchanged; the WASM build is a
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
| `wasm/build/evmlr.mjs` + `evmlr.wasm` | ES module factory (`createEvmlrModule`), for browsers, workers and Node |
| `wasm/evmlr.mjs` | JS wrapper, see below |
| `wasm/demo.mjs` | authority, voters, counter and auditor exchanging only bytes (`make -f Makefile.wasm demo`) |
| `wasm/build/evmlr_*.js` | the test/benchmark binaries; `node wasm/build/evmlr_voting.js test` |

The module uses neither WASM SIMD nor threads, so it runs on every engine and
needs no cross-origin isolation. Casting takes a while; in a page, run it in a
Web Worker.

## JavaScript API

```js
import { loadEvmlr } from './wasm/evmlr.mjs';
const evmlr = await loadEvmlr();

// Election authority
const authority = evmlr.setup(3);            // 3 mixing servers
const ppBytes = authority.exportParams();    // public: publish to voters
const skBytes = authority.exportKeys();      // secret: to the servers

// Voter page: the public parameters are all it needs
const election = evmlr.importParams(ppBytes);
const ballot = election.cast(voterId, [candidate]);
const ballotBytes = ballot.toBytes();        // send to the bulletin board

// Counting: the tally, and a transcript of every server's proof of shuffle
const counter = evmlr.importParams(ppBytes);
counter.importKeys(skBytes);
const { tally, transcript } =
  counter.countWithTranscript(ballots.map((b) => counter.importBallot(b)));

// Anyone: check the count from public data (parameters, ballots, transcript)
const auditor = evmlr.importParams(ppBytes);
const verified = auditor.verifyCounting(ballots.map((b) => auditor.importBallot(b)), transcript);
```

All bytes are `Uint8Array`s in the format below, so they move freely between
the WASM module and the native code. Import errors throw `EvmlrError`, whose
`status` is an `evmlr_ser_status_t`. Handles hold WASM memory; call `free()`
on elections and ballots when you are done with them.

## Serialization format

`evmlr_serialize.h` (also used natively) encodes the public parameters, the
secret keys, ballots with their proofs of correct encryption, proofs of
shuffle with their statements, and counting transcripts. Every encoding has:

- a 20-byte header: magic, format version, object type, and a fingerprint of
  the scheme parameters (`MOD_Q`, `DEGREE_N`, `K_LWE`, ...), so bytes from a
  build with different parameters are refused;
- fixed-width little-endian integers, the same on wasm32 and x86-64;
- polynomial coefficients bit-packed at 12 bits mod q (6 bits mod 2^zeta).

Decoding is strict: every dimension must be the one the protocol fixes,
coefficients must be reduced, and polynomials must be normalised. So every
object has exactly one encoding, and what is read back hashes into the
Fiat-Shamir transcripts like the original. A truncated, padded or tampered
ballot is refused before any proof is checked.

The public parameters carry only what is sampled: `H`, the random block
`H''` of `H' = (H'' | I)` and the servers' public keys. Everything else is
rebuilt on import, including the shuffle's commitment key (see below).

Sizes with the default parameters, for 3 servers: 28 KB of public
parameters, 100 KB per ballot, and a transcript of about 2.1 MB for 5 counted
ballots. The transcript grows linearly with the number of ballots.

### Verifiable counting

`evmlr_voting_counting_transcript` counts like `evmlr_voting_counting` and
also returns, for every server, the statement its shuffle proves (the
commitment `D`, the outputs `c_hat`, the inputs `c_star`) and the proof.
`evmlr_voting_verify_counting` checks such a transcript with public data
only:

- applying the revote policy to the ballots gives the number of ballots the
  transcript counted;
- server 0's inputs are the OTSE ciphertexts of those ballots' outer layers,
  in order;
- server j's inputs, for j > 0, are the ciphertexts server j − 1 output;
- every proof of shuffle verifies.

It returns the tally the last server's outputs hold. It does not show that
the servers decrypted correctly; that proof is not implemented yet (see the
top-level README).

The commitment key of the shuffle is expanded from a SHA-256 hash of the
public parameters and the number of ballots, with the AES-256-CTR stream.
Before, the counting servers sampled it, and they are the provers of the
shuffle. A prover who chooses the commitment key can choose one with a
trapdoor, which breaks binding and so the soundness of the proof. Now
nobody chooses the key, and every verifier recomputes the same one.

## What changes for WebAssembly

`evmlr_portable.h` defines `EVMLR_PORTABLE` for every Emscripten build. Without
Emscripten it is never set, so native builds compile the same code as before.

| Native | Problem in WASM | WASM replacement |
| --- | --- | --- |
| `getrandom(2)` (`evmlr_utils.c`, `evmlr_main.h`), `/dev/urandom` (`randombytes.c`) | no system calls | `getentropy()`, which Emscripten implements with `crypto.getRandomValues` (browser) or `crypto.randomFillSync` (Node), in chunks of 256 bytes |
| AES-256-CTR with AES-NI (`fastrandombytes.c`) | no AES-NI | `fastrandombytes_portable.c`: bitsliced, constant-time AES (BearSSL's `aes_ct` design) |
| FACCT Gaussian with SSE4.1 + VCL (`gaussian.cpp`) | [Emscripten has no SSE4.1 `_mm_round_pd`](https://emscripten.org/docs/porting/simd.html), and VCL needs x86 | `gaussian_portable.c`: the same FACCT algorithm in scalar C |
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

### Why not malb/dgs

[dgs](https://github.com/malb/dgs) is portable, but it would replace the
sampler with a different algorithm and its own randomness source, so it
would need patching to draw from `fastrandombytes()`. The scalar port keeps
the sampler the proofs were analysed with. For the same random stream it
returns the same samples as `gaussian.cpp`, except where libm and VCL round
`log`/`sin`/`cos` differently in the last bit and that flips a decision.
