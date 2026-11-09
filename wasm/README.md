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
| `wasm/evmlr.mjs` | JS wrapper: `loadEvmlr()` → `setup(k)`, `election.cast(voter, vote)`, `election.count(ballots)` |
| `wasm/build/evmlr_*.js` | the test/benchmark binaries; `node wasm/build/evmlr_voting.js test` |

The module uses neither WASM SIMD nor threads, so it runs on every engine and
needs no cross-origin isolation. Casting takes a while; in a page, run it in a
Web Worker.

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
