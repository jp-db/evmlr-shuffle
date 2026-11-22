# evmlr-shuffle

[![CI](https://github.com/jp-db/evmlr-shuffle/actions/workflows/ci.yml/badge.svg)](https://github.com/jp-db/evmlr-shuffle/actions/workflows/ci.yml)

Implementation of the Proof of Shuffle and primitives of the paper "Efficient Verifiable Mixnets from Lattices,
Revisited" by Jonathan Bootle, Vadim Lyubashevsky, and Antonio Merino-Gallardo (https://eprint.iacr.org/2025/658).

## Overview

This repository provides an implementation of the Proof of Shuffle protocol based on lattice-based cryptography. 
The implementation includes the necessary primitives and algorithms to create and verify shuffles of encrypted messages.

### Requirements

FLINT is required to run the code in this repository. 
We recommend installing FLINT via the package manager or directly from the [FLINT website](http://flintlib.org/).
The code has been tested with FLINT versions 3.3.1 and 3.6.0.

## Primitives

The implementation includes the following key primitives:

- **Commitment Scheme**: A lattice-based commitment scheme based on the instantiation by Ajtai.
- **MLPKE**: Module-Lattice Public Key Encryption scheme, this scheme also underlies Kyber's encryption.
- **OTSE**: A lattice-based One-Time Symmetric Encryption scheme.
- **HPKE**: A lattice-based Hybrid Public Key Encryption scheme, relies on **OTSE** and **MLPKE**.

## Running the Code

`make` compiles every binary; it no longer runs anything, so it is safe to use
as a quick syntax check. Parallel builds are supported: `make -j$(nproc)`.

| Command | What it does |
| --- | --- |
| `make` / `make build` | Compile every binary |
| `make test` | Run the test suites only — fast, this is the inner development loop |
| `make bench` | Run the benchmark suites only — slow |
| `make run` | Run tests and benchmarks for every binary (the old `make` behaviour) |
| `make <name>` | Build and run one binary, e.g. `make shuffle` or `make commit` |
| `make shared` | Build the shared library `libevmlr.so` (soname `libevmlr.so.0`), without the test/benchmark harness; link with `-levmlr -lflint -lgmp` |
| `make help` | List the available targets |
| `make clean` | Remove build artifacts |

The binaries themselves take an optional mode argument, so a single component
can be exercised directly:

```sh
./evmlr_shuffle.bin test    # tests only
./evmlr_shuffle.bin bench   # benchmarks only
./evmlr_shuffle.bin         # both (same as `all`)
```

Test failures are reported through the exit status, so `make test` fails loudly
and can be wired into CI.

### Continuous integration

`.github/workflows/ci.yml` builds and runs the test suites on every push and
pull request, in two configurations: the usual `-march=native` build and a
portable one (`ARCH=`) that checks the code still compiles without it.

The workflow pins `ubuntu-26.04` rather than `ubuntu-latest`, because this code
needs the FLINT >= 3.2 random-state API (`flint_rand_init`) and Ubuntu 24.04 —
which `ubuntu-latest` still resolves to — only packages FLINT 3.0.1.

Benchmarks are not run in CI; shared runners make the timings meaningless. To
treat compiler warnings as errors, add `EXTRA_CFLAGS=-Werror` to the build step.

Build flags can be overridden without editing the Makefile:

```sh
make ARCH=                      # portable build, without -march=native
make OPT=-O0                    # unoptimised, for debugging
make EXTRA_CFLAGS=-fsanitize=address
```

### WebAssembly

The voting protocol (encryption, proofs, shuffle) also builds to WebAssembly
with Emscripten, through a separate makefile that leaves the native build
untouched:

```sh
wasm/build-deps.sh          # once: GMP, MPFR and FLINT for wasm32
make -f Makefile.wasm test  # build, then run the test suites under Node
```

See [wasm/README.md](wasm/README.md) for the JS API and for what the WASM
build changes: the randomness source, and a portable AES in place of AES-NI.
The Gaussian sampler compiles unchanged on WASM SIMD, so the module needs an
engine with SIMD support (Chrome 91, Firefox 89, Safari 16.4, Node 16.4).

### Running an election

Each party runs only its own step, and only bytes pass between them:

1. every mixing server generates its own key pair and publishes the public
   key (`evmlr_voting_keygen`);
2. the public parameters are built from the election id and the k public
   keys (`evmlr_voting_pp_init`). The OTSE matrices and the shuffle's
   commitment key are expanded from hashes, so nobody chooses them;
3. voters cast ballots; the bulletin board accepts a ballot only if
   `evmlr_voting_check_ballot` holds. The proofs are bound to the election id
   and the voter id, so a ballot copied under another voter id, or replayed
   in another election, is refused;
4. server j = 0, ..., k − 1 in turn runs `evmlr_voting_mix_step` with its own
   key: it decrypts its layer, drops every ciphertext whose proof of
   encryption fails (publishing which and why), shuffles the rest and proves
   the shuffle. It hands a batch to server j + 1 and publishes a mix record;
5. anyone checks the k mix records against the board
   (`evmlr_voting_verify_counting`).

`evmlr_voting_setup` and `evmlr_voting_counting_transcript` still run all of
this in one process with every key, for tests and benchmarks only.

Voter authentication and eligibility are out of scope: they belong to the
application that uses the library. The library takes the voter id it is
given, binds it into the ballot's proofs and applies the revote policy to it
(the last valid ballot of each id counts), but it does not check who may
post under an id. The bulletin board must authenticate voters before it
accepts a ballot; anyone it lets post under an id can replace that voter's
ballot with a valid one of their own.

### Integration with an application

`evmlr_api.h` is the entry point for an external application: a C ABI on
byte strings, with no FLINT types in the header, one function per party and
step (`evmlr_keygen`, `evmlr_params`, `evmlr_cast`, `evmlr_check_ballot`,
`evmlr_first_batch`, `evmlr_mix`, `evmlr_mix_info`, `evmlr_transcript`,
`evmlr_verify`, `evmlr_results`). `make shared` builds it into `libevmlr.so`;
the WebAssembly module exports the same functions, and `wasm/evmlr.mjs` wraps
them for JavaScript. `make api` runs an election through it end to end.

`evmlr_serialize.h` defines the bytes: public parameters, a server's public
key and key pair, ballots, the batches servers hand on, mix records and
counting transcripts. They are platform-independent: a ballot cast in the
browser build verifies in the native one. The format is described in
[wasm/README.md](wasm/README.md#serialization-format).

To tweak the parameters of the schemes, mainly the number of messages, you can
modify `evmlr_params.h` in the root directory. Header dependencies are tracked,
so editing it rebuilds everything that depends on it.

## TODO

- [ ] **The inner-layer proofs of encryption deanonymise every ballot.** A
      ballot publishes a proof for each of its k layers, and each mix record
      publishes the server's outputs, which are the next layer's ciphertexts.
      The ML-PKE half of a ballot's layer-j proof verifies against exactly one
      of server j − 1's outputs, and its OTSE half against exactly one of
      server j's: anyone holding the board and the transcript can follow
      every ballot through every shuffle to its vote (checked: with 4 ballots
      and 2 servers, all 4 voters are linked to their votes). This predates
      the mix-step API, which carries the proofs the way counting always did.

      The paper does not have this problem, because it builds the proofs
      differently (Section 3.3 and the proof sizes in Section 5):
      - a voter gives *one* proof of knowledge that the outermost ciphertext
        is a correct layered encryption of some plaintext. The relations of
        every layer are inside that proof, and the inner ciphertexts are part
        of the witness, never public statements. Here each layer has its own
        proof, whose statement is that layer's ciphertext and, for the OTSE
        half, its plaintext (t = c − m);
      - the servers check all users' proofs before any decryption starts,
        and the proofs "can be deleted once the servers have verified them";
        nothing about an inner layer is checked, or carried, later;
      - each server publishes a proof of correct decryption (below), which
        together with the proof of shuffle is what makes the count
        verifiable. The voters' proofs only protect the servers' keys from
        malformed ciphertexts.

      So the fix is the paper's construction: a single proof over committed
      inner ciphertexts (the paper uses LaZer for it, 322 KB per voter with 4
      servers), checked in full by every server at the start, and inner
      proofs no longer in ballots, batches or rejections at inner servers.
      The paper's proof of correct encryption also proves the rounding
      d = ⌈H s̄⌋ of each layer, which the current one does not.
- [ ] Prove correct decryption (Section 3.3 of the paper): each server
      proves knowledge of its secret key, the decryption of every seed s̄_i,
      and the rounding of H s̄_i into d_i. The proof of shuffle shows that a
      server's outputs are a permutation of its inputs under the OTSE keys it
      committed to, but not that those are the keys its ML-PKE ciphertexts
      encrypt. For the same reason a rejection for a failed OTSE proof, which
      needs the decryption, cannot be checked by others (a rejection for the
      ML-PKE proof at server 0 can, and `evmlr_voting_verify_counting` does).
- [ ] Mix a single remaining ciphertext. The proof of shuffle needs N >= 2,
      so `evmlr_voting_mix_step` returns `EVMLR_VOTING_ERR_TOO_FEW` when one
      is left (which a shuffle could not hide anyway).
- [ ] Check all of a ballot before it goes on the board.
      `evmlr_voting_check_ballot` checks only the ML-PKE half of the
      outermost proof, which sits in the first ~2.3 KB of a ~100 KB ballot
      (k = 3). Flipping one bit almost anywhere else still passes it: about
      90% of single-bit changes in a sampled scan. The mix server whose
      layer was changed rejects the ciphertext later (HPKE or OTSE reason)
      and the vote is not counted, so integrity holds, but:
      - the board accepts ballots that will never count, and the voter is
        not told when casting;
      - since the last ballot of each voter counts, a corrupted revote
        replaces an earlier valid ballot, and the voter is left with no
        vote.

      The single proof of the first item, checked in full at casting time,
      would close this too.
- [ ] Verify a partial chain of mix records. `evmlr_verify` needs all k
      records and the board, so neither server j nor the ballot box can
      check servers 0..j − 1's proofs of shuffle before server j mixes.
      `evmlr_mix` only checks that its input batch is what the previous
      record says was output. A bad proof therefore shows up only in the
      final verification, after every server has mixed. A call that checks
      records 0..j − 1 against the board, without needing the rest, would
      let each server refuse a bad chain before it starts.
- [ ] Prove that a vote is well formed. A vote is any EVMLR_VOTE_LEN
      coefficients mod q, and nothing proves it encodes a valid answer (for
      an n-out-of-m question: coefficients in {0, 1} and at most n set).
      Applications must decode the votes after the count and leave out the
      malformed ones (the mixnet-voting ballot box does this and reports
      them). A proof of well-formedness, for n-out-of-m first and later for
      the polynomials ranked and transferable votes need, would let the
      board refuse such ballots instead.
- [ ] Optimize the implementation for performance. Multiplication in R_q now
      uses the CRT decomposition (`evmlr_crt.c`); the remaining headroom is in
      keeping operands in the split representation between operations instead
      of converting per multiplication.
- [ ] Fix the shift-exponent overflow in `gaussian.cpp` (the FACCT rejection
      check). When `exp(x)` is small enough that its biased exponent drops
      below 1003, `res_exponent` underflows and `1LL << res_exponent` shifts by
      2^64-1, which is undefined. UBSan reports it from the voting benchmark.
      Since a wrong rejection decision biases the sampler, this is worth
      treating as correctness-relevant rather than cosmetic. The sanitizer CI
      job reports but does not fail on it until this is resolved. The
      WebAssembly build compiles the same code with clang, which, like gcc,
      may assume the shift never overflows. Masking the count
      (`res_exponent & 63`, what x86 and WASM shifts do in hardware) fixes
      both builds.
- [ ] Free the keypairs allocated inside the benchmark timing loops (`bench()`
      in `evmlr_mlpke.c` and `evmlr_hpke.c` regenerate a keypair per iteration
      without clearing the previous one). Harness-only, but it is why leak
      detection is disabled for the benchmark phase in CI.
- [ ] Add more detailed documentation and comments in the code.
- [ ] Separate the tests and benchmarks into their own files (they can already be run independently via `make test` / `make bench`, but still live behind `#ifdef MAIN` in the implementation files).
- [x] Fix a bug where if a proof of shuffle is run multiple times it sometimes
      fails (this is likely due to a variable changing state unexpectedly).
      Not reproducible as of the memory-ownership fixes: 400 consecutive proofs
      (`./evmlr_shuffle.bin bench`, four runs of 100) completed without a
      failure, and the suites are clean under ASan, UBSan and valgrind. The
      ownership bugs fixed there -- outputs being re-initialised on top of a
      live allocation -- are a plausible cause, but this has not been confirmed,
      so the item stays open until someone reproduces it or agrees to close it. ([0cf34706](https://github.com/jp-db/evmlr-shuffle/commit/0cf34706b4b4cbb5a4c190d758a5be8de6653183))
- [ ] Implement the proof of correct decryption that each mixing server owes
      (knowledge of the secret key, the bound on the seed decryption error, and
      correct rounding). Without it the servers' decryption step is unverified,
      so the mixnet is not yet verifiable end to end.
- [x] Run the proof of shuffle in zero-knowledge, instead of only verifying the
      mathematical correctness. The shuffle now produces and verifies binary
      proofs for `D`, `P` and `W` plus a linear proof for `u`. ([ffee11bc](https://github.com/jp-db/evmlr-shuffle/commit/ffee11bc716529a62869a0d8e57ec2b44d0776d0))
- [x] Turn the proof of shuffle into a non-interactive proof using Fiat-Shamir
      with Aborts. Challenges are derived by hashing the transcript with
      SHA-256, and the linear and binary proofs use rejection sampling. ([99cbb089](https://github.com/jp-db/evmlr-shuffle/commit/99cbb0893fa74ece81fe7a9e353762b66d34cae0))

**WARNING**: This is a prototype implementation for research purposes only. It is not optimized for performance or security.
