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

To tweak the parameters of the schemes, mainly the number of messages, you can
modify `evmlr_params.h` in the root directory. Header dependencies are tracked,
so editing it rebuilds everything that depends on it.

## TODO

- [ ] Optimize the implementation for performance (consider using the chinese remainder theorem).
- [ ] Fix the shift-exponent overflow in `gaussian.cpp` (the FACCT rejection
      check). When `exp(x)` is small enough that its biased exponent drops
      below 1003, `res_exponent` underflows and `1LL << res_exponent` shifts by
      2^64-1, which is undefined. UBSan reports it from the voting benchmark.
      Since a wrong rejection decision biases the sampler, this is worth
      treating as correctness-relevant rather than cosmetic. The sanitizer CI
      job reports but does not fail on it until this is resolved.
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
