# evmlr-shuffle


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


## TODO

- [ ] Optimize the implementation for performance (consider using the chinese remainder theorem).
- [ ] Add more detailed documentation and comments in the code.
- [ ] Separate the tests and benchmarks into their own files (they can already be run independently via `make test` / `make bench`, but still live behind `#ifdef MAIN` in the implementation files).
- [ ] Fix a bug where if a proof of shuffle is run multiple times it sometimes
      fails (this is likely due to a variable changing state unexpectedly).
      Not reproducible as of the memory-ownership fixes: 400 consecutive proofs
      (`./evmlr_shuffle.bin bench`, four runs of 100) completed without a
      failure, and the suites are clean under ASan, UBSan and valgrind. The
      ownership bugs fixed there -- outputs being re-initialised on top of a
      live allocation -- are a plausible cause, but this has not been confirmed,
      so the item stays open until someone reproduces it or agrees to close it.
- [ ] Implement the proof of correct decryption that each mixing server owes
      (knowledge of the secret key, the bound on the seed decryption error, and
      correct rounding). Without it the servers' decryption step is unverified,
      so the mixnet is not yet verifiable end to end.
- [x] Run the proof of shuffle in zero-knowledge, instead of only verifying the
      mathematical correctness. The shuffle now produces and verifies binary
      proofs for `D`, `P` and `W` plus a linear proof for `u`.
- [x] Turn the proof of shuffle into a non-interactive proof using Fiat-Shamir
      with Aborts. Challenges are derived by hashing the transcript with
      SHA-256, and the linear and binary proofs use rejection sampling.

**WARNING**: This is a prototype implementation for research purposes only. It is not optimized for performance or security.
