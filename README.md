# evmlr-shuffle

Implementation of the Proof of Shuffle and primitives of the paper "Efficient Verifiable Mixnets from Lattices,
Revisited" by Jonathan Bootle, Vadim Lyubashevsky, and Antonio Merino-Gallardo (https://eprint.iacr.org/2025/658).

## Overview

This repository provides an implementation of the Proof of Shuffle protocol based on lattice-based cryptography. 
The implementation includes the necessary primitives and algorithms to create and verify shuffles of encrypted messages.

### Requirements

FLINT is required to run the code in this repository. 
We recommend installing FLINT via the package manager or directly from the [FLINT website](http://flintlib.org/).
The code has been tested with FLINT version 3.3.1.

## Primitives

The implementation includes the following key primitives:

- **Commitment Scheme**: A lattice-based commitment scheme based on the instantiation by Ajtai.
- **MLPKE**: Module-Lattice Public Key Encryption scheme, this scheme also underlies Kyber's encryption.
- **OTSE**: A lattice-based One-Time Symmetric Encryption scheme.
- **HPKE**: A lattice-based Hybrid Public Key Encryption scheme, relies on **OTSE** and **MLPKE**.

## Running the Code

To compile and execute the code just run `make` in the root directory. 
This will compile the source files and run the executables for the primitives and the Proof of Shuffle.

To run the tests of a single file, you can execute `make <primitive lowercase>` or `make shuffle`.

To tweak the parameters of the schemes, mainly the number of messages, you can modify the `evmlr_params.h` in the root directory.

## TODO

- [ ] Optimize the implementation for performance (consider using the chinese remainder theorem).
- [ ] Add more detailed documentation and comments in the code.
- [ ] Separate the tests and benchmarks from the main implementation.
- [ ] Fix a bug where if a proof of shuffle is run multiple times it sometimes fails (this is likely due to a variable changing state unexpectedly).
- [ ] Run in the proof of shuffle in zero-knowledge, instead of only verifying the mathematical correctness.
- [ ] Turn the proof of shuffle into a non-interactive proof using the Fiat-Shamir with Aborts.

**WARNING**: This is a prototype implementation for research purposes only. It is not optimized for performance or security.
