// Thin JavaScript wrapper over the WebAssembly build (wasm/build/evmlr.mjs).
//
//   import { loadEvmlr } from './evmlr.mjs';
//   const evmlr = await loadEvmlr();
//   const election = evmlr.setup(3);            // 3 mixing servers
//   const ballot = election.cast(1, [42]);      // voter 1 votes for 42
//   const tally = election.count([ballot]);     // [42]
//   ballot.free(); election.free();
//
// Casting encrypts the vote in one layer per server and proves every layer;
// it takes seconds, so in a page run it in a Web Worker. Handles own WASM
// memory and must be freed.
import createEvmlrModule from './build/evmlr.mjs';

export async function loadEvmlr(moduleArgs = {}) {
  const M = await createEvmlrModule(moduleArgs);
  const params = {
    q: M._evmlr_wasm_param(0),
    n: M._evmlr_wasm_param(1),
    kLWE: M._evmlr_wasm_param(2),
    kLWR: M._evmlr_wasm_param(3),
  };

  // Copies a JS array of integers into a fresh uint32 buffer on the heap.
  function u32(values) {
    const ptr = M._malloc(Math.max(1, values.length) * 4);
    M.HEAPU32.set(values.map((v) => v >>> 0), ptr >>> 2);
    return ptr;
  }

  class Ballot {
    constructor(election, ptr, voterId) {
      this.election = election;
      this.ptr = ptr;
      this.voterId = voterId;
    }
    free() {
      if (this.ptr) M._evmlr_wasm_ballot_free(this.election.ptr, this.ptr);
      this.ptr = 0;
    }
  }

  class Election {
    constructor(ptr, servers) {
      this.ptr = ptr;
      this.servers = servers;
    }

    // vote: coefficients of the vote polynomial in R_q (only the constant
    // term is reported back by count()).
    cast(voterId, vote) {
      if (!this.ptr) throw new Error('election freed');
      const coeffs = Array.isArray(vote) ? vote : [vote];
      if (coeffs.length > params.n) throw new RangeError(`vote has more than ${params.n} coefficients`);
      const buf = u32(coeffs);
      try {
        const ptr = M._evmlr_wasm_cast(this.ptr, voterId, buf, coeffs.length);
        if (!ptr) throw new Error('evmlr_wasm_cast failed');
        return new Ballot(this, ptr, voterId);
      } finally {
        M._free(buf);
      }
    }

    // Mixes and decrypts; the last ballot of each voter counts. Returns the
    // constant term of each vote, in shuffled order.
    count(ballots) {
      if (!this.ptr) throw new Error('election freed');
      if (ballots.length === 0) return [];
      const list = u32(ballots.map((b) => {
        if (b.election !== this || !b.ptr) throw new Error('ballot from another election, or freed');
        return b.ptr;
      }));
      const out = M._malloc(ballots.length * 4);
      try {
        const counted = M._evmlr_wasm_count(this.ptr, list, ballots.length, out);
        return Array.from(M.HEAPU32.subarray(out >>> 2, (out >>> 2) + counted));
      } finally {
        M._free(out);
        M._free(list);
      }
    }

    free() {
      if (this.ptr) M._evmlr_wasm_election_free(this.ptr);
      this.ptr = 0;
    }
  }

  return {
    params,
    setup(servers) {
      const ptr = M._evmlr_wasm_setup(servers);
      if (!ptr) throw new Error('evmlr_wasm_setup failed');
      return new Election(ptr, servers);
    },
    module: M,
  };
}
