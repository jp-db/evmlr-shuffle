// Thin JavaScript wrapper over the WebAssembly build (wasm/build/evmlr.mjs).
//
//   import { loadEvmlr } from './evmlr.mjs';
//   const evmlr = await loadEvmlr();
//
//   // Authority
//   const election = evmlr.setup(3);                 // 3 mixing servers
//   const ppBytes = election.exportParams();         // Uint8Array, public
//   const skBytes = election.exportKeys();           // Uint8Array, secret
//
//   // Voter page: needs only the public parameters
//   const voterView = evmlr.importParams(ppBytes);
//   const ballotBytes = voterView.cast(1, [42]).toBytes();
//
//   // Counting
//   const counter = evmlr.importParams(ppBytes);
//   counter.importKeys(skBytes);
//   const tally = counter.count([counter.importBallot(ballotBytes)]); // [42]
//
// The bytes are the format of evmlr_serialize.h, the same as the native
// build reads and writes. Casting encrypts the vote in one layer per server
// and proves every layer; it takes seconds, so in a page run it in a Web
// Worker. Handles own WASM memory: free() them (Election.free() does not free
// its ballots).
import createEvmlrModule from './build/evmlr.mjs';

export class EvmlrError extends Error {
  constructor(message, status) {
    super(message);
    this.name = 'EvmlrError';
    this.status = status;
  }
}

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

  // Calls fn(lenPtr) -> dataPtr, copies the bytes out and frees them.
  function takeBytes(fn, what) {
    const lenPtr = M._malloc(4);
    try {
      const ptr = fn(lenPtr);
      if (!ptr) throw new EvmlrError(`could not export ${what}`);
      const len = M.HEAPU32[lenPtr >>> 2];
      const out = M.HEAPU8.slice(ptr, ptr + len);
      M._free(ptr);
      return out;
    } finally {
      M._free(lenPtr);
    }
  }

  // Copies bytes into the heap and calls fn(ptr, len, statusPtr).
  function giveBytes(bytes, fn) {
    if (!(bytes instanceof Uint8Array)) throw new TypeError('expected a Uint8Array');
    const ptr = M._malloc(Math.max(1, bytes.length));
    const statusPtr = M._malloc(4);
    try {
      M.HEAPU8.set(bytes, ptr);
      M.HEAP32[statusPtr >>> 2] = 0;
      const result = fn(ptr, bytes.length, statusPtr);
      return { result, status: M.HEAP32[statusPtr >>> 2] };
    } finally {
      M._free(statusPtr);
      M._free(ptr);
    }
  }

  function fail(what, status) {
    throw new EvmlrError(`${what}: ${M.UTF8ToString(M._evmlr_wasm_strerror(status))}`, status);
  }

  class Ballot {
    constructor(election, ptr) {
      this.election = election;
      this.ptr = ptr;
    }
    get voterId() {
      return M._evmlr_wasm_ballot_voter(this.ptr);
    }
    toBytes() {
      if (!this.ptr) throw new Error('ballot freed');
      return takeBytes((lenPtr) => M._evmlr_wasm_export_ballot(this.election.ptr, this.ptr, lenPtr), 'ballot');
    }
    free() {
      if (this.ptr && this.election.ptr) M._evmlr_wasm_ballot_free(this.election.ptr, this.ptr);
      this.ptr = 0;
    }
  }

  class Election {
    constructor(ptr) {
      this.ptr = ptr;
    }

    get hasKeys() {
      return !!M._evmlr_wasm_has_sk(this.ptr);
    }

    exportParams() {
      return takeBytes((lenPtr) => M._evmlr_wasm_export_pp(this.ptr, lenPtr), 'public parameters');
    }

    exportKeys() {
      if (!this.hasKeys) throw new EvmlrError('this election holds no secret keys');
      return takeBytes((lenPtr) => M._evmlr_wasm_export_sk(this.ptr, lenPtr), 'secret keys');
    }

    importKeys(bytes) {
      const { result } = giveBytes(bytes, (ptr, len) => M._evmlr_wasm_import_sk(this.ptr, ptr, len));
      if (result !== 0) fail('importKeys', result);
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
        if (!ptr) throw new EvmlrError('cast failed');
        return new Ballot(this, ptr);
      } finally {
        M._free(buf);
      }
    }

    importBallot(bytes) {
      const { result, status } = giveBytes(bytes, (ptr, len, st) => M._evmlr_wasm_import_ballot(this.ptr, ptr, len, st));
      if (!result) fail('importBallot', status);
      return new Ballot(this, result);
    }

    // Mixes and decrypts; the last ballot of each voter counts. Needs the
    // secret keys. Returns the constant term of each vote, in shuffled order.
    count(ballots) {
      if (!this.ptr) throw new Error('election freed');
      if (!this.hasKeys) throw new EvmlrError('count needs the secret keys (importKeys)');
      if (ballots.length === 0) return [];
      const list = u32(ballots.map((b) => {
        if (b.election !== this || !b.ptr) throw new Error('ballot from another election, or freed');
        return b.ptr;
      }));
      const out = M._malloc(ballots.length * 4);
      try {
        const counted = M._evmlr_wasm_count(this.ptr, list, ballots.length, out);
        if (counted < 0) throw new EvmlrError('count failed');
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
      if (!ptr) throw new EvmlrError('setup failed');
      return new Election(ptr);
    },
    importParams(bytes) {
      const { result, status } = giveBytes(bytes, (ptr, len, st) => M._evmlr_wasm_import_pp(ptr, len, st));
      if (!result) fail('importParams', status);
      return new Election(result);
    },
    module: M,
  };
}
