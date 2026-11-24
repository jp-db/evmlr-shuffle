// JavaScript wrapper over the WebAssembly build of evmlr_api.h.
//
//   import { loadEvmlr } from './evmlr.mjs';
//   const evmlr = await loadEvmlr();
//
//   // Each mixing server, on its own machine
//   const { publicKey, keyPair } = evmlr.keygen();      // publish publicKey, keep keyPair
//
//   // Election authority: public keys in, public parameters out
//   const id = await evmlr.electionId('Board election 2026');
//   const pp = evmlr.params(id, [pk0, pk1, pk2]);
//
//   // Voter page
//   const ballot = evmlr.cast(pp, voterId, candidate);
//
//   // Bulletin board
//   evmlr.checkBallot(pp, ballot);                       // throws EvmlrError if invalid
//
//   // Servers, in turn
//   let batch = evmlr.firstBatch(pp, board);             // server 0's input
//   const mixes = [];
//   for (let j = 0; j < k; j++) {                        // server j, with its own key only
//     evmlr.verifyMixes(pp, board, mixes);               // servers 0..j-1 checked first
//     const step = evmlr.mix(pp, keyPairOfServer[j], batch, mixes[j - 1] ?? null);
//     mixes.push(step.mix);                              // publish
//     batch = step.batch;                                // hand to server j + 1
//   }
//
//   // Anyone
//   const votes = evmlr.verify(pp, board, evmlr.transcript(pp, mixes));
//
// Every value is a Uint8Array in the format of evmlr_serialize.h, the same
// bytes the native library reads and writes. A vote is a number (the first
// coefficient) or an array of up to evmlr.voteLen coefficients mod
// evmlr.modulus; votes come back as Uint16Arrays of voteLen coefficients.
// Voter ids are 64-bit signed integers, as in the native build: pass a
// BigInt or a safe-integer Number; they come back as BigInts.
// Casting and mixing take seconds: in a page, run them in a Web Worker.
import createEvmlrModule from './build/evmlr.mjs';

export class EvmlrError extends Error {
  constructor(message, status) {
    super(message);
    this.name = 'EvmlrError';
    this.status = status;
  }
}

/** Reasons a server gives for rejecting a ciphertext (evmlr_voting_reject_t). */
export const REJECT_REASONS = { 1: 'hpke-proof', 2: 'otse-proof' };

export async function loadEvmlr(moduleArgs = {}) {
  const M = await createEvmlrModule(moduleArgs);
  const modulus = M._evmlr_modulus();
  const voteLen = M._evmlr_vote_len();
  const ERR_INVALID = 12;
  const ERR_SPACE = 14;

  // Runs fn with an allocator whose blocks are all freed afterwards. The heap
  // can grow during a call, so M.HEAPU8 and friends are re-read every time.
  function scoped(fn) {
    const blocks = [];
    const alloc = (size) => {
      const p = M._malloc(Math.max(1, size));
      if (!p) throw new EvmlrError('out of memory', 9);
      blocks.push(p);
      return p;
    };
    try {
      return fn(alloc);
    } finally {
      for (const p of blocks) M._free(p);
    }
  }

  function bytesIn(alloc, bytes) {
    if (!(bytes instanceof Uint8Array)) throw new TypeError('expected a Uint8Array');
    const p = alloc(bytes.length);
    M.HEAPU8.set(bytes, p);
    return p;
  }

  // An array of byte strings: [pointer array, length array].
  function listIn(alloc, list) {
    const ptrs = alloc(4 * list.length);
    const lens = alloc(4 * list.length);
    list.forEach((b, i) => {
      const p = bytesIn(alloc, b);
      M.HEAPU32[(ptrs >>> 2) + i] = p;
      M.HEAPU32[(lens >>> 2) + i] = b.length;
    });
    return [ptrs, lens];
  }

  // An evmlr_bytes { uint8_t *data; size_t len; } (8 bytes on wasm32).
  function bytesOut(alloc) {
    const p = alloc(8);
    M.HEAPU32[p >>> 2] = 0;
    M.HEAPU32[(p >>> 2) + 1] = 0;
    return p;
  }

  function takeBytes(p) {
    const data = M.HEAPU32[p >>> 2];
    const len = M.HEAPU32[(p >>> 2) + 1];
    const out = M.HEAPU8.slice(data, data + len);
    M._evmlr_bytes_free(p);
    return out;
  }

  function check(what, status) {
    if (status !== 0) {
      throw new EvmlrError(`${what}: ${M.UTF8ToString(M._evmlr_strerror(status))}`, status);
    }
  }

  function voteIn(alloc, vote) {
    const coeffs = typeof vote === 'number' ? [vote] : Array.from(vote);
    if (coeffs.length > voteLen) throw new RangeError(`a vote has at most ${voteLen} coefficients`);
    for (const c of coeffs) {
      if (!Number.isInteger(c) || c < 0 || c >= modulus) throw new RangeError(`vote coefficients are integers in [0, ${modulus})`);
    }
    const p = alloc(2 * coeffs.length);
    M.HEAPU16.set(coeffs, p >>> 1);
    return [p, coeffs.length];
  }

  function votesOut(p, count) {
    const votes = [];
    for (let i = 0; i < count; i++) {
      votes.push(M.HEAPU16.slice((p >>> 1) + i * voteLen, (p >>> 1) + (i + 1) * voteLen));
    }
    return votes;
  }

  const readInt64 = (p) => (BigInt(M.HEAP32[(p >>> 2) + 1]) << 32n) | BigInt(M.HEAPU32[p >>> 2]);

  function voterIdIn(id) {
    if (typeof id === 'number' && !Number.isSafeInteger(id)) {
      throw new RangeError('a voter id Number must be a safe integer; pass a BigInt for larger ids');
    }
    const big = BigInt(id);
    if (big !== BigInt.asIntN(64, big)) throw new RangeError('a voter id is a 64-bit signed integer');
    return big;
  }

  return {
    modulus,
    voteLen,
    module: M,

    /** A 32-byte election id: SHA-256 of a description (string or bytes). */
    async electionId(description) {
      const data = typeof description === 'string' ? new TextEncoder().encode(description) : description;
      return new Uint8Array(await globalThis.crypto.subtle.digest('SHA-256', data));
    },

    /** A server's key pair: { publicKey (to publish), keyPair (to keep secret) }. */
    keygen() {
      return scoped((alloc) => {
        const pk = bytesOut(alloc), kp = bytesOut(alloc);
        check('keygen', M._evmlr_keygen(pk, kp));
        return { publicKey: takeBytes(pk), keyPair: takeBytes(kp) };
      });
    },

    /** Public parameters from the election id and the servers' public keys, in mixing order. */
    params(electionId, publicKeys) {
      if (!(electionId instanceof Uint8Array) || electionId.length !== 32) {
        throw new TypeError('the election id is 32 bytes (see electionId())');
      }
      return scoped((alloc) => {
        const out = bytesOut(alloc);
        const [ptrs, lens] = listIn(alloc, publicKeys);
        check('params', M._evmlr_params(out, bytesIn(alloc, electionId), publicKeys.length, ptrs, lens));
        return takeBytes(out);
      });
    },

    /** { servers, electionId } of public parameters. */
    paramsInfo(pp) {
      return scoped((alloc) => {
        const k = alloc(4), id = alloc(32);
        check('paramsInfo', M._evmlr_params_info(bytesIn(alloc, pp), pp.length, k, id));
        return { servers: M.HEAPU32[k >>> 2], electionId: M.HEAPU8.slice(id, id + 32) };
      });
    },

    /**
     * Encrypts and proves a vote for voterId (a BigInt, or a safe-integer
     * Number), and checks every layer before returning. Returns the ballot.
     */
    cast(pp, voterId, vote) {
      return scoped((alloc) => {
        const out = bytesOut(alloc);
        const [v, n] = voteIn(alloc, vote);
        check('cast', M._evmlr_cast(out, bytesIn(alloc, pp), pp.length, voterIdIn(voterId), v, n));
        return takeBytes(out);
      });
    },

    /**
     * The bulletin board's check: the ballot is for this election, unchanged
     * since it was cast, and bound to its voter id. Returns the voter id (a
     * BigInt); throws EvmlrError if invalid.
     */
    checkBallot(pp, ballot) {
      return scoped((alloc) => {
        const id = alloc(8);
        check('checkBallot', M._evmlr_check_ballot(bytesIn(alloc, pp), pp.length, bytesIn(alloc, ballot), ballot.length, id));
        return readInt64(id);
      });
    },

    /** checkBallot as a boolean, for filtering. */
    isValidBallot(pp, ballot) {
      try {
        this.checkBallot(pp, ballot);
        return true;
      } catch (e) {
        if (e instanceof EvmlrError && e.status === ERR_INVALID) return false;
        throw e;
      }
    },

    /** Server 0's input batch, from the ballots on the board (in board order). */
    firstBatch(pp, ballots) {
      return scoped((alloc) => {
        const out = bytesOut(alloc);
        const [ptrs, lens] = listIn(alloc, ballots);
        check('firstBatch', M._evmlr_first_batch(out, bytesIn(alloc, pp), pp.length, ptrs, lens, ballots.length));
        return takeBytes(out);
      });
    },

    /**
     * One server's step, with its own key pair. prevMix is the previous
     * server's mix record (null for server 0). Returns { batch, mix }: the
     * batch for the next server, and the record to publish.
     */
    mix(pp, keyPair, batch, prevMix = null) {
      return scoped((alloc) => {
        const outBatch = bytesOut(alloc), outMix = bytesOut(alloc);
        const prev = prevMix ? bytesIn(alloc, prevMix) : 0;
        check('mix', M._evmlr_mix(outBatch, outMix, bytesIn(alloc, pp), pp.length, bytesIn(alloc, keyPair), keyPair.length,
                                  bytesIn(alloc, batch), batch.length, prev, prevMix ? prevMix.length : 0));
        return { batch: takeBytes(outBatch), mix: takeBytes(outMix) };
      });
    },

    /** What a mix record says: { server, received, passedOn, rejected: [{ index, reason }] }. */
    mixInfo(pp, mix) {
      return scoped((alloc) => {
        const ppPtr = bytesIn(alloc, pp), mixPtr = bytesIn(alloc, mix);
        const f = alloc(16); // server, n_in, n_out, n_rejected
        const call = (cap, idx, why) =>
          M._evmlr_mix_info(ppPtr, pp.length, mixPtr, mix.length, f, f + 4, f + 8, f + 12, idx, why, cap);
        let st = call(0, 0, 0);
        const nRej = M.HEAPU32[(f >>> 2) + 3];
        let idx = 0, why = 0;
        if (st === ERR_SPACE) {
          idx = alloc(4 * nRej);
          why = alloc(nRej);
          st = call(nRej, idx, why);
        }
        check('mixInfo', st);
        const rejected = [];
        for (let r = 0; r < nRej; r++) {
          const code = M.HEAPU8[why + r];
          rejected.push({ index: M.HEAPU32[(idx >>> 2) + r], reason: REJECT_REASONS[code] ?? code });
        }
        return { server: M.HEAPU32[f >>> 2], received: M.HEAPU32[(f >>> 2) + 1], passedOn: M.HEAPU32[(f >>> 2) + 2], rejected };
      });
    },

    /**
     * Checks the first mixes.length records (server 0's first) against the
     * ballots on the board, as verify would. Server j calls it with the j
     * records before its own, before mixing. Returns how many ciphertexts
     * the last of them passed on (with no records: how many ballots count);
     * throws EvmlrError (status 12) if any fails.
     */
    verifyMixes(pp, ballots, mixes) {
      return scoped((alloc) => {
        const n = alloc(4);
        const [bp, bl] = listIn(alloc, ballots);
        const [mp, ml] = listIn(alloc, mixes);
        check('verifyMixes', M._evmlr_verify_mixes(n, bytesIn(alloc, pp), pp.length, bp, bl, ballots.length,
                                                   mp, ml, mixes.length));
        return M.HEAPU32[n >>> 2];
      });
    },

    /** The k mix records, in server order, as one transcript. */
    transcript(pp, mixes) {
      return scoped((alloc) => {
        const out = bytesOut(alloc);
        const [ptrs, lens] = listIn(alloc, mixes);
        check('transcript', M._evmlr_transcript(out, bytesIn(alloc, pp), pp.length, ptrs, lens, mixes.length));
        return takeBytes(out);
      });
    },

    /**
     * Checks a count from public data: every ballot on the board, in board
     * order, and the transcript. Returns the votes; throws EvmlrError (status
     * 12) if the count does not verify.
     */
    verify(pp, ballots, transcript) {
      return scoped((alloc) => {
        const cap = Math.max(1, ballots.length);
        const votes = alloc(2 * voteLen * cap), n = alloc(4);
        const [ptrs, lens] = listIn(alloc, ballots);
        check('verify', M._evmlr_verify(votes, cap, n, bytesIn(alloc, pp), pp.length, ptrs, lens, ballots.length,
                                        bytesIn(alloc, transcript), transcript.length));
        return votesOut(votes, M.HEAPU32[n >>> 2]);
      });
    },

    /** The votes in the last server's output batch. */
    results(pp, finalBatch) {
      return scoped((alloc) => {
        const ppPtr = bytesIn(alloc, pp), b = bytesIn(alloc, finalBatch), n = alloc(4);
        let st = M._evmlr_results(0, 0, n, ppPtr, pp.length, b, finalBatch.length);
        const count = M.HEAPU32[n >>> 2];
        if (st === ERR_SPACE) {
          const votes = alloc(2 * voteLen * count);
          check('results', M._evmlr_results(votes, count, n, ppPtr, pp.length, b, finalBatch.length));
          return votesOut(votes, count);
        }
        check('results', st);
        return [];
      });
    },
  };
}
