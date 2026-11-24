// An election under Node, through the WebAssembly library, with every party
// separate and only bytes between them:
//   make -f Makefile.wasm lib && node wasm/demo.mjs
//
//   - three mixing servers each generate their own key pair;
//   - the authority builds the public parameters from the public keys only;
//   - voters cast (voter 2 votes twice; the second vote counts);
//   - the bulletin board checks every ballot, and refuses one copied under
//     another voter id or changed anywhere;
//   - the servers mix one after the other, each with its own key, each
//     checking the records before its own and that its input is what the
//     previous server published;
//   - an auditor checks the count from public data.
import { loadEvmlr, EvmlrError } from './evmlr.mjs';

const evmlr = await loadEvmlr();
const ms = (t) => `${(performance.now() - t).toFixed(0)} ms`;
const expectError = (what, fn) => {
  try {
    fn();
  } catch (e) {
    if (!(e instanceof EvmlrError)) throw e;
    console.log(`${what} refused: ${e.message}`);
    return true;
  }
  console.log(`MISMATCH: ${what} accepted`);
  return false;
};
let ok = true;
console.log(`q = ${evmlr.modulus}, ${evmlr.voteLen} coefficients per vote`);

// Servers: each on its own, keeping its key pair.
const K = 3;
let t = performance.now();
const servers = Array.from({ length: K }, () => evmlr.keygen());
console.log(`keygen x${K}: ${ms(t)}; public key ${servers[0].publicKey.length} B, key pair ${servers[0].keyPair.length} B`);

// Authority: public keys in, public parameters out.
const electionId = await evmlr.electionId('evmlr demo election');
const pp = evmlr.params(electionId, servers.map((s) => s.publicKey));
console.log(`parameters: ${pp.length} B, ${evmlr.paramsInfo(pp).servers} servers`);

// Voters, and the bulletin board.
const board = [];
// Voter 5's id needs more than 32 bits: ids are 64-bit in every build.
for (const [voter, vote] of [[1, 10], [2, 99], [3, 30], [2, 20], [4, 40], [2n ** 40n + 5n, 50]]) {
  t = performance.now();
  const ballot = evmlr.cast(pp, voter, vote);
  const who = evmlr.checkBallot(pp, ballot);
  ok &&= who === BigInt(voter);
  board.push(ballot);
  console.log(`voter ${voter} casts ${vote}: ${ms(t)}, ${ballot.length} B, accepted`);
}

// Someone resubmits voter 1's ballot under voter 6's id (the voter id is the
// u64 right after the 20-byte header).
const copy = board[0].slice();
copy[20] = 6;
ok &&= expectError('ballot copied under voter 6', () => evmlr.checkBallot(pp, copy));
ok &&= expectError('truncated ballot', () => evmlr.checkBallot(pp, board[0].slice(0, -1)));
// One bit changed deep inside the ballot, in a proof of an inner layer.
const flipped = board[0].slice();
flipped[flipped.length - 1000] ^= 1;
ok &&= expectError('ballot changed in an inner proof', () => evmlr.checkBallot(pp, flipped));

// Mixing: server j takes the batch server j - 1 handed on, with its record.
let batch = evmlr.firstBatch(pp, board);
const mixes = [];
for (let j = 0; j < K; j++) {
  t = performance.now();
  ok &&= evmlr.verifyMixes(pp, board, mixes) === 5; // servers 0..j-1, before mixing
  const step = evmlr.mix(pp, servers[j].keyPair, batch, j > 0 ? mixes[j - 1] : null);
  const info = evmlr.mixInfo(pp, step.mix);
  console.log(`server ${j}: ${ms(t)}; ${info.received} in, ${info.passedOn} out, rejected ${JSON.stringify(info.rejected)}; record ${step.mix.length} B`);
  mixes.push(step.mix);
  batch = step.batch;
}
ok &&= expectError('mixing with another server\'s key', () => evmlr.mix(pp, servers[1].keyPair, evmlr.firstBatch(pp, board)));
const tally = evmlr.results(pp, batch).map((v) => v[0]);
console.log('tally (shuffled):', tally);

// Auditor: parameters, board and the records; no keys.
const transcript = evmlr.transcript(pp, mixes);
t = performance.now();
const verified = evmlr.verify(pp, board, transcript).map((v) => v[0]);
console.log(`audit: ${ms(t)}; transcript ${transcript.length} B; verified tally:`, verified);

const tampered = transcript.slice();
tampered[tampered.length - 100] ^= 1;
ok &&= expectError('tampered transcript', () => evmlr.verify(pp, board, tampered));
ok &&= expectError('count against a board missing a ballot', () => evmlr.verify(pp, board.slice(1), transcript));
const badMix = mixes[0].slice();
badMix[badMix.length - 100] ^= 1;
ok &&= expectError('a bad record of server 0, before server 1 mixes', () => evmlr.verifyMixes(pp, board, [badMix]));

// A one-ballot election is mixed and counted too.
let single = evmlr.firstBatch(pp, board.slice(0, 1));
const singleMixes = [];
for (let j = 0; j < K; j++) {
  const step = evmlr.mix(pp, servers[j].keyPair, single, j > 0 ? singleMixes[j - 1] : null);
  singleMixes.push(step.mix);
  single = step.batch;
}
const singleVotes = evmlr.verify(pp, board.slice(0, 1), evmlr.transcript(pp, singleMixes)).map((v) => v[0]);
console.log('one-ballot election, verified tally:', singleVotes);
ok &&= JSON.stringify(singleVotes) === '[10]';

const expected = [10, 20, 30, 40, 50];
const sorted = (xs) => JSON.stringify([...xs].sort((a, b) => a - b));
ok &&= sorted(tally) === JSON.stringify(expected) && JSON.stringify(verified) === JSON.stringify(tally);
console.log(ok ? 'OK' : `MISMATCH, expected ${expected}`);
process.exit(ok ? 0 : 1);
