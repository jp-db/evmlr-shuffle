// End-to-end check of the WebAssembly library under Node:
//   make -f Makefile.wasm lib && node wasm/demo.mjs
//
// Separate handles stand for the parties, and only bytes pass between them:
// the authority runs setup, the voters cast from the public parameters
// alone, the counter imports parameters, keys and ballots and publishes a
// transcript, and an auditor checks the count from public data only.
// Mirrors test_voting_protocol in evmlr_voting.c: voter 2 votes twice and
// only the second vote counts.
import { loadEvmlr, EvmlrError } from './evmlr.mjs';

const evmlr = await loadEvmlr();
const ms = (t) => `${(performance.now() - t).toFixed(0)} ms`;
console.log('params', evmlr.params);

// Authority
let t = performance.now();
const authority = evmlr.setup(3);
const ppBytes = authority.exportParams();
const skBytes = authority.exportKeys();
authority.free();
console.log(`setup (3 servers): ${ms(t)}; public parameters ${ppBytes.length} B, keys ${skBytes.length} B`);

// Voters: public parameters only
const voterView = evmlr.importParams(ppBytes);
console.log('voter page holds keys?', voterView.hasKeys);
const casts = [[1, 10], [2, 99], [3, 30], [2, 20], [4, 40], [5, 50]];
const ballotBytes = [];
for (const [voter, vote] of casts) {
  t = performance.now();
  const ballot = voterView.cast(voter, vote);
  ballotBytes.push(ballot.toBytes());
  ballot.free();
  console.log(`cast voter ${voter} -> ${vote}: ${ms(t)}, ${ballotBytes.at(-1).length} B`);
}
voterView.free();

// Counter
const counter = evmlr.importParams(ppBytes);
counter.importKeys(skBytes);

// A damaged ballot is refused at import, before any proof is checked.
const damaged = ballotBytes[0].slice(0, -1);
try {
  counter.importBallot(damaged);
  console.log('MISMATCH: truncated ballot accepted');
  process.exit(1);
} catch (e) {
  if (!(e instanceof EvmlrError)) throw e;
  console.log('truncated ballot refused:', e.message);
}

const ballots = ballotBytes.map((b) => counter.importBallot(b));
t = performance.now();
const { tally, transcript } = counter.countWithTranscript(ballots);
console.log(`count: ${ms(t)}; transcript ${transcript.length} B`);
console.log('tally (shuffled):', tally);
ballots.forEach((b) => b.free());
counter.free();

// Auditor: public parameters, ballots and transcript; no keys.
const auditor = evmlr.importParams(ppBytes);
const audited = ballotBytes.map((b) => auditor.importBallot(b));
t = performance.now();
const verified = auditor.verifyCounting(audited, transcript);
console.log(`audit: ${ms(t)}; verified tally:`, verified);

// A transcript with one byte of a proof changed does not verify.
const tampered = transcript.slice();
tampered[tampered.length - 100] ^= 1;
let refused = false;
try {
  auditor.verifyCounting(audited, tampered);
} catch (e) {
  if (!(e instanceof EvmlrError)) throw e;
  refused = true;
  console.log('tampered transcript refused:', e.message);
}
audited.forEach((b) => b.free());
auditor.free();

const expected = [10, 20, 30, 40, 50];
const sorted = (xs) => JSON.stringify([...xs].sort((a, b) => a - b));
const ok = sorted(tally) === JSON.stringify(expected)
  && JSON.stringify(verified) === JSON.stringify(tally) && refused;
console.log(ok ? 'OK' : `MISMATCH, expected ${expected}`);
process.exit(ok ? 0 : 1);
