// End-to-end check of the WebAssembly library under Node:
//   make -f Makefile.wasm lib && node wasm/demo.mjs
//
// Three separate handles stand for the three parties, and only bytes pass
// between them: the authority runs setup, the voters cast from the public
// parameters alone, and the counter imports parameters, keys and ballots.
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
const tally = counter.count(ballots);
console.log(`count: ${ms(t)}`);
console.log('tally (shuffled):', tally);

const expected = [10, 20, 30, 40, 50];
const ok = JSON.stringify([...tally].sort((a, b) => a - b)) === JSON.stringify(expected);
console.log(ok ? 'OK' : `MISMATCH, expected ${expected}`);

ballots.forEach((b) => b.free());
counter.free();
process.exit(ok ? 0 : 1);
