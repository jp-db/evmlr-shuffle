// End-to-end check of the WebAssembly library under Node:
//   make -f Makefile.wasm lib && node wasm/demo.mjs
// Mirrors test_voting_protocol in evmlr_voting.c: voter 2 votes twice and
// only the second vote counts.
import { loadEvmlr } from './evmlr.mjs';

const evmlr = await loadEvmlr();
console.log('params', evmlr.params);

const t0 = performance.now();
const election = evmlr.setup(3);
console.log(`setup (3 servers): ${(performance.now() - t0).toFixed(0)} ms`);

const casts = [[1, 10], [2, 99], [3, 30], [2, 20], [4, 40], [5, 50]];
const ballots = [];
for (const [voter, vote] of casts) {
  const t = performance.now();
  ballots.push(election.cast(voter, vote));
  console.log(`cast voter ${voter} -> ${vote}: ${(performance.now() - t).toFixed(0)} ms`);
}

const t1 = performance.now();
const tally = election.count(ballots);
console.log(`count: ${(performance.now() - t1).toFixed(0)} ms`);
console.log('tally (shuffled):', tally);

const expected = [10, 20, 30, 40, 50];
const ok = JSON.stringify([...tally].sort((a, b) => a - b)) === JSON.stringify(expected);
console.log(ok ? 'OK' : `MISMATCH, expected ${expected}`);

ballots.forEach((b) => b.free());
election.free();
process.exit(ok ? 0 : 1);
