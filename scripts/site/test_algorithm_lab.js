// Exact, bounded regression checks for the educational engines.
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const C = require('../../docs/algorithm-lab/core.js');
const multiples = [null, [5,1], [6,3], [10,6], [3,1], [9,16], [16,13], [0,6], [13,7], [7,6], [7,11], [13,10], [0,11], [16,4], [9,1], [3,16], [10,11], [6,14], [5,16]];
assert.equal(C.points.length, 18);
for (let i = 0; i < 19; i++) {
  assert.deepEqual(C.multiply(i), multiples[i]);
  for (let j = 0; j < 19; j++) assert.deepEqual(C.add(multiples[i], multiples[j]), multiples[(i + j) % 19]);
}
assert.equal(C.multiply(19), null);
for (let size = 2; size <= 9; size++) {
  const base = C.factorBase(size), pairs = C.pairs(base), full = C.relations(base);
  assert.equal(pairs.length, size * (size + 1) / 2);
  assert(full.every(r => C.checkRelation(base, r)));
  const solved = C.eliminate(full);
  assert.equal(solved.status, 'unique');
  const expected = base.map(p => multiples.findIndex(q => JSON.stringify(q) === JSON.stringify(p)));
  assert.deepEqual(solved.solution, expected);
  for (const step of solved.steps) for (const row of step.matrix) {
    assert.equal(C.mod(row.slice(0, -1).reduce((s, c, i) => s + c * expected[i], 0), 19), row.at(-1));
  }
  for (let k = 1; k < 19; k++) {
    const target = multiples[k], matches = C.decompose(base, target), recovery = C.recoverTarget(base, solved.solution, target);
    assert.equal(matches.length > 0, recovery !== null);
    if (recovery) { assert.equal(recovery.k, k); assert(recovery.verified); }
    for (const pair of matches) assert.deepEqual(C.add(base[pair.i], base[pair.j]), target);
  }
  const dependent = C.eliminate(C.relations(base, 'dependent'));
  assert.equal(dependent.status, 'underdetermined');
  assert.equal(dependent.rank, size - 1);
  assert.equal(dependent.solution, null);
  const badRows = C.relations(base, 'inconsistent');
  assert(!C.checkRelation(base, badRows.at(-1)));
  assert.equal(C.eliminate(badRows).status, 'inconsistent');
}
// Independent interpretation of each original teaching formula.
function original(kind, [x,y,z,w]) {
  if (z !== (x & y)) return false;
  if (kind === 'search') return Boolean((x || y) && (!x || y) && (x || !y)) && (z ^ w) === 0;
  if (kind === 'unsat') return (x ^ y) === 0 && (x ^ y) === 1 && (z ^ w) === 0;
  return (x ^ y ^ w) === 1 && (z ^ w) === 0;
}
const assignments = Array.from({length:16}, (_, mask) => Array.from({length:4}, (_, i) => (mask >> i) & 1));
for (const kind of ['gates', 'search', 'unsat']) {
  const example = C.satExample(kind), snapshot = JSON.stringify(example);
  const expectedModels = assignments.filter(a => original(kind, a));
  for (const mode of ['native', 'cnf']) {
    const encoded = C.encode(example, mode);
    for (const model of assignments) assert.equal(C.satisfies(encoded, model), original(kind, model));
    for (const first of [0, 1]) {
      const r = C.solveSAT(example, mode, first);
      assert.equal(r.status, expectedModels.length ? 'SAT' : 'UNSAT');
      assert(r.verified);
      if (r.model) assert(original(kind, r.model));
      assert.equal(r.counts.conflicts, r.steps.filter(s => s.kind === 'conflict').length);
      assert.equal(r.counts.decisions, r.steps.filter(s => s.kind === 'decide').length);
      assert.equal(r.counts.propagations, r.steps.filter(s => s.kind === 'propagate').length);
      if (kind === 'search' && first === 0) assert(r.steps.some(s => s.kind === 'backtrack'));
      assert.deepEqual(r.steps.at(-1).counts, r.counts);
    }
  }
  assert.equal(JSON.stringify(example), snapshot);
  const dimacs = C.dimacs(example).split('\n'), header = dimacs.find(l => l.startsWith('p ')).split(' ');
  const clauses = dimacs.filter(l => l && !/^[cp]/.test(l)).map(l => l.split(' ').map(Number));
  assert.equal(Number(header[3]), clauses.length);
  assert(clauses.every(c => c.pop() === 0));
  for (const model of assignments) assert.equal(clauses.every(c => c.some(lit => Boolean(model[Math.abs(lit)-1]) === (lit > 0))), original(kind, model));
}
// An empty clause must be detected at the root, without guessing.
assert.equal(C.solveSAT({clauses:[[]],xors:[]}, 'cnf').status, 'UNSAT');
const html = fs.readFileSync(path.join(__dirname, '../../docs/algorithm-lab.html'), 'utf8');
const ids = [...html.matchAll(/\bid="([^"]+)"/g)].map(m => m[1]);
assert.equal(new Set(ids).size, ids.length);
for (const anchor of html.matchAll(/href="#([^"]+)"/g)) assert(ids.includes(anchor[1]), anchor[1]);
console.log('PASS: 361 point sums, 8 factor-base sizes, 144 targets, modular row invariants, rank/contradiction cases, all SAT encodings and branches, DIMACS semantics, navigation.');
