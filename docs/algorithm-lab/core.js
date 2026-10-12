/* Small, exact teaching algorithms. No measured performance claims. */
(function (root) {
  'use strict';
  const FIELD = 17, ORDER = 19, G = [5, 1];
  const mod = (x, p) => ((x % p) + p) % p;
  const key = p => p ? p.join(',') : 'O';
  const equal = (p, q) => key(p) === key(q);
  function inverse(x, p) {
    x = mod(x, p);
    for (let n = 1; n < p; n++) if (mod(n * x, p) === 1) return n;
    throw new Error('This element has no multiplicative inverse.');
  }
  function onCurve(p) {
    return p === null || (p.length === 2 && p.every(x => Number.isInteger(x) && x >= 0 && x < FIELD) && mod(p[1] ** 2 - p[0] ** 3 - 2 * p[0] - 2, FIELD) === 0);
  }
  function add(p, q) {
    if (!p) return q && [...q];
    if (!q) return [...p];
    if (p[0] === q[0] && mod(p[1] + q[1], FIELD) === 0) return null;
    const slope = equal(p, q)
      ? mod((3 * p[0] ** 2 + 2) * inverse(2 * p[1], FIELD), FIELD)
      : mod((q[1] - p[1]) * inverse(q[0] - p[0], FIELD), FIELD);
    const x = mod(slope * slope - p[0] - q[0], FIELD);
    return [x, mod(slope * (p[0] - x) - p[1], FIELD)];
  }
  function multiply(k, p = G) {
    k = mod(k, ORDER);
    let result = null;
    for (let i = 0; i < k; i++) result = add(result, p);
    return result;
  }
  const points = [];
  for (let x = 0; x < FIELD; x++) for (let y = 0; y < FIELD; y++) if (onCurve([x, y])) points.push([x, y]);
  const representatives = points.filter(p => p[1] <= 8);
  function factorBase(size) {
    if (!Number.isInteger(size) || size < 2 || size > representatives.length) throw new Error('Choose 2–9 factor-base points.');
    return representatives.slice(0, size).map(p => [...p]);
  }
  function pairs(base) {
    const out = [];
    for (let i = 0; i < base.length; i++) for (let j = i; j < base.length; j++) {
      const coefficients = base.map((_, c) => Number(c === i) + Number(c === j));
      out.push({i, j, coefficients, point: add(base[i], base[j])});
    }
    return out;
  }
  function decompose(base, target) { return pairs(base).filter(p => equal(p.point, target)); }
  function eliminate(input, p = ORDER) {
    const matrix = input.map(row => row.map(x => mod(x, p)));
    const cols = matrix[0].length - 1, steps = [], pivots = [];
    const save = (message, rows = [], column = null) => steps.push({message, rows, column, matrix: matrix.map(row => [...row]), rank: pivots.length});
    save(`Start with ${matrix.length} equations in ${cols} unknowns, modulo ${p}.`);
    let row = 0;
    for (let col = 0; col < cols && row < matrix.length; col++) {
      const found = matrix.findIndex((r, i) => i >= row && r[col] !== 0);
      if (found < 0) continue;
      if (found !== row) {
        [matrix[row], matrix[found]] = [matrix[found], matrix[row]];
        save(`Swap R${row + 1} and R${found + 1}.`, [row, found], col);
      }
      const scale = inverse(matrix[row][col], p);
      matrix[row] = matrix[row].map(x => mod(x * scale, p));
      pivots.push(col);
      save(`R${row + 1} ← ${scale} · R${row + 1} (mod ${p}); pivot in column ℓ${col + 1}.`, [row], col);
      for (let r = 0; r < matrix.length; r++) {
        if (r === row || matrix[r][col] === 0) continue;
        const factor = matrix[r][col];
        matrix[r] = matrix[r].map((x, c) => mod(x - factor * matrix[row][c], p));
        save(`R${r + 1} ← R${r + 1} − ${factor} · R${row + 1} (mod ${p}).`, [r, row], col);
      }
      row++;
    }
    const inconsistent = matrix.some(r => r.slice(0, cols).every(x => x === 0) && r[cols] !== 0);
    const status = inconsistent ? 'inconsistent' : pivots.length < cols ? 'underdetermined' : 'unique';
    const solution = status === 'unique' ? Array(cols).fill(0) : null;
    if (solution) pivots.forEach((col, i) => { solution[col] = matrix[i][cols]; });
    save(status === 'unique' ? `Full rank ${cols}: each logarithm is determined.` : status === 'inconsistent' ? 'A row says 0 equals a nonzero residue: the system is inconsistent.' : `Rank ${pivots.length} < ${cols}: dependent equations leave free variables.`);
    return {steps, matrix, rank: pivots.length, status, solution};
  }
  function relations(base, mode = 'full') {
    // Enumerate known probes [a]G, then look up point decompositions. No logs
    // are supplied to elimination; enumeration is explicit teaching overhead.
    const probe = new Map(Array.from({length: ORDER}, (_, a) => [key(multiply(a)), a]));
    const candidates = pairs(base).sort((a, b) => Number(a.i === a.j) - Number(b.i === b.j));
    let rows = [];
    for (const relation of candidates) {
      const row = [...relation.coefficients, probe.get(key(relation.point))];
      if (eliminate([...rows, row]).rank > rows.length) rows.push(row);
      if (rows.length === base.length) break;
    }
    if (mode === 'dependent') rows = [...rows.slice(0, -1), [...rows[0]]];
    if (mode === 'inconsistent') rows.push([...rows[0].slice(0, -1), mod(rows[0].at(-1) + 1, ORDER)]);
    return rows;
  }
  function checkRelation(base, row) {
    const sum = base.reduce((s, point, i) => add(s, multiply(row[i], point)), null);
    return equal(sum, multiply(row.at(-1)));
  }
  function checkLogs(base, logs) { return logs !== null && logs.length === base.length && base.every((p, i) => equal(multiply(logs[i]), p)); }
  function recoverTarget(base, logs, target) {
    if (!checkLogs(base, logs)) return null;
    const relation = decompose(base, target)[0];
    if (!relation) return null;
    const k = mod(relation.coefficients.reduce((sum, c, i) => sum + c * logs[i], 0), ORDER);
    return {k, relation, verified: equal(multiply(k), target)};
  }
  const andClauses = [[-1, -2, 3], [1, -3], [2, -3]]; // z ↔ x ∧ y
  function satExample(kind) {
    if (kind === 'search') return {kind, name: 'Search and backtracking', clauses: [[1, 2], [-1, 2], [1, -2], ...andClauses.map(c => [...c])], xors: [{vars: [3, 4], rhs: 0}], description: '(x ∨ y) ∧ (¬x ∨ y) ∧ (x ∨ ¬y), z = x ∧ y, z ⊕ w = 0'};
    if (kind === 'unsat') return {kind, name: 'Contradictory parity', clauses: andClauses.map(c => [...c]), xors: [{vars: [1, 2], rhs: 0}, {vars: [1, 2], rhs: 1}, {vars: [3, 4], rhs: 0}], description: 'z = x ∧ y, x ⊕ y = 0, x ⊕ y = 1, z ⊕ w = 0'};
    return {kind: 'gates', name: 'AND and XOR', clauses: andClauses.map(c => [...c]), xors: [{vars: [1, 2, 4], rhs: 1}, {vars: [3, 4], rhs: 0}], description: 'z = x ∧ y, x ⊕ y ⊕ w = 1, z ⊕ w = 0'};
  }
  function xorToCNF(row) {
    const clauses = [];
    for (let mask = 0; mask < 2 ** row.vars.length; mask++) {
      const bits = row.vars.map((_, i) => (mask >> i) & 1);
      if (bits.reduce((a, b) => a ^ b, 0) !== row.rhs) clauses.push(row.vars.map((v, i) => bits[i] ? -v : v));
    }
    return clauses;
  }
  function encode(example, mode) {
    return {clauses: example.clauses.map(c => [...c]).concat(mode === 'cnf' ? example.xors.flatMap(xorToCNF) : []), xors: mode === 'cnf' ? [] : example.xors.map(r => ({vars: [...r.vars], rhs: r.rhs}))};
  }
  function clauseState(clause, model) {
    const values = clause.map(lit => model[Math.abs(lit) - 1] === null ? null : model[Math.abs(lit) - 1] === (lit > 0 ? 1 : 0));
    return values.includes(true) ? 'satisfied' : values.includes(null) ? 'open' : 'conflict';
  }
  function xorState(row, model) {
    if (row.vars.some(v => model[v - 1] === null)) return 'open';
    return row.vars.reduce((acc, v) => acc ^ model[v - 1], 0) === row.rhs ? 'satisfied' : 'conflict';
  }
  function satisfies(example, model) {
    return model.length === 4 && model.every(x => x === 0 || x === 1) && example.clauses.every(c => clauseState(c, model) === 'satisfied') && example.xors.every(r => xorState(r, model) === 'satisfied');
  }
  function exhaustiveModels(example) {
    return Array.from({length: 16}, (_, mask) => Array.from({length: 4}, (_, i) => (mask >> i) & 1)).filter(model => satisfies(example, model));
  }
  function solveSAT(example, mode = 'native', first = 0) {
    const constraints = encode(example, mode), steps = [], counts = {decisions: 0, propagations: 0, conflicts: 0};
    const names = ['x', 'y', 'z', 'w'];
    const save = (kind, message, assignment, depth, active = null) => steps.push({kind, message, assignment: [...assignment], depth, active, counts: {...counts}});
    const search = (assigned, depth) => {
      const model = [...assigned];
      while (true) {
        let forced = null;
        for (let i = 0; i < constraints.clauses.length; i++) {
          const c = constraints.clauses[i], state = clauseState(c, model);
          if (state === 'satisfied') continue;
          const open = c.filter(lit => model[Math.abs(lit) - 1] === null);
          if (!open.length) {
            counts.conflicts++; save('conflict', `Clause C${i + 1} is false. This branch cannot satisfy the formula.`, model, depth, ['clause', i]); return null;
          }
          if (open.length === 1) { forced = {v: Math.abs(open[0]) - 1, value: open[0] > 0 ? 1 : 0, active: ['clause', i], reason: `Unit clause C${i + 1}`}; break; }
        }
        if (!forced) for (let i = 0; i < constraints.xors.length; i++) {
          const r = constraints.xors[i], unknown = r.vars.filter(v => model[v - 1] === null);
          const assignedParity = r.vars.reduce((sum, v) => sum ^ (model[v - 1] ?? 0), 0);
          if (unknown.length === 0 && assignedParity !== r.rhs) {
            counts.conflicts++; save('conflict', `Parity row X${i + 1} has the wrong parity.`, model, depth, ['xor', i]); return null;
          }
          if (unknown.length === 1) { forced = {v: unknown[0] - 1, value: assignedParity ^ r.rhs, active: ['xor', i], reason: `Parity row X${i + 1}`}; break; }
        }
        if (!forced) break;
        model[forced.v] = forced.value; counts.propagations++;
        save('propagate', `${forced.reason} forces ${names[forced.v]} = ${forced.value}.`, model, depth, forced.active);
      }
      const next = model.indexOf(null);
      if (next < 0) { save('sat', 'Every original constraint is satisfied. Check the assignment in the truth table.', model, depth); return model; }
      for (const value of [first, 1 - first]) {
        counts.decisions++;
        const branch = [...model]; branch[next] = value;
        save('decide', `Try ${names[next]} = ${value} at decision level ${depth + 1}.`, branch, depth + 1);
        const result = search(branch, depth + 1);
        if (result) return result;
        save('backtrack', `Undo the branch ${names[next]} = ${value}; restore the assignments at level ${depth}.`, model, depth);
      }
      return null;
    };
    const initial = [null, null, null, null];
    save('start', 'Start with four unassigned Boolean variables.', initial, 0);
    const model = search(initial, 0);
    if (!model) save('unsat', 'Both branches are exhausted. The formula is UNSAT.', initial, 0);
    return {steps, model, status: model ? 'SAT' : 'UNSAT', constraints, counts, verified: model ? satisfies(example, model) : exhaustiveModels(example).length === 0};
  }
  function dimacs(example) {
    const {clauses} = encode(example, 'cnf');
    return `c Algorithm lab: ${example.name}; variables 1=x 2=y 3=z 4=w\np cnf 4 ${clauses.length}\n${clauses.map(c => c.join(' ') + ' 0').join('\n')}\n`;
  }
  const api = {FIELD, ORDER, G, mod, key, equal, inverse, onCurve, add, multiply, points, representatives, factorBase, pairs, decompose, eliminate, relations, checkRelation, checkLogs, recoverTarget, satExample, xorToCNF, encode, clauseState, xorState, satisfies, exhaustiveModels, solveSAT, dimacs};
  if (typeof module !== 'undefined' && module.exports) module.exports = api;
  else root.CryptoLab = api;
})(typeof globalThis !== 'undefined' ? globalThis : this);
