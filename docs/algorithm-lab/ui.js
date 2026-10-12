(function () {
  'use strict';
  const C = window.CryptoLab, $ = id => document.getElementById(id);
  const point = p => p ? `(${p[0]}, ${p[1]})` : 'O (point at infinity)';
  const names = ['x', 'y', 'z', 'w'];
  let base, target, matrixRun, matrixStep = 0, originalRows, satRun, satStep = 0, example;
  function stat(value, label) { return `<div><strong>${value}</strong><span>${label}</span></div>`; }
  function curvePlot(matches) {
    const scale = v => 48 + v * 22, height = v => 400 - v * 22;
    let svg = '<title id="curve-title">Curve points, selected factor base, and target</title><desc id="curve-desc">Discrete points on y squared equals x cubed plus 2x plus 2 modulo 17. Teal dots belong to the factor base; the amber ring marks the target.</desc>';
    for (let n = 0; n <= 16; n++) {
      svg += `<line x1="${scale(n)}" x2="${scale(n)}" y1="48" y2="400" stroke="#e3eae6"/><line x1="48" x2="400" y1="${height(n)}" y2="${height(n)}" stroke="#e3eae6"/>`;
      if (n % 4 === 0) svg += `<text x="${scale(n)}" y="423" text-anchor="middle" font-size="12">${n}</text><text x="32" y="${height(n) + 4}" text-anchor="end" font-size="12">${n}</text>`;
    }
    const used = new Set(matches.flatMap(m => [m.i, m.j]));
    for (const p of C.points) {
      const idx = base.findIndex(b => C.equal(b, p));
      svg += `<g><title>${idx >= 0 ? `F${idx + 1}: ` : ''}${point(p)}${C.equal(p, target) ? ' — target R' : ''}</title><circle cx="${scale(p[0])}" cy="${height(p[1])}" r="${idx >= 0 ? 6 : 4}" fill="${idx >= 0 ? '#167f73' : '#b7c4c8'}"/>`;
      if (used.has(idx)) svg += `<circle cx="${scale(p[0])}" cy="${height(p[1])}" r="10" fill="none" stroke="#167f73" stroke-dasharray="3 2"/>`;
      if (idx >= 0) svg += `<text x="${scale(p[0]) + 9}" y="${height(p[1]) - 9}" font-size="11">F${idx + 1}</text>`;
      svg += '</g>';
    }
    if (target) svg += `<circle cx="${scale(target[0])}" cy="${height(target[1])}" r="13" fill="none" stroke="#a26c27" stroke-width="3"/><text x="${scale(target[0]) + 16}" y="${height(target[1]) + 4}" font-size="12">R</text>`;
    svg += '<text x="223" y="447" font-size="12" text-anchor="middle">x (mod 17)</text><text x="18" y="26" font-size="12">y (mod 17)</text>';
    $('curve-plot').innerHTML = svg;
  }
  function drawFactor() {
    base = C.factorBase(Number($('base-size').value));
    const a = Number($('target-k').value);
    target = C.multiply(a);
    $('target-value').textContent = `a = ${a}`;
    const all = C.pairs(base), matches = C.decompose(base, target);
    const covered = new Set(all.filter(r => r.point !== null).map(r => C.key(r.point))).size;
    $('factor-stats').innerHTML = stat(base.length, 'factor-base points') + stat(all.length, 'unordered pairs') + stat(`${covered} / 18`, 'nonidentity targets covered');
    $('base-points').innerHTML = base.map((p, i) => `<span class="point-chip">F${i + 1} = ${point(p)}</span>`).join('');
    $('target-heading').textContent = `R = [${a}]G = ${point(target)}`;
    $('decompositions').className = 'result' + (matches.length ? '' : ' warning');
    $('decompositions').innerHTML = matches.length
      ? `<strong>${matches.length} verified decomposition${matches.length > 1 ? 's' : ''}</strong>` + matches.map(m => `<p>F${m.i + 1} + F${m.j + 1} = ${point(target)}<br><small>ℓ${m.i + 1} + ℓ${m.j + 1} ≡ ${a} (mod 19)</small></p>`).join('')
      : `<strong>No two-point decomposition in this base.</strong><p>All ${all.length} unordered pairs were checked. Move the target or enlarge the base.</p>`;
    curvePlot(matches);
  }
  function resetMatrix() {
    const mode = $('matrix-case').value;
    originalRows = C.relations(base, mode);
    matrixRun = C.eliminate(originalRows); matrixStep = 0;
    $('matrix-context').textContent = mode === 'inconsistent'
      ? 'Deliberately invalid input: the last right-hand side is increased by 1. Its point check fails; a real collector must reject this relation before solving. Continue to see the algebraic contradiction.'
      : mode === 'dependent' ? 'The last independent equation is replaced by a copy of the first. More rows do not guarantee more independent information.'
      : 'Every original row is verified by point addition. Probes are known multiples of G; only the factor-base logarithms are unknown to elimination.';
    $('relation-list').innerHTML = originalRows.map((r, i) => {
      const terms = r.slice(0, -1).flatMap((c, j) => c ? [`${c === 1 ? '' : c + '·'}F${j + 1}`] : []).join(' + ');
      const valid = C.checkRelation(base, r);
      return `<p class="${valid ? '' : 'invalid'}">R${i + 1}: ${terms} = [${r.at(-1)}]G — ${valid ? 'verified' : 'INVALID: point sum disagrees'}</p>`;
    }).join('');
    drawMatrix();
  }
  function drawMatrix() {
    const step = matrixRun.steps[matrixStep], done = matrixStep === matrixRun.steps.length - 1;
    $('relation-matrix').innerHTML = '<caption>Augmented matrix [A | a], all entries modulo 19</caption><thead><tr><th scope="col">Row</th>' + base.map((_, i) => `<th scope="col">ℓ${i + 1}</th>`).join('') + '<th class="rhs" scope="col">a</th></tr></thead><tbody>' + step.matrix.map((row, i) => `<tr class="${step.rows.includes(i) ? 'active-row' : ''}"><th scope="row">R${i + 1}</th>${row.map((v, j) => `<td class="${j === base.length ? 'rhs' : ''} ${j === step.column && step.rows.includes(i) ? 'pivot' : ''}">${v}</td>`).join('')}</tr>`).join('') + '</tbody>';
    $('matrix-position').textContent = `Step ${matrixStep + 1} of ${matrixRun.steps.length}`;
    $('matrix-step').textContent = step.message;
    $('matrix-rank').textContent = `Pivots found: ${step.rank} / ${base.length}`;
    $('matrix-back').disabled = matrixStep === 0;
    $('matrix-next').disabled = done; $('matrix-finish').disabled = done;
    const result = $('matrix-solution'); result.className = 'result';
    if (!done) { result.textContent = 'Step through elimination to recover the logarithms, then verify them on the curve.'; return; }
    if (matrixRun.status !== 'unique') {
      result.className += matrixRun.status === 'inconsistent' ? ' error' : ' warning';
      result.textContent = matrixRun.status === 'inconsistent' ? 'Rejected: a zero coefficient row has a nonzero right-hand side. No logarithm vector satisfies these equations.' : `Underdetermined: rank ${matrixRun.rank} with ${base.length} unknowns. An additional independent relation is needed.`;
      return;
    }
    const logs = matrixRun.solution, verified = C.checkLogs(base, logs), recovered = C.recoverTarget(base, logs, target);
    result.innerHTML = `<strong>${verified ? 'Every recovered logarithm passes [ℓᵢ]G = Fᵢ.' : 'Verification failed.'}</strong><p>${logs.map((v, i) => `ℓ${i + 1} = ${v}`).join(' · ')}</p>`;
    if (recovered) result.innerHTML += `<p>For the selected target: k = ℓ${recovered.relation.i + 1} + ℓ${recovered.relation.j + 1} ≡ <strong>${recovered.k}</strong> (mod 19).<br>[${recovered.k}]G = ${point(target)}: ${recovered.verified ? 'verified' : 'FAILED'}.</p>`;
    else result.innerHTML += '<p>The chosen target has no two-point decomposition in this base. Its scalar cannot be recovered by this decomposition route. Try another target above.</p>';
  }
  function resetSAT() {
    example = C.satExample($('sat-case').value);
    satRun = C.solveSAT(example, $('sat-encoding').value, Number($('sat-first').value)); satStep = 0;
    $('sat-formula').textContent = example.description;
    $('sat-size').textContent = `${satRun.constraints.clauses.length} CNF clauses + ${satRun.constraints.xors.length} native XOR rows; 4 Boolean variables. In the all-CNF view, each parity row is expanded into equivalent clauses.`;
    drawSAT();
  }
  function drawSAT() {
    const step = satRun.steps[satStep], done = satStep === satRun.steps.length - 1;
    $('sat-position').textContent = `Step ${satStep + 1} of ${satRun.steps.length} · level ${step.depth}`;
    $('sat-event').textContent = ({start: 'Unassigned', decide: 'Decision', propagate: 'Propagation', conflict: 'Conflict', backtrack: 'Backtrack', sat: 'SAT — model found', unsat: 'UNSAT — no model'})[step.kind];
    $('sat-step').textContent = step.message;
    $('sat-assignment').innerHTML = step.assignment.map((v, i) => `<div><span>${names[i]}</span><b>${v === null ? '?' : v}</b></div>`).join('');
    const constraint = (expression, state, kind, index) => `<div class="constraint ${state}${step.active && step.active[0] === kind && step.active[1] === index ? ' active' : ''}"><span>${kind === 'clause' ? 'C' : 'X'}${index + 1}: ${expression}</span><small>${state}</small></div>`;
    $('sat-constraints').innerHTML = satRun.constraints.clauses.map((clause, i) => constraint(clause.map(lit => (lit < 0 ? '¬' : '') + names[Math.abs(lit) - 1]).join(' ∨ '), C.clauseState(clause, step.assignment), 'clause', i)).join('') + satRun.constraints.xors.map((row, i) => constraint(`${row.vars.map(v => names[v - 1]).join(' ⊕ ')} = ${row.rhs}`, C.xorState(row, step.assignment), 'xor', i)).join('');
    $('sat-counts').innerHTML = stat(step.counts.decisions, 'decisions') + stat(step.counts.propagations, 'forced values') + stat(step.counts.conflicts, 'conflicts');
    $('sat-verification').textContent = done ? satRun.model ? 'Verified against the original formula and its exhaustive truth table.' : 'Verified: none of the 16 assignments satisfy the original formula.' : 'Counts accumulate over the visited branches, including assignments later undone.';
    $('sat-back').disabled = satStep === 0; $('sat-next').disabled = done; $('sat-finish').disabled = done;
    let rows = '';
    for (let mask = 0; mask < 16; mask++) {
      const assignment = Array.from({length: 4}, (_, i) => (mask >> i) & 1), valid = C.satisfies(example, assignment);
      const current = step.assignment.every((v, i) => v !== null && v === assignment[i]);
      rows += `<tr class="${valid ? 'satisfying' : ''} ${current ? 'current-model' : ''}">${assignment.map(v => `<td>${v}</td>`).join('')}<td>${valid ? 'Satisfies' : 'Fails'}</td></tr>`;
    }
    $('truth-table').innerHTML = '<caption>All 16 assignments checked against the original formula</caption><thead><tr>' + names.map(n => `<th scope="col">${n}</th>`).join('') + '<th scope="col">Original formula</th></tr></thead><tbody>' + rows + '</tbody>';
  }
  $('base-size').addEventListener('change', () => { drawFactor(); resetMatrix(); });
  $('target-k').addEventListener('input', () => { drawFactor(); drawMatrix(); });
  $('find-hit').addEventListener('click', () => {
    const a = Number($('target-k').value);
    for (let n = 1; n <= 18; n++) {
      const next = (a + n - 1) % 18 + 1;
      if (C.decompose(base, C.multiply(next)).length) { $('target-k').value = next; break; }
    }
    drawFactor(); drawMatrix();
  });
  $('matrix-case').addEventListener('change', resetMatrix);
  $('matrix-reset').addEventListener('click', resetMatrix);
  $('matrix-back').addEventListener('click', () => { matrixStep = Math.max(0, matrixStep - 1); drawMatrix(); });
  $('matrix-next').addEventListener('click', () => { matrixStep = Math.min(matrixRun.steps.length - 1, matrixStep + 1); drawMatrix(); });
  $('matrix-finish').addEventListener('click', () => { matrixStep = matrixRun.steps.length - 1; drawMatrix(); });
  ['sat-case', 'sat-encoding', 'sat-first'].forEach(id => $(id).addEventListener('change', resetSAT));
  $('sat-reset').addEventListener('click', resetSAT);
  $('sat-back').addEventListener('click', () => { satStep = Math.max(0, satStep - 1); drawSAT(); });
  $('sat-next').addEventListener('click', () => { satStep = Math.min(satRun.steps.length - 1, satStep + 1); drawSAT(); });
  $('sat-finish').addEventListener('click', () => { satStep = satRun.steps.length - 1; drawSAT(); });
  $('sat-download').addEventListener('click', () => {
    const url = URL.createObjectURL(new Blob([C.dimacs(example)], {type: 'text/plain;charset=utf-8'}));
    const link = document.createElement('a'); link.href = url; link.download = `algorithm-lab-${example.kind}.cnf`;
    document.body.append(link); link.click(); link.remove(); setTimeout(() => URL.revokeObjectURL(url), 1000);
  });
  drawFactor(); resetMatrix(); resetSAT();
})();
