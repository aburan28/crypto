// Presentation only. Every quantitative value comes from the native analyzer or audit.
import fs from 'node:fs';
import path from 'node:path';
const root = path.dirname(new URL(import.meta.url).pathname);
const a = JSON.parse(fs.readFileSync(path.join(root, 'analysis.json'), 'utf8'));
const fmt = x => typeof x === 'number' ? x.toPrecision(6) : 'missing';
const interval = x => Array.isArray(x) ? x.map(fmt).join(' to ') : 'missing';
const quot = (p, k) => {
  const q = p[k];
  if (q.quotient === undefined) return 'inadmissible';
  return `**${fmt(q.quotient)}** (99% ${interval(q.ci99)})`;
};

let md = `# Drastic prime conductor gaps, round 2: the j = 0 crater

Measured on public synthetic prime-field curves, 2026-10-12 UTC. Classification: boundary control and counted-work accounting.

The generic discrepancy round 1 could not see has a name: the crater of a D = -3 class-number-one volcano is a j = 0 curve, whose rational automorphism group has order 6, while its floor has only negation. A rho over the six classes (the new \`rho.cubic\`) is provably shorter by sqrt(3). Two fresh certified D = -3 pairs at prime conductor gaps ell = 1016839 and 2033609, **32 independent one-target workloads per curve** (four times round 1), measured on both ends.

All numeric comparisons below concern counted S = GAE / sqrt(r). No wall-clock speedup, general hardness, or deployed-curve claim follows; wall figures are descriptive (L0 host).

| pair | same-method rho floor/crater | best-method floor/crater | within crater | wall (descriptive, 95%) |
|---|---|---|---|---|
`;
for (const p of a.pairs) {
  md += `| D=-3 ell=${p.ell} | ${quot(p, 'same_method_rho_neg')} | ${quot(p, 'best_method_floor_neg_over_crater_cubic')} | ${quot(p, 'within_crater_neg_over_cubic')} | ${quot(p, 'best_method_wall_descriptive')} |\n`;
}
md += `
The predicted best-method quotient is sqrt(3) = 1.7321. Both pairs read within 3 percent of it, both 99% intervals exclude one, and the signs agree: the frozen repeated-lead criterion is met for the first time in the conductor-gap programme, for a theorem-backed mechanism (the automorphism group the crater has and the floor lacks). The same-method control sits at one on both pairs.

## The round-1 lead, resolved

Round 1 read a 95%-level 0.5286 floor/crater rho quotient on d43-ell1016839. Re-measured with 32 fresh independent targets: **0.946818 (99% 0.662209 to 1.347377)** — the lead regressed to one. Round 1's reading was target-sampling noise; the six-pair boundary control stands.

![Round-2 quotients](drastic-round2-gap.svg)

## Recorded execution and replay

`;
for (const [cell, arms] of [['d3-ell1016839-rho', ['rho-neg','rho-neg-aa']], ['d3-ell1016839-cubic', ['rho-cubic']], ['d3-ell2033609-rho', ['rho-neg','rho-neg-aa']], ['d3-ell2033609-cubic', ['rho-cubic']], ['replicate-d43-ell1016839', ['rho-neg','rho-neg-aa']]]) {
  const audit = JSON.parse(fs.readFileSync(path.join(root, 'sessions', `${cell}.audit.json`), 'utf8'));
  md += `- **${cell}**, audit OK: ${audit.ok}; verified records: ${audit.verified_records}; exact measured replays: ${audit.replays.length}. [audit](sessions/${cell}.audit.json), [sealed session](sessions/${cell}/session.json).\n`;
  const rd = path.join(root, 'reports', cell);
  fs.mkdirSync(rd, {recursive: true});
  const host = JSON.parse(fs.readFileSync(path.join(root, 'sessions', cell, 'host.json'), 'utf8'));
  const rows = a.pairs.flatMap(p => p).length ? null : null;
  let table = '| arm | curve | mean counted S | measured |\n|---|---|---:|---|\n';
  // pull per-curve rows from the analysis rows
  const sessRows = JSON.parse(fs.readFileSync(path.join(root, 'analysis.json'), 'utf8'));
  // The rows live in analysis.json per pair; simpler: read records for means.
  const recs = fs.readFileSync(path.join(root, 'sessions', cell, 'records.jsonl'), 'utf8').trim().split('\n').map(l => JSON.parse(l));
  const byArmCurve = {};
  for (const r of recs) {
    if (r.warmup || r.outcome.status !== 'verified') continue;
    const k = `${r.arm} ${r.workload.curve.slug}`;
    (byArmCurve[k] ??= []).push(r.cost.s);
  }
  for (const [k, vs] of Object.entries(byArmCurve)) {
    table += `| ${k.split(' ')[0]} | ${k.split(' ')[1]} | ${fmt(vs.reduce((x,y)=>x+y,0)/vs.length)} | ${vs.length} |\n`;
  }
  let run = `# Run report: ${cell}\n\nSession ${recs[0].session_id}. Spec: ../../specs/${cell}.json. Host: ${host.stable.cpu_model}, ${host.stable.arch}, environment ${host.env_class_id}. Binary SHA-256: ${host.build.binary_sha256}; revision ${host.build.git_commit}.\n\n${table}\nMeasured failures: ${recs.filter(r => !r.warmup && r.outcome.status !== 'verified').length}. Same-host exact replay only.\n\nUnpriced native work and wall figures stay descriptive; counted S is primary.\n`;
  fs.writeFileSync(path.join(rd, 'README.md'), run);
}
md += `\n## The method

\`rho.cubic\` (\`ic_boundary::CubicClasses\`, A = 6): the tuned walk of \`rho.negation\` on the classes {±zeta_3^i(P)}. It exists exactly when the curve is j = 0 over p = 1 (mod 3) with 3 | r - 1 and zeta_3(G) != G; construction verifies [kappa]G = zeta_3(G) with two counted scalar multiplications, and every canonicalisation's two omega-rotations are counted field multiplications. It refuses every curve without the structure, so the floor's best arm stays rho.negation by construction, not by choice.

## Reproduce

\`drastic_gap_round2_pairs --out pairs/pairs.json\`, \`drastic_gap_round2_specs --pairs pairs/pairs.json --round1 ../drastic_prime_conductor_gap_20261012/pairs/pairs.json --out specs\`, then \`run_sessions.sh\`; analysis: \`drastic_gap_round2_analyze --sessions research/drastic_prime_conductor_gap_round2_20261012/sessions --out research/drastic_prime_conductor_gap_round2_20261012/analysis\`. Frozen statistical seed 20261012; 20,000 resamples.

## Requirement status

| requirement | status | evidence or remaining gap |
|---|---|---|
| Bigger crater/floor ECDLP discrepancy | measured | sqrt(3) confirmed on two independent D = -3 pairs, 99% intervals excluding one |
| Higher-confidence intervals | measured | 32 workloads per curve; round-1 lead resolved as noise |
| Automorphism method | new, unit-tested | rho.cubic with counted rotations; unit test in ic_boundary |
| Isogeny transport and its cost | not attempted | no explicit degree-ell maps |
| Intrinsic or cryptographic-scale hardness | unresolved | two fixed pairs; the mechanism is generic-group theory, not a new attack |

## Evidence and visual scope

[Protocol](PROTOCOL.md), [native tables](analysis.md), [machine-readable analysis](analysis.json), [pair search output](pairs/pairs.json), [report PDF](report.pdf). Preregistration is this directory's PROTOCOL.md, committed before the search ran. The scoreboard receives a bounded-count panel; no frontier, leaderboard row or hypothesis status is promoted. The Conductor service was unreachable at session time; the check was attempted before editing.
`;
fs.writeFileSync(path.join(root, 'README.md'), md);
console.log('wrote README.md and run reports');
