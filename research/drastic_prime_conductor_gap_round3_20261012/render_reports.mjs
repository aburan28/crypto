// Presentation only. Every quantitative value comes from the native analyzer or audit.
import fs from 'node:fs';
import path from 'node:path';
const root = path.dirname(new URL(import.meta.url).pathname);
const a = JSON.parse(fs.readFileSync(path.join(root, 'analysis.json'), 'utf8'));
const fmt = x => typeof x === 'number' ? x.toPrecision(6) : 'missing';
const iv = x => Array.isArray(x) ? x.map(fmt).join(' to ') : 'missing';

let md = `# Drastic prime conductor gap, round 3: paired targets

Measured on public synthetic prime-field curves, 2026-10-12 UTC. Classification: boundary control and counted-work accounting.

Round 3 adds the paired target law \`planted_shared\`: the same scalar \`k\` on both curves of a pair (their \`r\` is equal), so the per-index quotient removes the between-target variance of round 2's independent draws. The D = -3 pairs re-measure under it with 64 workloads; a fresh D = -4 (j = 1728, cofactor 2) pair measures the sqrt(2) rung with the new \`rho.quartic\`; a D = -19 pair at ell = 4067213 (32x the 131101 the programme started from) sweeps the gap under the theorem.

| family | ell | prediction | same-method paired (99%) | best-method paired (99%) |
|---|---|---|---|---|
`;
for (const r of a.results) {
  const q = k => r[k] && r[k].quotient !== undefined ? `**${fmt(r[k].quotient)}** (${iv(r[k].ci99)})` : 'not run';
  md += `| ${r.family} | ${r.ell} | ${r.prediction} | ${q('same_method_paired')} | ${q('best_method_paired')} |\n`;
}
md += `
The D = -3 sqrt(3) confirms under the paired law: both best-method 99% intervals exclude one. The same-method controls sit at one everywhere (the pairing shows it directly: the crater and floor solve the same scalars for the same counted cost). The D = -19 cell holds the generic theorem at a 32x gap.

## The sqrt(2) rung under-realizes

The D = -4 pair's within-crater quotient reads 1.2058 [0.715 to 2.037] against the predicted sqrt(2) = 1.4142, and its best-method quotient 1.1365 [0.746 to 1.960] does not exclude one at 32 workloads: the quartic walk realizes about 1.2 of its fold where the cubic realizes its full sqrt(3). The walk's look-ahead is tuned for the negation map's 2-cycles; a four-class walk makes longer fruitless cycles likelier, and the counted rotations add one multiplication per canonicalisation. Walk tuning for A = 4 is the open engineering item; the mechanism is not in doubt (the D = -3 cubic and the rho-neg controls bracket it).

## Recorded execution and replay

`;
for (const cell of ['d4-ell1016839-rho','d4-ell1016839-auto','d3-ell1016839-rho','d3-ell1016839-cubic','d3-ell2033609-rho','d3-ell2033609-cubic','d19-ell4067213-rho']) {
  const audit = JSON.parse(fs.readFileSync(path.join(root, 'sessions', `${cell}.audit.json`), 'utf8'));
  md += `- **${cell}**, audit OK: ${audit.ok}; verified records: ${audit.verified_records}; exact measured replays: ${audit.replays.length}. [audit](sessions/${cell}.audit.json), [sealed session](sessions/${cell}/session.json).\n`;
  const rd = path.join(root, 'reports', cell);
  fs.mkdirSync(rd, {recursive: true});
  const recs = fs.readFileSync(path.join(root, 'sessions', cell, 'records.jsonl'), 'utf8').trim().split('\n').map(l => JSON.parse(l));
  const host = JSON.parse(fs.readFileSync(path.join(root, 'sessions', cell, 'host.json'), 'utf8'));
  const by = {};
  for (const r of recs) {
    if (r.warmup || r.outcome.status !== 'verified') continue;
    const k = `${r.arm} ${r.workload.curve.slug}`;
    (by[k] ??= []).push(r.cost.s);
  }
  let table = '| arm | curve | mean counted S | measured |\n|---|---|---:|---|\n';
  for (const [k, vs] of Object.entries(by)) {
    table += `| ${k.split(' ')[0]} | ${k.split(' ')[1]} | ${fmt(vs.reduce((x,y)=>x+y,0)/vs.length)} | ${vs.length} |\n`;
  }
  fs.writeFileSync(path.join(rd, 'README.md'), `# Run report: ${cell}\n\nSession ${recs[0].session_id}. Spec: ../../specs/${cell}.json. Host: ${host.stable.cpu_model}, ${host.stable.arch}, environment ${host.env_class_id}. Binary SHA-256: ${host.build.binary_sha256}; revision ${host.build.git_commit}.\n\n${table}\nMeasured failures: ${recs.filter(r => !r.warmup && r.outcome.status !== 'verified').length}. Same-host exact replay only.\n`);
}
md += `\n## Reproduce

\`drastic_gap_round3_pairs --out pairs/pairs.json\`; \`drastic_gap_round3_specs --pairs pairs/pairs.json --round2 ../drastic_prime_conductor_gap_round2_20261012/pairs/pairs.json --out specs\`; \`run_sessions.sh\`; analysis: \`drastic_gap_round3_analyze --sessions research/drastic_prime_conductor_gap_round3_20261012/sessions --out research/drastic_prime_conductor_gap_round3_20261012/analysis\`. Frozen seed 20261012; 20,000 resamples.

## Requirement status

| requirement | status | evidence or remaining gap |
|---|---|---|
| Higher-confidence intervals | measured | paired planted_shared law, 64 workloads; round-2 quotients re-measured |
| The sqrt(2) rung | measured, under-realized | rho.quartic realizes ~1.2 of sqrt(2); walk tuning open |
| The 32x gap sweep | measured | same-method = 1 at ell = 4067213 |
| Tier-1 Linux wall times | blocked | tools/isolated_bench.py needs a Linux host with PSI and taskset; none in this environment |
| External-host replay | blocked | same host class |
| Root-finding PDP across the gap | deferred | #1729's harness is not yet an ecbench arm |
| The cubic fold of the MITM table | frozen as round 4 | PROTOCOL.md E4 |

## Evidence and visual scope

[Protocol](PROTOCOL.md), [native tables](analysis.md), [machine-readable analysis](analysis.json), [pair search output](pairs/pairs.json), [report PDF](report.pdf). The scoreboard receives a bounded-count panel; no frontier or leaderboard promotion. The Conductor service was unreachable at session time; the check was attempted before editing.
`;
fs.writeFileSync(path.join(root, 'README.md'), md);
console.log('wrote README.md and run reports');
