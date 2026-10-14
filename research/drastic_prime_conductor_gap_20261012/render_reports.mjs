// Presentation only. Every quantitative value comes from the native analyzer or audit.
import fs from 'node:fs';
import path from 'node:path';
const root = path.dirname(new URL(import.meta.url).pathname);
const a = JSON.parse(fs.readFileSync(path.join(root, 'analysis.json'), 'utf8'));
const fmt = x => typeof x === 'number' ? x.toPrecision(6) : 'missing';
const interval = x => Array.isArray(x) ? x.map(fmt).join(' to ') : 'missing';
let md = `# Drastic prime conductor gaps, ECDLP measured on both ends

Measured on public synthetic prime-field curves, 2026-10-12 UTC. Classification: boundary control and counted-work accounting.

Six fresh certified CM pairs, two class-number-one discriminants (D = -19 and D = -43) at single-prime conductor gaps ell in {508433, 1016839, 2033609} (19 to 21 bits): 4 to 16 times the largest gap measured before (131101). Each curve has eight independent one-target workloads. The target law hashes the curve model: cross-curve samples are independent and no target was transported.

All numeric comparisons below concern counted S = group-addition equivalents / sqrt(r). The harness records unpriced work. A quotient of two lower bounds is a diagnostic, not a bound on complete attack cost. No wall-clock speedup, general hardness, or ECC2K-130 conclusion follows.

| pair | rho floor/crater | 95% interval | 99% interval | permutation p | IC floor/crater (99%) |
|---|---:|---|---|---:|---|
`;
for (const s of a.sessions) {
  const rho=s.arms['rho-neg'].cross_curve;
  const ic=s.arms['ic-mitm'] && s.arms['ic-mitm'].cross_curve;
  const icTxt = ic && ic.ratio_floor_over_crater !== undefined
    ? `${fmt(ic.ratio_floor_over_crater)} (${interval(ic.ci99)})` : 'not run at this size';
  md += `| ${s.cell} | ${fmt(rho.ratio_floor_over_crater)} | ${interval(rho.ci95)} | ${interval(rho.ci99)} | ${fmt(rho.permutation_p_two_sided)} | ${icTxt} |\n`;
}
md += `\nThe repeated-lead rule is applied exactly as frozen in the [protocol](PROTOCOL.md). This is not an equivalence proof.

![Independent-sample count comparisons](drastic-gap.svg)

Eight targets on a fixed curve are eight workloads on one curve, not eight curves. The floors were found by the prior round's method: random curve scan, one scalar multiplication per candidate, with the sample count recorded per pair.

## Recorded execution and replay

`;
for (const s of a.sessions) {
  const cell=path.basename(s.session);
  const audit=JSON.parse(fs.readFileSync(path.join(root,'sessions',`${cell}.audit.json`),'utf8'));
  md += `- **${cell}**, session ${s.session_id}: ${s.records} records, ${s.measured} measured; ${s.failed.length} measured failures. Audit OK: ${audit.ok}; verified records: ${audit.verified_records}; exact measured replays: ${audit.replays.length}. A/A pairs: ${s.aa.pairs}; mismatches: ${s.aa.mismatches.length}. [Run report](reports/${cell}/README.md), [audit](sessions/${cell}.audit.json), [sealed session](sessions/${cell}/session.json).\n`;
  const rd=path.join(root,'reports',cell);
  fs.mkdirSync(rd,{recursive:true});
  const host=JSON.parse(fs.readFileSync(path.join(root,'sessions',cell,'host.json'),'utf8'));
  let run=`# Run report: ${cell}\n\nSession: ${s.session_id}. Spec: ../../specs/${cell}.json. Exact launch/audit commands: ../../run_sessions.sh.\n\nHost: ${host.stable.cpu_model}, ${host.stable.arch}, environment ${host.env_class_id}. Binary SHA-256: ${host.build.binary_sha256}. Recorded revision: ${host.build.git_commit}; dirty: ${host.build.git_dirty}.\n\n${s.records} records: ${s.measured} measured, the rest retained warm-ups. Measured failures: ${s.failed.length}. Eight independent one-target workloads on each fixed curve; one measured round, one warm-up. All runs are cold. A/A matches: ${s.aa.pairs} pairs, ${s.aa.mismatches.length} mismatches.\n\nAudit OK: ${audit.ok}; ${audit.verified_records} verified records; ${audit.replays.length} exact measured replays. This is same-host replay, not independent-host validation.\n\n| curve | role | arm | mean counted S | observed range | lower bound | isolation levels |\n|---|---|---|---:|---|---|---|\n`;
  for(const arm of ['rho-neg','rho-neg-aa','ic-mitm']) {
    const armData = s.arms[arm];
    if (!armData || !armData.curves) continue;
    for(const r of armData.curves) {
      const levels=JSON.stringify(r.isolation_levels);
      run += `| ${r.slug} | ${r.role} | ${arm} | ${fmt(r.mean_s)} | ${fmt(r.min_s)} to ${fmt(r.max_s)} | ${r.lower_bound} | ${levels} |\n`;
    }
  }
  run += '\n![Recorded counts](counts.svg)\n\nDots are native-analyzer means; lines are observed ranges, not confidence intervals. No wall figure is interpreted. Raw records: ../../sessions/'+cell+'/records.jsonl; their SHA-256: '+s.records_sha256+'. All exact run IDs and phase means are in ../../analysis.json. No session file was edited.\n\n## Exact representation identities\n\n';
  for(const c of s.curve_manifests ?? []) {
    run += `- ${c.role}: ${c.curve_id}; ${c.curve_uid}. Full hash preimage and subgroup/generator are in ../../analysis.json. Metadata identity is not an additional mathematical certificate.\n`;
  }
  run += '\nUnpriced native work, memory and model-dependent primitive prices remain explicit in each record. Counted S and IC/rho quotients do not establish complete-cost advantages. No isogeny was computed or priced; no structural hardness or transport inference is made.\n';
  fs.writeFileSync(path.join(rd,'README.md'),run);
}
md += `\n## The pairs

\`\`\`json\n${JSON.stringify(JSON.parse(fs.readFileSync(path.join(root,'pairs','pairs.json'),'utf8')).pairs.map(p=>({d:p.d,k:p.k,ell:p.ell,t:p.t,p:p.p.toString(),N:p.group_order.toString(),floor_trials:p.floor_trials,j_crater:p.j_crater.toString(),j_floor:p.j_floor.toString()})), null, 1)}\n\`\`\`

Certification per pair: t^2-4q = -D*ell^2 with D fundamental class-number one fixes f_pi = ell; j = j(O_K) gives End = O_K (crater); j != j(O_K) with f_E | ell, f_E != 1 gives f_E = ell (floor). Orders verified by 64 random [N]P = O scalar multiplications per model before any session ran; every generator independently re-verified off-harness. Failures (wrong-trace twists, exhausted windows) are kept in [pairs/pairs.json](pairs/pairs.json).

## Reproduce

Build ecbench and the native examples from the pinned source. \`drastic_prime_gap_pairs --out pairs/pairs.json\`, then \`drastic_gap_specs --pairs pairs/pairs.json --out specs\`, then \`run_sessions.sh\` into fresh session directories; existing directories are never replaced. Run \`drastic_gap_analyze --sessions research/drastic_prime_conductor_gap_20261012/sessions --out research/drastic_prime_conductor_gap_20261012/analysis\`. The frozen statistical seed is 20261012; 20,000 independent-workload bootstrap resamples and label permutations.

## Requirement status

| requirement | status | evidence or remaining gap |
|---|---|---|
| Find drastic prime conductor gaps | measured | six certified pairs, ell up to 2033609, two discriminants; failures recorded |
| Measure ECDLP on both curves | measured for rho on all six pairs | counted one-target S with A/A controls and exact replays |
| Algebraic method coverage | ic-mitm on the two smallest-ell pairs | large-gap algebraic F4/SAT and MITM at the largest gaps remain open |
| Reproducible counts | audited on this host | external-host replay remains open |
| Endomorphism certification | class-number-one argument, no Sage | exact orders verified by scalar multiplication |
| Isogeny transport and its cost | not attempted | no explicit degree-ell maps; not assigned zero cost |
| Intrinsic or cryptographic-scale hardness | unresolved | six fixed pairs, two configurations; more independent classes required |

## Local validation

See [VALIDATION.md](VALIDATION.md).

## Evidence and visual scope

[Protocol](PROTOCOL.md), [native tables](analysis.md), [machine-readable analysis](analysis.json), [pair search output](pairs/pairs.json), [report PDF](report.pdf). Preregistration is this directory's PROTOCOL.md, committed before the search ran. The scoreboard receives a separate bounded-count panel for this measurement; previous fiber-round figures remain historical. No frontier, full-cost leaderboard or hypothesis status is promoted.
`;
fs.writeFileSync(path.join(root,'README.md'),md);
console.log('wrote README.md and run reports');
