// Presentation only: all estimates and validation decisions are native results.
import fs from 'node:fs';
import path from 'node:path';
const root = path.dirname(new URL(import.meta.url).pathname);
const a = JSON.parse(fs.readFileSync(path.join(root, 'analysis.json'), 'utf8'));
const f = n => typeof n === 'number' ? n.toPrecision(7) : 'missing';
const esc = s => String(s).replaceAll('&', '&amp;').replaceAll('<', '&lt;').replaceAll('>', '&gt;');
const cells = ['ell8209-q29', 'ell131101-q37'];
const svgStart = (w,h) => `<svg xmlns="http://www.w3.org/2000/svg" width="${w}" height="${h}" viewBox="0 0 ${w} ${h}"><rect width="100%" height="100%" fill="#fff"/><g font-family="DejaVu Sans,sans-serif" fill="#152338">`;
const text = (x,y,s,size=14) => `<text x="${x}" y="${y}" font-size="${size}">${esc(s)}</text>`;
let graph = svgStart(940,365) + text(25,30,'IC floor/crater: complete simulated instruction count',20);
graph += text(25,55,'Eight archived one-target workloads per fixed curve; independent 99% bootstrap intervals',13);
const crosses = a.comparisons.filter(c => c.kind === 'floor_over_crater_ic');
const max = Math.max(2, ...crosses.map(c => c.estimate.bootstrap_99[1] * 1.1));
const px = n => 260 + n/max*630;
graph += `<line x1="${px(1)}" y1="80" x2="${px(1)}" y2="260" stroke="#8a98a9" stroke-dasharray="5 5"/>`;
for (let i=0;i<crosses.length;i++) {
  const c=crosses[i], e=c.estimate, y=115+i*95;
  graph += text(25,y+5,c.cell,17);
  graph += `<line x1="${px(e.bootstrap_99[0])}" y1="${y}" x2="${px(e.bootstrap_99[1])}" y2="${y}" stroke="#1764a0" stroke-width="5"/><circle cx="${px(e.ratio_of_means)}" cy="${y}" r="7" fill="#1764a0"/>`;
  graph += text(260,y+34,`${f(e.ratio_of_means)}; 99% [${e.bootstrap_99.map(f).join(', ')}]`,13);
}
graph += `<line x1="260" y1="265" x2="890" y2="265" stroke="#53657c"/>`;
for (let tick=0;tick<=max;tick+=0.5) graph += text(px(tick)-8,286,f(tick),11);
graph += text(25,325,`Unit: Callgrind Ir ratio. Profiles verified: ${a.verified_profiles}/${a.expected_profiles}.`,13);
graph += text(25,347,'Implementation census; external verification excluded. No intrinsic-hardness or native wall-time claim.',12) + '</g></svg>';
fs.writeFileSync(path.join(root,'instruction-ratios.svg'),graph);

function typst(dir,title,body,svg) {
  const source = `#set page(paper: "a4", margin: 18mm)\n#set text(font: "DejaVu Sans", size: 9pt)\n#set par(leading: 5pt)\n#text(size: 18pt, weight: "bold", ${JSON.stringify(title)})\n#parbreak()\n#image(${JSON.stringify(svg)}, width: 100%)\n#parbreak()\n#text(${JSON.stringify(body)})\n`;
  fs.writeFileSync(path.join(dir,'report.typ'),source);
}
let md = '# Conductor-gap complete instruction census\n\n';
md += `Measured on public synthetic prime-field workloads, 2026-10-14 UTC. Classification: accounting and boundary control. Follow-up to [PR #1724](https://github.com/aburan28/crypto/pull/1724).\n\n`;
md += `The complete existing MITM IC and matched negation-rho implementations are measured in Callgrind simulated user-space instructions. The cold interval includes setup, queries, relation handling, elimination, recovery and internal verification. It excludes input loading, known-answer target construction, process launch, kernel instructions and the external independent point check. These exclusions remain unmeasured; they are not charged zero. No native wall-time speedup is claimed.\n\n`;
md += `Profiles verified against the archived scalars, points and deterministic phase counters: **${a.verified_profiles}/${a.expected_profiles}**. The 32 archived workloads remain separate one-target solves. No setup or target transport is shared. The original eight samples per fixed curve are reused; this is accounting on those same points, not a fresh-target replication.\n\n`;
md += `![Instruction ratios](instruction-ratios.svg)\n\n| Gap | IC floor/crater Ir | 99% interval |\n|---|---:|---|\n`;
for (const c of crosses) md += `| ${c.cell} | ${f(c.estimate.ratio_of_means)} | ${c.estimate.bootstrap_99.map(f).join(' to ')} |\n`;
md += `\nFrozen repeated instruction-lead criterion: **${a.repeated_instruction_lead ? 'met within this implementation census' : 'not met'}**. A/A within the frozen 2% gate: **${a.aa_within_two_percent}**. Maximum absolute relative deviation: ${f(a.aa_max_deviation)} (fraction). Intervals including one are not equivalence proofs. The estimates are implementation observations on two fixed pairs and cannot establish a conductor effect across isogeny classes.\n\n`;
md += `Comparison status: **${a.comparison_status}**. A failed A/A gate makes the instruction comparisons exploratory; failed or missing profiles make inference inadmissible.\n\n`;
md += `See [native per-curve and per-target tables](analysis.md), [machine-readable analysis](analysis.json), [frozen protocol](PROTOCOL.md), [frozen jobs and original records](jobs.json), [execution commands](COMMANDS.md), [validation](VALIDATION.md) and [PDF](report.pdf). The prior study's exact EC1/UID representation manifests and Sage evidence remain at [the archived input study](../ecbench_conductor_gap_20261012/README.md). No map is constructed and no transfer cost is assumed.\n\n`;
md += `## Run reports\n\n`;
for (const cell of cells) {
  const dir=path.join(root,'reports',cell);
  fs.mkdirSync(dir,{recursive:true});
  const summaries=a.summaries.filter(s=>s.cell===cell);
  const rows=a.rows.filter(s=>s.cell===cell);
  let run=`# Run report: ${cell}\n\nUnit: Callgrind simulated user-space instructions over the entire methods::solve call. Profiles requested: ${rows.length}; verified: ${rows.filter(r=>!r.error).length}. Each fixed curve retains its eight archived independent one-target workloads. No wall-clock value is interpreted.\n\n`;
  run += `Valgrind 3.24.0; original native ecbench SHA-256 429de171c022b71427759b47d411d8f8d6dd52e9f2daa819de23f4a67536194a. Source revision 7df956f97dad09c96dab12772b2b4cb244146f. CPU feature dispatch under Callgrind: unknown. Host and tool receipts: ../../host.txt and ../../tool-hashes.sha256. Exact commands and excluded costs: ../../COMMANDS.md and ../../PROTOCOL.md.\n\n`;
  run += '| ICV1 curve model | Arm | n | Mean Ir | Range |\n|---|---|---:|---:|---|\n';
  for (const s of summaries) run += `| ${s.slug} | ${s.arm} | ${s.n} | ${f(s.mean_ir)} | ${s.min_ir} to ${s.max_ir} |\n`;
  run += '\n![Recorded instruction means](instructions.svg)\n\nThe bars show native-analyzer means and observed ranges, not confidence intervals. All raw profile IDs, source run/record IDs, scalar and count validation errors are in ../../analysis.json. The four additional IC repeats are A/A diagnostics, excluded from primary means. Profile receipts and raw parts are retained losslessly in the study evidence archive; all original sealed sessions remain unchanged.\n';
  let plot=svgStart(1000,470)+text(20,28,`${cell}: recorded whole-solve instruction counts`,19);
  plot+=text(20,52,'Means and observed ranges; log10 axis in Callgrind Ir; n = 8 per primary arm and curve',12);
  const logs=summaries.flatMap(s=>[Math.log10(s.min_ir),Math.log10(s.max_ir)]);
  const lo=Math.floor(Math.min(...logs)), hi=Math.ceil(Math.max(...logs));
  const x=n=>365+(Math.log10(n)-lo)/(hi-lo)*590;
  for (let i=0;i<summaries.length;i++) {
    const s=summaries[i], y=90+i*47, color=s.arm==='ic-mitm'?'#1764a0':'#8b516b';
    plot+=text(20,y,s.slug,12)+text(250,y,s.arm,12);
    plot+=`<line x1="${x(s.min_ir)}" y1="${y-4}" x2="${x(s.max_ir)}" y2="${y-4}" stroke="${color}" stroke-width="4"/><circle cx="${x(s.mean_ir)}" cy="${y-4}" r="5" fill="${color}"/>`;
  }
  plot+='<line x1="365" y1="383" x2="955" y2="383" stroke="#53657c"/>';
  for (let t=lo;t<=hi;t++) plot+=text(x(10**t)-13,405,`10^${t}`,12);
  plot+=text(20,443,'Raw sources: frozen jobs.json and profile directories; all failures remain recorded.',12)+'</g></svg>';
  fs.writeFileSync(path.join(dir,'instructions.svg'),plot);
  fs.writeFileSync(path.join(dir,'README.md'),run);
  typst(dir,`Instruction census: ${cell}`,run,'instructions.svg');
  md+=`- [${cell} run report](reports/${cell}/README.md), [PDF](reports/${cell}/report.pdf).\n`;
}
md += '\n## Requirements and unresolved work\n\n| Requirement | Status and evidence |\n|---|---|\n| Continue the large-gap experiment | Whole-solve instruction census on all archived one-target workloads |\n| Previously unpriced user-space implementation work | Included in the Ir interval; kernel/hardware work remains outside it |\n| Scalar correctness and deterministic phase replay | Native validation in analysis.json; another-host replay remains open |\n| Prime-field F4/SAT comparison | Not implemented in the native prime pipeline; remains open |\n| Explicit isogeny and transport costs | Not attempted; remains open |\n| Representation controls and more independent classes | Not attempted; two fixed pairs cannot establish a structural effect |\n| Intrinsic or cryptographic-scale hardness | Unresolved; this census supplies no general security conclusion |\n\n';
md += '## Canonical visual scope\n\nThe scoreboard receives a separate instruction-unit panel. The earlier GAE panel and figures remain historical. No new curve, map or conductor is constructed, so conductor/isogeny graph topology is unchanged. No generic GAE floor, scaling exponent, native wall-time leaderboard, or hypothesis state is promoted.\n';
fs.writeFileSync(path.join(root,'README.md'),md);
typst(root,'Conductor-gap complete instruction census',md,'instruction-ratios.svg');
