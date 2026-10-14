// Presentation only: all estimates and validation decisions are native results.
import fs from 'node:fs';
import path from 'node:path';
const root = path.dirname(new URL(import.meta.url).pathname);
const a = JSON.parse(fs.readFileSync(path.join(root, 'analysis.json'), 'utf8'));
const config = JSON.parse(fs.readFileSync(path.join(root, 'profile-config.json'), 'utf8'));
const f = n => typeof n === 'number' ? n.toPrecision(7) : 'missing';
const esc = s => String(s).replaceAll('&', '&amp;').replaceAll('<', '&lt;').replaceAll('>', '&gt;');
const cells = ['ell8209-q29', 'ell131101-q37'];
const svgStart = (w,h) => `<svg xmlns="http://www.w3.org/2000/svg" width="${w}" height="${h}" viewBox="0 0 ${w} ${h}"><rect width="100%" height="100%" fill="#fff"/><g font-family="DejaVu Sans,sans-serif" fill="#152338">`;
const text = (x,y,s,size=14) => `<text x="${x}" y="${y}" font-size="${size}">${esc(s)}</text>`;
const incompleteGraph=(title,rows)=>svgStart(940,240)+text(25,35,title,20)
  +text(25,85,`Profiles requested: ${rows.length}; verified: ${rows.filter(r=>!r.error).length}.`,17)
  +text(25,125,'The full-panel verification gate failed; no primary comparison is admitted.',15)
  +text(25,163,'Every missing or failed profile remains in analysis.json with its reason.',14)
  +text(25,207,'A coverage diagram replaces numerical inference; no successful subset is analyzed.',13)+'</g></svg>';
let graph = svgStart(940,365) + text(25,30,'IC floor/crater: complete simulated instruction count',20);
graph += text(25,55,'Eight archived one-target workloads per fixed curve; independent 99% bootstrap intervals',16);
const crosses = a.comparisons.filter(c => c.kind === 'floor_over_crater_ic');
const max = Math.max(2, ...crosses.map(c => c.estimate.bootstrap_99[1] * 1.1));
const px = n => 260 + n/max*630;
graph += `<line x1="${px(1)}" y1="80" x2="${px(1)}" y2="260" stroke="#8a98a9" stroke-dasharray="5 5"/>`;
for (let i=0;i<crosses.length;i++) {
  const c=crosses[i], e=c.estimate, y=115+i*95;
  graph += text(25,y+5,c.cell,17);
  graph += `<line x1="${px(e.bootstrap_99[0])}" y1="${y}" x2="${px(e.bootstrap_99[1])}" y2="${y}" stroke="#1764a0" stroke-width="5"/><circle cx="${px(e.ratio_of_means)}" cy="${y}" r="7" fill="#1764a0"/>`;
  graph += text(260,y+34,`${f(e.ratio_of_means)}; 99% [${e.bootstrap_99.map(f).join(', ')}]`,16);
}
graph += `<line x1="260" y1="265" x2="890" y2="265" stroke="#53657c"/>`;
for (let tick=0;tick<=max;tick+=0.5) graph += text(px(tick)-8,286,tick.toFixed(1),15);
graph += text(25,325,`Unit: Callgrind Ir ratio. Profiles verified: ${a.verified_profiles}/${a.expected_profiles}.`,16);
graph += text(25,347,'Implementation census. External check excluded. No intrinsic-hardness or native wall-time claim.',15) + '</g></svg>';
if(!a.all_profiles_verified) graph=incompleteGraph('Instruction census: incomplete verification',a.rows);
fs.writeFileSync(path.join(root,'instruction-ratios.svg'),graph);

function typst(dir,title,body,svg) {
  const clean=s=>s.replace(/\[([^\]]+)\]\(([^)]+)\)/g,'$1 ($2)').replaceAll('**','').replaceAll('`','');
  const lines=body.split('\n');
  let content='';
  for(let i=0;i<lines.length;i++) {
    const line=lines[i];
    if(!line.trim() || line.startsWith('# ') || line.startsWith('![')) continue;
    if(line.startsWith('## ')) {
      content+=`#heading(level: 2, ${JSON.stringify(clean(line.slice(3)))})\n`;
    } else if(line.startsWith('|')) {
      const table=[];
      while(i<lines.length && lines[i].startsWith('|')) {
        const cells=lines[i].split('|').slice(1,-1).map(s=>clean(s.trim()));
        if(!cells.every(s=>/^:?-+:?$/.test(s))) table.push(cells);
        i++;
      }
      i--;
      content+=`#text(size: 7.5pt)[#table(columns: ${table[0].length}, inset: 4pt, ${table.flat().map(s=>`text(${JSON.stringify(s)})`).join(', ')})]\n#parbreak()\n`;
    } else {
      content+=`#text(${JSON.stringify(clean(line))})\n#parbreak()\n`;
    }
  }
  const source = `#set page(paper: "a4", margin: 18mm)\n#set text(font: "DejaVu Sans", size: 9pt)\n#set par(leading: 5pt)\n#text(size: 18pt, weight: "bold", ${JSON.stringify(title)})\n#parbreak()\n#image(${JSON.stringify(svg)}, width: 100%)\n#parbreak()\n${content}`;
  fs.writeFileSync(path.join(dir,'report.typ'),source);
}
let md = '# Conductor-gap complete instruction census\n\n';
md += `Measured on public synthetic prime-field workloads, 2026-10-14 UTC. Classification: accounting and boundary control. Follow-up to [PR #1724](https://github.com/aburan28/crypto/pull/1724).\n\n`;
md += `The existing MITM IC and matched negation-rho implementations are measured in Callgrind simulated user-space instructions over the complete methods::solve call. The interval includes cold setup inside the adapter, queries, relation handling, elimination, recovery and internal verification. It excludes input loading, method resolution, known-answer target construction, process launch, kernel instructions and the external independent point check. These exclusions remain unmeasured; they are not charged zero. No native wall-time speedup is claimed.\n\n`;
md += `Profiles verified against the archived scalars, points and deterministic phase counters: **${a.verified_profiles}/${a.expected_profiles}**. The 32 archived workloads remain separate one-target solves. No setup or target transport is shared. The original eight samples per fixed curve are reused; this is accounting on those same points, not a fresh-target replication.\n\n`;
md += `![Instruction ratios](instruction-ratios.svg)\n\n| Gap | IC floor/crater Ir | 99% interval |\n|---|---:|---|\n`;
for (const c of crosses) md += `| ${c.cell} | ${f(c.estimate.ratio_of_means)} | ${c.estimate.bootstrap_99.map(f).join(' to ')} |\n`;
md += `\nFrozen repeated instruction-lead criterion: **${a.repeated_instruction_lead ? 'met within this implementation census' : 'not met'}**. A/A within the frozen 2% gate: **${a.aa_within_two_percent}**. Maximum absolute relative deviation: ${f(a.aa_max_deviation)} (fraction). Intervals including one are not equivalence proofs. The estimates are implementation observations on two fixed pairs and cannot establish a conductor effect across isogeny classes.\n\n`;
md += `Comparison status: **${a.comparison_status}**. A failed A/A gate makes the instruction comparisons exploratory; failed or missing profiles make inference inadmissible.\n\n`;
md+='| ICV1 model | IC / matched rho Ir | Paired 99% interval |\n|---|---:|---|\n';
for(const c of a.comparisons.filter(c=>c.kind==='ic_over_rho')) md+=`| ${c.slug} | ${f(c.estimate.ratio_of_means)} | ${c.estimate.bootstrap_99.map(f).join(' to ')} |\n`;
md+='\nThese ratios compare the same one-target workloads within the frozen instruction interval; excluded work and native wall time remain unmeasured.\n\n';
if(!a.all_profiles_verified) {
  md+='No primary estimate is admitted. All missing or failed profiles remain below; successful subsets are not compared.\n\n| Profile | Recorded validation failure |\n|---|---|\n';
  for(const r of a.rows.filter(r=>r.error)) md+=`| ${r.id} | ${r.error.replaceAll('|','/')} |\n`;
  md+='\n';
}
md += `See [native per-curve and per-target tables](analysis.md), [machine-readable analysis](analysis.json), [frozen protocol](PROTOCOL.md), [amended protocol](PROTOCOL-v2.md), [frozen jobs and original records](jobs.json), [lossless raw evidence](raw-callgrind.tar.xz), [execution commands](COMMANDS.md), [validation](VALIDATION.md) and [PDF](report.pdf). The prior study's exact EC1/UID representation manifests and Sage evidence remain at [the archived input study](../ecbench_conductor_gap_20261012/README.md). No map is constructed and no transfer cost is assumed.\n\n`;
md += `## Run reports\n\n`;
for (const cell of cells) {
  const dir=path.join(root,'reports',cell);
  fs.mkdirSync(dir,{recursive:true});
  const summaries=a.summaries.filter(s=>s.cell===cell);
  const rows=a.rows.filter(s=>s.cell===cell);
  let run=`# Run report: ${cell}\n\nUnit: Callgrind simulated user-space instructions over the entire methods::solve call. Profiles requested: ${rows.length}; verified: ${rows.filter(r=>!r.error).length}. Each fixed curve retains its eight archived independent one-target workloads. No wall-clock value is interpreted.\n\n`;
  run += `Valgrind 3.24.0; v2 native ecbench SHA-256 ${config.binary_sha256}. Source snapshot ${config.source_commit}; phase dump suppression is enabled while complete collection and outer markers remain active. CPU feature dispatch under Callgrind: unknown. Host and tool receipts: ../../host.txt and ../../tool-hashes.sha256. Exact commands and excluded costs: ../../COMMANDS.md, ../../PROTOCOL.md and ../../PROTOCOL-v2.md.\n\n`;
  run += '| ICV1 curve model | Arm | n | Mean Ir | Range |\n|---|---|---:|---:|---|\n';
  for (const s of summaries) run += `| ${s.slug} | ${s.arm} | ${s.n} | ${f(s.mean_ir)} | ${s.min_ir} to ${s.max_ir} |\n`;
  run += '\n![Recorded instruction means](instructions.svg)\n\nThe marks show native-analyzer means and observed ranges, not confidence intervals. All raw profile IDs, source run/record IDs, scalar and count validation errors are in ../../analysis.json. The four additional IC repeats are A/A diagnostics, excluded from primary means. Profile receipts and raw parts are retained losslessly in ../../raw-callgrind.tar.xz; all original sealed sessions remain unchanged.\n';
  let plot=svgStart(1000,470)+text(20,28,`${cell}: recorded whole-solve instruction counts`,19);
  plot+=text(20,52,'Means and observed ranges; log10 axis in Callgrind Ir; n = 8 per primary arm and curve',16);
  const logs=summaries.flatMap(s=>[Math.log10(s.min_ir),Math.log10(s.max_ir)]);
  const lo=Math.floor(Math.min(...logs)), hi=Math.ceil(Math.max(...logs));
  const x=n=>420+(Math.log10(n)-lo)/(hi-lo)*535;
  for (let i=0;i<summaries.length;i++) {
    const s=summaries[i], y=90+i*47, color=s.arm==='ic-mitm'?'#1764a0':'#8b516b';
    plot+=text(20,y,s.slug,16)+text(290,y,s.arm,16);
    plot+=`<line x1="${x(s.min_ir)}" y1="${y-4}" x2="${x(s.max_ir)}" y2="${y-4}" stroke="${color}" stroke-width="4"/><circle cx="${x(s.mean_ir)}" cy="${y-4}" r="5" fill="${color}"/>`;
  }
  plot+='<line x1="420" y1="383" x2="955" y2="383" stroke="#53657c"/>';
  for (let t=lo;t<=hi;t++) plot+=text(x(10**t)-13,405,`10^${t}`,16);
  plot+=text(20,443,'Raw sources: frozen jobs.json and profile directories; all failures remain recorded.',16)+'</g></svg>';
  if(!a.all_profiles_verified) {
    plot=incompleteGraph(`${cell}: full-panel inference blocked`,rows);
    run+='\nThe complete-panel verification gate failed. No cell estimate is admitted, including for a cell whose individual profiles verified. This diagram records coverage only; detailed failure reasons remain in ../../analysis.json.\n';
  }
  fs.writeFileSync(path.join(dir,'instructions.svg'),plot);
  fs.writeFileSync(path.join(dir,'README.md'),run);
  typst(dir,`Instruction census: ${cell}`,run,'instructions.svg');
  md+=`- [${cell} run report](reports/${cell}/README.md), [PDF](reports/${cell}/report.pdf).\n`;
}
md += '\n## Requirements and unresolved work\n\n| Requirement | Status and evidence |\n|---|---|\n| Continue the large-gap experiment | Whole-solve instruction census on all archived one-target workloads |\n| Previously unpriced user-space implementation work | Included in the Ir interval; kernel/hardware work remains outside it |\n| Scalar correctness and deterministic phase replay | Native validation in analysis.json; another-host replay remains open |\n| Prime-field F4/SAT comparison | Not implemented in the native prime pipeline; remains open |\n| Explicit isogeny and transport costs | Not attempted; remains open |\n| Representation controls and more independent classes | Not attempted; two fixed pairs cannot establish a structural effect |\n| Intrinsic or cryptographic-scale hardness | Unresolved; this census supplies no general security conclusion |\n\n';
md += 'Interrupted attempts: [launch before preregistration](attempt-unpreregistered/README.md) and [phase-dump volume obstruction](attempt-phase-dumps/README.md). Both retain all raw files and supply no primary estimate. The [profiling obstruction](../../problems/conductor-gap-callgrind-phase-dump-volume.md) and [missing prime-field algebraic oracle](../../problems/prime-field-conductor-gap-algebraic-oracle.md) track their exact obligations.\n\n';
md += '## Canonical visual scope\n\nThe scoreboard receives a separate instruction-unit panel. The earlier GAE panel and figures remain historical. No new curve, map or conductor is constructed, so conductor/isogeny graph topology is unchanged. No generic GAE floor, scaling exponent, native wall-time leaderboard, or hypothesis state is promoted.\n';
fs.writeFileSync(path.join(root,'README.md'),md);
typst(root,'Conductor-gap complete instruction census',md,'instruction-ratios.svg');

const repo=path.resolve(root,'../..');
const galleryName='20261014-conductor-gap-instructions';
const gallery=path.join(repo,'benchmarks',galleryName);
fs.mkdirSync(gallery,{recursive:true});
const headline=crosses.length ? crosses.map(c=>`${c.cell}: ${f(c.estimate.ratio_of_means)} (99% ${c.estimate.bootstrap_99.map(f).join(' to ')})`).join('; ') : 'No primary estimate admitted';
fs.writeFileSync(path.join(gallery,'README.md'),`# Conductor-gap instruction census\n\nSource: [study and frozen evidence](../../research/ecbench_conductor_gap_ir_20261014/README.md). Unit: Callgrind simulated user-space instructions within methods::solve. ${a.verified_profiles}/${a.expected_profiles} profiles verified; comparison status: ${a.comparison_status}.\n\nIC floor/crater: ${headline}. No native wall-time or intrinsic-hardness conclusion. Method resolution, target construction, kernel work and external point verification remain unmeasured.\n\nPDFs, SVG graphs and 150 dpi PNGs are copied from the same report sources:\n\n- [Results report](report.pdf), [graph](report.svg), [PNG](report.png).\n- [ell8209-q29 run](ell8209-q29.pdf), [graph](ell8209-q29.svg), [PNG](ell8209-q29.png).\n- [ell131101-q37 run](ell131101-q37.pdf), [graph](ell131101-q37.svg), [PNG](ell131101-q37.png).\n- [Unpreregistered interrupted attempt](attempt-unpreregistered.pdf), [diagram](attempt-unpreregistered.svg), [PNG](attempt-unpreregistered.png).\n- [Phase-dump interrupted attempt](attempt-phase-dumps.pdf), [diagram](attempt-phase-dumps.svg), [PNG](attempt-phase-dumps.png).\n`);
const indexPath=path.join(repo,'benchmarks/README.md');
let index=fs.readFileSync(indexPath,'utf8');
const row=`| 2026-10-14 | [\`${galleryName}\`](${galleryName}/) | Complete methods::solve instruction census on two fixed conductor-gap pairs; native MITM IC and matched negation rho | IC floor/crater Ir: ${headline}; ${a.verified_profiles}/${a.expected_profiles} verified profiles; status ${a.comparison_status} | Xeon Platinum 8573C, shared virtual x86-64, instruction counts only; no native wall-time claim | [Study and reports](../research/ecbench_conductor_gap_ir_20261014/) |`;
if(index.includes(`](${galleryName}/)`)) {
  index=index.split('\n').map(line=>line.includes(`](${galleryName}/)`)?row:line).join('\n');
} else {
  const lines=index.split('\n');
  const first=lines.findIndex(line=>line.startsWith('| 20'));
  if(first<0) throw new Error('benchmark index table missing');
  lines.splice(first,0,row);
  index=lines.join('\n');
}
fs.writeFileSync(indexPath,index);
