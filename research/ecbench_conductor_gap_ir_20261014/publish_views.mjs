// Copies native observations into canonical presentation views; no inference.
import fs from 'node:fs';
import path from 'node:path';
const study=path.dirname(new URL(import.meta.url).pathname);
const repo=path.resolve(study,'../..');
const read=p=>fs.readFileSync(path.join(repo,p),'utf8');
const write=(p,s)=>fs.writeFileSync(path.join(repo,p),s);
const a=JSON.parse(read('research/ecbench_conductor_gap_ir_20261014/analysis.json'));
if(!a.all_profiles_verified) throw new Error('failure-bearing evidence requires a failure panel, not this complete-census view');
const anchor='conductor-gap-ir-20261014';
const source='https://github.com/aburan28/crypto/blob/main/research/ecbench_conductor_gap_ir_20261014/README.md';
const esc=s=>String(s).replaceAll('&','&amp;').replaceAll('<','&lt;').replaceAll('>','&gt;').replaceAll('"','&quot;');
const fmt=n=>Number.isFinite(n)?n.toPrecision(7):'missing';
const interval=a=>a.map(fmt).join('–');
const crosses=a.comparisons.filter(c=>c.kind==='floor_over_crater_ic');
const summary=crosses.map(c=>`${c.cell}: ${fmt(c.estimate.ratio_of_means)} [${interval(c.estimate.bootstrap_99)}]`).join('; ');
let table='';
for(const s of a.summaries) table+=`<tr><td>${esc(s.cell)}</td><td><code>${esc(s.slug)}</code></td><td>${esc(s.arm)}</td><td>${s.n}</td><td>${fmt(s.mean_ir)}</td><td>${fmt(s.mean_ir_over_sqrt_r)}</td></tr>\n`;
const panel=`<section class="panel" id="${anchor}">
<h2>Conductor gap: complete implementation instruction census <span class="chip">accounting</span></h2>
<p class="panel-summary">The archived prime-field IC and matched rho solves are counted in Callgrind user-space instructions. ${a.verified_profiles}/${a.expected_profiles} profiles passed native point/scalar and exact deterministic-counter replay.</p>
<p>IC floor/crater ratios and independent 99% intervals: ${esc(summary)}. Frozen repeated instruction lead: ${a.repeated_instruction_lead}. A/A within 2%: ${a.aa_within_two_percent}. Comparison status: ${esc(a.comparison_status)}.</p>
<div class="scroll"><table><caption>Cold methods::solve instruction census; eight archived one-target workloads per fixed curve and primary arm. Setup and all user-space work inside the adapter are included. Method resolution, target construction, kernel work and the external independent verification are excluded and remain unmeasured. Ir is a separate unit from the GAE generic floor.</caption><thead><tr><th>gap cell</th><th>ICV1 model</th><th>arm</th><th>n</th><th>mean Ir</th><th>mean Ir/sqrt(r)</th></tr></thead><tbody>${table}</tbody></table></div>
<p>The same archived targets are reused for accounting; no new map or target transport is constructed. These two fixed pairs cannot establish intrinsic hardness. Prime-field F4/SAT, representation controls and another-host validation remain open. <a href="${source}">Reports, PDFs, frozen jobs and lossless evidence</a>.</p>
</section>\n`;
let page=read('docs/index-calculus-scoreboard.html');
if(page.includes(`id="${anchor}"`)) throw new Error('panel already exists');
page=page.replace('</main>',panel+'</main>');
if(!page.includes(`id="${anchor}"`)) throw new Error('main marker missing');
const timeline=JSON.parse(read('docs/ic/progress-timeline.json'));
const series=timeline.series.find(s=>s.id==='ec-callgrind-diagnostics');
if(!series) throw new Error('instruction-unit series missing');
let rows='';
for(const c of a.comparisons.filter(c=>c.kind==='ic_over_rho')) {
  const point={date:'2026-10-14',value:c.estimate.ratio_of_means,reference:'matched',unit:'Callgrind Ir diagnostic',
    label:`Conductor gap ${c.cell}, ${c.slug}, archived eight-target IC/rho census`,
    detail:`Complete user-space implementation Ir; paired 95% interval ${interval(c.estimate.bootstrap_95)}. ${a.verified_profiles}/${a.expected_profiles} profiles verified. Comparison status: ${a.comparison_status}. Method resolution, target preparation, kernel work and external audit excluded; no native wall-time or intrinsic-hardness conclusion.`,
    class:'accounting',anchor:'#'+anchor,source};
  series.points.push(point);
  rows+=`<tr><td>${point.date}</td><td>${esc(series.label)}</td><td>${esc(point.label)}</td><td class="n">${fmt(point.value)}&times; <small>simulated instruction diagnostic</small></td><td>matched</td><td><span class="chip">accounting</span></td><td><a href="${point.anchor}">panel</a> · <a href="${source}">source</a></td></tr>\n`;
}
timeline.updated=[timeline.updated,'2026-10-14'].sort().at(-1);
const json=JSON.stringify(timeline,null,2)+'\n';
write('docs/ic/progress-timeline.json',json);
page=page.replace(/(<script type="application\/json" id="progress-data">)[\s\S]*?(<\/script>)/,`$1\n${json}$2`);
page=page.replace(/(<summary>Every plotted point as a table<\/summary>[\s\S]*?)(<\/tbody>)/,`$1${rows}$2`);
write('docs/index-calculus-scoreboard.html',page);
