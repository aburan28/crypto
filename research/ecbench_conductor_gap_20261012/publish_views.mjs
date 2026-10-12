// Static presentation of native-analyzer values; no solver or statistical computation.
import fs from 'node:fs';
import path from 'node:path';
const study=path.dirname(new URL(import.meta.url).pathname);
const repo=path.resolve(study,'../..');
const read=p=>fs.readFileSync(path.join(repo,p),'utf8');
const write=(p,s)=>fs.writeFileSync(path.join(repo,p),s);
const a=JSON.parse(read('research/ecbench_conductor_gap_20261012/analysis.json'));
const esc=s=>String(s).replaceAll('&','&amp;').replaceAll('<','&lt;').replaceAll('>','&gt;').replaceAll('"','&quot;');
const fmt=x=>Number.isFinite(x)?x.toPrecision(6):'missing';
const interval=x=>Array.isArray(x)?x.map(fmt).join('–'):'missing';
const anchor='conductor-gap-ic37-20261012';
const source='https://github.com/aburan28/crypto/blob/main/research/ecbench_conductor_gap_20261012/README.md';
let table='';
for(const s of a.sessions) {
  for(const arm of ['rho-neg','rho-neg-aa','ic-mitm']) {
    for(const r of s.arms[arm].curves) {
      table+=`<tr><td>${s.ell}</td><td><code>${esc(r.slug)}</code> (${esc(r.role)})</td><td>${esc(arm)}</td><td>${fmt(r.mean_s)}</td><td>${r.measured}</td><td>${r.lower_bound?'bounded':'priced'}</td></tr>\n`;
    }
  }
}
const comparisons=a.sessions.map(s=>`ell = ${s.ell}: IC floor/crater ${fmt(s.arms['ic-mitm'].cross_curve.ratio_floor_over_crater)}, 99% interval ${interval(s.arms['ic-mitm'].cross_curve.ci99)}`).join('; ');
const panel=`<section class="panel" id="${anchor}">
<h2>Large-prime conductor gap: IC measured at ell = 131101 <span class="chip">accounting</span></h2>
<p class="panel-summary">Native MITM-based IC on two fixed D = -19 CM pairs, eight independent targets per curve; both 99% floor/crater count intervals include one, and every measured execution replayed exactly.</p>
<p>${esc(comparisons)}. No repeated conductor lead under the frozen rule. Cross-curve targets are independent, with no isogeny transport. All 192 records verified, all 96 measured executions replayed identically; A/A counted costs match exactly.</p>
<div class="scroll"><table><caption>Cold counted S = GAE / sqrt(r), eight measured one-target workloads per row. Unpriced native work makes these lower bounds; no complete-cost or wall advantage is inferred.</caption><thead><tr><th>ell</th><th>curve / role</th><th>arm</th><th>mean counted S</th><th>workloads</th><th>scope</th></tr></thead><tbody>${table}</tbody></table></div>
<p>Coverage: negation rho and MITM-based IC. Large-gap algebraic F4/SAT comparisons, independent-host validation and representation controls remain open. <a href="${source}">Report, graphs, PDFs, protocol and evidence</a>.</p>
</section>\n`;
let page=read('docs/index-calculus-scoreboard.html');
if(page.includes(`id="${anchor}"`)) throw new Error('panel already exists');
page=page.replace('</main>',panel+'</main>');
if(!page.includes(`id="${anchor}"`)) throw new Error('main marker not found');
const timeline=JSON.parse(read('docs/ic/progress-timeline.json'));
const series=timeline.series.find(s=>s.id==='ec-bounded-counts');
if(!series) throw new Error('bounded count series missing');
let rows='';
for(const s of a.sessions) {
  const point={date:'2026-10-12',value:s.ic_vs_rho.floor.ic_over_rho_count_diagnostic,
    reference:'matched',unit:'counted GAE lower-bound diagnostic',
    label:`Large-prime conductor-gap IC, ell = ${s.ell}, fixed floor, eight one-target workloads`,
    detail:`Mean cold counted S quotient, with unpriced native work. ${s.session_id}; all measured runs exactly replayed. No full-cost, scaling, transport or intrinsic-hardness claim.`,
    class:'accounting',anchor:'#'+anchor,source};
  series.points.push(point);
  rows+=`<tr><td>${point.date}</td><td>${esc(series.label)}</td><td>${esc(point.label)}</td><td class="n">${fmt(point.value)}&times; <small>bounded count diagnostic</small></td><td>matched</td><td><span class="chip">accounting</span></td><td><a href="${point.anchor}">panel</a> · <a href="${source}">source</a></td></tr>\n`;
}
timeline.updated='2026-10-12';
const json=JSON.stringify(timeline,null,2)+'\n';
write('docs/ic/progress-timeline.json',json);
page=page.replace(/(<script type="application\/json" id="progress-data">)[\s\S]*?(<\/script>)/,`$1\n${json}$2`);
page=page.replace(/(<summary>Every plotted point as a table<\/summary>[\s\S]*?)(<\/tbody>)/,`$1${rows}$2`);
write('docs/index-calculus-scoreboard.html',page);
