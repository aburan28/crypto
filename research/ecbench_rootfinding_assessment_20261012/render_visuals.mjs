import fs from 'node:fs';
import path from 'node:path';

const here = path.dirname(new URL(import.meta.url).pathname);
const prime = JSON.parse(fs.readFileSync(path.join(here, 'summary-prime20.json'), 'utf8'));
const binary = JSON.parse(fs.readFileSync(path.join(here, 'summary-koblitz17.json'), 'utf8'));
function row(rows, arm) {
  const r = rows.find(x => x.arm === arm);
  if (!r || r.verified !== 24 || !Number.isFinite(r.mean_s) || r.mean_s <= 0) throw new Error(`missing verified ${arm}`);
  return r;
}
const bars = [
  {curve:'Prime 20-bit',arm:'rho',x:150,v:row(prime,'rho-neg').mean_s,color:'#667085',bounded:true},
  {curve:'Prime 20-bit',arm:'subtract',x:275,v:row(prime,'ic-subtract').mean_s,color:'#2563eb',bounded:true},
  {curve:'Prime 20-bit',arm:'pair table',x:400,v:row(prime,'ic-pair-table').mean_s,color:'#0f766e',bounded:true},
  {curve:'Koblitz n=17',arm:'rho',x:660,v:row(binary,'rho-strong').mean_s,color:'#667085',bounded:true},
  {curve:'Koblitz n=17',arm:'subtract',x:785,v:row(binary,'ic-subtract').mean_s,color:'#2563eb',bounded:false},
  {curve:'Koblitz n=17',arm:'pair table',x:910,v:row(binary,'ic-pair-table').mean_s,color:'#0f766e',bounded:false}
];
const esc = s => String(s).replaceAll('&','&amp;').replaceAll('<','&lt;').replaceAll('>','&gt;');
const y0=380, width=72;
const grid=[1,10,100,1000].map(v=>`<line x1="90" y1="${y0-85*Math.log10(v)}" x2="1015" y2="${y0-85*Math.log10(v)}" stroke="#d0d5dd" stroke-dasharray="4 5"/><text x="75" y="${y0-85*Math.log10(v)+5}" text-anchor="end" font-size="15" fill="#475467">${v}</text>`).join('');
const rects=bars.map(b=>{const h=85*Math.log10(b.v), y=y0-h;return `<rect x="${b.x-width/2}" y="${y.toFixed(2)}" width="${width}" height="${h.toFixed(2)}" rx="3" fill="${b.color}"/><text x="${b.x}" y="${(y-12).toFixed(2)}" text-anchor="middle" font-size="17" font-weight="700" fill="#101828">${b.v.toFixed(2)}${b.bounded?'†':''}</text><text x="${b.x}" y="409" text-anchor="middle" font-size="15" fill="#344054">${esc(b.arm)}</text>`}).join('');
const svg=`<svg xmlns="http://www.w3.org/2000/svg" width="1080" height="505" viewBox="0 0 1080 505"><rect width="1080" height="505" fill="white"/><text x="40" y="38" font-family="Arial, sans-serif" font-size="25" font-weight="700" fill="#101828">Cold one-target counted cost</text><text x="40" y="64" font-family="Arial, sans-serif" font-size="15" fill="#475467">Mean S = charged group-addition equivalents / sqrt(r); logarithmic vertical scale</text><g font-family="Arial, sans-serif">${grid}<line x1="90" y1="380" x2="1015" y2="380" stroke="#667085" stroke-width="2"/>${rects}<line x1="80" y1="430" x2="520" y2="430" stroke="#d0d5dd"/><line x1="560" y1="430" x2="1025" y2="430" stroke="#d0d5dd"/><text x="300" y="452" text-anchor="middle" font-size="17" font-weight="700" fill="#101828">Prime 20-bit, r=562427</text><text x="790" y="452" text-anchor="middle" font-size="17" font-weight="700" fill="#101828">Koblitz n=17, r=65587</text><text x="40" y="486" font-size="14" fill="#475467">24 verified runs per arm; † unpriced native work; all local wall grades L0.</text></g></svg>`;
fs.writeFileSync(path.join(here,'mean_s.svg'),svg);
const pipe=`<svg xmlns="http://www.w3.org/2000/svg" width="1120" height="300" viewBox="0 0 1120 300"><rect width="1120" height="300" fill="white"/><text x="28" y="38" font-family="Arial, sans-serif" font-size="24" font-weight="700" fill="#101828">Where the replacement acts</text><g font-family="Arial, sans-serif" font-size="16" text-anchor="middle"><rect x="25" y="90" width="160" height="74" rx="8" fill="#dbeafe" stroke="#2563eb"/><text x="105" y="120">Build factor base</text><text x="105" y="144" font-size="13">reusable setup</text><rect x="230" y="90" width="180" height="74" rx="8" fill="#dbeafe" stroke="#2563eb"/><text x="320" y="120">Build pair table</text><text x="320" y="144" font-size="13">optional, O(N²)</text><rect x="455" y="90" width="160" height="74" rx="8" fill="#ccfbf1" stroke="#0f766e"/><text x="535" y="120">Target queries</text><text x="535" y="144" font-size="13">relation trials</text><rect x="660" y="90" width="160" height="74" rx="8" fill="#ccfbf1" stroke="#0f766e"/><text x="740" y="120">Decompose</text><text x="740" y="144" font-size="13">subtract / probe</text><rect x="865" y="90" width="220" height="74" rx="8" fill="#ccfbf1" stroke="#0f766e"/><text x="975" y="120">Matrix + verification</text><text x="975" y="144" font-size="13">recovered scalar checked</text></g><g fill="none" stroke="#667085" stroke-width="3"><path d="M185 127H230M410 127H455M615 127H660M820 127H865"/></g><g fill="#667085"><path d="M226 121l8 6-8 6M451 121l8 6-8 6M656 121l8 6-8 6M861 121l8 6-8 6"/></g><text x="28" y="215" font-family="Arial, sans-serif" font-size="16" fill="#344054">Cold cost charges every box once per target. Reusable setup can be amortized only over a declared number of targets.</text><text x="28" y="248" font-family="Arial, sans-serif" font-size="16" fill="#344054">The paper changes the decomposition box; relation yield and matrix work do not disappear.</text></svg>`;
fs.writeFileSync(path.join(here,'pipeline.svg'),pipe);
