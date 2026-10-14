// Presentation only: render the sealed session's derived phase means.
import fs from 'node:fs';

const rows = JSON.parse(fs.readFileSync(new URL('./summary.json', import.meta.url), 'utf8'));
const order = ['pair-b32', 'pair-b64', 'rho-neg'];
const byArm = new Map(rows.map((r) => [r.arm, r]));
const color = { oracle_setup: '#5b7db9', relations: '#18a999', verify: '#ef9f35',
  factor_base: '#8b6bb3', linear_algebra: '#c36ba4', search: '#18a999', setup: '#5b7db9',
  internal_verification: '#ef9f35' };
const ymax = 25000;
const top = 55, bottom = 390, scale = (bottom - top) / ymax;
const esc = (s) => String(s).replaceAll('&', '&amp;').replaceAll('<', '&lt;').replaceAll('>', '&gt;');
const bars = order.map((arm, idx) => {
  const row = byArm.get(arm);
  const x = 155 + idx * 220;
  let y = bottom;
  const parts = row.phases.filter((p) => p.mean_gae > 0).map((p) => {
    const h = p.mean_gae * scale;
    y -= h;
    return `<rect x="${x}" y="${y.toFixed(2)}" width="110" height="${h.toFixed(2)}" fill="${color[p.name]}"/>`;
  }).join('');
  return `${parts}<text x="${x+55}" y="${(y-10).toFixed(2)}" text-anchor="middle" class="value">${Math.round(row.mean_gae).toLocaleString()}</text>
    <text x="${x+55}" y="416" text-anchor="middle" class="label">${esc(arm)}</text>`;
}).join('');
const grid = [0, 5000, 10000, 15000, 20000, 25000].map((v) => {
  const y = bottom - v * scale;
  return `<line x1="90" y1="${y}" x2="810" y2="${y}" class="grid"/>
    <text x="75" y="${y+5}" text-anchor="end" class="tick">${(v/1000).toFixed(0)}k</text>`;
}).join('');
const svg = `<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg" width="900" height="520" viewBox="0 0 900 520">
<style>text{font-family:Arial,Helvetica,sans-serif;fill:#1c2736}.title{font-size:22px;font-weight:700}.sub{font-size:14px}.tick{font-size:13px}.label{font-size:16px;font-weight:600}.value{font-size:16px;font-weight:700}.grid{stroke:#dbe2ea;stroke-width:1}</style>
<rect width="900" height="520" fill="white"/>
<text x="40" y="30" class="title">Cold counted cost by phase</text>
<text x="40" y="50" class="sub">20-bit prime subgroup, 8 public targets × 3 measured rounds per arm; ecbench.gae</text>
${grid}${bars}
<rect x="110" y="449" width="13" height="13" fill="${color.oracle_setup}"/><text x="130" y="460" class="tick">table or rho setup</text>
<rect x="300" y="449" width="13" height="13" fill="${color.relations}"/><text x="320" y="460" class="tick">relations or rho search</text>
<rect x="520" y="449" width="13" height="13" fill="${color.verify}"/><text x="540" y="460" class="tick">verification</text>
<text x="40" y="492" class="sub">ECBS1h382073f2a6ad · means from sealed records; prime native work remains unpriced; L0 wall</text>
</svg>`;
fs.writeFileSync(new URL('./cost.svg', import.meta.url), svg);
