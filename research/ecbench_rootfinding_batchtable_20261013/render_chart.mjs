// Presentation of sealed ecbench phase means and source-derived inversion
// counts. No benchmark or algorithm is run by this file.
import fs from 'node:fs';

const rows = JSON.parse(fs.readFileSync(new URL('./phase-metrics.json', import.meta.url), 'utf8'));
const scalar = rows.find((r) => r.arm === 'pair-scalar');
const batch = rows.find((r) => r.arm === 'pair-batch');
if (!scalar || !batch) throw new Error('both measured arms are required');
const esc = (x) => String(x).replaceAll('&', '&amp;').replaceAll('<', '&lt;').replaceAll('>', '&gt;');
const makePanel = ({x, title, unit, max, values, format}) => {
  const top = 105, bottom = 385, scale = (bottom - top) / max;
  const grid = [0, .25, .5, .75, 1].map((t) => {
    const y = bottom - t * (bottom - top);
    return `<line x1="${x+58}" y1="${y}" x2="${x+355}" y2="${y}" class="grid"/>
      <text x="${x+50}" y="${y+5}" text-anchor="end" class="tick">${format(max*t)}</text>`;
  }).join('');
  const bars = values.map((v, i) => {
    const bx = x + 125 + i * 125, h = v.value * scale;
    return `<rect x="${bx}" y="${bottom-h}" width="75" height="${h}" fill="${v.color}"/>
      <text x="${bx+37.5}" y="${Math.max(top+18,bottom-h-8)}" text-anchor="middle" class="value">${format(v.value)}</text>
      <text x="${bx+37.5}" y="413" text-anchor="middle" class="label">${esc(v.name)}</text>`;
  }).join('');
  return `<text x="${x+35}" y="65" class="panel">${esc(title)}</text>
    <text x="${x+35}" y="84" class="sub">${esc(unit)}</text>${grid}${bars}`;
};
const svg = `<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg" width="960" height="540" viewBox="0 0 960 540">
<style>text{font-family:Arial,Helvetica,sans-serif;fill:#1c2736}.title{font-size:22px;font-weight:700}.panel{font-size:18px;font-weight:700}.sub{font-size:14px}.tick{font-size:13px}.label{font-size:15px;font-weight:600}.value{font-size:16px;font-weight:700}.grid{stroke:#dbe2ea;stroke-width:1}</style>
<rect width="960" height="540" fill="white"/>
<text x="30" y="30" class="title">Prime pair-table setup: same keys, fewer inversions</text>
${makePanel({x:0, title:'Field inversions (derived)', unit:'B=64; counts from the two builder loops', max:4500,
  values:[{name:'scalar',value:scalar.setup_inversions_derived,color:'#5b7db9'},
          {name:'batch',value:batch.setup_inversions_derived,color:'#18a999'}],
  format:(v)=>Math.round(v).toLocaleString()})}
${makePanel({x:480, title:'Mean setup wall (L0)', unit:'microseconds; 24 measured runs per arm', max:280,
  values:[{name:'scalar',value:scalar.mean_setup_wall_ns/1000,color:'#5b7db9'},
          {name:'batch',value:batch.mean_setup_wall_ns/1000,color:'#18a999'}],
  format:(v)=>v.toFixed(0)})}
<text x="30" y="455" class="sub">Cold total: ${scalar.mean_total_gae.toFixed(3)} GAE in each arm; 32/32 paired outputs and counted phases identical.</text>
<text x="30" y="480" class="sub">ECBS1h4d0b9e4a02a7 · source-derived inversion counts; setup wall is descriptive only.</text>
<text x="30" y="505" class="sub">macOS L0; prime native work remains partly unpriced; this figure does not establish a whole-solve runtime gain.</text>
</svg>`;
fs.writeFileSync(new URL('./setup.svg', import.meta.url), svg);
