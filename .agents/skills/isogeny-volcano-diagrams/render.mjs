import fs from 'node:fs';
import path from 'node:path';
import crypto from 'node:crypto';
import {pathToFileURL} from 'node:url';

const need=(ok,msg)=>{if(!ok)throw Error(msg);};
const esc=s=>String(s).replaceAll('&','&amp;').replaceAll('<','&lt;').replaceAll('>','&gt;').replaceAll('"','&quot;').replaceAll("'",'&#39;');
const fmt=x=>Number(x).toLocaleString('en-US',{maximumSignificantDigits:7});
export function render(data){
  need(data&&Array.isArray(data.nodes)&&data.nodes.length>0,'nonempty nodes required');
  need(Array.isArray(data.edges),'edges required');
  need(typeof data.run_id==='string'&&data.run_id.length>0,'run_id required');
  const nodes=data.nodes,byId=new Map(),layers=new Map();
  for(const n of nodes){
    need(typeof n.id==='string'&&n.id.length>0&&!byId.has(n.id),'unique string node ids required');
    need(typeof n.curve_id==='string'&&n.curve_id.length>0,'full curve identity required');
    need(Number.isInteger(n.level)&&n.level>=0,'recorded nonnegative integer level required');
    if(data.metric)need((n.value===null||Number.isFinite(n.value))&&Number.isInteger(n.n)&&n.n>=0&&(n.value===null||n.n>0),'metric value and sample count invalid');
    byId.set(n.id,n);if(!layers.has(n.level))layers.set(n.level,[]);layers.get(n.level).push(n.id);
  }
  const neighbors=new Map(nodes.map(n=>[n.id,[]]));
  for(const e of data.edges){need(byId.has(e.source)&&byId.has(e.target),'edge endpoint missing');need(Number.isInteger(e.degree)&&e.degree>0,'positive degree required');need(typeof e.verified==='boolean','edge verification status required');neighbors.get(e.source).push(e.target);neighbors.get(e.target).push(e.source);}
  const levels=[...layers.keys()].sort((a,b)=>a-b),rank=new Map();
  const update=()=>{for(const l of levels)layers.get(l).forEach((id,i)=>rank.set(id,(i+.5)/layers.get(l).length));};
  for(const ids of layers.values())ids.sort();update();
  // Layered barycentric ordering handles surface cycles and arbitrary depth.
  // No recursive tree traversal: loops/cycles cannot drop or repeat vertices.
  for(let pass=0;pass<4;pass++)for(const direction of [1,-1])for(const l of direction===1?levels:[...levels].reverse()){
    const score=id=>{const ns=neighbors.get(id).filter(x=>direction*(byId.get(x).level-l)<0);return ns.length?ns.reduce((s,x)=>s+rank.get(x),0)/ns.length:rank.get(id);};
    layers.get(l).sort((a,b)=>score(a)-score(b)||a.localeCompare(b));update();
  }
  const width=Math.max(1100,Math.max(...[...layers.values()].map(x=>x.length))*75+120),height=245+levels.length*180;
  const pos=new Map();levels.forEach((l,j)=>layers.get(l).forEach(id=>pos.set(id,{x:60+rank.get(id)*(width-120),y:190+j*180})));
  let range=null;
  if(data.metric){need(typeof data.metric.unit==='string'&&typeof data.metric.boundary==='string'&&typeof data.metric.statistic==='string','metric units, statistic and boundary required');const vs=nodes.map(n=>n.value).filter(x=>x!==null);range=data.range??(vs.length?[Math.min(...vs),Math.max(...vs)]:null);if(range){need(range.length===2&&range.every(Number.isFinite)&&range[0]<=range[1],'invalid color range');need(vs.every(x=>x>=range[0]&&x<=range[1]),'measurement outside shared color range');}}
  const fill=n=>{if(!data.metric)return ['#b9dfed','#8bc8dd','#5bacc8','#3085aa','#155f87'][levels.indexOf(n.level)%5];if(n.value===null)return 'url(#missing)';const t=range[0]===range[1]?.5:(n.value-range[0])/(range[1]-range[0]);return `rgb(${Math.round(225-195*t)},${Math.round(241-135*t)},${Math.round(249-105*t)})`;};
  let svg=`<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${width} ${height}" role="img" aria-label="${esc(data.title??'Isogeny volcano')}"><defs><pattern id="missing" width="7" height="7" patternUnits="userSpaceOnUse"><rect width="7" height="7" fill="#eef2f5"/><path d="M0 7L7 0" stroke="#9aaab5"/></pattern><marker id="arrow" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="5" markerHeight="5" orient="auto"><path d="M0 0L10 5L0 10" fill="#647e90"/></marker></defs><rect width="${width}" height="${height}" fill="white"/><g font-family="sans-serif" fill="#183b50"><text x="25" y="34" font-size="23">${esc(data.title??'Isogeny volcano')}</text><text x="25" y="62" font-size="15">${esc(data.run_id)} · ${nodes.length} vertices · ${data.edges.length} directed maps · ${levels.length} recorded levels</text>`;
  const legend=data.metric?`${data.metric.label} · ${data.metric.statistic} · ${data.metric.unit} · ${range?fmt(range[0])+' to '+fmt(range[1]):'no measured values'}`:'Categorical level colors; no performance values';
  svg+=`<text x="25" y="89" font-size="15">${esc(legend)}</text><text x="25" y="113" font-size="13">${esc(data.metric?.boundary??'Exact topology from supplied records')} · Hatched = missing · Dashed map = unverified</text>`;
  for(const [j,l]of levels.entries())svg+=`<rect x="15" y="${140+j*180}" width="${width-30}" height="160" rx="10" fill="#f3f6f9"/><text x="30" y="${163+j*180}" font-size="14">Level ${l} · ${layers.get(l).length} vertices</text>`;
  const multiplicity=new Map();
  for(const [i,e]of data.edges.entries()){
    const a=pos.get(e.source),b=pos.get(e.target),key=[e.source,e.target].sort().join('\u0000'),index=multiplicity.get(key)??0;multiplicity.set(key,index+1);
    let d,lx,ly;
    if(e.source===e.target){const r=28+index*13;d=`M${a.x-10} ${a.y-14}C${a.x-r*2} ${a.y-r*2},${a.x+r*2} ${a.y-r*2},${a.x+10} ${a.y-14}`;lx=a.x;ly=a.y-r*1.5;}
    else{const dx=b.x-a.x,dy=b.y-a.y,len=Math.hypot(dx,dy),ux=dx/len,uy=dy/len,bend=8+Math.floor(index/2)*12;const cx=(a.x+b.x)/2-uy*bend,cy=(a.y+b.y)/2+ux*bend;d=`M${a.x+ux*15} ${a.y+uy*15}Q${cx} ${cy},${b.x-ux*17} ${b.y-uy*17}`;lx=cx+4;ly=cy-4;}
    svg+=`<g class="edge" data-edge="${esc(e.id??String(i))}"><path d="${d}" fill="none" stroke="#8ca4b5" stroke-width="1.1"${e.verified?'':' stroke-dasharray="4 3"'} marker-end="url(#arrow)"/><text x="${lx}" y="${ly}" font-size="10">${e.degree}</text><title>${esc(JSON.stringify(e))}</title></g>`;
  }
  for(const n of nodes){const p=pos.get(n.id);svg+=`<g class="node" data-node="${esc(n.id)}" data-level="${n.level}" data-curve-id="${esc(n.curve_id)}" tabindex="0"><circle cx="${p.x}" cy="${p.y}" r="14" fill="${fill(n)}" stroke="#254f66"/><text x="${p.x}" y="${p.y+32}" text-anchor="middle" font-size="20">${esc(n.id)}</text><text x="${p.x}" y="${p.y+47}" text-anchor="middle" font-family="monospace" font-size="9">${esc(n.curve_id.slice(0,12))}</text>${data.metric?`<text x="${p.x}" y="${p.y+64}" text-anchor="middle" font-size="13">${n.value===null?'missing':'n='+n.n}</text>`:''}<title>${esc(JSON.stringify(n))}</title></g>`;}
  svg+=`<text x="25" y="${height-22}" font-size="12">Full exact identities, map evidence and measurements: graph_data.json / interactive inspector.</text></g></svg>`;
  return {svg,range,width,height,nodes:nodes.length,edges:data.edges.length,levels};
}

export function write(data,out){need(!fs.existsSync(out),'refuse output overwrite');const result=render(data);fs.mkdirSync(out,{recursive:true});const graph=JSON.stringify(data,null,2)+'\n';fs.writeFileSync(path.join(out,'graph_data.json'),graph);fs.writeFileSync(path.join(out,'volcano.svg'),result.svg);
  const html=`<!doctype html><html><meta charset="utf-8"><title>${esc(data.title??'Isogeny volcano')}</title><style>@page{size:A4 landscape;margin:12mm}body{font:14px/1.5 system-ui;color:#183b50;margin:25px}h1{font-size:25px}svg{width:100%;height:auto}.figure{break-before:page}.node{cursor:pointer}pre{white-space:pre-wrap;overflow-wrap:anywhere;font-size:11px}@media print{.inspector{display:none}}</style><body><h1>${esc(data.title??'Isogeny volcano')}</h1><p>${esc(data.caption??'')}</p><p>Run: <code>${esc(data.run_id)}</code>. Exact data: graph_data.json. ${data.metric?esc(data.metric.label+'; '+data.metric.unit+'; '+data.metric.statistic+'; '+data.metric.boundary):'Colors represent topology levels only.'}</p><section class="figure">${result.svg}</section><section class="inspector"><h2>Exact curve and evidence</h2><pre id="inspect">Select a vertex.</pre></section><script>const DATA=${JSON.stringify(data).replaceAll('<','\\u003c')};for(const element of document.querySelectorAll('.node')){const show=()=>{document.getElementById('inspect').textContent=JSON.stringify(DATA.nodes.find(n=>n.id===element.dataset.node),null,2);};element.addEventListener('click',show);element.addEventListener('keydown',e=>{if(e.key==='Enter')show();});}</script></body></html>`;fs.writeFileSync(path.join(out,'index.html'),html);
  const manifest={schema:'isogeny-volcano-diagram/v1',run_id:data.run_id,node_count:result.nodes,directed_edge_count:result.edges,recorded_levels:result.levels,color_range:result.range,source_sha256:crypto.createHash('sha256').update(graph).digest('hex'),missing_nodes:data.metric?data.nodes.filter(n=>n.value===null).length:0,all_vertices_and_maps_retained:true,pdf_inspection_required:true};fs.writeFileSync(path.join(out,'render_manifest.json'),JSON.stringify(manifest,null,2)+'\n');return manifest;
}
if(process.argv[1]&&import.meta.url===pathToFileURL(process.argv[1]).href){need(process.argv.length===4,'usage: node render.mjs INPUT.json NEW_OUTPUT_DIRECTORY');console.log(JSON.stringify(write(JSON.parse(fs.readFileSync(process.argv[2])),process.argv[3])));}
