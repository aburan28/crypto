// Presentation-only regression checks against the actual self-contained export.
// No graph generation, ECDLP, curve arithmetic or benchmark aggregation here.
import assert from 'node:assert/strict';
import fs from 'node:fs';
import vm from 'node:vm';

const html=fs.readFileSync(process.argv[2],'utf8');
const encoded=html.match(/<script id="graph-data" type="application\/json">([\s\S]*?)<\/script>/)[1];
const code=html.match(/<script>\n([\s\S]*?)<\/script>/)[1];
const data=JSON.parse(encoded, (key,value)=>{
  if(typeof value==='number')assert.ok(Number.isSafeInteger(value),'Unsafe numeric metadata: '+key);
  return value;
});
assert.deepEqual(data.graphs.map(g=>[g.nodes.length,g.edges.length,Object.keys(g.strata).length]),[[8192,16384,14],[6561,13122,9],[10070,39794,28]]);
const expectedHashes=['8ec54bf16456fd0cf6a62f8758a68425e0151157a44c11e6c2e00ddff67a8f64','a72ec162a221914f12e9c8806652e865b99565b7ddbcba860af0bae923de76b2','621d236a1a604ef391978aabb8b419f10db6daa54977831e37fb515d9bdc9d96'];
assert.deepEqual(data.graphs.map(g=>g.sha256),expectedHashes);
const ctx=new Proxy({}, {get:(_,k)=>()=>{}});
class Element {
  constructor(id){this.id=id;this.children=[];this.events={};this.value='';this.textContent='';this.checked=['edges','audited'].includes(id);}
  appendChild(e){this.children.push(e);if(this.children.length===1&&this.id==='graph')this.value=e.value;return e;}
  replaceChildren(...e){this.children=e;if(this.id==='strata')this.value='';}
  addEventListener(name,fn){this.events[name]=fn;}
  getBoundingClientRect(){return {width:1200,height:750,left:0,top:0};}
  getContext(){return ctx;}
}
const elements=new Map();
const get=id=>{if(!elements.has(id))elements.set(id,new Element(id));return elements.get(id);};
get('graph-data').textContent=encoded;
const document={getElementById:get,createElement:()=>new Element(),documentElement:{dataset:{}}};
vm.runInNewContext(code,{document,window:{devicePixelRatio:1,addEventListener(){}},location:{search:'?graph=2&vertex=10069&conductor=1728'},URLSearchParams,console});
assert.equal(document.documentElement.dataset.ready,'true');
assert.equal(get('audited-label').textContent,'Mark '+data.checkpoint_inventory.vertices_audited.toLocaleString()+' audited vertices');
assert.match(get('selection').textContent,/Vertex 10069/);
assert.equal(get('strata').value,'1728');
for(let i=0;i<data.graphs.length;i++){
  get('graph').value=String(i);get('graph').events.change();
  const g=data.graphs[i];
  assert.match(get('status').textContent,new RegExp(g.nodes.length.toLocaleString()));
  for(const conductor of Object.keys(g.strata)){get('strata').value=conductor;get('strata').events.change();}
  for(const id of [0,g.nodes.length-1]){
    get('vertex').value=String(id);get('go').events.click();
    const model=JSON.parse(get('model').textContent),n=g.nodes[id];
    assert.equal(model.p,g.p);assert.equal(model.a,n[1]);assert.equal(model.b,n[2]);assert.equal(model.j,n[3]);
    assert.equal(typeof model.p,'string');assert.equal(typeof model.a,'string');
    const maps=g.edges.filter(e=>e[0]===id);assert.equal(get('neighbors').children.length,maps.length);
    if(maps.length){get('neighbors').children[0].events.click();const m=JSON.parse(get('map').textContent);assert.equal(m.kernel_json,maps[0][4]);assert.equal(m.scaling_u,maps[0][3]);assert.equal(typeof m.scaling_u,'string');}
  }
  const unseen=g.nodes.find(n=>!data.checkpoint_inventory.capsules.some(c=>c.graph===g.name&&c.vertex===n[0]));
  if(unseen){get('vertex').value=String(unseen[0]);get('go').events.click();assert.match(get('coverage').textContent,/No audited measurements in this frozen/);}
  get('vertex').value=String(g.nodes.length);get('go').events.click();assert.match(get('selection').textContent,/out of range/);
  get('vertex').value='0';get('go').events.click();assert.match(get('coverage').textContent,/Audited/);
}
console.log('PASS: 3 full graph populations/hashes, all 51 conductor filters, exact large field strings, root/floor inspectors, complete outgoing kernel lists, deep link, bounds and checkpoint labels.');
