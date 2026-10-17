import assert from 'node:assert/strict';
import {render} from './render.mjs';
const node=(id,level,value=null,n=0)=>({id,level,curve_id:'exact-curve-'+id,value,n});
const edge=(source,target,verified=true)=>({source,target,degree:2,verified});
const metric={label:'Recorded work',unit:'group operations',statistic:'mean',boundary:'online query only'};
let checks=0;
for(const depth of [0,1,4,8,20]){
  const nodes=Array.from({length:depth+1},(_,i)=>node('C'+i,i));
  const edges=Array.from({length:depth},(_,i)=>[edge('C'+i,'C'+(i+1)),edge('C'+(i+1),'C'+i)]).flat();
  const r=render({run_id:'depth-'+depth,nodes,edges});
  assert.equal((r.svg.match(/class="node"/g)??[]).length,depth+1);
  assert.equal((r.svg.match(/class="edge"/g)??[]).length,depth*2);
  assert.equal(r.levels.length,depth+1);checks+=3;
}
const graph={run_id:'cycle-loop-parallel',nodes:[node('a',0,0,1),node('b',0,4,2),node('c',1)],edges:[edge('a','b'),edge('b','a'),edge('a','a'),edge('a','b',false),edge('b','c')]};
const a=render({...graph,metric,range:[0,4]});assert.equal(a.nodes,3);assert.equal(a.edges,5);assert.match(a.svg,/stroke-dasharray/);assert.match(a.svg,/url\(#missing\)/);assert.match(a.svg,/n=1/);checks+=5;
const constant=render({run_id:'constant',metric,nodes:[node('x',0,13,100),node('y',6,13,100)],edges:[]});assert.deepEqual(constant.range,[13,13]);assert.equal((constant.svg.match(/fill="rgb\(128,174,197\)"/g)??[]).length,2);checks+=2;
assert.throws(()=>render({...graph,metric,range:[1,4]}),/outside/);assert.throws(()=>render({...graph,nodes:[...graph.nodes,graph.nodes[0]]}),/unique/);assert.throws(()=>render({...graph,edges:[edge('a','unknown')]}),/endpoint/);checks+=3;
const unsafe=render({run_id:'escape',nodes:[node('<script>',0)],edges:[]});assert(!unsafe.svg.includes('<script>'));checks++;
console.log(JSON.stringify({passed:true,checks,depths:[0,1,4,8,20],surface_cycles:true,loops:true,parallel_maps:true,missing_and_zero_distinct:true,constant_and_shared_ranges:true}));
