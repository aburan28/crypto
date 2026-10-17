// Aggregate existing immutable measurements; never reads scalar witnesses.
import fs from 'node:fs';
import readline from 'node:readline';
const root=process.argv[2],out=process.argv[3];
if(fs.existsSync(out))throw Error('write-once output');
fs.mkdirSync(out,{recursive:true});
const read=p=>JSON.parse(fs.readFileSync(p));
const stat=a=>a.length?{n:a.length,min:Math.min(...a),mean:a.reduce((x,y)=>x+y,0)/a.length,max:Math.max(...a)}:{n:0,min:null,mean:null,max:null};
const problem=Array.from({length:100},(_,instance)=>({instance,curves:[],metrics:{}}));
const icCurve=[];
for(let curve=0;curve<46;curve++){
  const rows=[];
  for await(const line of readline.createInterface({input:fs.createReadStream(root+'/runs/ic/curve_'+curve+'.jsonl'),crlfDelay:Infinity})){
    const r=JSON.parse(line);if(r.stage!=='target')continue;
    const metrics={gb_wall_ns:r.gb_wall_ns,gb_process_cpu_ns:r.gb_process_cpu_ns,construction_ns:r.construction_ns,root_and_relation_ns:r.root_and_relation_ns,whole_probe_wall_ns:r.whole_probe_wall_ns,dreg:r.dreg_top_forms_first_zero_hilbert_degree,...Object.fromEntries(Object.entries(r.gb_counters).map(([k,v])=>['gb.'+k,v]))};
    // Targets are emitted generator-major with the recorded key index.
    const id=r.generator*10+r.key_index;
    if(!Number.isInteger(id)||id<0||id>=100)throw Error('target identity');
    const p=problem[id];p.curves.push(curve);for(const [k,v]of Object.entries(metrics)){if(!Number.isFinite(v))throw Error('nonfinite metric');(p.metrics[k]??=[]).push(v);}
    rows.push(metrics);
  }
  if(rows.length!==100)throw Error('100 actual target probes required');
  icCurve.push({curve,metrics:Object.fromEntries(Object.keys(rows[0]).map(k=>[k,stat(rows.map(r=>r[k]))]))});
}
for(const p of problem){if(p.curves.length!==46||new Set(p.curves).size!==46)throw Error('paired curve coverage');for(const k of Object.keys(p.metrics))p.metrics[k]=stat(p.metrics[k]);}
fs.writeFileSync(out+'/ic_per_problem.json',JSON.stringify(problem,null,2)+'\n');
fs.writeFileSync(out+'/ic_per_curve.json',JSON.stringify(icCurve,null,2)+'\n');
const csv=['instance,metric,n,min,mean,max'];for(const p of problem)for(const[k,v]of Object.entries(p.metrics))csv.push([p.instance,k,v.n,v.min,v.mean,v.max].join(','));fs.writeFileSync(out+'/ic_per_problem.csv',csv.join('\n')+'\n');
const edgeRows=new Map(),pointRows=new Map();let maps=0;
for await(const line of readline.createInterface({input:fs.createReadStream(root+'/data/map_measurements.jsonl'),crlfDelay:Infinity})){
  const r=JSON.parse(line);const metrics={wall_ns:r.pair_wall_ns,...r.counters};
  for(const [map,key]of [[edgeRows,r.edge],[pointRows,r.instance]]){if(!map.has(key))map.set(key,{});const obj=map.get(key);for(const[k,v]of Object.entries(metrics)){if(!Number.isFinite(v))throw Error('map metric');(obj[k]??=[]).push(v);}}maps++;
}
const convert=m=>[...m].sort((a,b)=>a[0]-b[0]).map(([id,metrics])=>({id,metrics:Object.fromEntries(Object.entries(metrics).map(([k,v])=>[k,stat(v)]))}));
fs.writeFileSync(out+'/map_per_edge.json',JSON.stringify(convert(edgeRows),null,2)+'\n');fs.writeFileSync(out+'/map_per_problem.json',JSON.stringify(convert(pointRows),null,2)+'\n');
const summary={ic_target_probes:4600,ic_problem_summaries:100,ic_curves:46,map_records:maps,map_edges:edgeRows.size,map_instances:pointRows.size,scalar_witness_files_read:false,setup_cost_metadata_not_summed:true};
if(summary.map_records!==45000||summary.map_edges!==90||summary.map_instances!==100)throw Error('map coverage');fs.writeFileSync(out+'/summary.json',JSON.stringify(summary,null,2)+'\n');console.log(JSON.stringify(summary));
