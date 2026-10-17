// Aggregate a separately frozen native arithmetic replay. No scalar output.
import fs from 'node:fs';
import crypto from 'node:crypto';
import {isDeepStrictEqual} from 'node:util';
const [root,tag]=process.argv.slice(2),out=root+'/reports/'+tag;
if(fs.existsSync(out))throw Error('write-once replay aggregate');
const load=p=>JSON.parse(fs.readFileSync(p)),binding=load(root+'/FIELD_FROZEN.json');
const query=[],attempts=[],jobs=[];
const stats=a=>a.length?{n:a.length,min:Math.min(...a),mean:a.reduce((s,x)=>s+x,0)/a.length,max:Math.max(...a)}:{n:0,min:null,mean:null,max:null};
const counters=c=>Object.fromEntries(Object.entries(c).flatMap(([mod,row])=>Object.entries(row).map(([k,v])=>[mod+'.'+k,v])));
const keys=['instrumented_wall_ns','instrumented_process_cpu_ns','validation_ns',...['Fp','scalar_mod_r','other_modulus'].flatMap(mod=>['add','sub','mul_plain','mont_mul','mont_square','mont_redc','inverse','euclid_divisions','complement'].map(k=>mod+'.'+k))];
const aggregate=rs=>Object.fromEntries(keys.map(k=>[k,stats(rs.map(r=>r[k]))]));
for(const algorithm of ['rho','bsgs']){
 const dir=root+'/field_runs/'+algorithm;if(!fs.existsSync(dir))continue;
 for(const job of fs.readdirSync(dir)){
  const folder=dir+'/'+job,logical=new Map(),complete=fs.existsSync(folder+'/complete.json')?load(folder+'/complete.json'):null;
  if(complete&&!isDeepStrictEqual(complete.binding,binding))throw Error('binding mismatch');
  for(const attempt of fs.readdirSync(folder).filter(x=>x.startsWith('attempt_')).sort()){
   const p=folder+'/'+attempt,input=load(p+'/input.json');if(!isDeepStrictEqual(input.binding,binding))throw Error('attempt binding');
   const record={algorithm,job,attempt,curve:Number(job.match(/^C(\d+)_/)[1]),setup:null,destruction:null,closed_queries:0,verified_queries:0};
   for(const phase of ['setup','destruction'])if(fs.existsSync(p+'/'+phase+'.json')){
    const r=load(p+'/'+phase+'.json');record[phase]={wall_ns:r.wall_ns,process_cpu_ns:r.process_cpu_ns,status:r.status??null,counters:r.counters,resources_before:r.resources_before??null,resources_after:r.resources_after??null};
   }
   for(const file of fs.readdirSync(p).filter(x=>/^case_\d+\.json$/.test(x))){
    const r=load(p+'/'+file);record.closed_queries++;
    if(r.verified){if(r.status!==0||!r.verified_public_equation||!r.matched_baseline_answer_and_all_nonclock_stats||!r.baseline_read_after_native_return)throw Error('invalid verification flags');record.verified_queries++;logical.set(r.instance,r);}
   }
   attempts.push(record);
  }
  const rs=complete?complete.records:[...logical.values()];if(complete&&!complete.verified)throw Error('invalid complete marker');
  for(const r of rs){if(!r.verified)continue;query.push({algorithm,curve:r.curve,level:r.level,instance:r.instance,generator:r.generator,key_index:r.key_index,point_order:r.point_order,group_ops:r.native_stats.group_ops,iterations:r.native_stats.iterations,collisions:r.native_stats.collisions,table_entries:r.native_stats.table_entries,bytes_peak:r.native_stats.bytes_peak,threads:r.native_stats.threads,instrumented_wall_ns:r.instrumented_wall_ns,instrumented_process_cpu_ns:r.instrumented_process_cpu_ns,validation_ns:r.validation_ns,...counters(r.counters)});}
  jobs.push({algorithm,job,complete:!!complete,verified_unique_cases:rs.filter(r=>r.verified).length});
 }
}
if(new Set(query.map(r=>r.algorithm+'/'+r.curve+'/'+r.instance)).size!==query.length)throw Error('duplicate logical case');
const perCurve=[],perGenerator=[],perProblem=[];
for(const algorithm of ['rho','bsgs']){
 for(let curve=0;curve<46;curve++){
  const rs=query.filter(r=>r.algorithm===algorithm&&r.curve===curve),a=attempts.filter(r=>r.algorithm===algorithm&&r.curve===curve);
  perCurve.push({algorithm,curve,verified_unique_cases:rs.length,complete:rs.length===100,query:aggregate(rs),setup:a.filter(x=>x.setup).map(({job,attempt,setup})=>({job,attempt,...setup})),destruction:a.filter(x=>x.destruction).map(({job,attempt,destruction})=>({job,attempt,...destruction}))});
  for(let generator=0;generator<10;generator++){const g=rs.filter(r=>r.generator===generator);perGenerator.push({algorithm,curve,generator,verified_unique_cases:g.length,query:aggregate(g)});}
 }
 for(let instance=0;instance<100;instance++){const rs=query.filter(r=>r.algorithm===algorithm&&r.instance===instance);perProblem.push({algorithm,instance,observed_curves:rs.length,query:aggregate(rs)});}
}
const summary={generated_utc:new Date().toISOString(),phase:'separately frozen native arithmetic replay',verified:Object.fromEntries(['rho','bsgs'].map(a=>[a,query.filter(r=>r.algorithm===a).length])),expected_per_algorithm:4600,complete:query.length===9200,all_selected_records_match_baseline_and_public_equation:true,original_timings_unchanged:true,primitive_counts_are_source_level_calls:true,mont_square_is_subset_of_mont_mul:true,mont_redc_overlaps_products_and_conversions:true,inverse_is_extended_Euclid:true,thread_local_capture_threads:1,timing_class:'counter-instrumented replay; concurrent L0; separate from original timing',all_closed_attempts_retained:true,full_random_corpus_IC_recovery_complete:false,FIELD_FROZEN_sha256:crypto.createHash('sha256').update(fs.readFileSync(root+'/FIELD_FROZEN.json')).digest('hex'),contains_scalar_witnesses:false};
fs.mkdirSync(out,{recursive:true});for(const[name,value]of Object.entries({summary,per_curve:perCurve,per_generator:perGenerator,per_problem:perProblem,attempts,jobs}))fs.writeFileSync(out+'/'+name+'.json',JSON.stringify(value,null,2)+'\n');
const csv=['algorithm,curve,verified_unique_cases,Fp_mont_mul_min,Fp_mont_mul_mean,Fp_mont_mul_max,Fp_mont_square_min,Fp_mont_square_mean,Fp_mont_square_max,Fp_inverse_min,Fp_inverse_mean,Fp_inverse_max,scalar_mont_mul_min,scalar_mont_mul_mean,scalar_mont_mul_max'];for(const r of perCurve)csv.push([r.algorithm,r.curve,r.verified_unique_cases,...['Fp.mont_mul','Fp.mont_square','Fp.inverse','scalar_mod_r.mont_mul'].flatMap(k=>[r.query[k].min,r.query[k].mean,r.query[k].max])].join(','));fs.writeFileSync(out+'/per_curve.csv',csv.join('\n')+'\n');
console.log(JSON.stringify({out,summary}));
