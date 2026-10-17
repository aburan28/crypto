// Local observer only. Does not launch solvers or export data.
import fs from 'node:fs';
import {spawnSync} from 'node:child_process';
const root=process.argv[2],lock=root+'/field_report_observer.lock';
if(fs.existsSync(lock)){const old=JSON.parse(fs.readFileSync(lock));try{process.kill(old.pid,0);throw Error('observer already active');}catch(e){if(e.code!=='ESRCH')throw e;}}
fs.writeFileSync(lock,JSON.stringify({pid:process.pid,started_utc:new Date().toISOString()})+'\n');
const log=root+'/logs/field_report_observer_'+Date.now()+'_'+process.pid+'.jsonl';
const emit=r=>fs.appendFileSync(log,JSON.stringify({utc:new Date().toISOString(),...r})+'\n');
function run(script,args){const r=spawnSync(process.execPath,[root+'/source/'+script,...args],{encoding:'utf8',maxBuffer:16*1024*1024});emit({script,args,exit_status:r.status,stdout:r.stdout,stderr:r.stderr});if(r.status!==0)throw Error(script+' failed');}
function finish(){if(fs.existsSync(lock)&&JSON.parse(fs.readFileSync(lock)).pid===process.pid)fs.unlinkSync(lock);}
function tick(){try{
 const state=JSON.parse(fs.readFileSync(root+'/FIELD_LIVE_STATUS.json'));
 if(state.verified.rho===4600&&state.verified.bsgs===4600&&state.active_replay_worker_pids.length===0){
  const tag='field_complete_'+Date.now();run('field_aggregate.mjs',[root,tag]);run('phase_report.mjs',[root,tag+'_visuals','field',tag]);
  const result={local_aggregate:tag,local_visuals:tag+'_visuals',arithmetic_replay_complete:true,PDF_render_and_page_inspection_pending:true,PR_update_pending:true,exports_performed:false,full_random_corpus_IC_complete:false};
  fs.writeFileSync(root+'/FIELD_REPORT_STATUS.json',JSON.stringify(result,null,2)+'\n');emit(result);finish();return;
 }
 if(state.phase==='launcher_error'||state.phase==='workers_finished_with_failures'||state.phase==='workers_finished_inspect_exact_coverage')throw Error('controller finished without complete coverage');
 emit({phase:state.phase,verified:state.verified});setTimeout(tick,60000);
}catch(e){emit({observer_error:String(e),exports_performed:false});finish();process.exitCode=2;}}
tick();
