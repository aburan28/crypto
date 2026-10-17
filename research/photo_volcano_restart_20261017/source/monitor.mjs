// Local status/checkpoint observer only. Never exports files or starts solvers.
import fs from 'node:fs';
import path from 'node:path';
import {spawnSync} from 'node:child_process';
const root=path.resolve(process.argv[2]);
const lock=root+'/monitor.lock';
try{const old=JSON.parse(fs.readFileSync(lock));try{process.kill(old.pid,0);throw Error('monitor already active');}catch(e){if(e.code!=='ESRCH')throw e;}}catch(e){if(e.code!=='ENOENT'&&e.code!=='ESRCH')throw e;}
fs.writeFileSync(lock,JSON.stringify({pid:process.pid,started_utc:new Date().toISOString()})+'\n',{flag:fs.existsSync(lock)?'w':'wx'});
const journal=root+'/logs/progress_'+Date.now()+'_'+process.pid+'.jsonl';
let noWorkers=0;
function activeWorkers(){const ids=[];for(const p of fs.readdirSync('/proc').filter(x=>/^\d+$/.test(x))){try{const words=fs.readFileSync('/proc/'+p+'/cmdline','utf8').split('\0');if(words[0]===root+'/worker'&&words[1]===root)ids.push(Number(p));}catch(e){if(!['ENOENT','ESRCH','EACCES'].includes(e.code))throw e;}}return ids;}
function scan(){const totals={rho:0,bsgs:0},failed={rho:0,bsgs:0};for(const alg of ['rho','bsgs']){const dir=root+'/runs/'+alg;for(const job of fs.readdirSync(dir)){const at=dir+'/'+job,records=new Map();if(fs.existsSync(at+'/complete.json')){for(const r of JSON.parse(fs.readFileSync(at+'/complete.json')).records)records.set(r.instance,r);}else for(const a of fs.readdirSync(at).filter(x=>x.startsWith('attempt_')).sort())for(const f of fs.readdirSync(at+'/'+a).filter(x=>/^case_\d+\.json$/.test(x))){const r=JSON.parse(fs.readFileSync(at+'/'+a+'/'+f));records.set(r.instance,r);}for(const r of records.values())r.verified&&r.status===0?totals[alg]++:failed[alg]++;}}
  const pids=activeWorkers();const status={updated_utc:new Date().toISOString(),observer_pid:process.pid,active_native_worker_pids:pids,verified:totals,failed_cases_in_latest_attempts:failed,requested_per_algorithm:4600,native_campaign_complete:totals.rho===4600&&totals.bsgs===4600,full_random_IC_recoveries:0,full_requested_experiment_complete:false,raw_attempts_retained:true,exports_performed:false};
  fs.writeFileSync(root+'/LIVE_STATUS.json.tmp',JSON.stringify(status,null,2)+'\n');fs.renameSync(root+'/LIVE_STATUS.json.tmp',root+'/LIVE_STATUS.json');fs.appendFileSync(journal,JSON.stringify(status)+'\n');return status;
}
function finishReport(status){const tag='native_complete_'+Date.now();const tasks=[[process.execPath,[root+'/source/report.mjs',root,tag]],[process.execPath,[root+'/source/quantitative_report.mjs',root,tag+'_visuals',tag,'paired_diagnostics_0640']]];for(const[exe,args]of tasks){const result=spawnSync(exe,args,{encoding:'utf8',maxBuffer:16*1024*1024});fs.appendFileSync(journal,JSON.stringify({phase:'final_report',command:args,exit_status:result.status,stdout:result.stdout,stderr:result.stderr})+'\n');if(result.status!==0)return false;}
  fs.appendFileSync(journal,JSON.stringify({local_report:tag+'_visuals',pdf_render_and_page_inspection_pending:true,PR_update_pending:true,full_IC_recovery_still_unresolved:true})+'\n');return true;}
function tick(){try{const status=scan();if(status.native_campaign_complete){finishReport(status);fs.unlinkSync(lock);return;}noWorkers=status.active_native_worker_pids.length?0:noWorkers+1;if(noWorkers>=3){fs.appendFileSync(journal,JSON.stringify({stalled:true,reason:'No native workers observed in three consecutive polls; incomplete cases remain. No replacement solver started.'})+'\n');fs.unlinkSync(lock);process.exitCode=2;return;}setTimeout(tick,60000);}catch(e){fs.appendFileSync(journal,JSON.stringify({observer_error:String(e)})+'\n');process.exitCode=1;}}
tick();
