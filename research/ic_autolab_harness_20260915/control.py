#!/usr/bin/env python3
"""Launch and supervise one real AutoLab/Harbor trial on existing local compute."""
import argparse
import datetime as dt
import json
import os
import re
import resource
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import time

HERE=Path(__file__).resolve().parent
REPO=HERE.parents[1]
UPSTREAM=Path('/home/ubuntu/crypto-work/autolab-upstream-20260915')
HARBOR=UPSTREAM/'.venv/bin/harbor'

def read(path):return json.loads(path.read_text())
def write(path,value):
    temporary=path.with_suffix('.tmp');temporary.write_text(json.dumps(value,indent=2,sort_keys=True)+'\n');temporary.replace(path)
def command(args,**kwargs):return subprocess.run(args,text=True,capture_output=True,check=True,**kwargs).stdout.strip()
def containers(run):
    found=[]
    for trial in (run/'jobs'/run.name).glob('*'):
        if (trial/'config.json').exists():
            # Harbor's Compose project name is the trial directory basename.
            project=trial.name.lower()
            if not re.match(r'^[a-z0-9]',project):project='0'+project
            project=re.sub(r'[^a-z0-9_-]','-',project)
            ids=command(['docker','ps','-aq','--filter','label=com.docker.compose.project='+project])
            found.extend(ids.splitlines())
    return sorted(set(found))

def launch(image_name="ic-autolab:20260915", agent="free"):
    stamp=dt.datetime.now(dt.timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    run=HERE/'runs'/('autolab-'+stamp);run.mkdir(parents=True)
    shutil.copytree(HERE/'task',run/'task')
    (run/'adapter').mkdir();shutil.copy2(HERE/'free_agent.py',run/'adapter/free_agent.py')
    image=command(['docker','image','inspect',image_name,'--format','{{.Id}}'])
    task=run/'task/task.toml';task.write_text(task.read_text().replace('ic-autolab:20260915',image))
    script=run/'task/tests/test.sh'
    script.write_text(script.read_text().replace('export CARGO_HOME=/opt/cargo', 'export CARGO_HOME=/opt/cargo\nexport IC_AUTOLAB_ROUND_NAME=round-'+stamp.lower()))
    config={'job_name':run.name,'jobs_dir':str(run/'jobs'),'n_attempts':1,'n_concurrent_trials':1,
            'quiet':True,'retry':{'max_retries':0},'environment':{'type':'docker','delete':False,'kwargs':{'keep_containers':True}},
            'agents':[{'import_path':'free_agent:FreeOpenCode','model_name':'opencode/mimo-v2.5-free',
                       'override_setup_timeout_sec':120}], 'tasks':[{'path':str(run/'task')}],
            'verifier':{'max_timeout_sec':3300}}
    if agent=='nop':config['agents']=[{'name':'nop'}]
    write(run/'job.json',config)
    unit='ic-autolab-'+stamp.lower()
    write(run/'launch.json',{'unit':unit,'created_utc':stamp,'image_id':image,
        'autolab_commit':command(['git','rev-parse','HEAD'],cwd=UPSTREAM),
        'harbor_version':command([str(UPSTREAM/'.venv/bin/python'),'-c','import importlib.metadata; print(importlib.metadata.version("harbor"))']),
        'baseline_round':str(REPO/'research/ic_candidate_tournament_20260915/runs/round-0002'),
        'wall_budget_seconds':7200,'agent_budget_seconds':0 if agent=='nop' else 3600,'verifier_budget_seconds':3300,
        'cpus':2,'memory_gib':8,'maximum_allocated_container_cpu_hours':4,'maximum_allocated_gib_hours':16,
        'inference_policy':('No model calls; frozen source-candidate verification' if agent=='nop' else 'Active zero-price OpenCode models only, live checked before every coding session; no credentials or paid fallback'),
        'agent_mode':agent,
        'local_compute_price_usd_per_hour':None,'cloud_resources_created':False,
        'promotion_gate':'<=0.8 complete Ir cost ratio; paired 95% upper <1; max cell ratio <=1.1; verified fresh holdouts and replay',
        'merge_policy':'Preserve isolated candidate and evidence; no automatic library merge'})
    write(run/'status.json',{'status':'starting','launch':str(run/'launch.json')})
    write(HERE/'latest.json',{'run':str(run),'unit':unit})
    result=command(['systemd-run','--user','--unit',unit,'--property=RuntimeMaxSec=7200',
                    '--property=TimeoutStopSec=60','--property=KillMode=mixed',
                    '/usr/bin/python3',str(HERE/'control.py'),'supervise','--run',str(run)])
    print(json.dumps({'run':str(run),'unit':unit,'systemd':result},indent=2))

def supervise(run):
    started=time.time();deadline=started+7140
    write(run/'status.json',{'status':'starting','started_unix':started,'deadline_unix':started+7200})
    env={key:os.environ[key] for key in ('PATH','HOME','XDG_RUNTIME_DIR','DBUS_SESSION_BUS_ADDRESS','LANG') if key in os.environ}
    env['PYTHONPATH']=str(run/'adapter');env['PYTHONUNBUFFERED']='1'
    process=None;outcome='inconclusive';error=None
    def stop_signal(signum,frame):raise TimeoutError('Campaign service deadline reached')
    signal.signal(signal.SIGTERM,stop_signal)
    try:
        with (run/'harbor.log').open('w') as log:
            process=subprocess.Popen([str(HARBOR),'run','--config',str(run/'job.json')],cwd=UPSTREAM,env=env,
                                     stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
            while process.poll() is None:
                if time.time()>=deadline:raise TimeoutError('Two-hour campaign budget exhausted')
                state={'status':'running','started_unix':started,'deadline_unix':started+7200,
                       'elapsed_seconds':time.time()-started,'harbor_pid':process.pid}
                for p in (run/'jobs'/run.name).glob('*/agent/status.json'):
                    state['agent']=read(p)
                for p in (run/'jobs'/run.name).glob('*/verifier/round-*/state.json'):
                    state['verification']=read(p)
                ids=containers(run);state['containers']=ids
                for cid in ids:
                    raw=subprocess.run(['docker','exec',cid,'cat','/sys/fs/cgroup/cpu.stat'],capture_output=True,text=True,timeout=10)
                    if raw.returncode==0:
                        with (run/'resources.jsonl').open('a') as ledger:
                            ledger.write(json.dumps({'unix':time.time(),'container':cid,'cpu_stat':raw.stdout})+'\n')
                write(run/'status.json',state)
                time.sleep(15)
        outcomes=list((run/'jobs'/run.name).glob('*/verifier/outcome.json'))
        if outcomes:outcome=read(outcomes[0])['status']
        if process.returncode and not outcomes:raise RuntimeError('Harbor exited '+str(process.returncode))
        rounds=list((run/'jobs'/run.name).glob('*/verifier/round-*/decision.json'))
        for decision in rounds:
            with (run/'publication.log').open('w') as log:
                subprocess.run(['/usr/bin/python3',str(REPO/'research/ic_candidate_tournament_20260915/report.py'),
                                '--round',str(decision.parent),'--scoreboard',str(REPO/'docs/index-calculus-scoreboard.html')],
                               check=True,stdout=log,stderr=subprocess.STDOUT,timeout=max(1,deadline-time.time()))
        write(run/'completion.json',{'status':outcome,'harbor_returncode':process.returncode,
              'wall_seconds':time.time()-started,'published_reports':[str(p.parent/'REPORT.md') for p in rounds]})
    except Exception as exc:
        error=str(exc)
        write(run/'failure.json',{'status':'inconclusive','error':error,'elapsed_seconds':time.time()-started})
    finally:
        if process and process.poll() is None:
            os.killpg(process.pid,signal.SIGTERM)
            try:process.wait(timeout=10)
            except subprocess.TimeoutExpired:os.killpg(process.pid,signal.SIGKILL);process.wait()
        cleanup(run)
        write(run/'status.json',{'status':outcome if error is None else 'inconclusive','error':error,
              'finished_unix':time.time(),'wall_seconds':time.time()-started,'active':False,
              'host_child_cpu_seconds':resource.getrusage(resource.RUSAGE_CHILDREN).ru_utime+resource.getrusage(resource.RUSAGE_CHILDREN).ru_stime})

def cleanup(run):
    # Scope cleanup to this trial's Compose project only. Preserve unfinished work.
    for cid in containers(run):
        destination=run/'interrupted-artifacts'/cid;destination.mkdir(parents=True,exist_ok=True)
        subprocess.run(['docker','stop','--time','5',cid],capture_output=True,timeout=15)
        names=['best','NOTES.md','candidate.json','src','Cargo.toml','Cargo.lock','build.rs','examples']
        outcomes=list((run/'jobs'/run.name).glob('*/verifier/outcome.json'))
        if not outcomes or read(outcomes[0]).get('status')=='inconclusive':names.append('probes')
        copies=[]
        for name in names:
            result=subprocess.run(['docker','cp',cid+':/app/'+name,str(destination/name)],capture_output=True,timeout=25)
            copies.append({'path':name,'returncode':result.returncode})
        write(destination/'copy-status.json',copies)
        subprocess.run(['docker','rm',cid],capture_output=True,timeout=15)
    # Networks/images are retained for reproducibility; no unrelated service is touched.

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=['launch','supervise','cleanup']);parser.add_argument('--run',type=Path)
    parser.add_argument('--image',default='ic-autolab:20260915')
    parser.add_argument('--agent',choices=['free','nop'],default='free')
    args=parser.parse_args()
    if args.action=='launch':launch(args.image,args.agent)
    elif args.action=='supervise':supervise(args.run.resolve())
    else:cleanup(args.run.resolve())
