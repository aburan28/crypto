#!/usr/bin/env python3
"""Bounded coding sessions, verified zero pricing, no host/provider credentials."""
import argparse
import datetime as dt
import json
import os
from pathlib import Path
import signal
import subprocess
import time
import urllib.request

PREFERRED=['mimo-v2.5-free','big-pickle','nemotron-3.5-lightning-free']
LOG=Path('/logs/agent')

def write(path,value):
    temporary=path.with_suffix('.tmp');temporary.write_text(json.dumps(value,indent=2)+'\n');temporary.replace(path)

def fetch(url):
    req=urllib.request.Request(url,headers={'User-Agent':'opencode/1.18.30'})
    with urllib.request.urlopen(req,timeout=30) as response:return json.load(response)

def free_models(catalog,active):
    available={row['id'] for row in active['data']};allowed=[]
    for name in PREFERRED:
        item=catalog['opencode']['models'].get(name,{})
        cost=item.get('cost',{})
        if (name in available and item.get('status') not in ('deprecated','retired')
            and item.get('tool_call') is True and cost.get('input')==0 and cost.get('output')==0
            and all(type(value) in (int,float) and value==0 for value in cost.values())):
            allowed.append(name)
    if not allowed:raise RuntimeError('No active approved model has verified zero prices')
    return allowed

def config(model):
    permissions={'*':'deny','read':'allow','glob':'allow','grep':'allow','edit':'allow',
                 'bash':'allow','external_directory':{'*':'deny','/opt/*':'allow'},
                 'task':'deny','question':'deny','skill':'deny','webfetch':'deny','websearch':'deny'}
    return {'$schema':'https://opencode.ai/config.json','model':'opencode/'+model,
            'small_model':'opencode/'+model,'enabled_providers':['opencode'],
            'provider':{'opencode':{'whitelist':[model]}},'default_agent':'research',
            'share':'disabled','snapshot':False,'autoupdate':False,'permission':permissions,
            'agent':{'research':{'description':'Improve verified complete index-calculus cost',
                       'mode':'primary','model':'opencode/'+model,'steps':40,'permission':permissions}}}

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--seconds',type=int,default=3600)
    parser.add_argument('--instruction',required=True);args=parser.parse_args()
    started=time.time();deadline=started+args.seconds;LOG.mkdir(parents=True,exist_ok=True)
    total_cost=0;total_tools=0;finished=[]
    for cycle in range(4):
        remaining=deadline-time.time()
        if remaining<60:break
        run=LOG/f'cycle-{cycle+1:02d}';run.mkdir()
        try:
            catalog=fetch('https://models.dev/api.json');active=fetch('https://opencode.ai/zen/v1/models')
            allowed=free_models(catalog,active);model=allowed[cycle%len(allowed)]
            write(run/'models.json',catalog)
            write(run/'pricing.json',{'checked_at':dt.datetime.now(dt.timezone.utc).isoformat(),
                'model':model,'allowed':allowed,'cost':catalog['opencode']['models'][model]['cost'],'paid_fallback':False})
            env={key:os.environ[key] for key in ('PATH','LANG','LC_ALL','CARGO_HOME') if key in os.environ}
            env.update(HOME='/home/researcher',USER='researcher',LOGNAME='researcher',
                XDG_CONFIG_HOME='/home/researcher/.config',XDG_DATA_HOME='/home/researcher/.local/share',
                XDG_CACHE_HOME='/home/researcher/.cache',XDG_STATE_HOME='/home/researcher/.local/state',
                OPENCODE_CONFIG_CONTENT=json.dumps(config(model)),OPENCODE_MODELS_PATH=str(run/'models.json'),
                OPENCODE_DISABLE_MODELS_FETCH='1',OPENCODE_DISABLE_AUTOUPDATE='1',
                OPENCODE_DISABLE_CLAUDE_CODE='1',OPENCODE_DISABLE_LSP_DOWNLOAD='1',
                OPENCODE_DISABLE_DEFAULT_PLUGINS='1',OPENCODE_EXPERIMENTAL_OUTPUT_TOKEN_MAX='8192',
                OPENCODE_EXPERIMENTAL_BASH_DEFAULT_TIMEOUT_MS='720000')
            prompt=args.instruction+f'\nThis is coding session {cycle+1}/4. You have at most {int(min(900,remaining))} seconds. '
            prompt+='First inspect existing probe results and NOTES.md, then make and measure a concrete candidate. Do not simply repeat the same configuration. Run at least one probe and leave notes grounded in its output.'
            write(LOG/'status.json',{'status':'coding','cycle':cycle+1,'model':model,'started_unix':started,
                  'deadline_unix':deadline,'reported_inference_cost_usd':total_cost,'tools_completed':total_tools})
            with (run/'events.jsonl').open('w') as stdout,(run/'stderr.txt').open('w') as stderr:
                process=subprocess.Popen(['/usr/local/bin/opencode','run','--pure','--format','json','--agent','research',
                    '--model','opencode/'+model,'--dir','/app',prompt],env=env,cwd='/app',stdout=stdout,stderr=stderr,start_new_session=True)
                timed_out=False;observed_cost=0;observed_lines=0
                try:
                    end=min(deadline,time.time()+900)
                    while process.poll() is None and time.time()<end:
                        time.sleep(3)
                        lines=(run/'events.jsonl').read_text().splitlines()
                        for line in lines[observed_lines:]:
                            try:event=json.loads(line)
                            except json.JSONDecodeError:continue
                            if event.get('type')=='step_finish':
                                cost=event.get('part',{}).get('cost')
                                if type(cost) not in (int,float) or cost!=0:
                                    raise RuntimeError('Missing or nonzero model cost; stopping free-only campaign')
                                observed_cost+=cost
                        observed_lines=len(lines)
                    if process.poll() is None:timed_out=True
                finally:
                    if process.poll() is None:
                        os.killpg(process.pid,signal.SIGTERM)
                        try:process.wait(timeout=10)
                        except subprocess.TimeoutExpired:os.killpg(process.pid,signal.SIGKILL);process.wait()
            events=[]
            for line in (run/'events.jsonl').read_text().splitlines():
                try:events.append(json.loads(line))
                except json.JSONDecodeError:pass
            costs=[e.get('part',{}).get('cost') for e in events if e.get('type')=='step_finish']
            if any(type(c) not in (int,float) or c!=0 for c in costs):raise RuntimeError('Unverified free-model cost')
            tools=sum(e.get('type')=='tool_use' for e in events)
            total_cost+=sum(costs);total_tools+=tools
            result={'cycle':cycle+1,'model':model,'returncode':process.returncode,'timed_out':timed_out,
                    'tool_events':tools,'steps_with_cost':len(costs),'reported_cost_usd':sum(costs),
                    'errors':[e for e in events if e.get('type')=='error']}
            write(run/'result.json',result);finished.append(result)
            with (LOG/'opencode.txt').open('a') as combined:combined.write((run/'events.jsonl').read_text())
            if len(list(Path('/app/probes').glob('*/result.json')))>=12:break
        except Exception as error:
            write(run/'failure.json',{'error':str(error)});finished.append({'cycle':cycle+1,'error':str(error)})
            if 'cost' in str(error).lower():break
    write(LOG/'status.json',{'status':'finished','wall_seconds':time.time()-started,
          'reported_inference_cost_usd':total_cost,'tools_completed':total_tools,'cycles':finished,
          'cost_status':'provider zero-price catalog and completed event costs; unfinished requests not independently invoiced'})
    if total_tools==0:raise RuntimeError('No model tool work completed; campaign is not a successful launch')

if __name__=='__main__':main()
