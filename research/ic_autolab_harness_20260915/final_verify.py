#!/usr/bin/env python3
"""Lock the submitted candidate, then run the unchanged fresh-holdout tournament."""
import json
import os
import re
from pathlib import Path
import secrets
import shutil
import subprocess
import sys
import time
sys.path.insert(0,'/opt/evaluator')
import tournament as t
from probe import APP,BASE,BASE_BUILT,CONFIG,validate

out=Path('/logs/verifier')
out.mkdir(parents=True,exist_ok=True)
started=time.time()
try:
    # Stop all remaining processes of the coding user before reading the candidate.
    subprocess.run(['pkill','-KILL','-u','researcher'],check=False)
    selected=APP/'best/source' if (APP/'best/source').is_dir() else APP
    config_path=APP/'best/candidate-config.json' if selected!=APP else APP/'candidate.json'
    config=t.read(config_path)
    changed=validate(selected,config)
    t.write(out/'submission.json',{'selected_source':str(selected),'config':config,'changed_files':changed,
                                 'locked_unix':time.time(),'merge_status':'not merged; source changes require review'})
    # Always preserve notes and development evidence, including rejected/failed probes.
    for name in ('NOTES.md','candidate.json'):
        if (APP/name).is_file():shutil.copy2(APP/name,out/name)
    if (APP/'probes').is_dir():shutil.copytree(APP/'probes',out/'development')
    if not changed and config==CONFIG:
        t.write(out/'outcome.json',{'status':'incumbent_retained','reason':'No changed source/config submitted','promotion':False})
        (out/'reward.txt').write_text('0.0\n')
    else:
        registry=[{'id':'incumbent','hypothesis':'Previous audited batch16 winner','config':CONFIG},
                  {'id':'autolab','hypothesis':'Agent-selected development checkpoint; see NOTES.md',
                   'config':config,'source_root':str(selected)}]
        t.write(out/'registry.json',registry,exclusive=True)
        seed=secrets.randbits(63)
        t.write(out/'holdout-seed.json',{'seed':seed,'generated_after_candidate_locked':True},exclusive=True)
        round_name=os.environ.get('IC_AUTOLAB_ROUND_NAME','round-autolab-20260915')
        t.require(re.fullmatch(r'round-[a-z0-9_-]+',round_name) is not None,'invalid round name')
        round_dir=out/round_name
        def run(args,name):
            with (out/name).open('w') as log:
                subprocess.run(args,check=True,stdout=log,stderr=subprocess.STDOUT,env=t.child_env())
        run(['python3','/opt/evaluator/tournament.py','prepare','--out',str(round_dir),
             '--source-root',str(BASE),'--candidates',str(out/'registry.json'),
             '--profile','pilot','--seed',str(seed),'--timeout','30','--max-processes','1400'],'prepare.log')
        runner=round_dir/'evaluator/tournament.py'
        run(['python3',str(runner),'run','--round',str(round_dir)],'tournament.log')
        run(['python3',str(runner),'verify','--round',str(round_dir)],'audit.log')
        decision=t.read(round_dir/'decision.json')
        t.write(out/'outcome.json',{'status':decision['status'],'promotion':decision['status']=='promoted',
                                  'decision':decision,'audit':'passed','wall_seconds':time.time()-started})
        (out/'reward.txt').write_text('1.0\n' if decision['status']=='promoted' else '0.0\n')
except Exception as error:
    t.write(out/'outcome.json',{'status':'inconclusive','promotion':False,'error':str(error),'wall_seconds':time.time()-started})
    (out/'reward.txt').write_text('0.0\n')
    raise
finally:
    usage={}
    for name in ('cpu.stat','memory.peak','memory.events'):
        p=Path('/sys/fs/cgroup')/name
        if p.exists():usage[name]=p.read_text()
    t.write(out/'resources.json',{'wall_seconds_final_stage':time.time()-started,'cgroup':usage,
        'compute_pricing':'Existing local machine; no cloud instances purchased; dollar CPU/electricity price unspecified',
        'budget':'2 CPUs, 8 GiB, 2 hours total wall; <=4 allocated CPU-hours and <=16 GiB-hours',
        'cost_unit':t.UNIT,'excluded_from_solver_Ir':'kernel/device, profiler execution, build and independent audit overhead; retained as research resource cost'})

    subprocess.run(['chmod','-R','a+rwX',str(out)],check=False)
