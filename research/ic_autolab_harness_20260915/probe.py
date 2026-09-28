#!/usr/bin/env python3
"""Development-only complete-cost probe, using the frozen tournament checker."""
import argparse
import copy
import fcntl
import json
import os
from pathlib import Path
import random
import shutil
import subprocess
import sys
import time
import uuid
sys.path.insert(0, '/opt/evaluator')
import tournament as t

APP = Path('/app')
BASE = Path('/opt/baseline')
BASE_BUILT = Path('/opt/baseline-built')
CONFIG = {'solver':'pair_table','linear_algebra':'sparse','batch_trials':16,'max_trials':4096,'summands':3}
ALLOWED = {'src/binary_ecc.rs','src/cryptanalysis/koblitz_index_calculus.rs',
           'src/cryptanalysis/koblitz_sparse_la.rs','src/cryptanalysis/koblitz_factor_base_search.rs'}

def validate(source, config):
    t.require(isinstance(config,dict), 'configuration must be an object')
    t.require(set(config) <= set(CONFIG)|{'collection_window','sparse'}, 'unknown configuration key')
    t.require(config.get('summands') == 3 and config.get('max_trials') == 4096, 'changed fixed m or trial cap')
    t.require(config.get('solver') in ('pair_table','enumerate'), 'unsupported solver')
    t.require(config.get('linear_algebra') in ('dense','sparse'), 'unsupported linear algebra')
    batch=config.get('batch_trials')
    t.require(type(batch) is int and 1 <= batch <= 4096,'invalid batch size')
    baseline=t.read(BASE_BUILT/'source-manifest.json')
    changed=[]
    for rel, expected in baseline.items():
        path=source/rel
        t.require(path.is_file() and not path.is_symlink(), 'missing or symlink source: '+rel)
        if t.digest(path)!=expected:
            t.require(rel in ALLOWED, 'changed frozen source: '+rel)
            changed.append(rel)
            # Reject added counter controls or environment/fixture plumbing. This is
            # a guardrail; an audit of generated code remains required before merging.
            import difflib,re
            added='\n'.join(line[1:] for line in difflib.unified_diff((BASE/rel).read_text().splitlines(),path.read_text().splitlines())
                              if line.startswith('+') and not line.startswith('+++'))
            t.require(not re.search(r'(?i)valgrind|callgrind|asm!|std::process|std::env|include!|target_seeds|algorithm_seed',added),
                      'added instrumentation, external execution or fixture plumbing')
    actual={str(p.relative_to(source)) for p in (source/'src').rglob('*') if p.is_file()}
    expected={r for r in baseline if r.startswith('src/')}
    t.require(actual==expected,'added or deleted source files')
    return changed

def probe():
    started=time.time()
    probes=APP/'probes';probes.mkdir(exist_ok=True)
    with (APP/'probe.lock').open('a+') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        t.require(len(list(probes.iterdir())) < 12, '12-probe campaign cap reached')
        out=probes/(time.strftime('%Y%m%dT%H%M%S',time.gmtime())+'-'+uuid.uuid4().hex[:6])
        out.mkdir()
        try:
            config=t.read(APP/'candidate.json')
            changed=validate(APP,config)
            t.write(out/'candidate-config.json',config,exclusive=True)
            candidate_binary, manifest=t.snapshot_build(APP,out)
            shutil.rmtree(out/'build')  # Frozen source/binary suffice; discard reproducible compiler intermediates.
            shutil.copy2(BASE_BUILT/'worker',out/'baseline-worker')
            arms=[{'id':'incumbent','config':CONFIG,'binary_relative':'baseline-worker'},
                  {'id':'candidate','config':config,'binary_relative':'worker'},
                  {'id':'rho','config':CONFIG,'binary_relative':'baseline-worker'}]
            for arm in arms:
                arm['binary_sha256']=t.digest(out/arm['binary_relative'])
            limits={'timeout_seconds':30,'memory_bytes':8*1024**3,
                    'cpu':max(os.sched_getaffinity(0)),'worker_threads':1,'max_profiled_jobs':72}
            contract={'unit':t.UNIT,'limits':limits,'repetitions':3,'development_only':True,
                      'source_manifest_sha256':t.objhash(manifest),'arms':arms}
            t.write(out/'contract.json',contract,exclusive=True)
            rng=random.Random(916_202_609_151)
            cases=[]
            for degree,a in [(13,0),(17,1),(19,0),(23,0)]:
                for index in range(2):
                    job={'mode':'fixture','degree':degree,'curve_a':a,'target_seeds':[rng.getrandbits(64)],
                         'algorithm_seed':rng.getrandbits(64),
                         'factor_base':{'kind':'subgroup_orbits','seed':43,'points':6*degree},'config':CONFIG}
                    raw=subprocess.run([str(out/'baseline-worker')],input=json.dumps(job),capture_output=True,
                                       text=True,check=True,timeout=30,env=t.child_env())
                    fixture=json.loads(raw.stdout)['fixture'];t.Curve(fixture)
                    cases.append({'id':f'n{degree}a{a}-{index:03d}','cell':f'n{degree}a{a}',
                                  'job':job,'fixture':fixture,'fixture_sha256':t.objhash(fixture)})
            t.write(out/'fixtures.json',cases,exclusive=True)
            schedule=[(case,rep) for case in cases for rep in range(3)];rng.shuffle(schedule)
            rows=[]
            for case,rep in schedule:
                order=arms.copy();rng.shuffle(order)
                for arm in order:
                    row=t.run_trial(out,contract,'development',case,arm,rep);rows.append(row)
                    t.write(out/'state.json',{'completed_pairs':len(rows),'status':row['status'],'updated_unix':time.time()})
            comparison=t.comparison(rows,'candidate')
            result={'status':'complete','scope':'development only; no promotion', 'comparison':comparison,
                    'rho_over_incumbent':t.comparison(rows,'rho'),'changed_files':changed,
                    'candidate_config':config,'unit':t.UNIT,'total_profile_native_pairs':len(rows),
                    'failures':[{k:r.get(k) for k in ('case','arm','repetition','status','reason')} for r in rows if r['status']!='VERIFIED'],
                    'phase_totals':{arm['id']:{phase:sum(r['phase_costs'].get(phase,0) for r in rows if r['arm']==arm['id'] and r['status']=='VERIFIED')
                                            for phase in sorted(t.PHASES)} for arm in arms},
                    'wall_seconds':time.time()-started}
            t.write(out/'result.json',result,exclusive=True)
            best=APP/'best';previous=t.read(best/'result.json') if (best/'result.json').exists() else None
            if comparison.get('eligible') and not result['failures'] and comparison['candidate_over_baseline'] < (previous['comparison']['candidate_over_baseline'] if previous else 1.0):
                if best.exists(): shutil.rmtree(best)
                best.mkdir();shutil.copytree(out/'source',best/'source')
                for name in ('result.json','candidate-config.json','source-manifest.json'):shutil.copy2(out/name,best/name)
                t.write(best/'origin.json',{'probe':str(out),'source_manifest_sha256':t.objhash(manifest)})
            print(json.dumps({'out':str(out),**result},indent=2),flush=True)
        except Exception as error:
            t.write(out/'failure.json',{'status':'failed','error':str(error),'wall_seconds':time.time()-started})
            raise

if __name__=='__main__':
    probe()
