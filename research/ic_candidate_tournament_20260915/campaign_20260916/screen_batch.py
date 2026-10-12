"""Bounded development screen; never promotes or claims confirmed rho parity."""
import copy
import json
import os
from pathlib import Path
import random
import subprocess
import sys

WORK = Path(__file__).resolve().parent
ROOT = WORK.parent
sys.path.insert(0, str(ROOT))
import tournament as t


def main():
    out = WORK / 'batch16-development-screen'
    out.mkdir(exist_ok=False)
    parent = ROOT / 'runs/round-0004'
    arms = t.read(parent / 'candidates.json')
    selected = [copy.deepcopy(next(a for a in arms if a['id'] == name))
                for name in ['incumbent', 'folded_lift']]
    selected += [t.synthetic_arm('rho', selected[0])]
    for arm in selected:
        arm['binary_relative'] = str(parent / arm['binary_relative'])
        arm['binary_sha256'] = t.digest(Path(arm['binary_relative']))
    siblings = set()
    for part in Path('/sys/devices/system/cpu/cpu1/topology/thread_siblings_list').read_text().strip().split(','):
        bounds = list(map(int,part.split('-')))
        siblings.update(range(bounds[0],bounds[-1]+1))
    available = sorted(set(os.sched_getaffinity(0))-siblings)
    profile_cpu = available[-1]
    controller_cpu = available[0]
    contract = {'unit': t.UNIT, 'targets': 16, 'repetitions': 3,
                'limits': {'timeout_seconds': 60, 'memory_bytes': 8*1024**3, 'cpu': profile_cpu},
                'scope': 'Development screen only: four cells, no promotion or confirmed parity claim.',
                'boundary': 'Complete cold 16-target job; all setup and all target solves charged. Rho uses the existing per-target solver API on the same constructed curve.',
                'floor': 'K relation-producing trials/instructions for K full-rank columns; weak implementation floor.',
                'hypothesis': 'Sharing IC setup across 16 targets may close the instruction/time gap.',
                'followup_gate': 'Only a fresh full tournament with confirmation and replay may establish parity.',
                'seed': 2026091690, 'parent_contract_sha256': t.digest(parent/'contract.json')}
    t.write(out/'contract.json', contract, exclusive=True)
    t.write(out/'candidates.json', selected, exclusive=True)
    rng = random.Random(contract['seed'])
    cases = []
    for n,a in [(13,0),(17,1),(19,0),(23,0)]:
        job = {'mode':'fixture','degree':n,'curve_a':a,
               'target_seeds':[rng.getrandbits(64) for _ in range(16)],
               'algorithm_seed':rng.getrandbits(64),
               'factor_base':{'kind':'subgroup_orbits','seed':43,'points':6*n},
               'config':selected[0]['config']}
        process = subprocess.run([str(parent/'worker')],input=json.dumps(job),
                                 capture_output=True,text=True,env=t.child_env(),check=True,timeout=60)
        fixture = json.loads(process.stdout)['fixture']
        t.Curve(fixture)
        cases.append({'id':f'n{n}a{a}-screen','cell':f'n{n}a{a}','job':job,
                      'fixture':fixture,'fixture_sha256':t.objhash(fixture)})
    t.write(out/'fixtures.json', cases, exclusive=True)
    # Keep controller/audit work away from the single-target measurement CPU.
    os.sched_setaffinity(0,{controller_cpu})
    rows=[]
    schedule=[(case,rep) for case in cases for rep in range(3)]
    rng.shuffle(schedule)
    for case,rep in schedule:
        order=list(selected);rng.shuffle(order)
        for arm in order:
            assert t.digest(Path(arm['binary_relative'])) == arm['binary_sha256']
            row=t.run_trial(out,contract,'development',case,arm,rep)
            rows.append(row)
            print(json.dumps({'case':case['id'],'arm':arm['id'],'rep':rep,'status':row['status']}),flush=True)
    for case in cases:
        for arm in selected:
            for rep in range(3):
                t.verify_receipt(t.trial_path(out,'development',case,arm['id'],rep),contract,case,arm)
    result={'scope':contract['scope'],'verified':sum(r['status']=='VERIFIED' for r in rows),
            'runs':len(rows),'candidate_over_incumbent':t.comparison(rows,'folded_lift'),
            'candidate_over_rho':t.comparison(rows,'folded_lift',baseline='rho')}
    t.write(out/'result.json', result, exclusive=True)
    print(json.dumps(result),flush=True)


if __name__ == '__main__':
    main()
