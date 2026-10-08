"""Confirm a batch-parity candidate after the single-target parent is promoted."""
import copy
import difflib
import fcntl
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

from descent_mechanisms import lazy_descent, fast_descent_check

WORK = Path(__file__).resolve().parent
ROOT = WORK.parent
PARENT = ROOT / 'runs/round-0004'
ROUND = ROOT / 'runs/round-0005-batch16'
REL = Path('src/cryptanalysis/koblitz_index_calculus.rs')


def read(path):
    return json.loads(path.read_text())


def main():
    deadline = time.monotonic() + 3*3600
    while not (PARENT/'summaries/selection.json').exists():
        if time.monotonic()>deadline: raise TimeoutError('parent selection missing')
        time.sleep(5)
    selected=read(PARENT/'summaries/selection.json')['provisional_challenger']
    arm=next(a for a in read(PARENT/'candidates.json') if a['id']==selected)
    source=PARENT/arm['source_directory']
    config=arm['config']
    original=(source/REL).read_text()
    registry=[{'id':'incumbent','config':config,'parent':selected,
               'parent_round':str(PARENT),'hypothesis':'Previous single-target winner, freshly measured on complete cold 16-target jobs.'}]
    for name,transform in [('lazy_descent',lazy_descent),('fast_descent_check',fast_descent_check),
                            ('combined_descent',lambda s:fast_descent_check(lazy_descent(s)))]:
        destination=WORK/'round5-sources'/name
        shutil.copytree(source,destination)
        changed=transform(original)
        (destination/REL).write_text(changed)
        (WORK/('round5-'+name+'.patch')).write_text(''.join(difflib.unified_diff(
            original.splitlines(True),changed.splitlines(True),fromfile='a/'+str(REL),tofile='b/'+str(REL))))
        registry.append({'id':name,'source_root':str(destination),'config':copy.deepcopy(config),
                         'parent':'incumbent','hypothesis':{
                             'lazy_descent':'Try the first existing probe before initializing the other 63 walks, preserving probe order and checks.',
                             'fast_descent_check':'Use the cached exact single-word curve for each recovered-log check; retain general final worker verification.',
                             'combined_descent':'Measure lazy walk initialization and the faster internal scalar check together.'}[name],
                         'falsification':'Any incorrect/missing target rejects; promotion and rho parity require both metric gates on confirmation and replay.'})
    registry_path=WORK/'round-0005-candidates.json'
    registry_path.write_text(json.dumps(registry,indent=2)+'\n')
    siblings=set()
    for part in Path('/sys/devices/system/cpu/cpu1/topology/thread_siblings_list').read_text().strip().split(','):
        bounds=list(map(int,part.split('-')));siblings.update(range(bounds[0],bounds[-1]+1))
    available=sorted(set(os.sched_getaffinity(0))-siblings)
    cpu=available[-1]
    (WORK/'ROUND5.md').write_text(f'''# Round 0005 pre-registration: complete cold batches of 16

This is a separate workload panel. It does not establish single-target parity.
Each job constructs one curve and 16 independent public-hash targets, charges all
factor-base/table/precomputation work once, solves every target, performs every
internal check and general final scalar check, and verifies all outputs with the
independent Python checker. Rho solves the same 16 targets through the existing
per-target signed-Frobenius solver API on the same constructed curve. No cache or
precomputation outside the timed/profiled job is admitted for either method.

Parent: round-0004 selection `{selected}`, conditional on full parent promotion.
Source: `{source}`. Fresh seed 2026091605; target count 16; five confirmation
cells, 60 fresh job fixtures (960 targets per repetition and arm), three
repetitions, independent replay. Budget 1800 paired jobs. Single logical CPU {cpu},
8 GiB address-space cap, 60-second watchdog per process. Builds avoid both threads
of the parent's measurement core, and measured execution waits for its completion.

Three source candidates: lazy first descent probe, cached fast verification of
recovered logs, and both. Base support, m=3, source accounting boundary and all
verification obligations stay fixed. The first-probe candidate keeps the same
probe order, avoids repeating the first failed probe, and respects the trial cap.

Reference is matched rho. The floor remains K instructions for K required
independent relation columns; it is weak and cannot establish a non-generic
advance. Class: engineering. Within this panel, promotion requires >=20% lower
instructions and native wall, upper paired 95% limits <1 and every cell <=1.10,
on confirmation and replay. Parity requires candidate/rho upper paired 95% limits
and every cell ratio <=1.10 in BOTH metrics on BOTH final stages.

The preceding 36-trial development screen completed correctly but did not meet
the parity rule: ratio estimates near one hid a regression in the smallest cell.
That screen selects hypotheses only; its cases are not reused for confirmation.
No mathematical or family-wide speedup is claimed. Comparisons use the shipped
rho implementation; additional cross-target rho optimizations are outside this
measured comparison and must be measured if introduced.
''')
    os.sched_setaffinity(0,set(available))
    with (WORK/'prepare-round-0005.log').open('w') as log:
        subprocess.run([sys.executable,str(ROOT/'tournament.py'),'prepare','--out',str(ROUND),
                        '--source-root',str(source),'--candidates',str(registry_path),
                        '--seed','2026091605','--targets','16','--cpu',str(cpu),
                        '--require-native-progress'],stdout=log,stderr=subprocess.STDOUT,check=True)
    print('Batch16 tournament prepared; waiting for parent promotion.',flush=True)
    with (PARENT/'operation.lock').open('a+') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        decision=read(PARENT/'decision.json')
        if decision['status']!='promoted' or decision['winner']!=selected:
            raise RuntimeError('Parent did not promote selected candidate; preserve prepared proposal.')
    with (ROUND/'operation.jsonl').open('w') as log:
        subprocess.run([sys.executable,str(ROUND/'evaluator/tournament.py'),'run','--round',str(ROUND)],
                       stdout=log,stderr=subprocess.STDOUT,check=True)
    print('Batch16 tournament completed.',flush=True)


if __name__=='__main__': main()
