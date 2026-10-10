"""Follow the exact solver_07 panel with linear support and root recovery."""
import hashlib
import importlib.util
import json
from pathlib import Path
import platform
import subprocess
import sys
import time
import traceback
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2];CODE=ROOT/'ecc2k130/codegen'
sys.path.insert(0,str(CODE));sys.path.insert(0,str(HERE.parent/'solver_07'))
import image_solver
import field
import curves
import nagaocompare
spec=importlib.util.spec_from_file_location('campaign07',HERE.parent/'solver_07/run.py');previous=importlib.util.module_from_spec(spec);spec.loader.exec_module(previous)


def validate():
    oracle=previous.PairOracle(5,3);f=oracle.f;c=oracle.curve;checked=0;models=0
    for target in nagaocompare.affinePoints(f,c):
        search=image_solver.Search(f,c,3,target,time.perf_counter()+30)
        reference=previous.optimized.Search(f,c,3,target,time.perf_counter()+30)
        # For every candidate, the image test must return exactly the roots
        # admitted by the previous full cubic support test, including negatives.
        for a,b,invB,z in search.candidates():
            want=reference.recover(a,b,invB);got=search.recover(a,b,invB,z)
            if want!=got:raise ArithmeticError('image/cubic support mismatch')
            models+=1
        result=image_solver.cell(5,3,[f.toCoords(v) for v in target],'enumerate',30)
        if result['status']!='complete' or {tuple(x) for x in result['solutions']}!=oracle.expected(target):raise ArithmeticError('all-target oracle mismatch')
        checked+=1
    return {'all_affine_five_bit_targets':checked,'candidate_support_comparisons':models,'image_equals_cubic_support':True,'complete_oracle_equality':True}


def main():
    raw=HERE/'raw.jsonl'
    if raw.exists():raise SystemExit('Evidence exists')
    contract=json.loads((HERE/'contract.json').read_text())
    sources=list(CODE.glob('*.py'))+[HERE/'run.py',HERE/'image_solver.py',HERE/'contract.json']
    for folder in ('solver_02','solver_05','solver_06','solver_07'):sources+=list((HERE.parent/folder).glob('*.py'))
    priorPath=HERE.parent/'solver_07/raw.jsonl';sources.append(priorPath)
    provenance={'commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),'python':platform.python_version(),'command':[sys.executable]+sys.argv,'sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}}
    rows=[]
    with raw.open('x') as out:
        def record(r):rows.append(r);out.write(json.dumps(r)+'\n');out.flush()
        record({'kind':'provenance',**provenance})
        try:record({'kind':'validation',**validate()});print('validation passed',flush=True)
        except Exception:record({'kind':'validation_failure','traceback':traceback.format_exc()});raise
        prior=[json.loads(x) for x in priorPath.read_text().splitlines()]
        for panel in contract['panels']:
            n,d=panel['n'],panel['d'];t=time.perf_counter();oracle=previous.PairOracle(n,d);f=oracle.f
            record({'kind':'oracle','n':n,'d':d,'validation_setup_seconds':time.perf_counter()-t})
            for old in prior:
                if old.get('kind')!='trial' or (old['n'],old['d'],old['variant'])!=(n,d,'quadratic-optimized'):continue
                try:
                    r=image_solver.cell(n,d,old['target'],old['mode'],contract['seconds_per_instance'])
                    expected=oracle.expected(tuple(f.fromCoords(x) for x in old['target']));got={tuple(x) for x in r['solutions']}
                    if not got<=expected or (r['status']=='complete' and got!=expected):raise ArithmeticError('benchmark mismatch')
                    record({'kind':'trial','stratum':old['stratum'],'expected_count':len(expected),**r})
                    print(n,old['stratum'],old['mode'],r['status'],len(got),round(r['all_phase_seconds'],3),flush=True)
                except Exception:record({'kind':'error','prior':old,'traceback':traceback.format_exc()})
    groups=[]
    for panel in contract['panels']:
        for mode in ('first','enumerate'):
            for st in ('uniform','known_decomposable'):
                rs=[r for r in rows if r['kind']=='trial' and (r['n'],r['mode'],r['stratum'])==(panel['n'],mode,st)]
                secs=sum(r['all_phase_seconds'] for r in rs);rels=sum(r['verified_unique_relations'] for r in rs)
                groups.append({**panel,'variant':'quadratic-image','mode':mode,'stratum':st,'attempted':len(rs),'resolved':sum(r['status']!='timeout' for r in rs),'resolved_within_budget':sum(r['status']!='timeout' and r['within_budget'] for r in rs),'verified_relations':rels,'all_phase_seconds':secs,'diagnostic_seconds_per_relation':secs/rels if rels else None,'field_api_operations':sum(r['field_api_counts']['totals']['fieldOperations'] for r in rs),'full_dlp_S':None,'rho_ratio':None})
    (HERE/'summary.json').write_text(json.dumps({'provenance':provenance,'groups':groups,'errors':[r for r in rows if r['kind']=='error'],'limits':contract['limits']},indent=2)+'\n')

if __name__=='__main__':main()
