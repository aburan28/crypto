"""Deterministic quadratic function search on solver_04's exact fresh panel."""
import hashlib
import json
from pathlib import Path
import platform
import subprocess
import sys
import time
import traceback
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2];CODE=ROOT/'ecc2k130/codegen'
sys.path.insert(0,str(CODE));sys.path.insert(0,str(HERE.parent/'solver_02'))
import run as baseline
import quadratic
import curves
import field
import nagaocompare
import nagaoannihilator as scalar


def cell(n,d,targetCoords,mode,budget):
    start=time.perf_counter();f=scalar.CountedField(n);curve=curves.Curve(f)
    target=tuple(f.fromCoords(x) for x in targetCoords)
    solutions=set();seen=set();duplicates=0;candidates=0;first=None;complete=True
    f.phase='search'
    for a,b in quadratic.candidates(f,curve,d,target):
        if time.perf_counter()-start>=budget:complete=False;break
        candidates+=1
        if (a,b) in seen:duplicates+=1;continue
        seen.add((a,b));f.phase='support_extract_verify'
        xs=quadratic.recover(f,curve,d,target,a,b)
        f.phase='search'
        if xs is not None:
            solutions.add(xs)
            if first is None:first=time.perf_counter()-start
            if mode=='first':complete=False;break
    status='first' if mode=='first' and solutions else ('complete' if complete else 'timeout')
    return {'n':n,'d':d,'target':targetCoords,'mode':mode,'variant':'quadratic-function','status':status,
        'solutions':[list(x) for x in sorted(solutions)],'verified_unique_relations':len(solutions),
        'candidate_functions':candidates,'duplicate_functions':duplicates,
        'all_phase_seconds':time.perf_counter()-start,'first_verified_seconds':first,
        'field_api_counts':f.report(),'full_dlp_S':None,'rho_ratio':None}


def main():
    raw=HERE/'raw.jsonl'
    if raw.exists():raise SystemExit('Evidence exists')
    contract=json.loads((HERE/'contract.json').read_text())
    source=[HERE/'quadratic.py',HERE/'run.py',HERE/'contract.json',HERE.parent/'solver_04/raw.jsonl']+list(CODE.glob('*.py'))+[HERE.parent/'solver_02/run.py',HERE.parent/'solver_02/support.py']
    prov={'commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),'python':platform.python_version(),'command':[sys.executable]+sys.argv,'sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in source}}
    rows=[]
    with raw.open('x') as out:
        def record(r):rows.append(r);out.write(json.dumps(r)+'\n');out.flush()
        record({'kind':'provenance',**prov})
        try:
            f=field.Onb(5);curve=curves.Curve(f);truth,_=baseline.oracle(f,curve,3);checked=0
            for target in nagaocompare.affinePoints(f,curve):
                r=cell(5,3,[f.toCoords(x) for x in target],'enumerate',30)
                if r['status']!='complete' or {tuple(x) for x in r['solutions']}!=truth.get(target,set()):raise ArithmeticError('all-target validation mismatch')
                checked+=1
            record({'kind':'validation','all_affine_targets':checked,'complete_projected_oracle_equality':True})
        except Exception:record({'kind':'validation_failure','traceback':traceback.format_exc()});raise
        old=[json.loads(x) for x in (HERE.parent/'solver_04/raw.jsonl').read_text().splitlines()]
        for n,d in [(5,3),(9,4),(11,5)]:
            f=field.Onb(n);curve=curves.Curve(f);t=time.perf_counter();truth,_=baseline.oracle(f,curve,d)
            record({'kind':'oracle','n':n,'validation_setup_seconds':time.perf_counter()-t})
            for prior in old:
                if prior.get('kind')!='trial' or prior['n']!=n or prior['variant']!='rr-norm':continue
                try:
                    r=cell(n,d,prior['target'],prior['mode'],contract['seconds_per_instance'])
                    target=tuple(f.fromCoords(x) for x in prior['target']);expected=truth.get(target,set());got={tuple(x) for x in r['solutions']}
                    if not got<=expected or (r['status']=='complete' and got!=expected):raise ArithmeticError('oracle mismatch')
                    record({'kind':'trial','stratum':prior['stratum'],'expected_count':len(expected),**r})
                    print(n,prior['stratum'],prior['mode'],r['status'],len(got),flush=True)
                except Exception:record({'kind':'error','prior':prior,'traceback':traceback.format_exc()})
    groups=[]
    for n in (5,9,11):
        for mode in ('first','enumerate'):
            for st in ('uniform','known_decomposable'):
                rs=[r for r in rows if r['kind']=='trial' and (r['n'],r['mode'],r['stratum'])==(n,mode,st)]
                secs=sum(r['all_phase_seconds'] for r in rs);rels=sum(r['verified_unique_relations'] for r in rs)
                groups.append({'n':n,'mode':mode,'stratum':st,'attempted':len(rs),'resolved':sum(r['status']!='timeout' for r in rs),'verified_relations':rels,'all_phase_seconds':secs,'diagnostic_seconds_per_verified_relation':secs/rels if rels else None})
    (HERE/'summary.json').write_text(json.dumps({'provenance':prov,'groups':groups,'errors':[r for r in rows if r['kind']=='error']},indent=2)+'\n')

if __name__=='__main__':main()
