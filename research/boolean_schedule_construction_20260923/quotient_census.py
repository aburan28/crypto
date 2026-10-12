#!/usr/bin/env python3
"""Exact discovery census of necessary affine fibers after equation projection."""
import argparse
import datetime
import hashlib
import json
from pathlib import Path
import shutil

HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def basis(columns):
    rows=[]
    for value in columns:
        for pivot,row in rows:
            if value&pivot:value^=row
        if value:rows.append((value&-value,value))
    return rows
def project(value,rows):
    for pivot,row in rows:
        if value&pivot:value^=row
    return value
def coefficients(polys):
    out={}
    for e,poly in enumerate(polys):
        for monomial in poly:out[monomial]=out.get(monomial,0)^(1<<e)
    return {m:v for m,v in out.items() if v}
def analyze(polys,n,k):
    assert 0<=k<=n<=24 and len(polys)<=32
    original=coefficients(polys);low=(1<<k)-1
    quadratic=sorted(m for m in original if (m&low).bit_count()==2)
    rows=basis([original[m] for m in quadratic])
    transformed={m:project(value,rows) for m,value in original.items()}
    transformed={m:v for m,v in transformed.items() if v}
    assert all((m&low).bit_count()<=1 for m in transformed)
    return {'n':n,'low_variables':k,'equations':len(polys),'low_quadratic_columns':len(quadratic),
        'annihilated_rank':len(rows),'quotient_dimension':len(polys)-len(rows),
        'dimension_minus_low_variables':len(polys)-len(rows)-k,
        'elimination_basis':rows,'original_coefficients':original,'projected_coefficients':transformed,
        'all_low_quadratics_zero':all(project(original[m],rows)==0 for m in quadratic),
        'performance_ratio':None,'complete_solver_implemented':False}
def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--out',type=Path,required=True)
    out=parser.parse_args().out.resolve();out.mkdir(parents=True,exist_ok=False)
    source=HERE/'run_01'
    for name,digest in json.loads((source/'manifest.json').read_text())['files'].items():assert sha(source/name)==digest,name
    shutil.copyfile(Path(__file__),out/'quotient_census.py');cells=[]
    for n in [12,16,20,24]:
        for seed in [17,937]:
            for family in ['planted','cross_planted','unplanted']:
                path=source/f'n{n}-{seed}-{family}.jsonl'
                with path.open() as f:fixture=json.loads(next(f))
                for k in [4,5,6]:
                    cells.append({'seed':seed,'family':family,'input_sha256':sha(path),**analyze(fixture['polys'],n,k)})
    result={'schema_version':1,'recorded_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'scope':'Discovery structure only. Projected consistency is necessary, not sufficient; every original equation remains required.',
        'input_manifest_sha256':sha(source/'manifest.json'),'cells':cells}
    (out/'results.json').write_text(json.dumps(result,indent=2)+'\n')
    (out/'manifest.json').write_text(json.dumps({'files':{p.name:sha(p) for p in sorted(out.iterdir()) if p.is_file()}},indent=2)+'\n')
    print(json.dumps([{'n':n,'k':k,'ranks':sorted(set(c['annihilated_rank'] for c in cells if c['n']==n and c['low_variables']==k)),
        'quotient_dimensions':sorted(set(c['quotient_dimension'] for c in cells if c['n']==n and c['low_variables']==k))} for n in [12,16,20,24] for k in [4,5,6]]))
if __name__=='__main__':main()
