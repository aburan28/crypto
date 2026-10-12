"""Held-out beyond-eleven-bit paired campaign, with an independent pair oracle."""
import hashlib
import importlib.util
import json
from pathlib import Path
import platform
import random
import subprocess
import sys
import time
import traceback
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2];CODE=ROOT/'ecc2k130/codegen'
sys.path.insert(0,str(CODE));sys.path.insert(0,str(HERE.parent/'solver_02'))
import curves
import field
import nagaocompare
import nagaoannihilator as scalar
import optimized


def load(name,path):
    sys.path.insert(0,str(path.parent))
    spec=importlib.util.spec_from_file_location(name,path);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module
old=load('old_quadratic',HERE.parent/'solver_05/run.py')
s4=load('matched_s4',HERE.parent/'solver_06/run.py')


class PairOracle:
    def __init__(self,n,d):
        self.f=field.Onb(n);self.curve=curves.Curve(self.f);f=self.f;c=self.curve
        self.base={}
        for xBits in range(1,1<<d):
            point=c.pointFromX(f.fromCoords(xBits))
            if point is not None:self.base[xBits]=[point,c.neg(point)]
        self.pairs={};xs=sorted(self.base)
        for i,x in enumerate(xs):
            for y in xs[i+1:]:
                for p in self.base[x]:
                    for q in self.base[y]:
                        self.pairs.setdefault(c.add(p,q),[]).append((x,y))

    def expected(self,target):
        f=self.f;c=self.curve;r=f.toCoords(target[0]);out=set()
        for z,points in self.base.items():
            if z==r:continue
            for point in points:
                residual=c.add(target,c.neg(point))
                for x,y in self.pairs.get(residual,[]):
                    if z not in (x,y) and r not in (x,y):out.add(tuple(sorted((x,y,z))))
        return out

    def randomUniform(self,rng):
        f=self.f;c=self.curve
        while True:
            x=rng.randrange(1<<f.m);sign=rng.randrange(2)
            if x==0:
                if sign==0:return (0,f.one())
                continue
            point=c.pointFromX(f.fromCoords(x))
            if point is not None:return c.neg(point) if sign else point

    def randomSupported(self,rng):
        f=self.f;c=self.curve
        while True:
            xs=rng.sample(sorted(self.base),3);total=None
            for x in xs:total=c.add(total,self.base[x][rng.randrange(2)])
            if total is not None and f.toCoords(total[0]) not in xs:return total


def validate():
    oracle=PairOracle(5,3);f=oracle.f;c=oracle.curve;checked=0
    # Independently compare pair lookup with the earlier signed-triple enumeration.
    truth,_=old.baseline.oracle(f,c,3)
    for target in nagaocompare.affinePoints(f,c):
        expected=oracle.expected(target)
        if expected!=truth.get(target,set()):raise ArithmeticError('pair/triple oracle mismatch')
        r=optimized.cell(5,3,[f.toCoords(x) for x in target],'enumerate',30)
        if r['status']!='complete' or {tuple(x) for x in r['solutions']}!=expected:raise ArithmeticError('optimized exhaustive mismatch')
        # Candidate sequence equality catches cancellation and batch-inversion mistakes,
        # including coefficients rejected by the final support predicate.
        search=optimized.Search(f,c,3,target,time.perf_counter()+30)
        new=[(a,b) for a,b,_ in search.candidates()]
        prior=list(old.quadratic.candidates(f,c,3,target))
        if new!=prior:raise ArithmeticError('candidate sequence changed')
        checked+=1
    return {'all_affine_five_bit_targets':checked,'pair_oracle_equals_triple_oracle':True,'complete_sets_and_candidate_sequences_match':True}


def main():
    raw=HERE/'raw.jsonl';summary=HERE/'summary.json'
    if raw.exists() or summary.exists():raise SystemExit('Evidence exists; use a new folder')
    contract=json.loads((HERE/'contract.json').read_text())
    source=list(CODE.glob('*.py'))+[HERE/'run.py',HERE/'optimized.py',HERE/'contract.json']
    for folder in ('solver_02','solver_05','solver_06'):source+=list((HERE.parent/folder).glob('*.py'))
    source+=[ROOT/'src/cryptanalysis/koblitz_symmetrised.rs']
    excluded={}
    historical=[]
    for folder in ('solver_02','solver_04'):
        path=HERE.parent/folder/'raw.jsonl';source.append(path)
        for line in path.read_text().splitlines():
            row=json.loads(line)
            if row.get('kind')=='trial':excluded.setdefault(row['n'],set()).add(tuple(row['target']))
    prov={'commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),'python':platform.python_version(),'command':[sys.executable]+sys.argv,'sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in source}}
    rows=[]
    with raw.open('x') as out:
        def record(r):rows.append(r);out.write(json.dumps(r)+'\n');out.flush()
        record({'kind':'provenance',**prov})
        try:record({'kind':'validation',**validate()});print('validation passed',flush=True)
        except Exception:record({'kind':'validation_failure','traceback':traceback.format_exc()});raise
        for panel in contract['panels']:
            n,d=panel['n'],panel['d'];t=time.perf_counter();oracle=PairOracle(n,d);f=oracle.f;c=oracle.curve
            record({'kind':'oracle','n':n,'d':d,'signed_base_size':2*len(oracle.base),'pair_entries':sum(len(v) for v in oracle.pairs.values()),'validation_setup_seconds':time.perf_counter()-t})
            rng=random.Random(contract['seed']+n*100+d)
            for stratum in ('uniform','known_decomposable'):
                targets=[]
                while len(targets)<contract['targets_per_stratum']:
                    target=oracle.randomUniform(rng) if stratum=='uniform' else oracle.randomSupported(rng)
                    coords=tuple(f.toCoords(v) for v in target)
                    if coords not in excluded.get(n,set()) and coords not in targets:targets.append(coords)
                record({'kind':'targets','n':n,'d':d,'stratum':stratum,'targets':targets})
                for coords in targets:
                    target=tuple(f.fromCoords(v) for v in coords);expected=oracle.expected(target)
                    for mode in ('first','enumerate'):
                        # Rotate order deterministically to reduce fixed-order bias.
                        variants=contract['variants'][:];rng.shuffle(variants)
                        for variant in variants:
                            try:
                                if variant=='quadratic-optimized':r=optimized.cell(n,d,coords,mode,contract['seconds_per_instance'])
                                elif variant=='quadratic-original':r=old.cell(n,d,coords,mode,contract['seconds_per_instance']);r['variant']=variant
                                else:r=s4.runCell(f,c,d,target,expected,'s4-symmetric',mode,contract['seconds_per_instance'])
                                got={tuple(x) for x in r['solutions']}
                                if not got<=expected or (r['status']=='complete' and got!=expected):raise ArithmeticError('benchmark oracle mismatch')
                                r['within_budget']=r['all_phase_seconds']<=contract['seconds_per_instance']
                                record({'kind':'trial','stratum':stratum,'expected_count':len(expected),**r})
                                print(n,d,stratum,mode,variant,r['status'],len(got),round(r['all_phase_seconds'],3),flush=True)
                            except Exception:record({'kind':'error','n':n,'d':d,'stratum':stratum,'mode':mode,'variant':variant,'target':coords,'traceback':traceback.format_exc()})
    groups=[]
    for panel in contract['panels']:
        for variant in contract['variants']:
            for mode in ('first','enumerate'):
                for st in ('uniform','known_decomposable'):
                    rs=[r for r in rows if r['kind']=='trial' and (r['n'],r['d'],r['variant'],r['mode'],r['stratum'])==(panel['n'],panel['d'],variant,mode,st)]
                    seconds=sum(r['all_phase_seconds'] for r in rs);relations=sum(r['verified_unique_relations'] for r in rs)
                    ops=sum(r['field_api_counts']['totals']['fieldOperations'] for r in rs) if variant.startswith('quadratic') else None
                    groups.append({**panel,'variant':variant,'mode':mode,'stratum':st,'attempted':len(rs),'resolved':sum(r['status']!='timeout' for r in rs),'resolved_within_budget':sum(r['status']!='timeout' and r['within_budget'] for r in rs),'verified_relations':relations,'all_phase_seconds':seconds,'diagnostic_seconds_per_relation':seconds/relations if relations else None,'field_api_operations':ops,'full_dlp_S':None,'rho_ratio':None})
    summary.write_text(json.dumps({'provenance':prov,'groups':groups,'errors':[r for r in rows if r['kind']=='error'],'limits':contract['limits']},indent=2)+'\n')

if __name__=='__main__':main()
