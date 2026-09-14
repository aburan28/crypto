"""Fresh-target first-relation and complete-enumeration tests; one all-phase budget."""
import hashlib
import importlib.metadata
import json
from pathlib import Path
import platform
import random
import subprocess
import sys
import time
import traceback

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
CODE=ROOT/'ecc2k130/codegen'
sys.path.insert(0,str(CODE))
sys.path.insert(0,str(HERE.parent/'solver_02'))
import run as previous
import fixedu
import cnf
import curves
import decomp
import field
import indexcalc
import nagaocompare
import nagaodecomp
import nagaoannihilator as scalar
import support


def validate():
    f=field.Onb(5);d=3;target=nagaocompare.affinePoints(f,curves.Curve(f))[10]
    p,out=fixedu.build(5);checked=0;accepted=0
    ls=scalar.subspacePolynomial(f,d)
    for aa in range(32):
        for bb in range(1,32):
            h,_=scalar.residualNorm(f,target,f.fromCoords(aa),f.fromCoords(bb))
            cubic=not any(scalar.annihilatorRemainder(f,ls,h)) and bool(h[0]) and bool(scalar.polyEval(f,h,target[0]))
            anyAccept=False
            for uu in range(1,8):
                u=f.fromCoords(uu);x=f.add(h[2],u);v=f.add(h[1],f.mul(x,u))
                vals={'a':aa,'b':bb,'r':f.toCoords(target[0]),'s':f.toCoords(target[1]),'u':uu}
                bits={(name,i):val>>i&1 for name,val in vals.items() for i in range(5)}
                evaluated=p.evaluate(bits,out)
                got=[sum((evaluated[k*5+i] or 0)<<i for i in range(5)) for k in range(6)]
                want=[x,v,scalar.polyEval(f,h,x),h[0],scalar.polyEval(f,h,target[0]),scalar.polyEval(f,[v,u,f.fromCoords(31)],x)]
                if got!=[f.toCoords(t) for t in want]:raise ArithmeticError('fixed-u vector mismatch')
                membership=fixedu.reduceBits(got[1],fixedu.imageBasis(f,d,uu))==0
                accept=got[0]<8 and got[2]==0 and all(got[3:]) and membership
                anyAccept|=accept;accepted+=bool(accept);checked+=1
            if anyAccept!=cubic:raise ArithmeticError('fixed-u support equivalence mismatch')
    return {'coefficient_branch_checks':checked,'admitted_branch_models':accepted,'vectors_and_existential_support_match':True}


def runCell(f,curve,d,target,expected,variant,mode,budget):
    start=time.perf_counter();deadline=start+budget;n=f.m
    phases={'construct_load':0.0,'solve':0.0,'extract_verify':0.0}
    solutions=set();calls=0;branches=0;rejected=0;duplicates=0;allDone=True
    firstAt=None;cnfTotals={'vars':0,'clauses':0,'xors':0};branchRows=[]
    us=list(range(1,1<<d)) if variant=='fixed-u' else [None]
    for j,u in enumerate(us):
        if time.perf_counter()>=deadline:allDone=False;break
        branchStart=time.perf_counter()
        branchDeadline=min(deadline,branchStart+(deadline-branchStart)/(len(us)-j))
        c=cnf.Cnf();xr=f.toCoords(target[0]);yr=f.toCoords(target[1])
        if variant=='fixed-u':
            p,outputs=fixedu.build(n);variables=fixedu.encode(p,outputs,f,d,target,u,c)
        elif variant=='rr-norm':
            p,outputs=nagaodecomp.buildSystem(n,f.n,min(n,12),'norm')
            variables=nagaodecomp.encode(p,outputs,n,n,xr,yr,c,'norm')
        else:
            p,outputs=decomp.buildSystem(n,f.n,3,min(n,12))
            xs=decomp.encode(p,outputs,n,3,n,xr,c,orderPoints=False)
            nagaodecomp.addRestrictedDomain(c,xs,xr);variables={'pvars':xs}
        if variant!='fixed-u':
            for xs in variables['pvars']:
                for bit in xs[d:]:c.addClause([-bit])
        solver=previous.loadSolver(c);built=time.perf_counter()
        phases['construct_load']+=built-branchStart;branches+=1
        for key,value in c.stats().items():cnfTotals[key]+=value
        done=False;branchCalls=0
        while time.perf_counter()<branchDeadline:
            before=time.perf_counter()
            ok,model=solver.solve(time_limit=max(.001,branchDeadline-before))
            solved=time.perf_counter();phases['solve']+=solved-before;calls+=1;branchCalls+=1
            if ok is None:break
            if not ok:done=True;break
            vals={i:v for i,v in enumerate(model) if i}
            if not nagaocompare.cnfSatisfied(c,vals):raise ArithmeticError('CNF model invalid')
            if variant=='fixed-u':
                dec=lambda key:sum(int(vals[v])<<i for i,v in enumerate(variables[key]))
                xs=support.recover(f,curve,d,target,dec('a'),dec('b'))
                block=[v for key in ('a','b') for v in variables[key]]
            else:
                xs=sorted(sum(int(vals[v])<<i for i,v in enumerate(x)) for x in variables['pvars'])
                block=[v for x in variables['pvars'] for v in x]
                if variant=='rr-norm' and not nagaodecomp.reconstructWitness(f,curve,nagaodecomp.decode(vals,variables,c),target):raise ArithmeticError('invalid norm witness')
            if any(x==0 or x>=1<<d or x==xr for x in xs) or len(set(xs))!=3:raise ArithmeticError('invalid domain')
            valid=indexcalc.liftAndCheck(f,curve,xs,target) is not None
            if valid:
                if tuple(xs) not in expected:raise ArithmeticError('oracle mismatch')
                if tuple(xs) in solutions:duplicates+=1
                solutions.add(tuple(xs))
                if firstAt is None:firstAt=time.perf_counter()-start
            else:
                if variant!='chained-s3':raise ArithmeticError('invalid RR relation')
                rejected+=1
            solver.add_clause([-v if vals[v] else v for v in block])
            phases['extract_verify']+=time.perf_counter()-solved
            if valid and mode=='first':break
        branchRows.append({'u':u,'calls':branchCalls,'exhausted':done,'elapsed_seconds':time.perf_counter()-branchStart,'cnf':c.stats()})
        if not done:allDone=False
        if solutions and mode=='first':break
    if allDone and solutions!=expected:raise ArithmeticError('complete set differs from exhaustive oracle')
    status='first' if mode=='first' and solutions else ('complete' if allDone else 'timeout')
    return {'variant':variant,'mode':mode,'n':n,'d':d,'target':[xr,yr],'status':status,
        'expected_count':len(expected),'solutions':[list(x) for x in sorted(solutions)],
        'verified_unique_relations':len(solutions),'first_verified_seconds':firstAt,
        'all_phase_seconds':time.perf_counter()-start,'phase_seconds':phases,
        'branches':branches,'calls':calls,'duplicate_models':duplicates,'rejected_candidates':rejected,
        'branch_records':branchRows,'constructed_cnf_totals':cnfTotals,
        'solver_operations':None,'full_dlp_S':None,'rho_ratio':None}


def main():
    raw=HERE/'raw.jsonl';summary=HERE/'summary.json'
    if raw.exists() or summary.exists():raise SystemExit('Evidence exists; choose a new sibling folder')
    contract=json.loads((HERE/'contract.json').read_text())
    sources=[HERE/'run.py',HERE/'fixedu.py',HERE/'contract.json',HERE.parent/'solver_02/run.py',HERE.parent/'solver_02/support.py']+list(CODE.glob('*.py'))
    provenance={'commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        'command':[sys.executable]+sys.argv,'python':platform.python_version(),'pycryptosat':importlib.metadata.version('pycryptosat'),
        'sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}}
    old=[json.loads(x) for x in (HERE.parent/'solver_02/raw.jsonl').read_text().splitlines()]
    rows=[]
    with raw.open('x') as out:
        def record(row):rows.append(row);out.write(json.dumps(row)+'\n');out.flush()
        record({'kind':'provenance',**provenance})
        try:check=validate();record({'kind':'validation',**check});print(check,flush=True)
        except Exception:record({'kind':'validation_failure','traceback':traceback.format_exc()});raise
        for panel in contract['panels']:
            n,d=panel['n'],panel['d'];f=field.Onb(n);curve=curves.Curve(f)
            t=time.perf_counter();truth,base=previous.oracle(f,curve,d)
            excluded={tuple(r['target']) for r in old if r.get('kind')=='trial' and r['n']==n}
            coords=lambda t:tuple(f.toCoords(v) for v in t)
            rng=random.Random(contract['seed']+n)
            uniform=rng.sample([t for t in nagaocompare.affinePoints(f,curve) if coords(t) not in excluded],4)
            supported=rng.sample([t for t in truth if coords(t) not in excluded],4)
            record({'kind':'oracle','n':n,'d':d,'setup_validation_seconds':time.perf_counter()-t,'base_abscissas':base,'supported_targets':len(truth),'excluded_old_targets':sorted(excluded)})
            for stratum,targets in [('uniform',uniform),('known_decomposable',supported)]:
                for target in targets:
                    for mode in ('first','enumerate'):
                        for variant in contract['variants']:
                            try:
                                r=runCell(f,curve,d,target,truth.get(target,set()),variant,mode,contract['seconds_per_instance'])
                                record({'kind':'trial','stratum':stratum,**r});print(n,stratum,mode,variant,r['status'],r['verified_unique_relations'],flush=True)
                            except Exception:record({'kind':'error','n':n,'stratum':stratum,'mode':mode,'variant':variant,'target':coords(target),'traceback':traceback.format_exc()})
    groups=[]
    for panel in contract['panels']:
        for variant in contract['variants']:
            for mode in ('first','enumerate'):
                for st in ('uniform','known_decomposable'):
                    rr=[r for r in rows if r['kind']=='trial' and (r['n'],r['variant'],r['mode'],r['stratum'])==(panel['n'],variant,mode,st)]
                    seconds=sum(r['all_phase_seconds'] for r in rr);relations=sum(r['verified_unique_relations'] for r in rr)
                    groups.append({'n':panel['n'],'d':panel['d'],'variant':variant,'mode':mode,'stratum':st,'attempted':len(rr),
                        'resolved':sum(r['status']!='timeout' for r in rr),'timeouts':sum(r['status']=='timeout' for r in rr),
                        'verified_relations':relations,'all_phase_seconds':seconds,'diagnostic_seconds_per_verified_relation':seconds/relations if relations else None})
    summary.write_text(json.dumps({'provenance':provenance,'validation':check,'groups':groups,'errors':[r for r in rows if r['kind']=='error'],'limits':contract['metrics']},indent=2)+'\n')


if __name__=='__main__':main()
