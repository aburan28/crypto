"""Frozen matched-subspace first-relation experiment, with exact toy oracles."""
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
import cnf as cnfmod
import curves
import decomp
import field
import indexcalc
import nagaocompare
import nagaodecomp
import nagaoannihilator as scalar
import pycryptosat
import support


def oracle(f,curve,d):
    base={x:curve.pointFromX(f.fromCoords(x)) for x in range(1,1<<d)}
    base={x:p for x,p in base.items() if p is not None}
    xs=sorted(base)
    out={}
    for i in range(len(xs)):
        for j in range(i+1,len(xs)):
            for k in range(j+1,len(xs)):
                coords=(xs[i],xs[j],xs[k])
                for mask in range(8):
                    total=None
                    for bit,x in enumerate(coords):
                        p=base[x]
                        total=curve.add(total,curve.neg(p) if mask>>bit&1 else p)
                    if total is not None and f.toCoords(total[0]) not in coords:
                        out.setdefault(total,set()).add(coords)
    return out,len(base)


def validateCircuit():
    f=field.Onb(5);curve=curves.Curve(f);d=3
    target=nagaocompare.affinePoints(f,curve)[10]
    p,roots,nz=support.build(5,d)
    ls=scalar.subspacePolynomial(f,d)
    checked=0
    for a in range(32):
        for b in range(1,32):
            vals=support.inputs(f,d,target,a,b)
            bits={(name,i):value>>i&1 for name,value in vals.items() for i in range(5)}
            got=p.evaluate(bits,roots+nz)
            vectors=[sum((got[k*5+i] or 0)<<i for i in range(5)) for k in range(5)]
            h,_=scalar.residualNorm(f,target,f.fromCoords(a),f.fromCoords(b))
            remainder=scalar.annihilatorRemainder(f,ls,h)
            want=[f.toCoords(v) for v in remainder+[0]*(3-len(remainder))]
            want.extend([f.toCoords(h[0]),f.toCoords(scalar.polyEval(f,h,target[0]))])
            if vectors!=want:
                raise ArithmeticError('circuit/scalar mismatch at a=%d b=%d'%(a,b))
            checked+=1
    c=cnfmod.Cnf();v=support.encode(p,roots,nz,f,d,target,c)
    s=loadSolver(c)
    result,_=s.solve(assumptions=[-x for x in v['b']])
    if result is not False:
        raise ArithmeticError('zero-b chart was not excluded')
    return {'checked_coefficient_pairs':checked,'target':[f.toCoords(v) for v in target],
            'all_remainder_and_nonzero_vectors_match':True,'zero_b_rejected':True}


def loadSolver(c):
    s=pycryptosat.Solver(threads=1)
    for row in c.clauses:s.add_clause(row)
    for row,rhs in c.xors:s.add_xor_clause(row,rhs)
    return s


def cell(f,curve,d,target,expected,variant,budget):
    start=time.perf_counter();n=f.m;c=cnfmod.Cnf()
    xr,yr=[f.toCoords(v) for v in target]
    if variant.startswith('support'):
        p,roots,nz=support.build(n,d)
        variables=support.encode(p,roots,nz,f,d,target,c)
        measuredRoots=roots+nz
    else:
        if variant=='chained-s3':
            p,roots=decomp.buildSystem(n,f.n,3,min(n,12))
            xs=decomp.encode(p,roots,n,3,n,xr,c,orderPoints=False)
            nagaodecomp.addRestrictedDomain(c,xs,xr)
            variables={'pvars':xs}
        else:
            p,roots=nagaodecomp.buildSystem(n,f.n,min(n,12),'norm')
            variables=nagaodecomp.encode(p,roots,n,n,xr,yr,c,'norm')
        for xs in variables['pvars']:
            for bit in xs[d:]:c.addClause([-bit])
        measuredRoots=roots
    built=time.perf_counter()
    s=loadSolver(c);loaded=time.perf_counter()
    deadline=loaded+budget
    status='timeout';solution=None;calls=0;rejected=0;branches=0
    assumptions=[[]]
    if variant=='support-conditioned':
        assumptions=[[v if value>>i&1 else -v for i,v in enumerate(variables['b'][:2])]
                     for value in range(4)]
    exhausted=True
    for assume in assumptions:
        branches+=1
        branchDone=False
        while time.perf_counter()<deadline:
            ok,model=s.solve(assumptions=assume,time_limit=max(0.001,deadline-time.perf_counter()))
            calls+=1
            if ok is None:break
            if not ok:branchDone=True;break
            vals={i:v for i,v in enumerate(model) if i}
            if not nagaocompare.cnfSatisfied(c,vals):raise ArithmeticError('invalid CNF model')
            if variant.startswith('support'):
                decode=lambda name:sum(int(vals[v])<<i for i,v in enumerate(variables[name]))
                xs=support.recover(f,curve,d,target,decode('a'),decode('b'))
                block=[v for name in ('a','b') for v in variables[name]]
            else:
                xs=sorted(sum(int(vals[v])<<i for i,v in enumerate(x)) for x in variables['pvars'])
                if variant=='rr-norm':
                    cert=nagaodecomp.decode(vals,variables,c)
                    if not nagaodecomp.reconstructWitness(f,curve,cert,target):
                        raise ArithmeticError('norm certificate invalid')
                block=[v for x in variables['pvars'] for v in x]
            if any(x==0 or x>=1<<d or x==xr for x in xs) or len(set(xs))!=3:
                raise ArithmeticError('factor-base restriction violated')
            if indexcalc.liftAndCheck(f,curve,xs,target) is not None:
                if tuple(xs) not in expected:raise ArithmeticError('oracle does not contain relation')
                solution=xs;status='sat';break
            if variant!='chained-s3':raise ArithmeticError('non-liftable RR model')
            rejected+=1;s.add_clause([-v if vals[v] else v for v in block])
        if status=='sat':break
        if not branchDone:exhausted=False;break
    if status!='sat' and exhausted:
        status='unsat'
        if expected:raise ArithmeticError('false UNSAT: oracle has a decomposition')
    ended=time.perf_counter()
    return {'n':n,'d':d,'variant':variant,'target':[xr,yr],'expected_sat':bool(expected),
            'expected_projected_solutions':len(expected),'status':status,'verified':status!='timeout',
            'solution':solution,'calls':calls,'branches':branches,'rejected_candidates':rejected,
            'cnf':c.stats(),'ir_bit_operations_per_evaluation':p.bitOpCount(measuredRoots),
            'diagnostic_seconds':{'construct':built-start,'load':loaded-built,'solve_verify':ended-loaded},
            'solver_operations':None,'full_dlp_S':None,'rho_ratio':None}


def main():
    raw=HERE/'raw.jsonl';summary=HERE/'summary.json'
    if raw.exists() or summary.exists():raise SystemExit('Evidence exists; use a new sibling directory')
    contract=json.loads((HERE/'contract.json').read_text())
    sources=[HERE/'run.py',HERE/'support.py',HERE/'contract.json']+[CODE/x for x in
        ('decomp.py','cnf.py','field.py','curves.py','ir.py','build.py','nagaodecomp.py',
         'nagaocompare.py','nagaoannihilator.py','indexcalc.py')]
    provenance={'command':[sys.executable]+sys.argv,'python':platform.python_version(),
        'pycryptosat':importlib.metadata.version('pycryptosat'),
        'commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        'dirty':subprocess.check_output(['git','status','--porcelain'],cwd=ROOT,text=True),
        'sha256':{str(x.relative_to(ROOT)):hashlib.sha256(x.read_bytes()).hexdigest() for x in sources}}
    rows=[]
    with raw.open('x') as out:
        def record(r):rows.append(r);out.write(json.dumps(r)+'\n');out.flush()
        record({'kind':'provenance',**provenance})
        try:
            check=validateCircuit();record({'kind':'circuit_validation',**check});print(check,flush=True)
        except Exception:
            record({'kind':'validation_failure','traceback':traceback.format_exc()});raise
        for panel in contract['panels']:
            f=field.Onb(panel['n']);c=curves.Curve(f);d=panel['d']
            started=time.perf_counter();truth,baseSize=oracle(f,c,d)
            rng=random.Random(contract['seed']+panel['n'])
            uniform=rng.sample(nagaocompare.affinePoints(f,c),4)
            planted=rng.sample(list(truth),min(4,len(truth)))
            record({'kind':'oracle','panel':panel,'base_abscissas':baseSize,'supported_targets':len(truth),
                    'diagnostic_seconds':time.perf_counter()-started})
            for label,targets in [('uniform',uniform),('known_decomposable',planted)]:
                for target in targets:
                    for variant in contract['variants']:
                        try:
                            r=cell(f,c,d,target,truth.get(target,set()),variant,contract['seconds_per_instance'])
                            record({'kind':'trial','stratum':label,**r})
                            print(panel,label,variant,r['status'],flush=True)
                        except Exception:
                            record({'kind':'error','panel':panel,'stratum':label,'variant':variant,
                                    'target':[f.toCoords(v) for v in target],'traceback':traceback.format_exc()})
    groups=[]
    for panel in contract['panels']:
        for variant in contract['variants']:
            for label in ('uniform','known_decomposable'):
                rs=[r for r in rows if r['kind']=='trial' and r['n']==panel['n'] and r['variant']==variant and r['stratum']==label]
                groups.append({'n':panel['n'],'d':panel['d'],'variant':variant,'stratum':label,'attempted':len(rs),
                    **{status:sum(r['status']==status for r in rs) for status in ('sat','unsat','timeout')},
                    'diagnostic_total_seconds':sum(sum(r['diagnostic_seconds'].values()) for r in rs),
                    'cnf_variables':rs[0]['cnf']['vars'] if rs else None})
    with summary.open('x') as out:
        json.dump({'provenance':provenance,'validation':check,'groups':groups,
            'errors':[r for r in rows if r['kind']=='error'],'limits':contract['metrics']},out,indent=2)


if __name__=='__main__':main()
