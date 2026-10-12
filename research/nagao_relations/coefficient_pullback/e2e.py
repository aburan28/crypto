"""Tiny cold ECDLP pipeline exercising the novel oracle without exposing secrets."""
import argparse
import hashlib
import json
from pathlib import Path
import random
import sys
import time
import traceback

import experiment as stage
import nagaoannihilator as scalar
import curves
import indexcalc
import image_solver
import pullback

HERE = Path(__file__).resolve().parent


class ScalarOps:
    def __init__(self, order):
        self.order = order
        self.counts = {'additions':0, 'multiplications':0, 'inversions':0}

    def add(self, a, b):
        self.counts['additions'] += 1
        return (a+b) % self.order

    def mul(self, a, b):
        self.counts['multiplications'] += 1
        return a*b % self.order

    def inv(self, a):
        self.counts['inversions'] += 1
        return pow(a, -1, self.order)


class Rows:
    def __init__(self, width, ops):
        self.width = width; self.ops = ops; self.pivots = {}

    def add(self, row):
        m = self.ops
        for pivot, basis in sorted(self.pivots.items()):
            if row[pivot]:
                scale = row[pivot]
                row = [m.add(x, -m.mul(scale,y)) for x,y in zip(row,basis)]
        pivot = next((i for i in range(self.width) if row[i]), None)
        if pivot is None:
            if row[-1]: raise ArithmeticError('inconsistent verified relations')
            return False
        inverse = m.inv(row[pivot])
        self.pivots[pivot] = [m.mul(x,inverse) for x in row]
        return True

    def solve(self):
        if len(self.pivots) != self.width: raise ValueError('rank deficient')
        m = self.ops; answer = [0]*self.width
        for pivot, row in sorted(self.pivots.items(),reverse=True):
            value = row[-1]
            for i in range(pivot+1,self.width):
                value = m.add(value,-m.mul(row[i],answer[i]))
            answer[pivot] = value
        return answer


def setup(panel, seed):
    f = scalar.CountedField(panel['n']); curve = scalar.CountedCurve(f)
    full = panel['curve_order']; order = panel['subgroup_order']; h = full//order
    points = stage.nagaocompare.affinePoints(f,curve)
    if len(points)+1 != full: raise ArithmeticError('wrong curve order')
    A = next(p for p in points if curve.mul(p,full) is None
             and curve.mul(p,full//2) is not None and curve.mul(p,full//order) is not None)
    G = curve.mul(A,h)
    if G is None or curve.mul(G,order) is not None: raise ArithmeticError('wrong subgroup')
    secret = random.Random(seed+999).randrange(1,order)
    Q = curve.mul(G,secret)
    return f,curve,A,G,Q,secret,h


def rho(curve, f, G, Q, ops, seed, deadline):
    rng = random.Random(seed); iterations = 0
    def step(state):
        p,a,b = state
        part = f.toCoords(p[0])%3 if p is not None else 0
        if part == 0: return curve.add(p,G),ops.add(a,1),b
        if part == 1: return curve.dbl(p),ops.mul(a,2),ops.mul(b,2)
        return curve.add(p,Q),a,ops.add(b,1)
    while time.perf_counter() < deadline:
        a,b = rng.randrange(ops.order),rng.randrange(ops.order)
        start = (curve.add(curve.mul(G,a),curve.mul(Q,b)),a,b)
        left = right = start
        for _ in range(4*ops.order):
            if time.perf_counter() >= deadline: raise TimeoutError('rho deadline')
            left = step(left); right = step(step(right)); iterations += 3
            if left[0] == right[0]:
                den = ops.add(right[2],-left[2])
                if den:
                    k = ops.mul(ops.add(left[1],-right[1]),ops.inv(den))
                    if curve.mul(G,k) == Q: return k,iterations
                break
    raise TimeoutError('rho deadline')


def run(panel, seed, variant, contract):
    started = time.perf_counter(); deadline = started+contract['run_seconds']
    f,curve,A,G,Q,secret,h = setup(panel,seed)
    order = panel['subgroup_order']; ops = ScalarOps(order)
    setupTime = time.perf_counter()-started
    attempts = []; phases = {'setup':setupTime,'factor_base':0.,'relation_collection':0.,'matrix':0.,'individual':0.,'rho_walk':0.,'final_verification':0.}
    output = {'panel':panel,'seed':seed,'variant':variant,'status':'incomplete','attempts':attempts,
              'generator':[f.toCoords(v) for v in G],'target':[f.toCoords(v) for v in Q],
              'subgroup_order':order,'full_dlp_S':None,'rho_ratio':None}
    try:
        if variant == 'rho':
            f.phase = 'rho_walk'; t=time.perf_counter()
            recovered, iterations = rho(curve,f,G,Q,ops,seed,deadline)
            phases['rho_walk'] += time.perf_counter()-t; output['rho_iterations']=iterations
        else:
            t=time.perf_counter(); f.phase='factor_base'
            mapping={}; columns={}
            for x in range(1,1<<panel['d']):
                p=curve.pointFromX(f.fromCoords(x))
                if p is None: continue
                projected=curve.mul(p,h)
                if projected is None: mapping[x]=None;continue
                negative=curve.neg(projected); key=min(projected,negative)
                if key not in columns:columns[key]=len(columns)
                mapping[x]=(columns[key],1 if projected==key else -1)
            phases['factor_base'] += time.perf_counter()-t
            matrix=Rows(len(columns),ops); rng=random.Random(seed)
            def decompose(target, phase, coefficient):
                if time.perf_counter()>=deadline:raise TimeoutError('DLP deadline')
                coords=[f.toCoords(x) for x in target]
                solver=pullback.cell if variant=='coefficient-pullback' else image_solver.cell
                result=solver(panel['n'],panel['d'],coords,'first',min(contract['decomposition_seconds'],max(.001,deadline-time.perf_counter())))
                attempts.append({'phase':phase,'coefficient':coefficient,'result':result})
                if not result['solutions']:return None
                xs=result['solutions'][0];f.phase=phase+'_verification'
                lifted=indexcalc.liftAndCheck(f,curve,xs,target)
                if lifted is None:raise ArithmeticError('independent relation check failed')
                _,signs=lifted;row=[0]*len(columns)
                for x,sign in zip(xs,signs):
                    if mapping[x] is not None:
                        i,orientation=mapping[x];row[i]=ops.add(row[i],sign*orientation)
                return row
            for _ in range(contract['target_attempt_cap']):
                t=time.perf_counter();f.phase='relation_generation'
                s=rng.randrange(1,panel['curve_order']);target=curve.mul(A,s)
                row=decompose(target,'relations',s)
                phases['relation_collection'] += time.perf_counter()-t
                if row is None:continue
                t=time.perf_counter();matrix.add(row+[s%order]);phases['matrix']+=time.perf_counter()-t
                if len(matrix.pivots)==len(columns):break
            output.update(columns=len(columns),rank=len(matrix.pivots))
            t=time.perf_counter();logs=matrix.solve();phases['matrix']+=time.perf_counter()-t
            rng=random.Random(seed+100000);recovered=None
            t=time.perf_counter()
            for _ in range(contract['target_attempt_cap']):
                f.phase='individual_generation';s=rng.randrange(panel['curve_order'])
                target=curve.add(Q,curve.mul(A,s))
                if target is None:continue
                row=decompose(target,'individual',s)
                if row is None:continue
                value=0
                for coefficient,log in zip(row,logs):value=ops.add(value,ops.mul(coefficient,log))
                recovered=ops.mul(ops.add(value,-s),ops.inv(h%order));break
            phases['individual']+=time.perf_counter()-t
            if recovered is None:raise TimeoutError('individual attempt cap')
        f.phase='final_verification';t=time.perf_counter()
        if recovered != secret or curve.mul(G,recovered) != Q:raise ArithmeticError('scalar recovery failed')
        phases['final_verification']+=time.perf_counter()-t
        output.update(status='verified',recovered_scalar=recovered,planted_scalar=secret)
    except (TimeoutError,ValueError) as exc:
        output.update(status='incomplete',error=str(exc))
    output.update(cold_wall_seconds=time.perf_counter()-started,phase_seconds=phases,
                  orchestration_field_api_counts=f.report(),scalar_modular_counts=ops.counts)
    totals={k:output['orchestration_field_api_counts']['totals'][k] for k in ['additions','multiplications','squarings','fieldOperations']}
    for attempt in attempts:
        for key in totals:totals[key]+=attempt['result']['field_api_counts']['totals'][key]
    output['all_phase_field_api_counts']=totals
    return output


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=False)
    contract=json.loads((HERE/'e2e_contract.json').read_text());rows=[]
    sources=list(HERE.glob('*.py'))+[HERE/'e2e_contract.json',HERE.parent/'solver_08/image_solver.py']+list(stage.CODE.glob('*.py'))
    provenance={'source_sha256':{str(p.relative_to(stage.ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},'contract':contract}
    (args.output/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    with (args.output/'raw.jsonl').open('x') as out:
        for panel in contract['panels']:
            for seed in contract['seeds']:
                variants=list(contract['variants']);random.Random(seed).shuffle(variants)
                for variant in variants:
                    try:row=run(panel,seed,variant,contract)
                    except Exception:row={'panel':panel,'seed':seed,'variant':variant,'status':'error','traceback':traceback.format_exc()}
                    rows.append(row);out.write(json.dumps(row)+'\n');out.flush()
                    print(panel,seed,variant,row['status'],row.get('cold_wall_seconds'),flush=True)
    compact=[{k:v for k,v in row.items() if k!='attempts'} for row in rows]
    (args.output/'summary.json').write_text(json.dumps({'runs':compact,'contract':contract},indent=2)+'\n')
    if any(row['status']!='verified' for row in rows):raise SystemExit(1)


if __name__=='__main__':main()
