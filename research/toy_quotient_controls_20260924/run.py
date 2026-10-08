"""Frozen bounded experiment. No external curve/target or size parameters."""
import argparse
import gzip
import hashlib
import json
import platform
import random
import sys
import time
from collections import defaultdict
from itertools import product
from pathlib import Path
import algebra as alg
from curves import certified_neighbors, O
from matrix import measure
from independent import audit

HERE = Path(__file__).resolve().parent
PHASES = ('setup','encoding','solve_extract','reconstruct','verify')


def digest(data):
    return hashlib.sha256(json.dumps(data,sort_keys=True,separators=(',',':')).encode()).hexdigest()


def opcode_call(function, *args):
    """One exclusive phase, counting every Python callee, never the observer."""
    count = 0
    previous = sys.gettrace()
    if previous is not None:
        raise RuntimeError('opcode profiler requires no existing trace hook')
    def observer(frame,event,arg):
        nonlocal count
        if event=='call':
            frame.f_trace_opcodes = True
            return observer
        if event=='opcode':
            count += 1
        return observer
    frame = sys._getframe()
    old_opcode_flag = frame.f_trace_opcodes
    frame.f_trace_opcodes = True  # required before settrace on CPython 3.12
    try:
        sys.settrace(observer)
        result = function(*args)
    finally:
        sys.settrace(None)
        frame.f_trace_opcodes = old_opcode_flag
    return result,count


def pipeline(support,target,b,m,variant,contract,instrument=False):
    phase_counts, durations = {}, {}
    def phase(name,function,*args):
        start = time.perf_counter()
        if instrument:
            value,phase_counts[name] = opcode_call(function,*args)
        else:
            value = function(*args)
        durations[name] = time.perf_counter()-start
        return value
    prepared = phase('setup',alg.setup,support,m,variant,contract)
    ring,generators = phase('encoding',alg.encode,support,target,b,m,variant,prepared,contract)
    solved = phase('solve_extract',alg.solve,ring,generators)
    slots = phase('reconstruct',alg.reconstruct,solved['roots'],support,m,variant,prepared)
    verified,lift_checks = phase('verify',alg.verify,slots,support,target)
    payload = dict(n=ring.n, generators=generators, input_sha256=digest(generators),
        quotient_domain=prepared.get('domain'), solutions=verified,
        verified_solution_sha256=digest(verified), lift_checks=lift_checks, **solved)
    return payload,phase_counts,durations


def oracle(support,target,m):
    points = [point for pair in support for point in pair]
    return sorted({tuple(sorted(ps)) for ps in product(points,repeat=m)
                   if alg.Curve().total(ps)==target})


def run(output):
    contract = json.loads((HERE/'contract.json').read_text())
    output.mkdir(parents=True,exist_ok=True)
    if (output/'raw.json.gz').exists() or (output/'partial.jsonl').exists():
        raise FileExistsError('use a new results directory; evidence is immutable')
    source,h,neighbors = certified_neighbors()
    models = [('source',source,lambda p:p)] + [(f'neighbor_{phi.target.b}',phi.target,phi) for phi,_ in neighbors]
    choices = sorted({min(p,source.neg(p)) for p in h if p!=O and p[0]!=0})
    cases = []
    with (output/'partial.jsonl').open('w') as journal:
        for split,seed in contract['seeds'].items():
            reps = random.Random(seed).sample(choices,4)
            pairs = [(p,source.neg(p)) for p in reps]
            for model,curve,mapping in models:
                support = [tuple(mapping(p) for p in pair) for pair in pairs]
                for m in contract['summands']:
                    for ti in contract['target_indices']:
                        target = mapping(h[ti])
                        truth = oracle(support,target,m)
                        for variant in contract['variants']:
                            identity = dict(split=split,seed=seed,model=model,b=curve.b,
                                            summands=m,target_index=ti,target=target,support=support,variant=variant)
                            try:
                                replays = [pipeline(support,target,curve.b,m,variant,contract)
                                           for _ in range(contract['repetitions'])]
                                result = dict(replays[0][0])
                                if any(replay[0]!=result for replay in replays[1:]):
                                    raise AssertionError('repetition mismatch')
                                if result['solutions']!=truth:
                                    raise AssertionError('verified signed multiset coverage differs from exhaustive oracle')
                                profiled,opcodes,_ = pipeline(support,target,curve.b,m,variant,contract,True)
                                if profiled!=result:
                                    raise AssertionError('profiling changed the result')
                                reference = measure(alg.Ring(result['n']),result['generators'])
                                if reference['roots']!=result['roots'] or reference['completion_degree']!=result['completion_degree']:
                                    raise AssertionError('independent F4 matrix path mismatch')
                                independent = None
                                spec = contract['independent_audit']
                                if seed in spec['seeds'] and ti in spec['target_indices']:
                                    independent = audit(result['n'],result['generators'],result['roots'],spec['timeout_seconds_per_case'],spec['f5b_python_call_cap'])
                                result.pop('solutions')
                                case = dict(**identity,**result,status='VERIFIED',independent=independent,
                                    phase_opcodes=opcodes,total_opcodes=sum(opcodes.values()),
                                    uninstrumented_phase_seconds=[r[2] for r in replays],
                                    repetition_sha256=[digest(r[0]) for r in replays],
                                    oracle_solution_sha256=digest(truth),oracle_solution_count=len(truth))
                            except Exception as error:
                                case = dict(**identity,status='FAILED',error=f'{type(error).__name__}: {error}')
                            cases.append(case)
                            journal.write(json.dumps(case,separators=(',',':'))+'\n')
                            journal.flush()
            print(f'completed split {split}: {len(cases)} cases',flush=True)
    sources = [HERE/'algebra.py',HERE/'run.py',HERE/'independent.py',
               alg.LEGACY/'curves.py',alg.LEGACY/'matrix.py']
    result = dict(contract=contract,contract_sha256=digest(contract),
        python=platform.python_version(),python_implementation=platform.python_implementation(),
        source_sha256={str(p.relative_to(HERE.parent)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
        certificates=[cert for _,cert in neighbors],cases=cases,
        full_dlp_total_operations=None,S=None,rho_ratio=None,floor_ratio=None,full_dlp_speedup=None)
    (output/'raw.json.gz').write_bytes(gzip.compress(json.dumps(result,separators=(',',':')).encode(),mtime=0))
    summary = summarize(result)
    (output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    (output/'partial.jsonl').unlink()
    print(json.dumps({k:v for k,v in summary.items() if k not in ('rows','paired_controls')},indent=2))
    return result


def summarize(data):
    cases = data['cases']
    groups = defaultdict(list)
    for case in cases:
        groups[case['split'],case['model'],case['summands'],case['variant']].append(case)
    rows = []
    for (split,model,m,variant),items in sorted(groups.items()):
        verified = [c for c in items if c['status']=='VERIFIED']
        rows.append(dict(split=split,model=model,summands=m,variant=variant,cases=len(items),
            verified=len(verified),total_opcodes=sum(c['total_opcodes'] for c in verified) if len(verified)==len(items) else None,
            phase_opcodes={p:sum(c['phase_opcodes'][p] for c in verified) for p in PHASES},
            matrix_xors=sum(t['total_matrix_xors'] for c in verified for t in c['traces']),
            degree_min=min((c['completion_degree'] for c in verified),default=None),
            degree_max=max((c['completion_degree'] for c in verified),default=None),
            lift_checks=sum(c['lift_checks'] for c in verified)))
    baseline = {(r['split'],r['model'],r['summands']):r for r in rows if r['variant']=='ordered'}
    for row in rows:
        base = baseline[row['split'],row['model'],row['summands']]
        row['opcode_ratio_to_ordered'] = row['total_opcodes']/base['total_opcodes'] if row['total_opcodes'] is not None and base['total_opcodes'] else None
    indexed = {(c['seed'],c['model'],c['summands'],c['target_index'],c['variant']):c for c in cases}
    paired = defaultdict(lambda:dict(lower=0,equal=0,higher=0,opcode_lower=0,opcode_equal=0,opcode_higher=0))
    neighbors = defaultdict(lambda:dict(lower=0,equal=0,higher=0))
    for c in cases:
        if c['status']!='VERIFIED':
            continue
        base = indexed[c['seed'],c['model'],c['summands'],c['target_index'],'ordered']
        if base['status']=='VERIFIED':
            for name,key in [('completion_degree',''),('total_opcodes','opcode_')]:
                direction = 'lower' if c[name]<base[name] else 'higher' if c[name]>base[name] else 'equal'
                paired[c['variant']][key+direction] += 1
        if c['model']!='source':
            base = indexed[c['seed'],'source',c['summands'],c['target_index'],c['variant']]
            if base['status']=='VERIFIED':
                direction = 'lower' if c['completion_degree']<base['completion_degree'] else 'higher' if c['completion_degree']>base['completion_degree'] else 'equal'
                neighbors[c['variant']][direction] += 1
    audits = [c['independent'] for c in cases if c.get('independent') is not None]
    correctness = all(c['status']=='VERIFIED' for c in cases) and all(a['status']=='VERIFIED' for a in audits)
    fresh = [r for r in rows if r['variant']=='invariant' and r['split'].startswith('fresh')]
    return dict(systems=len(cases),verified=sum(c['status']=='VERIFIED' for c in cases),
        independent_audits=len(audits),independent_verified=sum(a['status']=='VERIFIED' for a in audits),
        f5b_verified=sum(a.get('f5b_status')=='VERIFIED' for a in audits),
        f5b_budget_exceeded=sum(a.get('f5b_status')=='BUDGET_EXCEEDED' for a in audits),
        rows=rows,paired_controls=dict(paired),neighbor_degree_pairs_by_variant=dict(neighbors),
        frozen_success=correctness and len(fresh)==20 and all(r['opcode_ratio_to_ordered']<=.8 for r in fresh),
        classification='accounting',metric='complete decomposition-call Python opcode proxy',full_dlp_speedup=None)


if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output',type=Path,default=HERE/'results/run-001')
    args = parser.parse_args()
    run(args.output)
