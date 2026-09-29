"""Run the frozen, exhaustive, tiny-field neighbor comparison."""
import argparse
import gzip
import hashlib
import json
import random
import time
from collections import defaultdict
from itertools import product
from pathlib import Path
from curves import O, certified_neighbors, decomposition_value
from matrix import Ring, measure

HERE = Path(__file__).resolve().parent


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def encode(curve, support, target, summands, canonical):
    ring = Ring(2*summands)
    values, indices = [], []
    for assignment in range(1 << ring.n):
        slots = tuple((assignment >> (2*j)) & 3 for j in range(summands))
        indices.append(slots)
        values.append(decomposition_value([support[k][0][0] for k in slots], target, curve.b))
    generators = [ring.anf([(v >> bit) & 1 for v in values]) for bit in range(8)]
    if canonical:
        generators.append(ring.anf([int(tuple(sorted(s)) != s) for s in indices]))
    return ring, generators, indices


def ground_truth(curve, support, target, indices):
    """Independent point-addition oracle, including every y lift and repetition."""
    answer = {}
    checks = 0
    for a, slots in enumerate(indices):
        valid = []
        for signs in product(range(2), repeat=len(slots)):
            points = [support[k][s] for k, s in zip(slots, signs)]
            checks += 1
            if curve.total(points) == target:
                valid.append(signs)
        if valid:
            answer[a] = valid
    return answer, checks


def run():
    contract = json.loads((HERE/'contract.json').read_text())
    started = time.monotonic()
    source, subgroup, neighbors = certified_neighbors()
    choices = sorted({min(p, source.neg(p)) for p in subgroup if p != O and p[0] != 0})
    cases = []
    for split, seed in contract['seeds'].items():
        representatives = random.Random(seed).sample(choices, 4)
        source_support = [(p, source.neg(p)) for p in representatives]
        models = [('source', source, lambda p: p)] + [
            (f'neighbor_{phi.target.b}', phi.target, phi) for phi, _ in neighbors]
        for summands in contract['summands']:
            for target_index, original_target in enumerate(subgroup):
                paired_truth = {}
                for name, curve, mapping in models:
                    support = [tuple(mapping(p) for p in pair) for pair in source_support]
                    target = mapping(original_target)
                    for encoding in contract['encodings']:
                        ring, generators, indices = encode(curve, support, target, summands, encoding == 'canonical')
                        truth, lift_checks = ground_truth(curve, support, target, indices)
                        if encoding == 'canonical':
                            truth = {a: lifts for a, lifts in truth.items()
                                     if indices[a] == tuple(sorted(indices[a]))}
                        signature = digest(truth)
                        if name == 'source':
                            paired_truth[encoding] = signature
                        elif signature != paired_truth[encoding]:
                            raise AssertionError('transport changed the actual decomposition workload')
                        repetitions = [measure(ring, generators) for _ in range(contract['repetitions'])]
                        if any(r != repetitions[0] for r in repetitions[1:]):
                            raise AssertionError('nondeterministic degree/counter result')
                        measured = repetitions[0]
                        if set(truth) - set(measured['roots']):
                            raise AssertionError('Semaev encoding lost a genuine decomposition')
                        lifted = {a: truth[a] for a in measured['roots'] if a in truth}
                        if lifted != truth:
                            raise AssertionError('solver extraction/lifting mismatch')
                        cases.append(dict(split=split, seed=seed, summands=summands,
                            target_index=target_index, model=name, b=curve.b, encoding=encoding,
                            target=target, support=support, input_sha256=digest(generators),
                            generators=generators, truth_sha256=signature,
                            verified_x_assignments=len(truth), verified_signed_tuples=sum(map(len, truth.values())),
                            rejected_algebraic_roots=len(set(measured['roots'])-set(truth)),
                            equation_value_evaluations=1 << ring.n,
                            mobius_xors=(8+(encoding == 'canonical'))*ring.n*(1 << (ring.n-1)),
                            exhaustive_lift_checks=lift_checks,
                            repetition_sha256=[digest(r) for r in repetitions], **measured))
    return dict(contract_sha256=digest(contract), contract=contract,
        source_sha256={p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                       for p in [HERE/'curves.py', HERE/'matrix.py', HERE/'experiment.py']},
        certificates=[cert for _, cert in neighbors], cases=cases,
        elapsed_seconds=time.monotonic()-started,
        accounting=dict(field_table_entries=65536, field_table_setup_included_in_matrix_xors=False,
            matrix_xors_are_stage_only=True, native_rust_execution=False,
            full_dlp_total_operations=None, S=None, rho_ratio=None, floor_ratio=None, speedup=None))


def summarize(result):
    grouped = defaultdict(list)
    for case in result['cases']:
        grouped[(case['split'], case['encoding'], case['model'])].append(case)
    rows = []
    for (split, encoding, model), cases in sorted(grouped.items()):
        for solver in ['f4', 'f5']:
            traces = [t for c in cases for t in c['traces'][solver]]
            rows.append(dict(split=split, encoding=encoding, model=model, solver=solver,
                cases=len(cases), degree_min=min(c['completion_degree'] for c in cases),
                degree_max=max(c['completion_degree'] for c in cases),
                degree_sum=sum(c['completion_degree'] for c in cases),
                matrix_xors=sum(t['total_matrix_xors'] for t in traces),
                criterion_xors=sum(t['criterion_xors'] for t in traces),
                matrix_rows=sum(t['matrix_rows'] for t in traces),
                pruned_rows=sum(t['pruned_rows'] for t in traces),
                all_verified=all(c['status']=='VERIFIED' for c in cases)))
    source = {(r['split'],r['encoding'],r['solver']):r for r in rows if r['model']=='source'}
    for row in rows:
        base = source[row['split'], row['encoding'], row['solver']]
        row['matrix_xor_ratio_to_source'] = row['matrix_xors']/base['matrix_xors']
    indexed = {(c['split'],c['encoding'],c['summands'],c['target_index'],c['model']): c for c in result['cases']}
    degree_pairs = defaultdict(lambda: dict(lower=0, equal=0, higher=0))
    for c in result['cases']:
        if c['model'] != 'source':
            base = indexed[c['split'],c['encoding'],c['summands'],c['target_index'],'source']
            key = 'lower' if c['completion_degree'] < base['completion_degree'] else 'higher' if c['completion_degree'] > base['completion_degree'] else 'equal'
            degree_pairs[c['model']][key] += 1
    success = all(r['matrix_xor_ratio_to_source'] <= .9 for r in rows if r['model']!='source' and r['solver']=='f5') and all(v['higher']==0 for v in degree_pairs.values())
    return dict(rows=rows, degree_pairs=dict(degree_pairs), success=success,
        systems=len(result['cases']), solver_runs=len(result['cases'])*2*result['contract']['repetitions'],
        rejected_algebraic_roots=sum(c['rejected_algebraic_roots'] for c in result['cases']),
        classification='accounting', full_dlp_speedup=None)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, default=HERE/'results')
    args = parser.parse_args()
    result = run()
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output/'raw.json.gz').write_bytes(gzip.compress((json.dumps(result, separators=(',', ':'))+'\n').encode(), mtime=0))
    summary = summarize(result)
    (args.output/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps({k:v for k,v in summary.items() if k!='rows'}, indent=2))
