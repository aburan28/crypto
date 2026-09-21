"""Rerun frozen and fresh matched benchmarks; never overwrite evidence."""
import argparse
import functools
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
import types

HERE = Path(__file__).resolve().parent
PRIOR = HERE.parent/'coefficient_pullback'
sys.path.insert(0, str(PRIOR))
import experiment as previous
import pullback

candidateSpec = importlib.util.spec_from_file_location('pullback_incremental_solver', HERE/'optimized.py')
optimized = importlib.util.module_from_spec(candidateSpec)
candidateSpec.loader.exec_module(optimized)

spec = importlib.util.spec_from_file_location('incremental_dlp_pipeline', PRIOR/'e2e.py')
pipeline = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pipeline)
ROOT = previous.ROOT


def readRows(path):
    return [json.loads(line) for line in path.read_text().splitlines()]


def freeze(contract):
    frozen = []
    prior = readRows(PRIOR/'results/raw.jsonl')
    for row in prior:
        if row['kind'] == 'targets':
            frozen.extend({'n':row['n'], 'd':row['d'], 'stratum':row['stratum'],
                           'target':target, 'corpus':'frozen'} for target in row['targets'])
    historical = prior[:]
    for folder in sorted(HERE.parent.glob('solver_*')):
        if (folder/'raw.jsonl').exists():
            historical += readRows(folder/'raw.jsonl')
    excluded = {}
    for row in historical:
        if row.get('kind') == 'trial':
            excluded.setdefault(row['n'], set()).add(tuple(row['target']))
    fresh = []
    for panel in contract['panels']:
        n,d = panel['n'],panel['d']
        oracle = previous.previous.PairOracle(n,d)
        rng = random.Random(contract['fresh_seed']+n*100+d)
        selected = excluded.get(n,set()).copy()
        for stratum in contract['strata']:
            count = 0
            while count < contract['targets_per_stratum']:
                point = oracle.randomUniform(rng) if stratum == 'uniform' else oracle.randomSupported(rng)
                coords = tuple(oracle.f.toCoords(x) for x in point)
                if coords in selected:
                    continue
                selected.add(coords)
                fresh.append({**panel,'stratum':stratum,'target':coords,'corpus':'fresh'})
                count += 1
    return {'contract_sha256':hashlib.sha256((HERE/'contract.json').read_bytes()).hexdigest(),
            'targets':frozen+fresh}


def validate():
    oracle = previous.previous.PairOracle(5,3)
    f,curve = oracle.f,oracle.curve
    points = previous.nagaocompare.affinePoints(f,curve)
    occurrences = []
    for target in points:
        sets = []
        for cls in (pullback.Search,optimized.NormalizedSearch,optimized.CachedSearch):
            search = cls(f,curve,3,target,time.perf_counter()+60)
            coefficients = set()
            for a,b,invB,z in search.candidates():
                coefficients.add((a,b,z))
            sets.append(coefficients)
        if sets[0] != sets[1] or sets[1] != sets[2]:
            raise ArithmeticError(('branch coefficient set mismatch',target))
        occurrences.append(len(sets[0]))
        expected = oracle.expected(target)
        counts = []
        for variant in ('normalized-pullback','cached-pullback'):
            row = optimized.cell(5,3,[f.toCoords(x) for x in target],'enumerate',60,variant)
            if row['status'] != 'complete' or {tuple(x) for x in row['solutions']} != expected:
                raise ArithmeticError(('complete oracle mismatch',target,variant))
            counts.append(row['field_api_counts'])
        if counts[0] != counts[1]:
            raise ArithmeticError('cache changed field-operation accounting')
    return {'all_affine_five_bit_targets':len(points),'per_branch_coefficients_equal':True,
            'coefficient_occurrences':sum(occurrences),'complete_sets_equal':True,
            'normalized_cached_field_counts_equal':True}


def runDlp(panel,seed,variant,contract):
    # This module instance has its own namespace. Only its oracle binding is
    # substituted; the frozen pipeline function and all other modules stay intact.
    original = pipeline.pullback
    if variant in ('normalized-pullback','cached-pullback'):
        pipeline.pullback = types.SimpleNamespace(cell=functools.partial(optimized.cell,variant=variant))
    try:
        route = 'rho' if variant == 'rho' else 'coefficient-pullback'
        row = pipeline.run(panel,seed,route,contract)
        row['variant'] = variant
        return row
    finally:
        pipeline.pullback = original


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--freeze',type=Path)
    parser.add_argument('--targets',type=Path)
    parser.add_argument('--output',type=Path)
    parser.add_argument('--validate-only',action='store_true')
    args = parser.parse_args()
    contract = json.loads((HERE/'contract.json').read_text())
    if args.freeze:
        frozen = freeze(contract)
        with args.freeze.open('x') as out:
            json.dump(frozen,out,indent=2);out.write('\n')
        print('Frozen',len(frozen['targets']),'matched targets.')
        return
    if args.output is None or args.targets is None:
        parser.error('--targets and --output required for a run')
    args.output.mkdir(parents=True,exist_ok=False)
    targets = json.loads(args.targets.read_text())
    if targets['contract_sha256'] != hashlib.sha256((HERE/'contract.json').read_bytes()).hexdigest():
        raise ValueError('target/contract mismatch')
    sources = list(previous.CODE.glob('*.py'))+list(PRIOR.glob('*.py'))+list(HERE.glob('*.py'))+[HERE/'contract.json',args.targets]
    for folder in ['solver_02','solver_04','solver_05','solver_06','solver_07','solver_08']:
        sources += list((HERE.parent/folder).glob('*.py'))
    sources += [PRIOR/'results/raw.jsonl']+list(HERE.parent.glob('solver_*/raw.jsonl'))
    provenance = {'contract':contract,'targets':targets,'python':platform.python_version(),
                  'platform':platform.platform(),'pycryptosat':previous.importlib.metadata.version('pycryptosat'),
                  'cpu_model':next((s.split(':',1)[1].strip() for s in Path('/proc/cpuinfo').read_text().splitlines() if s.startswith('model name')),None),
                  'base_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
                  'command':[sys.executable]+sys.argv,
                  'source_sha256':{str(p.resolve().relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}}
    (args.output/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    rows = [];errors = []
    with (args.output/'raw.jsonl').open('x') as out:
        def record(row):
            rows.append(row);out.write(json.dumps(row)+'\n');out.flush()
        try:
            validation = validate();record({'kind':'validation',**validation});print(validation,flush=True)
        except Exception:
            record({'kind':'validation_error','traceback':traceback.format_exc()});raise
        if not args.validate_only:
            oracles = {};rng = random.Random(contract['fresh_seed'])
            for target in targets['targets']:
                n,d = target['n'],target['d']
                if (n,d) not in oracles:
                    oracles[n,d] = previous.previous.PairOracle(n,d)
                oracle = oracles[n,d];f,curve = oracle.f,oracle.curve
                point = tuple(f.fromCoords(x) for x in target['target'])
                expected = oracle.expected(point)
                for mode in contract['modes']:
                    variants = list(contract['variants']);rng.shuffle(variants)
                    for variant in variants:
                        try:
                            if variant == 'coefficient-pullback':
                                result = pullback.cell(n,d,target['target'],mode,contract['seconds_per_instance'])
                            elif variant in ('normalized-pullback','cached-pullback'):
                                result = optimized.cell(n,d,target['target'],mode,contract['seconds_per_instance'],variant)
                            else:
                                result = previous.previous.s4.runCell(f,curve,d,point,expected,variant,mode,contract['seconds_per_instance'])
                            got = {tuple(x) for x in result['solutions']}
                            if not got <= expected or (result['status']=='complete' and got!=expected):
                                raise ArithmeticError('benchmark oracle mismatch')
                            record({'kind':'stage',**target,'expected_count':len(expected),**result})
                            print(target['corpus'],n,d,mode,variant,result['status'],len(got),round(result['all_phase_seconds'],4),flush=True)
                        except Exception:
                            error = {'kind':'error',**target,'mode':mode,'variant':variant,'traceback':traceback.format_exc()}
                            errors.append(error);record(error)
            for panel in contract['full_dlp']['panels']:
                for seed in contract['full_dlp']['seeds']:
                    variants = list(contract['full_dlp']['variants']);rng.shuffle(variants)
                    for variant in variants:
                        try:
                            result = runDlp(panel,seed,variant,contract['full_dlp'])
                            record({'kind':'e2e','corpus':'frozen' if seed in contract['full_dlp']['frozen_seeds'] else 'fresh',**result})
                            print('e2e',panel['n'],seed,variant,result['status'],flush=True)
                        except Exception:
                            error = {'kind':'error','panel':panel,'seed':seed,'variant':variant,'traceback':traceback.format_exc()}
                            errors.append(error);record(error)
    (args.output/'summary.json').write_text(json.dumps({'validation':validation,'errors':errors,
        'stage_runs':sum(r['kind']=='stage' for r in rows),'e2e_runs':sum(r['kind']=='e2e' for r in rows),
        'verified_scalars':sum(r['kind']=='e2e' and r['status']=='verified' for r in rows)},indent=2)+'\n')
    if errors or any(r['kind']=='e2e' and r['status']!='verified' for r in rows):
        raise SystemExit(1)


if __name__ == '__main__':
    main()
