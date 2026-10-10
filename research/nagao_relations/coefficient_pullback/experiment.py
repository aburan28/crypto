"""Frozen paired exploration: coefficient pullback, image baseline, S4 and S3."""
import argparse
import hashlib
import importlib.metadata
import importlib.util
import json
from pathlib import Path
import platform
import random
import subprocess
import sys
import time
import traceback

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
CODE = ROOT / 'ecc2k130/codegen'
for folder in [CODE, HERE.parent / 'solver_08', HERE.parent / 'solver_07']:
    sys.path.insert(0, str(folder))
import image_solver
import pullback
import nagaocompare

spec = importlib.util.spec_from_file_location('campaign07', HERE.parent / 'solver_07/run.py')
previous = importlib.util.module_from_spec(spec); spec.loader.exec_module(previous)


def validate():
    oracle = previous.PairOracle(5, 3); f = oracle.f; curve = oracle.curve
    tripleTruth, _ = previous.old.baseline.oracle(f, curve, 3)
    targets = nagaocompare.affinePoints(f, curve)
    priorCount = newCount = 0
    for target in targets:
        expected = oracle.expected(target)
        if expected != tripleTruth.get(target, set()):
            raise ArithmeticError('pair/triple ground truth mismatch')
        validSets = []
        for cls in [image_solver.Search, pullback.Search]:
            search = cls(f, curve, 3, target, time.perf_counter()+60)
            valid = set()
            for a, b, invB, z in search.candidates():
                if cls is image_solver.Search: priorCount += 1
                else: newCount += 1
                got = search.recover(a, b, invB, z)
                if got is not None:
                    valid.add((a, b))
            validSets.append(valid)
        if validSets[0] != validSets[1]:
            raise ArithmeticError(('valid coefficient sets changed', target, validSets))
        result = pullback.cell(5, 3, [f.toCoords(x) for x in target], 'enumerate', 60)
        if result['status'] != 'complete' or {tuple(x) for x in result['solutions']} != expected:
            raise ArithmeticError(('pullback/oracle mismatch', target))
        for variant in ['s4-symmetric', 'chained-s3']:
            control = previous.s4.runCell(f, curve, 3, target, expected, variant, 'enumerate', 30)
            if control['status'] != 'complete' or {tuple(x) for x in control['solutions']} != expected:
                raise ArithmeticError(('Semaev control/oracle mismatch', target, variant))
    return {'all_affine_five_bit_targets': len(targets), 'original_candidate_occurrences': priorCount,
            'pullback_candidate_occurrences': newCount, 'valid_coefficient_sets_equal': True,
            'complete_projected_sets_equal': True, 'pair_oracle_matches_triple_oracle': True,
            'semaev_control_target_checks': 2*len(targets)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--validate-only', action='store_true')
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    contract = json.loads((HERE/'contract.json').read_text())
    sources = list(CODE.glob('*.py')) + list(HERE.glob('*.py')) + [HERE/'contract.json']
    for folder in ['solver_02', 'solver_04', 'solver_05', 'solver_06', 'solver_07', 'solver_08']:
        sources += list((HERE.parent/folder).glob('*.py'))
    historical = []
    for folder in sorted(HERE.parent.glob('solver_*')):
        path = folder/'raw.jsonl'
        if path.exists():
            sources.append(path)
            historical.extend(json.loads(line) for line in path.read_text().splitlines())
    excluded = {}
    for row in historical:
        if row.get('kind') == 'trial':
            excluded.setdefault(row['n'], set()).add(tuple(row['target']))
    provenance = {'commit': subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
                  'python': platform.python_version(), 'pycryptosat': importlib.metadata.version('pycryptosat'),
                  'command': [sys.executable]+sys.argv,
                  'sha256': {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}}
    (args.output/'contract.json').write_text(json.dumps(contract,indent=2)+'\n')
    rows = []
    with (args.output/'raw.jsonl').open('x') as out:
        def record(row):
            rows.append(row); out.write(json.dumps(row)+'\n'); out.flush()
        record({'kind':'provenance', **provenance})
        try:
            checks = validate(); record({'kind':'validation', **checks}); print(checks,flush=True)
        except Exception:
            record({'kind':'validation_failure','traceback':traceback.format_exc()}); raise
        if not args.validate_only:
            for panel in contract['panels']:
                n, d = panel['n'], panel['d']; t = time.perf_counter()
                oracle = previous.PairOracle(n,d); f = oracle.f; curve = oracle.curve
                record({'kind':'oracle','n':n,'d':d,'signed_base_size':2*len(oracle.base),
                        'pair_entries':sum(len(v) for v in oracle.pairs.values()),
                        'validation_setup_seconds':time.perf_counter()-t})
                rng = random.Random(contract['seed']+n*100+d)
                selected = set(excluded.get(n,set()))
                for stratum in contract['strata']:
                    targets = []
                    while len(targets) < contract['targets_per_stratum']:
                        target = oracle.randomUniform(rng) if stratum == 'uniform' else oracle.randomSupported(rng)
                        coords = tuple(f.toCoords(v) for v in target)
                        if coords not in selected:
                            selected.add(coords); targets.append(coords)
                    record({'kind':'targets','n':n,'d':d,'stratum':stratum,'targets':targets})
                    for coords in targets:
                        target = tuple(f.fromCoords(v) for v in coords); expected = oracle.expected(target)
                        for mode in contract['modes']:
                            variants = list(contract['variants']); rng.shuffle(variants)
                            for variant in variants:
                                try:
                                    if variant == 'quadratic-image':
                                        result = image_solver.cell(n,d,coords,mode,contract['seconds_per_instance'])
                                    elif variant == 'coefficient-pullback':
                                        result = pullback.cell(n,d,coords,mode,contract['seconds_per_instance'])
                                    else:
                                        result = previous.s4.runCell(f,curve,d,target,expected,variant,mode,contract['seconds_per_instance'])
                                    got = {tuple(x) for x in result['solutions']}
                                    if not got <= expected or (result['status'] == 'complete' and got != expected):
                                        raise ArithmeticError('candidate/oracle mismatch')
                                    result['within_budget'] = result['all_phase_seconds'] <= contract['seconds_per_instance']
                                    record({'kind':'trial','stratum':stratum,'expected_count':len(expected),**result})
                                    print(n,d,stratum,mode,variant,result['status'],len(got),round(result['all_phase_seconds'],4),flush=True)
                                except Exception:
                                    record({'kind':'error','n':n,'d':d,'stratum':stratum,'mode':mode,'variant':variant,
                                            'target':coords,'traceback':traceback.format_exc()})
        groups = []
        for panel in contract['panels']:
            for variant in contract['variants']:
                for mode in contract['modes']:
                    for stratum in contract['strata']:
                        rr = [r for r in rows if r['kind']=='trial' and (r['n'],r['d'],r['variant'],r['mode'],r['stratum']) == (panel['n'],panel['d'],variant,mode,stratum)]
                        if not rr: continue
                        ops = {name:sum(r['field_api_counts']['totals'][name] for r in rr) for name in ['additions','multiplications','squarings','fieldOperations']} if variant in ['coefficient-pullback','quadratic-image'] else None
                        groups.append({**panel,'variant':variant,'mode':mode,'stratum':stratum,'attempted':len(rr),
                                       'resolved':sum(r['status']!='timeout' for r in rr),
                                       'resolved_within_budget':sum(r['status']!='timeout' and r['within_budget'] for r in rr),
                                       'verified_relations':sum(r['verified_unique_relations'] for r in rr),
                                       'all_phase_seconds':sum(r['all_phase_seconds'] for r in rr),
                                       'field_api_counts':ops,'full_dlp_S':None,'rho_ratio':None,'cost_over_floor':None})
        result = {'provenance':provenance,'validation':checks,'groups':groups,
                  'errors':[r for r in rows if r['kind']=='error'],'contract':contract}
        (args.output/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
        if result['errors']: raise SystemExit(1)


if __name__ == '__main__':
    main()
