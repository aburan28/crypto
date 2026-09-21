"""Replay frozen subfield certificates and audit matched coverage/aggregates."""
import hashlib
import importlib.util
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
FOLDER = HERE / 'subfield_01'
sys.path.insert(0, str(FOLDER))
spec = importlib.util.spec_from_file_location('subfield_runner', FOLDER / 'run.py')
runner = importlib.util.module_from_spec(spec)
spec.loader.exec_module(runner)


def audit():
    contract = json.loads((FOLDER / 'contract.json').read_text())
    raw = [json.loads(s) for s in (FOLDER / 'raw.jsonl').read_text().splitlines()]
    summary = json.loads((FOLDER / 'summary.json').read_text())
    provenance = raw[0]
    for path, digest in provenance['sha256'].items():
        if hashlib.sha256((ROOT / path).read_bytes()).hexdigest() != digest:
            raise ArithmeticError('source/configuration hash mismatch: ' + path)
    expectedKeys = {'provenance', 'validation', 'target_generation', 'oracle',
                    'instance', 'trial', 'supplemental'}
    if set(r['kind'] for r in raw) != expectedKeys:
        raise ArithmeticError('missing or unexpected record kinds')
    if len([r for r in raw if r['kind'] == 'trial']) != 144:
        raise ArithmeticError('incomplete matched panel')
    if len([r for r in raw if r['kind'] == 'supplemental']) != 8:
        raise ArithmeticError('incomplete supplemental panel')
    oracles = {}
    instances = {}
    certificates = 0
    for row in raw:
        if row['kind'] != 'instance':
            continue
        body = {k: v for k, v in row.items() if k not in ('kind', 'sha256')}
        if hashlib.sha256(json.dumps(body, sort_keys=True).encode()).hexdigest() != row['sha256']:
            raise ArithmeticError('instance hash mismatch')
        key = row['n'], row['d']
        if key not in oracles:
            f, c = runner.makeCurve(row['n'], row['k'], row['coefficients'])
            for coefficient in (c.a2, c.a6):
                if coefficient in (0, f.one()) or f.frob(coefficient, row['k']) != coefficient:
                    raise ArithmeticError('non-F2 subfield coefficient invariant failed')
            oracles[key] = runner.PairOracle(f, c, row['d'])
        oracle = oracles[key]
        target = tuple(oracle.f.fromCoords(v) for v in row['target'])
        expected = oracle.expected(target)
        if expected != {tuple(x) for x in row['expected']}:
            raise ArithmeticError('recorded expected set differs from replayed group oracle')
        instances[row['sha256']] = row, expected, oracle
        trials = [r for r in raw if r['kind'] == 'trial' and r['instance_sha256'] == row['sha256']]
        if sorted((r['variant'], r['mode']) for r in trials) != sorted(
                (v, m) for v in contract['variants'] for m in ('first', 'enumerate')):
            raise ArithmeticError('matched coverage mismatch')
    # The same exact target is used across all three dimensions.
    for n in contract['field_degrees']:
        for cohort in ('development', 'holdout'):
            for stratum in ('uniform', 'known_decomposable'):
                rr = [r for r, _, _ in instances.values() if
                      (r['n'], r['cohort'], r['stratum']) == (n, cohort, stratum)]
                if sorted(r['d'] for r in rr) != contract['dimensions'] or len(
                        {tuple(r['target']) for r in rr}) != 1:
                    raise ArithmeticError('targets changed across dimensions')
    for row in raw:
        if row['kind'] not in ('trial', 'supplemental'):
            continue
        instance, expected, oracle = instances[row['instance_sha256']]
        for key in ('n', 'k', 'd', 'coefficients', 'target', 'cohort', 'stratum'):
            if row[key] != instance[key]:
                raise ArithmeticError('trial/instance mismatch: ' + key)
        runner.checkRow(row, oracle.f, oracle.c, expected)
        if len(row['certificates']) != row['verified_unique_relations']:
            raise ArithmeticError('missing signed-point certificates')
        for cert in row['certificates']:
            if sorted(p[0] for p in cert['points']) != cert['xs']:
                raise ArithmeticError('certificate projection mismatch')
        certificates += len(row['certificates'])
        if abs(sum(row['phase_seconds'].values()) - row['all_phase_seconds']) > 1e-7:
            raise ArithmeticError('phase timers are not exclusive/exhaustive')
        if any(row[k] is not None for k in ('common_operations', 'speedup', 'full_dlp_S', 'rho_ratio', 'floor_ratio')):
            raise ArithmeticError('unmeasured cost was populated')
    for group in summary['groups']:
        keys = ('kind', 'n', 'd', 'variant', 'mode', 'stratum')
        rr = [r for r in raw if all(r.get(k) == group[k] for k in keys)]
        checks = {'attempts': len(rr), 'resolved': sum(r['status'] != 'timeout' for r in rr),
                  'resolved_within_budget': sum(r['status'] != 'timeout' and r['within_budget'] for r in rr),
                  'verified_relations': sum(r['verified_unique_relations'] for r in rr),
                  'all_phase_seconds': sum(r['all_phase_seconds'] for r in rr),
                  'field_api_operations': sum(r['field_api_counts']['totals']['fieldOperations'] for r in rr)
                  if group['variant'] == 'hybrid-image' else None}
        if any(group[k] != v for k, v in checks.items()):
            raise ArithmeticError('summary aggregate mismatch')
    return {'matched_trials': 144, 'supplemental_trials': 8,
            'independent_oracle_instances': len(instances), 'signed_certificates_replayed': certificates,
            'source_and_instance_hashes_match': True, 'dimension_target_identity': True,
            'summary_aggregates_match': True, 'validation_failures': 0,
            'curve_addition_call_count': None, 'curve_doubling_call_count': None,
            'counter_note': 'Only field additions, multiplications, squarings and inversionCalls are instrumented. Inherited curve-call fields in the raw CountedField report are zero placeholders, not measurements; they are excluded from all aggregates.',
            'raw_sha256': hashlib.sha256((FOLDER / 'raw.jsonl').read_bytes()).hexdigest()}


if __name__ == '__main__':
    print(json.dumps(audit(), indent=2))
