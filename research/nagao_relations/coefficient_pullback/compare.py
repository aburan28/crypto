"""Audit frozen runs and emit diagnostics without treating censored work as a gain."""
import argparse
import hashlib
import json
from pathlib import Path

import experiment


def read_rows(path):
    return [json.loads(line) for line in path.read_text().splitlines()]


def totals(rows, key):
    return {name: sum(row[key][name] for row in rows) for name in rows[0][key]}


def audit(stage_dir, e2e_dir):
    raw = read_rows(stage_dir / 'raw.jsonl')
    assert not any(r['kind'] in ('error', 'validation_failure') for r in raw)
    trials = [r for r in raw if r['kind'] == 'trial']
    complete_pairs = []
    key = lambda r: (r['n'], r['d'], tuple(r['target']), r['mode'])
    baseline = {key(r): r for r in trials if r['variant'] == 'quadratic-image'}
    for candidate in trials:
        if candidate['variant'] != 'coefficient-pullback':
            continue
        before = baseline[key(candidate)]
        if before['status'] == 'timeout' or candidate['status'] == 'timeout':
            continue
        # Enumeration must produce the identical set. First mode compares one
        # independently verified relation for the identical target.
        if candidate['mode'] == 'enumerate':
            assert candidate['solutions'] == before['solutions']
        assert candidate['verified_unique_relations'] == before['verified_unique_relations']
        b = before['field_api_counts']['totals']
        c = candidate['field_api_counts']['totals']
        complete_pairs.append({'key': key(candidate), 'baseline': b, 'candidate': c,
                               'field_api_sum_ratio': b['fieldOperations']/c['fieldOperations']})
    runs = read_rows(e2e_dir / 'raw.jsonl')
    oracles = {}; checked = 0; distinct = set()
    for row in runs:
        assert row['status'] == 'verified', row
        n, d = row['panel']['n'], row['panel']['d']
        oracle = oracles.setdefault((n, d), experiment.previous.PairOracle(n, d))
        f, curve = oracle.f, oracle.curve
        g = tuple(f.fromCoords(x) for x in row['generator'])
        q = tuple(f.fromCoords(x) for x in row['target'])
        assert curve.mul(g, row['recovered_scalar']) == q
        for attempt in row['attempts']:
            result = attempt['result']
            target = tuple(f.fromCoords(x) for x in result['target'])
            expected = oracle.expected(target)
            got = {tuple(xs) for xs in result['solutions']}
            assert got <= expected
            if result['status'] == 'complete':
                assert got == expected
            checked += 1
            distinct.add((n, d, tuple(result['target'])))
    grouped = []
    for n in sorted({row['panel']['n'] for row in runs}):
        for variant in ['quadratic-image', 'coefficient-pullback', 'rho']:
            rr = [r for r in runs if r['panel']['n'] == n and r['variant'] == variant]
            grouped.append({'n': n, 'variant': variant, 'verified': len(rr),
                            'field_api_counts': totals(rr, 'all_phase_field_api_counts'),
                            'scalar_modular_counts': totals(rr, 'scalar_modular_counts'),
                            'phase_seconds': totals(rr, 'phase_seconds'),
                            'cold_wall_seconds': sum(r['cold_wall_seconds'] for r in rr),
                            'full_dlp_S': None, 'rho_ratio': None, 'cost_over_floor': None})
    for row in runs:
        if row['variant'] != 'coefficient-pullback':
            continue
        before = next(r for r in runs if r['panel'] == row['panel'] and r['seed'] == row['seed']
                      and r['variant'] == 'quadratic-image')
        assert before['generator'] == row['generator'] and before['target'] == row['target']
        path = lambda r: [(a['phase'], a['coefficient'], a['result']['target'], a['result']['solutions'])
                          for a in r['attempts']]
        assert path(before) == path(row), 'first-solution path changed; compare workloads explicitly'
    return {'stage_trial_count': len(trials), 'complete_stage_pairs': complete_pairs,
            'e2e_attempt_oracle_checks': checked, 'e2e_distinct_attempt_targets': len(distinct),
            'e2e_matched_attempt_paths': True, 'e2e_groups': grouped,
            'interpretation': 'Field API vectors and wall times are diagnostics; no calibrated common-operation speedup or normalized S is established.'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage', type=Path, required=True)
    parser.add_argument('--e2e', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = audit(args.stage, args.e2e)
    result['input_sha256'] = {str(p): hashlib.sha256(p.read_bytes()).hexdigest()
                              for p in [args.stage/'raw.jsonl', args.e2e/'raw.jsonl', Path(__file__)]}
    with args.output.open('x') as out:
        json.dump(result, out, indent=2)
        out.write('\n')
    print('Validated', result['stage_trial_count'], 'stage runs and',
          result['e2e_attempt_oracle_checks'], 'full-pipeline oracle calls.')
