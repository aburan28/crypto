"""Compare a complete candidate run against immutable baseline evidence."""
import argparse
import hashlib
import json
import math
from pathlib import Path

from run import aggregate


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_run(path):
    summary = json.loads((path / 'summary.json').read_text())
    metadata = summary['metadata']
    for name, key in [('raw.jsonl', 'raw_sha256')]:
        if sha256(path / name) != summary[key]:
            raise ValueError(f'{path}: {name} integrity mismatch')
    for name, key in [('contract.json', 'contract_sha256'), ('corpus_manifest.json', 'corpus_sha256')]:
        if sha256(path / name) != metadata[key]:
            raise ValueError(f'{path}: {name} integrity mismatch')
    contract = json.loads((path / 'contract.json').read_text())
    records = [json.loads(line) for line in (path / 'raw.jsonl').read_text().splitlines()]
    expected = {(c, v['id'], r) for c in contract['cases'] for v in contract['variants']
                for r in range(contract['repetitions'])}
    slots = [(r['case'], r['variant'], r['repetition']) for r in records]
    if len(slots) != len(expected) or set(slots) != expected:
        raise ValueError(f'{path}: missing, duplicate, or unexpected trial slots')
    for row in records:
        if not row['verified'] or row['status'] not in ('SAT', 'UNSAT', 'SAT_ALGEBRAIC_ONLY'):
            raise ValueError(f'{path}: unverified or incomplete result {row["case"]}')
        if row['status'] == 'SAT' and not row['point_relation_verified']:
            raise ValueError(f'{path}: SAT lacks a point relation')
        if row['status'] == 'SAT_ALGEBRAIC_ONLY' and (
                row['case'] not in contract['algebraic_only_controls'] or row['point_relation_verified']):
            raise ValueError(f'{path}: algebraic-only control mislabeled as a relation')
        if type(row['conflicts']) is not int or row['conflicts'] < 0:
            raise ValueError(f'{path}: invalid conflict counter')
        if row['returncode'] != 0:
            raise ValueError(f'{path}: nonzero solver exit')
        try:
            printed_conflicts = int(row['stdout'].strip().splitlines()[-1])
        except (ValueError, IndexError) as exc:
            raise ValueError(f'{path}: missing raw conflict counter') from exc
        if printed_conflicts != row['conflicts']:
            raise ValueError(f'{path}: conflict counter does not match stdout')
        if not math.isfinite(row['process_wall_s']) or row['process_wall_s'] < 0:
            raise ValueError(f'{path}: invalid process time')
    cells, groups = aggregate(records, contract)
    if summary['status'] != 'complete' or not all(c['repeatable'] for c in cells):
        raise ValueError(f'{path}: incomplete or non-repeatable run')
    if cells != summary['cells'] or groups != summary['groups']:
        raise ValueError(f'{path}: saved aggregates do not match raw records')
    # All configurations must decide the same instance, regardless of filename label.
    for case in contract['cases']:
        if len({r['status'] for r in records if r['case'] == case}) != 1:
            raise ValueError(f'{path}: configurations disagree on {case}')
    return summary, records


def ratio(candidate, baseline):
    return candidate / baseline if baseline else (1.0 if candidate == 0 else None)


def compare(baseline_path, candidate_path):
    baseline, br = load_run(baseline_path)
    candidate, cr = load_run(candidate_path)
    for key in ('suite_id', 'contract_sha256', 'corpus_sha256', 'runner_sha256',
                'checker_sha256', 'counter_unit'):
        if baseline['metadata'][key] != candidate['metadata'][key]:
            raise ValueError(f'incompatible {key}; version the suite instead of moving its target')
    old = {(r['case'], r['variant'], r['repetition']): r for r in br}
    for row in cr:
        reference = old[(row['case'], row['variant'], row['repetition'])]
        if row['normalized_sha256'] != reference['normalized_sha256']:
            raise ValueError(f'input fingerprint changed: {row["case"]}')
        if row['status'] != reference['status']:
            raise ValueError(f'outcome changed: {row["case"]}')
    b_cells = {(c['case'], c['variant']): c for c in baseline['cells']}
    pairs = []
    for c in candidate['cells']:
        b = b_cells[(c['case'], c['variant'])]
        pairs.append({'case': c['case'], 'variant': c['variant'], 'status': c['status'],
                      'baseline_median_conflicts': b['median_conflicts'],
                      'candidate_median_conflicts': c['median_conflicts'],
                      'conflict_ratio_candidate_over_baseline': ratio(c['median_conflicts'], b['median_conflicts']),
                      'counter_regression_over_10_percent': c['median_conflicts'] > 1.1*b['median_conflicts'],
                      'secondary_process_time_ratio_candidate_over_baseline': ratio(c['median_process_wall_s'], b['median_process_wall_s'])})
    variants = []
    for v in sorted({p['variant'] for p in pairs}):
        ps = [p for p in pairs if p['variant'] == v]
        bc = sum(p['baseline_median_conflicts'] for p in ps)
        cc = sum(p['candidate_median_conflicts'] for p in ps)
        regressions = [p['case'] for p in ps if p['counter_regression_over_10_percent']]
        variants.append({'variant': v, 'cases': len(ps), 'baseline_sum_median_conflicts': bc,
                         'candidate_sum_median_conflicts': cc, 'conflict_ratio_candidate_over_baseline': ratio(cc, bc),
                         'cases_regressing_over_10_percent': regressions,
                         'diagnostic_20_percent_counter_target_met': cc <= .8*bc and not regressions,
                         'full_dlp_S': None, 'cost_over_rho': None, 'cost_over_floor': None})
    matched_host = all(baseline['metadata'][k] == candidate['metadata'][k]
                       for k in ('cpu', 'platform', 'cpu_affinity'))
    return {'status': 'correctness_and_compatibility_pass', 'classification': 'accounting',
            'baseline_solver_sha256': baseline['metadata']['solver_sha256'],
            'candidate_solver_sha256': candidate['metadata']['solver_sha256'],
            'baseline_raw_sha256': baseline['raw_sha256'], 'candidate_raw_sha256': candidate['raw_sha256'],
            'counter_unit': baseline['metadata']['counter_unit'],
            'timing_host_metadata_matches': matched_host,
            'timing_warning': 'Secondary observations only; host metadata cannot prove equal load. Rebenchmark both binaries on the same host for a timing claim.',
            'variants': variants, 'per_case': pairs,
            'full_dlp_claim': False,
            'note': 'Counter targets apply only when solver counter meaning is unchanged. A conflict reduction is not a calibrated total-operation improvement.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('baseline', type=Path)
    parser.add_argument('candidate', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--require-counter-target', action='store_true')
    args = parser.parse_args()
    if args.output.exists():
        parser.error('output must not already exist')
    try:
        report = compare(args.baseline, args.candidate)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with args.output.open('x') as stream:
            json.dump(report, stream, indent=2, allow_nan=False); stream.write('\n')
        print(json.dumps({'status': report['status'], 'variants': report['variants']}, indent=2))
        if args.require_counter_target and not all(v['diagnostic_20_percent_counter_target_met'] for v in report['variants']):
            raise SystemExit(1)
    except (ValueError, KeyError, OSError) as exc:
        parser.error(str(exc))


if __name__ == '__main__':
    main()
