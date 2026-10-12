"""Matched ephemeral, durable and resumed fixed-parameter correctness campaign."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import time

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'ecc2k130/codegen'))
import indexcalc_fixed as fixed
from test_indexcalc_fixed import fixture, transcript


def one(document, mode, directory):
    start = time.perf_counter_ns()
    reports = []
    if mode == 'resumed':
        with fixed.Campaign(document, directory) as campaign:
            reports.append(campaign.run(stage='pairs', pairBudget=7, seconds=2))
        with fixed.Campaign(document, directory) as campaign:
            reports.append(campaign.run(stage='collect', attempts=1, pairBudget=10000, seconds=2))
        with fixed.Campaign(document, directory) as campaign:
            reports.append(campaign.run(stage='collect', attempts=255, pairBudget=0, seconds=2))
        with fixed.Campaign(document, directory) as campaign:
            reports.append(campaign.run(stage='solve', attempts=256, pairBudget=0, seconds=2))
            rows = transcript(campaign)
    else:
        with fixed.Campaign(document, directory, memory=mode == 'ephemeral') as campaign:
            reports.append(campaign.run(attempts=256, pairBudget=10000, seconds=2))
            rows = transcript(campaign)
    wall = time.perf_counter_ns() - start
    final = reports[-1]
    return {'mode': mode, 'whole_process_scope_ns': wall, 'reports': reports,
            'status': final['status'], 'targets': final['targets'],
            'transcript_sha256': fixed.digest(rows), 'relations': rows}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    contract = fixed.readJson(Path(__file__).with_name('contract.json'))
    raw = []
    summary = {'classification': 'engineering capability; accounting diagnostics only',
               'matched_cases': 0, 'complete_runs': {}, 'whole_process_scope_ns': {},
               'verified_targets': {}, 'disagreements': 0, 'warm_replays': 0,
               'warm_new_pair_candidates': 0, 'warm_new_attempts': 0,
               'total_common_operations': None, 'S': None, 'rho_ratio': None, 'floor_ratio': None,
               'timing_scope': 'Python API including campaign construction, report archival and connection close; no timing speedup claim',
               'contract_sha256': fixed.digest(contract)}
    with (args.out / 'raw.jsonl').open('x') as stream:
        for m in contract['suite']['degrees']:
            for summands in contract['suite']['summands']:
                for seed in contract['suite']['seeds']:
                    document, expected = fixture(m, seed, summands)
                    for repetition in range(contract['suite']['repetitions']):
                        matched = []
                        for mode in ('ephemeral', 'durable', 'resumed'):
                            with tempfile.TemporaryDirectory() as directory:
                                row = one(document, mode, directory)
                                row.update(degree=m, summands=summands, seed=seed, repetition=repetition,
                                           split='holdout' if seed in contract['suite']['holdout_seeds'] else 'training',
                                           parameters=document, expected=expected)
                                if mode == 'durable':
                                    start = time.perf_counter_ns()
                                    with fixed.Campaign(document, directory) as campaign:
                                        warm = campaign.run(attempts=0, pairBudget=0, seconds=2)
                                    row['warm_wall_ns'] = time.perf_counter_ns() - start
                                    row['warm'] = warm
                                    fixed.require(warm['status'] == 'complete', 'warm replay failed')
                                    summary['warm_replays'] += 1
                                    summary['warm_new_pair_candidates'] += warm['reuse'].get('pair_candidates_built', 0)
                                    summary['warm_new_attempts'] += warm['reuse'].get('new_attempts', 0)
                            fixed.require(row['status'] == 'complete', 'incomplete matched small run')
                            for name, scalar in expected.items():
                                fixed.require(row['targets'][name]['verified'] and int(row['targets'][name]['scalar']) == scalar,
                                              'incorrect scalar')
                            summary['complete_runs'][mode] = summary['complete_runs'].get(mode, 0) + 1
                            summary['verified_targets'][mode] = summary['verified_targets'].get(mode, 0) + len(expected)
                            summary['whole_process_scope_ns'][mode] = summary['whole_process_scope_ns'].get(mode, 0) + row['whole_process_scope_ns']
                            stream.write(json.dumps(row, sort_keys=True) + '\n')
                            stream.flush()
                            matched.append(row)
                        fixed.require(len({r['transcript_sha256'] for r in matched}) == 1,
                                      'durable/resumed relation transcript diverged')
                        summary['matched_cases'] += 1
        profile = fixed.readJson(ROOT / 'docs/ic/params/ecc2k130-fixed.json')
        with tempfile.TemporaryDirectory() as directory:
            wide = []
            with fixed.Campaign(profile, directory) as campaign:
                wide.append(campaign.run(stage='pairs', pairBudget=11))
            with fixed.Campaign(profile, directory) as campaign:
                wide.append(campaign.run(stage='pairs', pairBudget=13))
            with fixed.Campaign(profile, directory, solver='sat') as campaign:
                wide.append(campaign.run(stage='collect', attempts=1))
            summary['degree131'] = {'orbit_columns': wide[0]['orbit_columns'],
                'signed_base_points': wide[0]['signed_base_points'],
                'coverage_ceiling': wide[0]['coverage_ceiling'],
                'pair_candidates_committed': 24,
                'relations_saved': wide[-1]['relations_saved'],
                'status': wide[-1]['status'], 'complete_dlp': False}
            (args.out / 'degree131.json').write_text(json.dumps(wide, indent=2) + '\n')
    fixed.require(summary['warm_new_pair_candidates'] == 0 and summary['warm_new_attempts'] == 0,
                  'warm replay recomputed completed work')
    sources = ['ecc2k130/codegen/indexcalc_fixed.py', 'ecc2k130/codegen/indexcalc_e2e.py',
               'ecc2k130/codegen/curves.py', 'ecc2k130/codegen/field.py',
               'ecc2k130/codegen/test_indexcalc_fixed.py', 'src/bin/ic/fixed.rs',
               'docs/ic/params/ecc2k130-fixed.json', 'research/fixed_parameters_20260915/run.py']
    summary['source_sha256'] = {name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest() for name in sources}
    (args.out / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
