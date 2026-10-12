#!/usr/bin/env python3
"""Fail-closed correctness checks and secondary, paired runtime diagnostics."""
import argparse
from collections import Counter, defaultdict
import hashlib
import json
import math
from pathlib import Path
import random
import statistics


def interval(ratios):
    if not ratios:
        return None
    rng = random.Random(20260914)
    logs = [math.log(r) for r in ratios]
    samples = sorted(math.exp(statistics.mean(rng.choices(logs, k=len(logs)))) for _ in range(4000))
    return dict(candidate_over_reference_geomean=math.exp(statistics.mean(logs)),
                paired_bootstrap_95=[samples[100], samples[3899]], pairs=len(ratios),
                limitation='descriptive shared-host timing; repeated seeds are not independent algorithm instances')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('directory', type=Path)
    args = ap.parse_args()
    records = [json.loads(line) for line in (args.directory / 'processes.jsonl').read_text().splitlines()]
    contract = json.loads((args.directory / 'contract.json').read_text())
    assert json.loads((args.directory / 'complete.json').read_text()) == dict(complete=True, processes=len(records))
    assert Counter(r['kind'] for r in records) == dict(native=84, stage=240, dlp=480, rho=60)
    for r in records:
        r['rows'] = [json.loads(line) for line in r['stdout'].splitlines() if line.startswith('{')]
        assert r['status'] in ('finished', 'timeout'), (r['command'], r['stderr'])
        if r['status'] == 'finished':
            assert r['returncode'] == 0
        for row in r['rows']:
            if row['phase'] == 'target':
                assert row['correct']
            if row['phase'] == 'complete_roots' and not row['stats']['exhausted']:
                assert row['correct'] and row['root_hash'] == row['truth_hash']
            for check in row.get('attempts', []):
                assert check['correct']
    native = defaultdict(dict)
    for r in records:
        if r['kind'] == 'native':
            assert r['status'] == 'finished'
            native[(r['case'], r['split'], r['rep'])][r['variant']] = r
    native_checks = 0
    for pair in native.values():
        a = [x for x in pair['reference']['rows'] if x['phase'] == 'target']
        b = [x for x in pair['candidate']['rows'] if x['phase'] == 'target']
        assert len(a) == len(b)
        for x, y in zip(a, b):
            for k in ['target', 'input_hash', 'output_hash', 'root_word_xors', 'solutions', 'reductions', 'exhausted', 'correct']:
                assert x[k] == y[k], (k, pair)
            native_checks += 1
    groups = defaultdict(list)
    for r in records:
        if r['kind'] != 'native':
            groups[(r['kind'], r['case'], r['variant'])].append(r)
    summaries = []
    census = {}
    keyed = {}
    for (kind, ci, variant), runs in sorted(groups.items()):
        row = dict(kind=kind, case=ci, configuration=contract['stage_cases'][ci], variant=variant,
                   processes=len(runs), statuses=dict(Counter(r['status'] for r in runs)),
                   verified_targets=0, completed_processes=0, cold_ns_completed=[], reductions=0,
                   total_calibrated_operations=None, S=None, rho_ratio=None, floor_ratio=None,
                   classification='engineering experiment; no calibrated advance established')
        stage_yield = []
        for r in runs:
            by_phase = {x['phase']: x for x in r['rows']}
            targets = [x for x in r['rows'] if x['phase'] == 'target']
            if 'census' in by_phase:
                key = f'{ci}:{r["seed"]}'
                c = {k: v for k, v in by_phase['census'].items() if k != 'validation_ns'}
                if key in census:
                    assert census[key]['census'] == c
                census[key] = dict(census=c, setup={k: v for k, v in by_phase['setup'].items() if not k.endswith('_ns')})
            if kind == 'stage':
                wins = sum(x['witness'] is not None for x in targets)
                row['verified_targets'] += wins
                row.setdefault('targets_observed', 0)
                row['targets_observed'] += len(targets)
                complete = r['status'] == 'finished' and len(targets) == 8 and not any(x['exhausted'] for x in targets)
                if complete:
                    ns = by_phase['setup']['setup_ns'] + sum(by_phase['target_setup'][k] for k in ['field_setup_ns', 'target_generation_ns'])
                    ns += sum(x['solve_ns'] + x['verification_ns'] for x in targets)
                    row['completed_processes'] += 1
                    row['cold_ns_completed'].append(ns)
                    stage_yield.append(wins)
                    keyed[(kind, ci, variant, r['seed'], r['rep'])] = dict(ns=ns, wins=wins,
                        targets=[(x['target'], x['truth']) for x in targets], input=by_phase['setup']['domain_hash'])
                row['reductions'] += sum(x.get('reductions', 0) for x in targets)
            else:
                terminal = by_phase.get(kind)
                if terminal and r['status'] == 'finished':
                    row['completed_processes'] += 1
                    if terminal['verified']:
                        row['verified_targets'] += 1
                        row['cold_ns_completed'].append(terminal['cold_ns'])
                        keyed[(kind, ci, variant, r['seed'], r['rep'])] = dict(ns=terminal['cold_ns'], wins=1,
                            input=by_phase.get('setup', {}).get('domain_hash'), target=by_phase.get('setup', {}).get('target_hash'),
                            known_log=r['known_log'])
                    row['reductions'] += terminal.get('reductions', 0)
                    row.setdefault('relations', 0)
                    row['relations'] += terminal.get('relations', 0)
                    row.setdefault('trials', 0)
                    row['trials'] += terminal.get('trials', 0)
        row['median_cold_ns_completed'] = statistics.median(row['cold_ns_completed']) if row['cold_ns_completed'] else None
        row['cold_ns_per_witness_complete_stage_processes'] = sum(row['cold_ns_completed']) / sum(stage_yield) if sum(stage_yield) else None
        summaries.append(row)
    comparisons = []
    for ci in range(len(contract['stage_cases'])):
        for kind in ['stage', 'dlp']:
            seeds = contract['stage']['seeds'] if kind == 'stage' else [s[0] for s in contract['dlp']['seeds_and_logs']]
            variants = ['charts_f4', 'charts_linear', 'enumerate'] if kind == 'stage' else ['charts_linear', 'enumerate', 'pair_table']
            for ref in (['ambient_f4', 'enumerate'] if kind == 'stage' else ['ambient_f4', 'enumerate', 'pair_table', 'rho']):
                for variant in variants:
                    if ref == variant:
                        continue
                    for split in ['all', 'fresh']:
                        pairs, cost_ratios = [], []
                        for seed in seeds:
                            if split == 'fresh' and seed not in ([937] if kind == 'stage' else [503, 607]):
                                continue
                            for rep in range(3):
                                a = keyed.get((kind, ci, ref, seed, rep))
                                if ref == 'rho':
                                    curve = contract['stage_cases'][ci]
                                    rho_ci = next(i for i, c in enumerate(contract['stage_cases'])
                                                  if all(c[k] == curve[k] for k in ['n', 'k', 'a', 'b']))
                                    a = keyed.get(('rho', rho_ci, 'rho', seed, rep))
                                b = keyed.get((kind, ci, variant, seed, rep))
                                if a and b:
                                    if ref != 'rho':
                                        assert a['input'] == b['input']
                                    if kind == 'stage':
                                        assert a['targets'] == b['targets'] and a['wins'] == b['wins']
                                    else:
                                        assert a['known_log'] == b['known_log']
                                        if ref != 'rho':
                                            assert a['target'] == b['target']
                                    pairs.append(b['ns'] / a['ns'])
                                    if a['wins'] and b['wins']:
                                        cost_ratios.append((b['ns'] / b['wins']) / (a['ns'] / a['wins']))
                        comparisons.append(dict(kind=kind, case=ci, reference=ref, candidate=variant, split=split,
                            timing=interval(pairs), cost_per_witness=interval(cost_ratios),
                            qualification='only matched completed workloads; missing pairs retained in variant summaries and block an unqualified claim'))
    result = dict(native_root_set_comparisons=native_checks,
        native_processes=84, complete_root_set_checks=sum(x['phase']=='complete_roots' and x['correct'] for r in records for x in r['rows']),
        full_dlp_oracle_checks=sum(len(x.get('attempts',[])) for r in records for x in r['rows']),
        statuses=dict(Counter(r['status'] for r in records)), summaries=summaries, census=census, comparisons=comparisons,
        analysis_source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        analysis_note='Comparator additionally reports the already frozen pair-table and rho references; numerical acceptance and measured sources unchanged.',
        headline='structural and stage diagnostics only; no calibrated total-cost or exponent advance',
        total_calibrated_operations=None, S=None, rho_ratio=None, floor_ratio=None)
    (args.directory / 'comparison.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('summaries','census','comparisons')}, indent=2))


if __name__ == '__main__':
    main()
