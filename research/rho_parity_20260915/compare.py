#!/usr/bin/env python3
"""Check every saved oracle input and compare paired complete workloads."""
import argparse
from collections import Counter, defaultdict
import gzip
import itertools
import json
import math
from pathlib import Path
import random
import statistics


def interval(by_seed):
    if not by_seed:
        return None
    # Replicates share an algorithm instance. Resample seeds, not repetitions.
    values = [statistics.mean(v) for v in by_seed.values()]
    rng = random.Random(20260915)
    samples = sorted(math.exp(statistics.mean(rng.choices(values, k=len(values)))) for _ in range(4000))
    return dict(ratio=math.exp(statistics.mean(values)), ci95=[samples[100], samples[3899]],
                seeds=len(values), pairs=sum(map(len, by_seed.values())),
                unit='cold runtime ratio; paired seed-cluster bootstrap; shared-host diagnostic')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('directory', type=Path)
    ap.add_argument('--allow-witness-order-change', action='store_true')
    args = ap.parse_args()
    records = [json.loads(x) for x in gzip.open(args.directory/'processes.jsonl.gz', 'rt')]
    contract = json.loads((args.directory/'contract.json').read_text())
    assert json.loads((args.directory/'complete.json').read_text()) == dict(processes=len(records), complete=True)
    assert Counter(r['kind'] for r in records) == dict(native=84, stage=240, dlp=1440, rho=240)
    pairs, runs, ledger = defaultdict(dict), {}, defaultdict(list)
    checks = dict(native_targets=0, complete_root_sets=0, oracle_attempts=0, matched_dlp_inputs=0)
    for r in records:
        assert r['status'] in ('finished', 'timeout'), (r['command'], r['stderr'])
        rows = [json.loads(x) for x in r['stdout'].splitlines() if x.startswith('{')]
        r['rows'] = rows
        for d in rows:
            if d['phase'] == 'target':
                assert d['correct']
            if d['phase'] == 'complete_roots':
                assert not d['stats']['exhausted'] and d['correct'] and d['root_hash'] == d['truth_hash']
                checks['complete_root_sets'] += 1
            for a in d.get('attempts', []):
                assert a['correct']
                checks['oracle_attempts'] += 1
        key = (r['kind'], r['case'], r['variant'], r.get('seed', r.get('split')), r['rep'])
        assert r['revision'] not in pairs[key]
        pairs[key][r['revision']] = r
        if r['kind'] in ('dlp', 'rho'):
            terminal = next((d for d in rows if d['phase'] == r['kind']), None)
            r['terminal'] = terminal
            runs[key+(r['revision'],)] = r
            ledger[(r['case'], r['variant'], r['revision'])].append(r)
    for key, p in pairs.items():
        assert set(p) == {'reference', 'candidate'}
        a, b = p['reference'], p['candidate']
        if key[0] == 'native':
            assert a['status'] == b['status'] == 'finished'
            xs = [x for x in a['rows'] if x['phase'] == 'target']
            ys = [x for x in b['rows'] if x['phase'] == 'target']
            assert len(xs) == len(ys)
            for x, y in zip(xs, ys):
                for k in ['target', 'input_hash', 'output_hash', 'solutions', 'reductions', 'exhausted']:
                    assert x[k] == y[k], (key, k)
                checks['native_targets'] += 1
        if key[0] in ('stage', 'dlp') and a['status'] == b['status'] == 'finished':
            sa = next(d for d in a['rows'] if d['phase'] == 'setup')
            sb = next(d for d in b['rows'] if d['phase'] == 'setup')
            for k in ['basis', 'domain_hash', 'points', 'signed_orbits', 'target_hash', 'subgroup', 'cofactor']:
                assert sa[k] == sb[k], (key, k)
        if key[0] == 'dlp' and a['status'] == b['status'] == 'finished':
            x, y = a['terminal'], b['terminal']
            assert x['verified'] == y['verified'], key
            for aa, bb in zip(x.get('attempts', []), y.get('attempts', [])):
                for k in ['trial', 'target', 'a', 'b']:
                    assert aa[k] == bb[k], (key, k)
            if not args.allow_witness_order_change:
                assert x.get('attempts') == y.get('attempts'), key
                for k in ['trials', 'relations', 'matrix_rank', 'matrix_columns']:
                    assert x.get(k) == y.get(k), (key, k)
            checks['matched_dlp_inputs'] += 1
    summary = []
    for (case, variant, revision), rs in sorted(ledger.items()):
        verified = [r['terminal'] for r in rs if r['status'] == 'finished' and r['terminal'] and r['terminal']['verified']]
        summary.append(dict(case=case, variant=variant, revision=revision, attempts=len(rs),
                            verified=len(verified), statuses=dict(Counter(r['status'] for r in rs)),
                            cold_ns=statistics.median(d['cold_ns'] for d in verified) if verified else None,
                            total_operations=None, S=None, rho_cost_ratio=None, floor_ratio=None,
                            classification='engineering runtime diagnostic'))
    comparisons = []
    fresh = {p[0] for p in contract['fresh']}
    for ci, variant, split, reference in itertools.product(range(10), contract['dlp_variants'], ['all', 'fresh'], ['reference', 'rho', 'pair_table']):
        if reference == 'pair_table' and variant == 'pair_table':
            continue
        by_seed = defaultdict(list)
        missing = 0
        for seed, log in contract['frozen']+contract['fresh']:
            if split == 'fresh' and seed not in fresh:
                continue
            for rep in range(3):
                b = runs.get(('dlp', ci, variant, seed, rep, 'candidate'))
                if reference == 'rho':
                    c = contract['stage_cases'][ci]
                    rho_ci = next(j for j, cc in enumerate(contract['stage_cases']) if all(c[k] == cc[k] for k in ['n','k','a','b']))
                    a = runs.get(('rho', rho_ci, 'rho', seed, rep, 'candidate'))
                elif reference == 'pair_table':
                    a = runs.get(('dlp', ci, 'pair_table', seed, rep, 'candidate'))
                else:
                    a = runs.get(('dlp', ci, variant, seed, rep, 'reference'))
                if all(r and r['status'] == 'finished' and r['terminal'] and r['terminal']['verified'] for r in [a,b]):
                    by_seed[seed].append(math.log(b['terminal']['cold_ns']/a['terminal']['cold_ns']))
                else:
                    missing += 1
        comparisons.append(dict(case=ci, variant=variant, split=split, reference=reference,
                                missing=missing, runtime=interval(by_seed)))
    out = dict(processes=len(records), checks=checks, summaries=summary, comparisons=comparisons,
               outcome='runtime diagnostics; total calibrated costs unmeasured',
               witness_order_change=args.allow_witness_order_change)
    (args.directory/'comparison.json').write_text(json.dumps(out, indent=2)+'\n')
    print(json.dumps(checks))
    for r in comparisons:
        if r['variant'] == 'charts_linear' and r['split'] == 'fresh' and r['reference'] in ('reference', 'rho'):
            print(r['case'], r['reference'], r['runtime'], 'missing', r['missing'])


if __name__ == '__main__':
    main()
