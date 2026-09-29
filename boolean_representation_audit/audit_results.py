"""Recompute toy outputs and reconcile every saved timing receipt."""
import gzip
import hashlib
import json
import math
from pathlib import Path
import statistics

from measure import digest, median_interval
from backends import compute, semantic_result, verify

ROOT = Path(__file__).resolve().parent


def main():
    for line in (ROOT / 'PRE_RUN_SHA256SUMS').read_text().splitlines():
        expected, path = line.split('  ', 1)
        assert hashlib.sha256((ROOT / path).read_bytes()).hexdigest() == expected, path
    run = ROOT / 'runs' / '20260929'
    manifest = json.loads((run / 'SAMPLE_MANIFEST.json').read_text())
    raw = (run / manifest['path']).read_bytes()
    assert len(raw) == manifest['bytes']
    assert hashlib.sha256(raw).hexdigest() == manifest['sha256']
    report = json.loads((run / 'results.json').read_text())
    inputs = json.loads((ROOT / 'inputs.json').read_text())
    with gzip.open(run / 'samples.jsonl.gz', 'rt') as source:
        records = [json.loads(line) for line in source]
    samples = {r['filename']: r['data'] for r in records}
    assert len(samples) == len(records) == 504
    assert len(report['cases']) == len(inputs) == 12
    assert not report['failures'] and 'after' in report
    assert report['end_to_end_cryptographic_speedup'] is None
    used = set()
    computations = 0
    def get(name, case, arm, traced=False):
        nonlocal computations
        assert name not in used
        used.add(name)
        s = samples[name]
        assert s['case'] == case['name'] and s['arm'] == arm
        assert len(s['runs']) == (1 if traced else 3)
        computations += len(s['runs'])
        assert s['affinity']['affinity'] == [s['affinity']['cpu']]
        for r in s['runs']:
            assert r['compute_s'] > 0 and r['verify_s'] > 0
            assert math.isclose(r['compute_plus_verify_s'], r['compute_s'] + r['verify_s'])
        return s
    def mean(s, metric):
        return statistics.mean(r[metric] for r in s['runs'])
    for case, frozen in zip(report['cases'], inputs):
        assert all(case[k] == v for k, v in frozen.items())
        expected = {}
        for strategy, schedule in case['schedules'].items():
            arm_samples = {}
            for backend, arm in schedule['arms'].items():
                result = compute(case['n'], case['generators'], backend=backend, strategy=strategy)
                assert verify(case['n'], case['generators'], result)['verified'] == 'exact_ideal'
                fingerprint = digest(semantic_result(result))
                assert fingerprint == expected.setdefault(strategy, fingerprint)
                assert fingerprint == arm['semantic_sha256'] and result['stats'] == arm['stats']
                raw = [get(name, case, f'{backend}/{strategy}') for name in arm['receipts']]
                memory = get(arm['traced_receipt'], case, f'{backend}/{strategy}', True)
                for sample in raw + [memory]:
                    assert sample['semantic_sha256'] == fingerprint
                    assert sample['stats'] == result['stats']
                    assert sample['verified'] == {'rank': result['stats']['rank'], 'verified': 'exact_ideal'}
                assert arm['peak_compute_python_bytes'] == memory['peak_compute_python_bytes']
                assert arm['process_peak_rss_kib'] == [s['process_peak_rss_kib'] for s in raw]
                for metric in ('compute_s', 'compute_plus_verify_s'):
                    values = [mean(s, metric) for s in raw]
                    assert arm[metric] == {'median': statistics.median(values), 'minimum': min(values)}
                arm_samples[backend] = raw
            for metric in ('compute_s', 'compute_plus_verify_s'):
                ratios = [mean(b, metric) / mean(a, metric)
                          for a, b in zip(arm_samples['sparse'], arm_samples['packed'])]
                assert ratios == [p[metric] for p in schedule['pairs']]
                assert schedule[metric + '_packed_over_sparse'] == {
                    'median': statistics.median(ratios), 'bootstrap_95pct': median_interval(ratios)}
        for pair in case['aa']:
            a, b = [get(pair[k], case, 'sparse/frontier') for k in ('a', 'b')]
            assert a['semantic_sha256'] == b['semantic_sha256'] == expected['frontier']
            assert pair['b_over_a'] == mean(b, 'compute_s') / mean(a, 'compute_s')
    assert used == samples.keys()
    print(json.dumps({'verified': True, 'cases': len(inputs), 'samples': len(samples),
                      'verified_computations': computations, 'fresh_output_recomputations': 48}))


if __name__ == '__main__':
    main()
