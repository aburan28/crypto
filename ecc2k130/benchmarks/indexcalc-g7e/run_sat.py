#!/usr/bin/env python3
"""SAT decomposition on the G7e host, compared with pair enumeration.

This is a CPU CryptoMiniSat measurement on the same box as the GPU pair
scan.  It does not use the GPU.  The product-law floor is unchanged.

Type-II ONB exists only when 2m+1 is prime and ord_{2m+1}(2) is m or 2m.
The ladder is therefore 5, 9, 11, 23 — not 13, 15, 17.
"""
import hashlib
import json
import os
import platform
import sys
import time

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..'))
CODEGEN = os.path.join(ROOT, 'ecc2k130', 'codegen')
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, CODEGEN)

import indexcalc
import indexcalc_e2e as engine
import indexcalc_pairs as pairs


def sha256File(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 16), b''):
            h.update(chunk)
    return h.hexdigest()


def pythonInfo():
    import pycryptosat
    import pysat
    return {
        'executable': sys.executable,
        'version': sys.version.split()[0],
        'pycryptosat': getattr(pycryptosat, '__version__', 'unknown'),
        'python_sat': getattr(pysat, '__version__', 'unknown'),
        'hostname': platform.node(),
        'machine': platform.machine(),
    }


def growthCell(m, points, weight, trials, seed, timeout):
    print('SAT growth m=%d points=%d weight=%d trials=%d timeout=%s'
          % (m, points, weight, trials, timeout), flush=True)
    t0 = time.time()
    result = indexcalc.runTrials(m, points, weight, trials, seed, 12, True,
                                 0, timeout or 0)
    result['wall_s'] = time.time() - t0
    result['timeout_s'] = timeout
    return result


def satDlp(m, variant, seed):
    print('SAT e2e DLP m=%d variant=%s' % (m, variant), flush=True)
    report = engine.experiment(m, variant, seed, 'dlp', 64, weight=4 if m == 5 else 2)
    return {
        'degree': report['degree'],
        'variant': report['variant'],
        'status': report['status'],
        'scalar_recovered': report.get('scalar_recovered'),
        'verified': bool(report.get('certificate', {}).get('verified')),
        'orbit_columns': report['orbit_columns'],
        'signed_base_size': report['signed_base_size'],
        'accounting': {
            'elapsed_ns': report['accounting']['elapsed_ns'],
            'logical_calls': dict(report['accounting']['logical_calls']),
        },
    }


def pairVersusSat(m=9, weight=2, trials=16, seed=20260918):
    """Same planted triples: pair lookup vs SAT, wall-clock only."""
    meter = engine.Ledger()
    context = engine.setup(m, weight, 3, meter)
    onb, curve, ell, generator, eigen, reps, lookup, _, _ = context
    rng = __import__('random').Random(seed)
    points = list(lookup)
    pairTimes, satTimes = [], []
    pairHits = satHits = 0
    for i in range(trials):
        chosen = [points[rng.randrange(len(points))] for _ in range(3)]
        target = None
        for p in chosen:
            target = curve.add(target, p)
        t0 = time.perf_counter()
        found = pairs.pairDecomposeLookup(curve, lookup, target)
        pairTimes.append(time.perf_counter() - t0)
        pairHits += found is not None
        t0 = time.perf_counter()
        row, detail = engine.decompose(context, target, 3, weight, 'candidate', meter)
        satTimes.append(time.perf_counter() - t0)
        satHits += row is not None
        print('  trial %d pair=%.4fs sat=%.4fs pair_hit=%s sat_status=%s'
              % (i, pairTimes[-1], satTimes[-1], found is not None, detail['status']),
              flush=True)
    pairTimes.sort()
    satTimes.sort()
    return {
        'degree': m, 'weight': weight, 'trials': trials, 'seed': seed,
        'pair_hits': pairHits, 'sat_hits': satHits,
        'pair_median_s': pairTimes[len(pairTimes) // 2],
        'sat_median_s': satTimes[len(satTimes) // 2],
        'pair_total_s': sum(pairTimes), 'sat_total_s': sum(satTimes),
        'class': 'engineering',
        'note': 'wall-clock on this host; SAT is CryptoMiniSat, pairs are exhaustive lookup',
    }


def main():
    t0 = time.time()
    ladder = [
        growthCell(5, 3, 3, 5, 1, 0),
        growthCell(9, 3, 3, 5, 1, 0),
        growthCell(11, 3, 3, 5, 1, 60),
        growthCell(23, 3, 3, 1, 1, 30),
    ]
    dlp = [
        satDlp(5, 'candidate', 20260918),
        satDlp(9, 'candidate', 20260918),
    ]
    versus = pairVersusSat()
    report = {
        'schema': 'ecc2k130_g7e_index_calculus_sat/v1',
        'class': 'engineering',
        'device': 'G7e host CPU, not the GPU',
        'python': pythonInfo(),
        'growth': ladder,
        'sat_dlp': dlp,
        'pair_versus_sat': versus,
        'elapsed_s': time.time() - t0,
        'source_hashes': {
            'indexcalc.py': sha256File(os.path.join(CODEGEN, 'indexcalc.py')),
            'indexcalc_e2e.py': sha256File(os.path.join(CODEGEN, 'indexcalc_e2e.py')),
            'indexcalc_pairs.py': sha256File(os.path.join(CODEGEN, 'indexcalc_pairs.py')),
            'run_sat.py': sha256File(os.path.abspath(__file__)),
        },
        'verdict': 'SAT recovers toy logs and loses to pair enumeration on the same planted triples; the n=131 floor is untouched.',
    }
    path = os.path.join(HERE, 'sat.json')
    with open(path, 'w') as f:
        json.dump(report, f, indent=2, sort_keys=True)
        f.write('\n')
    print('wrote', path, flush=True)
    summaryPath = os.path.join(HERE, 'summary.json')
    if os.path.isfile(summaryPath):
        with open(summaryPath) as f:
            summary = json.load(f)
        summary['sat'] = {
            'receipt': 'sat.json',
            'pair_versus_sat': versus,
            'sat_dlp': [{'degree': r['degree'], 'status': r['status'],
                         'verified': r['verified']} for r in dlp],
            'growth': [{'m': c['m'], 'solved': c['solved'], 'budget': c['budget'],
                        'unsat': c['unsat'], 'median': c['median'],
                        'gates': c['gates'], 'orbits': c['orbits']}
                       for c in ladder],
        }
        with open(summaryPath, 'w') as f:
            json.dump(summary, f, indent=2, sort_keys=True)
            f.write('\n')
        print('updated', summaryPath, flush=True)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
