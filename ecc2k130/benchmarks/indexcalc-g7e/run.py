#!/usr/bin/env python3
"""Run the G7e pair-enumeration IC measurements and freeze summary.json.

Invoked from this directory:

    python3 run.py
    python3 run.py --skip-gpu          # CPU toy DLP and accounting only
    python3 run.py --cuda ../../build/indexcalc-cuda
"""
import argparse
import hashlib
import json
import os
import platform
import subprocess
import sys
import time

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..'))
CODEGEN = os.path.join(ROOT, 'ecc2k130', 'codegen')
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, CODEGEN)

import indexcalc_pairs as pairs  # noqa: E402


def sha256File(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 16), b''):
            h.update(chunk)
    return h.hexdigest()


def runToyDlp():
    reports = []
    for m, weight, seed in ((5, 4, 20260918), (9, 2, 20260918)):
        report = pairs.recoverLog(m, weight, seed, attempts=256)
        reports.append(report)
        if report['status'] != 'complete' or not report['verified']:
            raise SystemExit('toy DLP failed at m=%d: %s' % (m, report['status']))
    return reports


def runGpu(binary, rawPath):
    cmd = [binary, '--weight', '2', '--planted', '8', '--search-generator',
           '--cpu', '--bench-points', '2048', '--bench-steps', '64',
           '--json', rawPath]
    print('+', ' '.join(cmd), flush=True)
    subprocess.check_call(cmd)
    with open(rawPath) as f:
        return json.load(f)


def merge(toy, gpu):
    base = gpu.get('factor_base_points', 0) if gpu else 0
    relations = base / 131.0 if base else 0
    projection = pairs.projectAttack(base if base else 4000, relations if relations else 4000 / 131.0,
                                     pairs.R_131, pairs.ORDER_131)
    host = {
        'hostname': platform.node(),
        'machine': platform.machine(),
        'python': sys.version.split()[0],
        'platform': platform.platform(),
    }
    out = {
        'schema': 'ecc2k130_g7e_index_calculus/v1',
        'question': 'Can a G7e GPU make Hamming-weight index calculus on ECC2K-130 cheaper than Pollard rho?',
        'verdict': 'NO. Pair enumeration is the cheapest oracle measured here; the GPU moves wall-clock, not S.',
        'unit': 'field products; S = products / sqrt(r)',
        'class': 'engineering',
        'host': host,
        'boundaries': {
            'rho_log2_operations': pairs.RHO_LOG2,
            'rho_S': pairs.sScore(2 ** pairs.RHO_LOG2, pairs.R_131),
            'product_law_floor_log2': pairs.log2(3 * pairs.ORDER_131),
            'product_law_note': 'm·#E field-point candidates for m=3, independent of |F| at leading order',
            'falsification': 'an oracle whose streaming product count at some |F| falls a factor 2^70 below exhaustive pair search, with every hit verified on the curve',
        },
        'toy_dlp': toy,
        'gpu': gpu,
        'projection': projection,
        'source_hashes': {
            'indexcalccuda.cu': sha256File(os.path.join(ROOT, 'ecc2k130', 'src', 'indexcalccuda.cu')),
            'indexcalc_pairs.py': sha256File(os.path.join(CODEGEN, 'indexcalc_pairs.py')),
        },
    }
    if gpu:
        rate = gpu.get('bench_affine_adds_per_second') or gpu.get('generator_gpu_pairs_per_second')
        products = projection['streaming_field_products']
        if rate and products < float('inf'):
            # Each pair is two affine adds.
            seconds = (products / pairs.PRODUCTS_PER_AFFINE_ADD) / rate
            out['practicality'] = {
                'affine_adds_per_second': rate,
                'projected_streaming_seconds': seconds,
                'projected_streaming_log2_seconds': pairs.log2(seconds) if seconds > 0 else None,
                'note': 'wall-clock at the measured add rate; not the metric',
            }
    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--skip-gpu', action='store_true')
    parser.add_argument('--cuda', default=os.path.join(ROOT, 'ecc2k130', 'build', 'indexcalc-cuda'))
    args = parser.parse_args()
    t0 = time.time()
    print('toy pair-enumeration DLP', flush=True)
    toy = runToyDlp()
    gpu = None
    rawPath = os.path.join(HERE, 'raw-gpu.json')
    if not args.skip_gpu:
        if not os.path.isfile(args.cuda):
            raise SystemExit('missing CUDA binary %s; run make indexcalc-cuda' % args.cuda)
        gpu = runGpu(args.cuda, rawPath)
    summary = merge(toy, gpu)
    summary['elapsed_s'] = time.time() - t0
    out = os.path.join(HERE, 'summary.json')
    with open(out, 'w') as f:
        json.dump(summary, f, indent=2, sort_keys=True)
        f.write('\n')
    print('wrote', out, flush=True)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
