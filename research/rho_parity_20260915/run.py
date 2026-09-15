#!/usr/bin/env python3
"""Frozen matched revision comparisons; retain raw output including failures."""
import argparse
import gzip
import hashlib
import itertools
import json
import os
from pathlib import Path
import platform
import resource
import random
import subprocess
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
FROZEN = [(101, 5), (202, 5), (503, 7), (607, 11)]
FRESH = [(1201, 13), (1429, 17), (1601, 19), (1873, 23)]
CONFIRM_SEEDS = [2203, 2213, 2221, 2237]
ORDERS = [71, 127, 127, 127, 991, 991, 2003, 661, 661, 661]


def make_workloads(cases, confirmation):
    if not confirmation:
        return [FROZEN+FRESH for _ in cases]
    workloads = []
    for ci, c in enumerate(cases):
        canonical = next(j for j, cc in enumerate(cases) if all(c[k] == cc[k] for k in ['n','k','a','b']))
        fresh = [(seed, random.Random(f'rho-confirmation-20260915:{canonical}:{seed}').randrange(1, ORDERS[ci]))
                 for seed in CONFIRM_SEEDS]
        workloads.append(FROZEN+fresh)
    return workloads


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--reference', type=Path, required=True)
    ap.add_argument('--candidate', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--reference-revision', default='a2963f0d6fb5cb74a20e678cf4def0c7cbafa1d9')
    ap.add_argument('--confirm', action='store_true')
    args = ap.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    bins = {'reference': args.reference.resolve(), 'candidate': args.candidate.resolve()}
    contract = json.loads((HERE.parent / 'weil_factor_composition_20260914/contract.json').read_text())
    workloads = make_workloads(contract['stage_cases'], args.confirm)
    contract.update(reference_revision=args.reference_revision,
                    frozen=FROZEN, fresh=workloads[0][4:],
                    fresh_note='Seed list; exact scalar labels for every curve are in dlp_workloads.',
                    confirmation=args.confirm, dlp_workloads=workloads,
                    stage_seeds=[17, 937, 2203, 2237] if args.confirm else [17, 937, 1201, 1873],
                    holdout_status='previously unused confirmation seeds and uniform nonzero scalars' if args.confirm else
                        'introduced in iteration 1; reused as validation inputs in subsequent iterations',
                    dlp_variants=['charts_linear', 'enumerate', 'pair_table'])
    (args.output / 'contract.json').write_text(json.dumps(contract, indent=2)+'\n')
    (args.output / 'protocol.md').write_bytes((HERE / 'README.md').read_bytes())
    sources = ['src/cryptanalysis/weil_charts.rs', 'src/cryptanalysis/koblitz_fast.rs',
               'src/cryptanalysis/koblitz_index_calculus.rs', 'src/cryptanalysis/koblitz_groebner.rs',
               'examples/weil_factor_composition.rs', 'examples/f4_linear_algebra_bench.rs', 'Cargo.lock']
    sources += [str(p.relative_to(ROOT)) for p in HERE.glob('*.py')]
    hashes = {p: sha(ROOT / p) for p in sources}
    for rev, path in bins.items():
        for name in ['weil_factor_composition', 'f4_linear_algebra_bench']:
            hashes[rev+'/'+name] = sha(path / name)
    (args.output / 'hashes.json').write_text(json.dumps(hashes, indent=2)+'\n')
    cpu = min(os.sched_getaffinity(0))
    os.sched_setaffinity(0, {cpu})
    env = {k: v for k, v in os.environ.items() if not k.startswith('IC_')}
    env.update(RAYON_NUM_THREADS='1', IC_PREPROCESS_CACHE='off', IC_REDUCTION_CACHE='off')
    hardware = dict(platform=platform.platform(), cpu_affinity=[cpu], threads=1,
                    cpu_model=next(s.split(':', 1)[1].strip() for s in Path('/proc/cpuinfo').read_text().splitlines() if s.startswith('model name')),
                    rust='1.90.0 release; redis-cache feature; no RUSTFLAGS', caches='off', address_space_mib=2048)
    (args.output / 'hardware.json').write_text(json.dumps(hardware, indent=2)+'\n')

    def limits():
        resource.setrlimit(resource.RLIMIT_AS, (2048*1024**2,)*2)

    archive = gzip.open(args.output / 'processes.jsonl.gz', 'wt')
    count = 0

    def run(kind, rev, command, timeout=10, **key):
        nonlocal count
        started = time.monotonic()
        command = list(map(str, command))
        status, code = 'finished', None
        try:
            r = subprocess.run(command, env=env, capture_output=True, timeout=timeout, preexec_fn=limits)
            stdout, stderr, code = r.stdout, r.stderr, r.returncode
            if code:
                status = 'error'
        except subprocess.TimeoutExpired as e:
            status, stdout, stderr = 'timeout', e.stdout or b'', e.stderr or b''
        record = dict(kind=kind, revision=rev, command=command, status=status, returncode=code,
                      wall_seconds=time.monotonic()-started, stdout=stdout.decode(), stderr=stderr.decode(), **key)
        archive.write(json.dumps(record, separators=(',', ':'))+'\n')
        archive.flush()
        count += 1
        if count % 100 == 0 or status != 'finished':
            print(count, kind, rev, key, status, flush=True)

    def revisions(rep):
        return ['reference', 'candidate'] if rep % 2 == 0 else ['candidate', 'reference']

    for ci, (n, ell, m) in enumerate(contract['frozen_regression']['native_cases']):
        for split, rep in itertools.product(['frozen', 'fresh'], range(3)):
            for rev in revisions(rep):
                run('native', rev, [bins[rev]/'f4_linear_algebra_bench', n, ell, m, split], 60,
                    case=ci, split=split, rep=rep, variant='default')
    for ci, seed, rep in itertools.product(range(10), contract['stage_seeds'], range(3)):
        for rev in revisions(rep):
            run('stage', rev, [bins[rev]/'weil_factor_composition', 'stage', ci, 'charts_linear', seed],
                case=ci, seed=seed, rep=rep, variant='charts_linear')
    for ci, wi, rep in itertools.product(range(10), range(8), range(3)):
        seed, log = workloads[ci][wi]
        for variant in contract['dlp_variants']:
            for rev in revisions(rep):
                run('dlp', rev, [bins[rev]/'weil_factor_composition', 'dlp', ci, variant, seed, log],
                    case=ci, seed=seed, known_log=log, rep=rep, variant=variant)
    seen = set()
    for ci, c in enumerate(contract['stage_cases']):
        curve = tuple(c[k] for k in ['n', 'k', 'a', 'b'])
        if curve in seen:
            continue
        seen.add(curve)
        for (seed, log), rep in itertools.product(workloads[ci], range(3)):
            for rev in revisions(rep):
                run('rho', rev, [bins[rev]/'weil_factor_composition', 'rho', ci, 'rho', seed, log],
                    case=ci, seed=seed, known_log=log, rep=rep, variant='rho')
    archive.close()
    (args.output / 'complete.json').write_text(json.dumps(dict(processes=count, complete=True))+'\n')
    print('complete', count, flush=True)


if __name__ == '__main__':
    main()
