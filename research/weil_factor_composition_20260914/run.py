#!/usr/bin/env python3
"""Run the frozen suite sequentially on one CPU, retaining every outcome."""
import argparse
import hashlib
import itertools
import json
import os
from pathlib import Path
import platform
import resource
import subprocess
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--reference-native', type=Path, required=True)
    ap.add_argument('--candidate-native', type=Path, required=True)
    ap.add_argument('--experiment', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    contract = json.loads((HERE / 'contract.json').read_text())
    (args.output / 'contract.json').write_bytes((HERE / 'contract.json').read_bytes())
    binaries = {k: getattr(args, k).resolve() for k in ['reference_native', 'candidate_native', 'experiment']}
    hashes = {k: sha(p) for k, p in binaries.items()}
    assert hashes['reference_native'] != hashes['candidate_native'], 'require independent revision builds'
    sources = ['src/cryptanalysis/weil_charts.rs', 'src/cryptanalysis/koblitz_index_calculus.rs',
               'src/cryptanalysis/koblitz_groebner.rs', 'examples/weil_factor_composition.rs',
               'examples/f4_linear_algebra_bench.rs', 'Cargo.lock']
    sources += [str(p.relative_to(ROOT)) for p in HERE.glob('*.py')] + [str((HERE / 'contract.json').relative_to(ROOT))]
    hashes.update({p: sha(ROOT / p) for p in sources})
    (args.output / 'hashes.json').write_text(json.dumps(hashes, indent=2) + '\n')
    cpu = min(os.sched_getaffinity(0))
    os.sched_setaffinity(0, {cpu})
    env = {k: v for k, v in os.environ.items() if not k.startswith('IC_')}
    env.update(RAYON_NUM_THREADS='1', IC_PREPROCESS_CACHE='off', IC_REDUCTION_CACHE='off')
    (args.output / 'hardware.json').write_text(json.dumps(dict(platform=platform.platform(), cpu_affinity=[cpu],
        cpu_model=next(s.split(':', 1)[1].strip() for s in Path('/proc/cpuinfo').read_text().splitlines() if s.startswith('model name')),
        rust='1.90.0 release; redis-cache feature; no RUSTFLAGS', caches='off', threads=1,
        address_space_mib=2048, comparison='ambient baseline and new charts in same experiment binary; unchanged default additionally compared to reference revision binary'), indent=2) + '\n')

    def limits():
        cap = 2048 * 1024**2
        resource.setrlimit(resource.RLIMIT_AS, (cap, cap))

    archive = (args.output / 'processes.jsonl').open('x')
    count = 0

    def run(kind, command, timeout, **key):
        nonlocal count
        command = list(map(str, command))
        started = time.monotonic()
        status, code = 'finished', None
        try:
            r = subprocess.run(command, env=env, capture_output=True, timeout=timeout, preexec_fn=limits)
            stdout, stderr, code = r.stdout, r.stderr, r.returncode
            if code:
                status = 'error'
        except subprocess.TimeoutExpired as e:
            status, stdout, stderr = 'timeout', e.stdout or b'', e.stderr or b''
        record = dict(kind=kind, command=command, status=status, returncode=code,
                      wall_seconds=time.monotonic() - started, stdout=stdout.decode(), stderr=stderr.decode(), **key)
        archive.write(json.dumps(record, separators=(',', ':')) + '\n')
        archive.flush()
        count += 1
        print(count, kind, key, status, code, flush=True)

    for ci, (n, ell, m) in enumerate(contract['frozen_regression']['native_cases']):
        for split, rep in itertools.product(['frozen', 'fresh'], range(3)):
            for variant in (['reference', 'candidate'] if rep % 2 == 0 else ['candidate', 'reference']):
                run('native', [binaries[variant + '_native'], n, ell, m, split], 60,
                    case=ci, split=split, rep=rep, variant=variant)
    for ci, seed, rep in itertools.product(range(len(contract['stage_cases'])), contract['stage']['seeds'], range(3)):
        vs = contract['stage']['variants']
        for variant in (vs if rep % 2 == 0 else list(reversed(vs))):
            run('stage', [binaries['experiment'], 'stage', ci, variant, seed], contract['stage']['timeout_seconds'],
                case=ci, seed=seed, rep=rep, variant=variant)
    for ci, (seed, log), rep in itertools.product(range(len(contract['stage_cases'])), contract['dlp']['seeds_and_logs'], range(3)):
        vs = contract['dlp']['variants']
        for variant in (vs if rep % 2 == 0 else list(reversed(vs))):
            run('dlp', [binaries['experiment'], 'dlp', ci, variant, seed, log], contract['dlp']['timeout_seconds'],
                case=ci, seed=seed, known_log=log, rep=rep, variant=variant)
    # One matched rho workload per distinct curve, not per factor-base family.
    seen = set()
    for ci, c in enumerate(contract['stage_cases']):
        curve = tuple(c[k] for k in ['n', 'k', 'a', 'b'])
        if curve in seen:
            continue
        seen.add(curve)
        for (seed, log), rep in itertools.product(contract['dlp']['seeds_and_logs'], range(3)):
            run('rho', [binaries['experiment'], 'rho', ci, 'rho', seed, log], 10,
                case=ci, seed=seed, known_log=log, rep=rep, variant='rho')
    archive.close()
    (args.output / 'complete.json').write_text(json.dumps({'processes': count, 'complete': True}) + '\n')


if __name__ == '__main__':
    main()
