#!/usr/bin/env python3
"""Preserve matched cold native F4/DLP runs, including failures and holdouts."""
import argparse
import hashlib
import itertools
import json
import os
import pathlib
import platform
import resource
import subprocess
import time

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--reference-dir', type=pathlib.Path, required=True)
    parser.add_argument('--candidate-dir', type=pathlib.Path, required=True)
    parser.add_argument('--kernel', type=pathlib.Path, required=True)
    parser.add_argument('--output', type=pathlib.Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    contract = json.loads((HERE / 'contract.json').read_text())
    (args.output / 'contract.json').write_text(json.dumps(contract, indent=2) + '\n')
    env = {k: v for k, v in os.environ.items() if not k.startswith('IC_')}
    env.update(RAYON_NUM_THREADS='1', IC_PREPROCESS_CACHE='off', IC_REDUCTION_CACHE='off')
    cpu = min(os.sched_getaffinity(0))
    os.sched_setaffinity(0, {cpu})
    binaries = {'kernel': args.kernel.resolve()}
    for variant, directory in [('reference', args.reference_dir), ('candidate', args.candidate_dir)]:
        for name in ['ic', 'f4_linear_algebra_bench']:
            binaries[variant + '-' + name] = (directory / name).resolve()
    sources = [ROOT / 'src/cryptanalysis/koblitz_groebner.rs',
               ROOT / 'src/cryptanalysis/f4_rref_tests.rs',
               ROOT / 'examples/f4_linear_algebra_bench.rs', HERE / 'run.py', HERE / 'compare.py',
               HERE / 'contract.json', ROOT / 'Cargo.lock']
    hashes = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
    reference_source = subprocess.check_output(['git', 'show',
        contract['base_revision'] + ':src/cryptanalysis/koblitz_groebner.rs'], cwd=ROOT)
    hashes['reference_koblitz_groebner.rs'] = hashlib.sha256(reference_source).hexdigest()
    hashes.update({name: hashlib.sha256(p.read_bytes()).hexdigest() for name, p in binaries.items()})
    (args.output / 'hashes.json').write_text(json.dumps(hashes, indent=2) + '\n')
    hardware = {'platform': platform.platform(), 'cpu_affinity': [cpu],
                'cpu_model': next((s.split(':', 1)[1].strip() for s in pathlib.Path('/proc/cpuinfo').read_text().splitlines() if s.startswith('model name')), None),
                'rust': '1.90.0', 'profile': 'release, redis-cache feature, no RUSTFLAGS',
                'cache': 'off', 'rayon_threads': 1}
    (args.output / 'hardware.json').write_text(json.dumps(hardware, indent=2) + '\n')
    records = []

    def limits():
        cap = contract['address_space_mib'] * 1024**2
        resource.setrlimit(resource.RLIMIT_AS, (cap, cap))

    def run(name, command):
        command = list(map(str, command))
        start = time.monotonic()
        status, code = 'finished', None
        try:
            result = subprocess.run(command, env=env, capture_output=True,
                                    timeout=contract['process_timeout_seconds'], preexec_fn=limits)
            stdout, stderr, code = result.stdout, result.stderr, result.returncode
        except subprocess.TimeoutExpired as error:
            status, stdout, stderr = 'timeout', error.stdout or b'', error.stderr or b''
        (args.output / (name + '.stdout')).write_bytes(stdout)
        (args.output / (name + '.stderr')).write_bytes(stderr)
        records.append({'name': name, 'command': command, 'status': status, 'returncode': code,
                        'wall_seconds': time.monotonic() - start})
        (args.output / 'processes.json').write_text(json.dumps(records, indent=2) + '\n')
        print(name, status, code, flush=True)

    run('kernel', [binaries['kernel'], 'cryptanalysis::koblitz_groebner::rref_tests::paired_kernel_benchmark',
                   '--exact', '--ignored', '--nocapture', '--test-threads=1'])
    for (n, ell, m), split, rep in itertools.product(contract['stages']['cases'],
            contract['stages']['splits'], range(contract['stages']['paired_repetitions'])):
        for variant in (['reference', 'candidate'] if rep % 2 == 0 else ['candidate', 'reference']):
            run(f'stage-n{n}-l{ell}-m{m}-{split}-r{rep}-{variant}',
                [binaries[variant + '-f4_linear_algebra_bench'], n, ell, m, split])
    for n, solver, (seed, log, split), rep in itertools.product(contract['dlp']['degrees'],
            contract['dlp']['solvers'], contract['dlp']['fixtures'], range(contract['dlp']['paired_repetitions'])):
        for variant in (['reference', 'candidate'] if rep % 2 == 0 else ['candidate', 'reference']):
            run(f'dlp-n{n}-{solver}-s{seed}-k{log}-{split}-r{rep}-{variant}',
                [binaries[variant + '-ic'], 'run', '--degree', n, '--solver', solver,
                 '--seed', seed, '--known-log', log, '--batch', contract['dlp']['batch'],
                 '--max-trials', contract['dlp']['max_trials'], '--json'])
    print('Saved', len(records), 'process records', flush=True)


if __name__ == '__main__':
    main()
