#!/usr/bin/env python3
"""Portable IC development loop; the frozen tournament owns promotion.

Native screens run complete cold ECDLP jobs with the existing independent
certificate checker. They do not fabricate operation counts or promote winners.
"""
import argparse
import copy
import fcntl
import json
import math
import os
from pathlib import Path
import platform
import random
import re
import shutil
import signal
import statistics
import subprocess
import sys
import threading
import time

from oracle import Curve, InvalidEvidence, require, verify
from portfolio import factorial_candidates, recombine
from tournament import BASE_CONFIG, child_env, digest, objhash, parse_cells, read, write

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
EVALUATOR = ('autolab.py', 'portfolio.py', 'tournament.py', 'oracle.py')


def doctor(source):
    required = ['Cargo.toml', 'Cargo.lock', 'examples/ic_tournament_worker.rs',
                'src/cryptanalysis/koblitz_index_calculus.rs']
    checks = {p: (source/p).is_file() for p in required}
    vg = shutil.which('valgrind')
    version = subprocess.check_output([vg, '--version'], text=True).strip() if vg else None
    return {'source_root': str(source), 'source_files': checks, 'host': platform.uname()._asdict(),
            'native_build_ready': all(checks.values()) and bool(shutil.which('cargo')),
            'instruction_protocol_ready': platform.system() == 'Linux' and platform.machine() == 'x86_64'
                and version == 'valgrind-3.22.0', 'valgrind': version,
            'baseline_quality': 'Requires fresh rho sensitivity screen and arithmetic equivalence tests; no global-best certification.',
            'native_scope': 'Development diagnostic; no operation-count or hardware-independent speedup claim.',
            'stages': {
                'factor_base': ['subgroup_orbits', 'factor', 'divisor', 'frobenius_union', 'pruned', 'two_torsion_saturated'],
                'pdp': ['pair_table', 'enumerate', 'f4', 'f5', 'inherited_f4', 'sat_xor', 'sat_cnf'],
                'relation_processing': ['verify', 'deduplicate', 'rank over subgroup scalar field'],
                'linear_algebra': ['dense', 'sparse filtering and block Wiedemann'],
                'descent': ['verified factor-base relation per target'],
                'reference': ['signed Frobenius rho, packed arithmetic, batched inversions'],
                'external_adapters': 'WDSat/GPU/other solvers require executable pinning, equivalent certificates and child/device accounting before admission.'}}


def run_native(binary, job, directory, timeout):
    """Time the entire child, with a watchdog that kills its process group.

    On macOS, no fictitious affinity or memory cap is reported. These timings
    are development diagnostics, not substitutes for the Linux promotion gate.
    """
    directory.mkdir(parents=True, exist_ok=False)
    write(directory/'job.json', job, exclusive=True)
    expired = threading.Event()
    with (directory/'stdout.json').open('wb') as stdout, (directory/'stderr.txt').open('wb') as stderr:
        start = time.monotonic()
        process = subprocess.Popen([str(binary)], stdin=subprocess.PIPE, stdout=stdout,
            stderr=stderr, env=child_env(), start_new_session=True)
        def kill():
            try:
                os.killpg(process.pid, signal.SIGKILL)
                expired.set()
            except ProcessLookupError:
                pass
        timer = threading.Timer(timeout, kill)
        timer.daemon = True
        timer.start()
        try:
            process.communicate(json.dumps(job, allow_nan=False).encode())
        finally:
            timer.cancel()
            timer.join()
            if process.poll() is None:
                kill()
                process.wait()
        elapsed = time.monotonic()-start
    return {'exit_code': process.returncode, 'status': 'TIMEOUT' if expired.is_set() else 'EXITED',
            'whole_process_wall_seconds': elapsed, 'timeout_seconds': timeout,
            'memory_cap_bytes': None, 'cpu_affinity': None, 'peak_rss_bytes': None}


def source_manifest(source):
    paths = [source/p for p in ('Cargo.toml', 'Cargo.lock', 'build.rs', '.cargo/config.toml',
                              'examples/ic_tournament_worker.rs') if (source/p).is_file()]
    paths += [p for p in (source/'src').rglob('*') if p.is_file()]
    # Match the existing snapshot builder's explicit compile-time dependencies.
    for rust in (source/'src').rglob('*.rs'):
        for relative in re.findall(r'include(?:_str|_bytes)?!\s*\(\s*"([^"]+)"', rust.read_text()):
            dep = (rust.parent/relative).resolve()
            require(dep.is_relative_to(source), 'source dependency escapes checkout')
            require(dep.is_file(), 'missing compile-time dependency: '+str(dep))
            paths.append(dep)
    return {str(p.relative_to(source)): digest(p) for p in sorted(set(paths))}


def check_arms(arms):
    require(isinstance(arms, list) and 2 <= len(arms) <= 16, 'candidate count must be 2..16')
    require(arms[0]['id'] == 'incumbent', 'first candidate must be incumbent')
    ids = [a['id'] for a in arms]
    require(len(set(ids)) == len(ids), 'duplicate candidate id')
    require(all(re.fullmatch(r'[a-z][a-z0-9_-]{0,31}', a) for a in ids), 'invalid candidate id')
    require(not set(ids) & {'rho', 'aa_control'} and not any(x.startswith('rho_') for x in ids), 'reserved candidate id')
    require(all(not any(k in a for k in ('source_root', 'binary_relative', 'source_directory'))
                for a in arms), 'native screen uses one frozen source; use tournament for source candidates')
    require(len({objhash(a['config']) for a in arms}) == len(arms), 'duplicate configuration')


def paired_wall_ratio(contract, rows, candidate, baseline='incumbent'):
    """Equal-cell mean of paired case log ratios; repetitions aren't targets."""
    cells = {}
    for case in contract['cases']:
        pair = [[r for r in rows if r['case']==case['id'] and r['arm']==arm]
                for arm in (baseline,candidate)]
        for group in pair:
            if (len(group)!=contract['repetitions'] or
                {r['repetition'] for r in group}!=set(range(contract['repetitions'])) or
                any(r['status']!='VERIFIED' for r in group)):
                return None
        values = [statistics.median(r['process']['whole_process_wall_seconds'] for r in g) for g in pair]
        require(all(math.isfinite(v) and v>0 for v in values), 'invalid paired time')
        job = case.get('job',{})
        cell = (job.get('degree'),job.get('curve_a'))
        cells.setdefault(cell,[]).append(math.log(values[1]/values[0]))
    return math.exp(statistics.mean(statistics.mean(v) for v in cells.values())) if cells else None


def summarize(contract, rows):
    # Check this on both first execution and independent replay of a screen.
    for case in contract['cases']:
        valid = [r for r in rows if r['case'] == case['id'] and r['status'] == 'VERIFIED'
                 and not r['arm'].startswith('rho_')]
        for arm in contract['arms']:
            require(len({r['certificate']['factor_base_sha256'] for r in valid
                         if r['arm'] == arm['id']}) <= 1, 'unstable factor-base support within arm')
        if contract['comparison_kind'] == 'fixed-support':
            require(len({r['certificate']['factor_base_sha256'] for r in valid}) <= 1,
                    'changed factor base in fixed-support screen')
    table = []
    for arm in contract['arms']:
        records = sorted([r for r in rows if r['arm'] == arm['id']],
                         key=lambda r: (r['case'], r['repetition']))
        expected = len(contract['cases'])*contract['repetitions']
        complete = len(records) == expected and all(r['status'] == 'VERIFIED' for r in records)
        table.append({'arm': arm['id'], 'mode': arm['mode'], 'scheduled': expected,
            'verified': sum(r['status'] == 'VERIFIED' for r in records),
            'complete': complete,
            'median_cold_wall_seconds': statistics.median(r['process']['whole_process_wall_seconds']
                for r in records) if complete else None,
            'total_operations': None, 'S': None, 'ratio_to_floor': None,
            'rank_range': [min(r['certificate']['rank'] for r in records),
                           max(r['certificate']['rank'] for r in records)] if complete and arm['mode']=='ic' else None,
            'base_sizes': sorted({r['certificate']['signed_base_size'] for r in records}) if complete and arm['mode']=='ic' else [],
            'failures': [{'case': r['case'], 'repetition': r['repetition'], 'reason': r.get('reason')}
                         for r in records if r['status'] != 'VERIFIED']})
    for row in table:
        row['diagnostic_time_over_incumbent'] = paired_wall_ratio(contract, rows, row['arm'])
    aa = paired_wall_ratio(contract,rows,'aa_control')
    # This preserves every candidate; development results never lock a winner.
    return {'status': 'DEVELOPMENT_DIAGNOSTIC', 'promotion_eligible': False,
            'reason': 'No instruction accounting, held-out confirmation, or baseline quality certification.',
            'unit': 'complete cold child wall seconds', 'table': table,
            'ratio_definition': 'Equal-cell geometric mean of paired case median time ratios; diagnostic only.',
            'aa_diagnostic': {'paired_ratio':aa, 'within_five_percent':aa is not None and .95<=aa<=1.05,
                              'runtime_claim_eligible':False},
            'all_jobs_verified': all(r['complete'] for r in table),
            'retained_candidates': [a['id'] for a in contract['arms'] if a['mode'] == 'ic' and a['id'] != 'aa_control'],
            'next_action': 'Retain distinct mechanisms; run actual combinations on fresh development cases. Freeze one challenger before confirmation.',
            'global_optimum_claim': False}


def verify_trial(root, contract, case, arm, repetition):
    directory = root/'trials'/case['id']/arm['id']/f'rep-{repetition}'
    receipt = read(directory/'receipt.json')
    require(receipt['case'] == case['id'] and receipt['arm'] == arm['id']
        and receipt['repetition'] == repetition, 'changed trial identity')
    require(receipt['case_sha256'] == objhash(case) and receipt['arm_sha256'] == objhash(arm), 'changed trial binding')
    require(set(receipt['artifacts']) == {'job.json', 'stdout.json', 'stderr.txt', 'process.json'}, 'missing raw evidence')
    for name, expected in receipt['artifacts'].items():
        require(digest(directory/name) == expected, 'changed trial artifact: '+name)
    job = dict(copy.deepcopy(case['job']), config=arm['config'], mode=arm['mode'])
    require(read(directory/'job.json') == job, 'wrong executed job')
    require(read(directory/'process.json') == receipt['process'], 'changed process accounting')
    p = receipt['process']
    require(math.isfinite(p['whole_process_wall_seconds']) and p['whole_process_wall_seconds'] > 0, 'invalid wall measurement')
    if receipt['status'] == 'VERIFIED':
        require(p['status'] == 'EXITED' and p['exit_code'] == 0, 'failed process marked verified')
        proof = verify(read(directory/'stdout.json'), case['fixture'], expected_mode=arm['mode'], summands=arm['config']['summands'])
        require(proof == receipt['certificate'], 'changed certificate')
        require(not proof.get('degenerate_descents', 0), 'generic direct relation admitted as IC')
    else:
        require(receipt['status'] in ('TIMEOUT', 'INVALID_OR_INCOMPLETE') and receipt['certificate'] is None,
                'invalid failure receipt')
    return receipt


def audit(root, *, complete=True):
    contract = read(root/'contract.json')
    require(objhash(contract) == read(root/'seal.json')['sha256'], 'changed contract')
    for name, expected in contract['pinned'].items():
        require(digest(root/name) == expected, 'changed pinned file: '+name)
    for name, expected in contract['source_manifest'].items():
        require(digest(root/'source'/name) == expected, 'changed source snapshot: '+name)
    require(all(digest(HERE/n) == contract['pinned']['evaluator/'+n] for n in EVALUATOR),
            'use the frozen evaluator')
    rows = []
    for case in contract['cases']:
        for arm in contract['arms']:
            for rep in range(contract['repetitions']):
                path = root/'trials'/case['id']/arm['id']/f'rep-{rep}'
                if not path.exists() and not complete:
                    continue
                rows.append(verify_trial(root, contract, case, arm, rep))
    if complete:
        require(read(root/'summary.json') == summarize(contract, rows), 'changed summary')
    return contract, rows


def prepare(args):
    source, root = args.source_root.resolve(), args.out.resolve()
    arms = read(args.candidates)
    check_arms(arms)
    require(not root.exists(), 'output already exists')
    require(1 <= args.cases <= 30 and 1 <= args.repetitions <= 5 and 0 < args.timeout <= 300,
            'bounded cases/repetitions/timeout required')
    cells = parse_cells(args.cells)
    expected = len(cells)*args.cases*args.repetitions*(len(arms)+4)
    require(expected <= args.max_processes, 'planned screen exceeds process budget')
    prior = None
    if args.source_screen:
        previous = args.source_screen.resolve()
        prior = read(previous/'contract.json')
        require(objhash(prior) == read(previous/'seal.json')['sha256'], 'changed source-screen contract')
        for name, h in prior['pinned'].items():
            require(digest(previous/name) == h, 'changed source-screen pinned artifact')
        manifest = prior['source_manifest']
        source = previous/'source'
        for name, h in manifest.items():
            require(digest(source/name) == h, 'changed source-screen source')
    else:
        manifest = source_manifest(source)
    root.mkdir(parents=True)
    for name in manifest:
        dest = root/'source'/name
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source/name, dest)
    # Snapshot first, then compile that exact snapshot. Source may be edited by
    # another task without changing what this screen builds or records.
    for name, h in manifest.items():
        require(digest(root/'source'/name) == h, 'source changed during snapshot')
    if prior is not None:
        shutil.copy2(previous/'worker', root/'worker')
        shutil.copy2(previous/'build.log', root/'build.log')
    else:
        with (root/'build.log').open('w') as log:
            subprocess.run(['cargo', 'build', '--release', '--offline', '--locked', '--jobs', '2',
                '--example', 'ic_tournament_worker', '--target-dir', str(root/'build')],
                cwd=root/'source', env=child_env(), stdout=log, stderr=subprocess.STDOUT, check=True)
        from tournament import built_worker
        shutil.copy2(built_worker(root/'build'), root/'worker')
    for name in EVALUATOR:
        (root/'evaluator').mkdir(exist_ok=True)
        shutil.copy2(HERE/name, root/'evaluator'/name)
    cases = []
    rng = random.Random(args.seed)
    for n, a in cells:
        for i in range(args.cases):
            job = {'mode': 'fixture', 'degree': n, 'curve_a': a, 'target_seeds': [rng.getrandbits(64)],
                'algorithm_seed': rng.getrandbits(64), 'factor_base': {'kind':'subgroup_orbits','seed':43,'points':6*n},
                'config': BASE_CONFIG}
            case_id = f'n{n}a{a}-{i:03d}'
            proc = run_native(root/'worker', job, root/'fixture_generation'/case_id, args.timeout)
            write(root/'fixture_generation'/case_id/'process.json', proc, exclusive=True)
            require(proc['exit_code'] == 0, 'fixture generation failed')
            fixture = read(root/'fixture_generation'/case_id/'stdout.json')['fixture']
            c = Curve(fixture)
            require(fixture['target_scalar_constructed'] is False, 'fixture includes planted scalar')
            require(all(c.decode(p) is not None for p in fixture['targets']), 'invalid public targets')
            cases.append({'id': case_id, 'job': job, 'fixture': fixture})
    for arm in arms:
        arm['mode'] = 'ic'
    arms.append(dict(copy.deepcopy(arms[0]), id='aa_control'))
    for walks in (1, 8, 32):
        arms.append({'id': f'rho_w{walks}', 'mode': 'rho', 'config': dict(BASE_CONFIG, rho_parallel_walks=walks),
                     'hypothesis': 'Same-target reference sensitivity; packed signed-Frobenius walk.'})
    pinned = {'worker': digest(root/'worker'), 'build.log': digest(root/'build.log')}
    pinned.update({'evaluator/'+n: digest(root/'evaluator'/n) for n in EVALUATOR})
    pinned.update({str(p.relative_to(root)):digest(p) for p in (root/'fixture_generation').rglob('*') if p.is_file()})
    contract = {'schema_version': 1, 'purpose': 'native development screen only',
        'seed': args.seed, 'arms': arms, 'cases': cases, 'repetitions': args.repetitions,
        'timeout_seconds': args.timeout, 'max_processes': args.max_processes,
        'comparison_kind': args.comparison_kind, 'pinned': pinned, 'source_manifest': manifest,
        'compiler': prior['compiler'] if prior else subprocess.check_output(['rustc', '--version'], text=True).strip(),
        'build_reused_from': {'path':str(previous), 'contract_sha256':objhash(prior)} if prior else None,
        'host': platform.uname()._asdict(), 'promotion_eligible': False,
        'operation_counts': None, 'memory_cap_bytes': None, 'cpu_affinity': None,
        'boundary': 'Complete child, input/targets/base/failed attempts/solve/descent/verification/reporting/exit; independent Python audit outside.',
        'base_floor': 'Coverage ceiling min(1, binomial(B+m-1,m)/(r-1)); log rank requires K independent scalar-field rows. No wall-time conversion.'}
    write(root/'contract.json', contract, exclusive=True)
    write(root/'seal.json', {'sha256': objhash(contract)}, exclusive=True)
    return root


def run(root):
    with (root/'operation.lock').open('a+') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        contract, rows = audit(root, complete=False)
        blocks = [(c, r) for c in contract['cases'] for r in range(contract['repetitions'])]
        rng = random.Random(contract['seed'])
        rng.shuffle(blocks)
        schedule = []
        for case, rep in blocks:
            order = list(contract['arms'])
            rng.shuffle(order)
            schedule.extend((case, arm, rep) for arm in order)
        for case, arm, rep in schedule:
            directory = root/'trials'/case['id']/arm['id']/f'rep-{rep}'
            if directory.exists():
                continue  # audit above checked every existing directory
            job = dict(copy.deepcopy(case['job']), mode=arm['mode'], config=arm['config'])
            p = run_native(root/'worker', job, directory, contract['timeout_seconds'])
            write(directory/'process.json', p, exclusive=True)
            row = {'arm': arm['id'], 'case': case['id'], 'repetition': rep,
                'case_sha256': objhash(case), 'arm_sha256': objhash(arm), 'process': p,
                'status': 'INVALID_OR_INCOMPLETE', 'certificate': None}
            try:
                require(p['status'] != 'TIMEOUT', 'timeout')
                require(p['exit_code'] == 0, 'worker failed or incomplete')
                proof = verify(read(directory/'stdout.json'), case['fixture'], expected_mode=arm['mode'],
                               summands=arm['config']['summands'])
                require(not proof.get('degenerate_descents', 0), 'generic direct relation admitted as IC')
                row.update(status='VERIFIED', certificate=proof)
            except (ValueError, KeyError, TypeError) as exc:
                row.update(reason=str(exc), status='TIMEOUT' if p['status'] == 'TIMEOUT' else 'INVALID_OR_INCOMPLETE')
            row['artifacts'] = {p.name: digest(p) for p in directory.iterdir() if p.is_file()}
            write(directory/'receipt.json', row, exclusive=True)
            rows.append(row)
        write(root/'summary.json', summarize(contract, rows))
        audit(root)
        print(json.dumps({'status': 'AUDITED_DEVELOPMENT_SCREEN', 'trials': len(rows), 'summary': str(root/'summary.json')}))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    p = commands.add_parser('doctor'); p.add_argument('--source-root', type=Path, default=ROOT)
    p = commands.add_parser('propose')
    p.add_argument('--panel', choices=['implementation','algebra','factor-base'], default='implementation')
    p.add_argument('--base-config', type=Path); p.add_argument('--out', type=Path, required=True)
    p = commands.add_parser('recombine')
    p.add_argument('--candidates', type=Path, required=True); p.add_argument('--retain', required=True)
    p.add_argument('--out', type=Path, required=True)
    p = commands.add_parser('screen')
    p.add_argument('--source-root', type=Path, default=ROOT); p.add_argument('--out', type=Path, required=True)
    p.add_argument('--source-screen', type=Path, help='Reuse this screen\'s exact hashed source and binary, ignoring current source edits.')
    p.add_argument('--candidates', type=Path, required=True); p.add_argument('--cells', default='13a0,17a1')
    p.add_argument('--cases', type=int, default=3); p.add_argument('--repetitions', type=int, default=3)
    p.add_argument('--seed', type=int, required=True); p.add_argument('--timeout', type=float, default=15)
    p.add_argument('--max-processes', type=int, default=360)
    p.add_argument('--comparison-kind', choices=['fixed-support','factor-base-policy'], default='fixed-support')
    for cmd in ('run', 'verify'):
        p = commands.add_parser(cmd); p.add_argument('--round', type=Path, required=True)
    args = parser.parse_args()
    if args.command == 'doctor':
        print(json.dumps(doctor(args.source_root.resolve()), indent=2))
    elif args.command == 'propose':
        write(args.out, factorial_candidates(read(args.base_config) if args.base_config else BASE_CONFIG,
                                            panel=args.panel), exclusive=True)
    elif args.command == 'recombine':
        write(args.out, recombine(read(args.candidates), args.retain.split(',')), exclusive=True)
    elif args.command == 'screen':
        root = prepare(args)
        subprocess.run([sys.executable, str(root/'evaluator/autolab.py'), 'run', '--round', str(root)], check=True)
    elif args.command == 'run':
        run(args.round.resolve())
    else:
        _, rows = audit(args.round.resolve())
        print(json.dumps({'status': 'VERIFIED', 'trial_receipts': len(rows), 'promotion_eligible': False}))


if __name__ == '__main__':
    try:
        main()
    except (InvalidEvidence, OSError, ValueError, KeyError, subprocess.SubprocessError) as exc:
        sys.exit(str(exc))
