#!/usr/bin/env python3
"""Fixed archived toy vectors for candidate wiring; no performance selection."""
import argparse
import copy
import json
import os
from pathlib import Path
import platform
import subprocess
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from autolab import run_native
from driver_admission import make_admission, producer_metadata, run_record
from identity import sha256
from oracle import InvalidEvidence, require
from tournament import digest, execute, parse_profiles, read, write


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prepared', type=Path, required=True)
    parser.add_argument('--worker', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--profile', action='store_true')
    args = parser.parse_args()
    prepared, worker, out = args.prepared.resolve(), args.worker.resolve(), args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    panel = read(HERE.parent/'goal_20260924/improvement/round1.json')
    manifest = read(prepared/'source-manifest.json')
    require(sha256(manifest) == panel['candidate_source_sha256'], 'unregistered source')
    for name, value in manifest.items():
        require(digest(prepared/'source'/name) == value, 'changed source '+name)
    metadata = producer_metadata(prepared/'source', manifest,
        subprocess.check_output(['rustc', '--version'], text=True).strip())
    if args.profile:
        require(platform.system() == 'Linux' and platform.machine() == 'x86_64' and
                subprocess.check_output(['valgrind', '--version'], text=True).strip() == 'valgrind-3.22.0',
                'profile controls require calibrated Linux amd64')
    cpu = sorted(os.sched_getaffinity(0))[-1] if args.profile else None
    resources = dict(cpu=cpu, worker_threads=1, memory_bytes=8*1024**3 if args.profile else None,
                     timeout_seconds='180')
    host_id = sha256(platform.uname()._asdict())
    outcomes = []
    for vector in ('n13-public.json', 'n23.json'):
        stored = read(HERE/'testdata'/vector)
        stored = stored['ic'] if 'ic' in stored else stored
        fixture = stored['report']['fixture']
        for arm in panel['candidates'][1:]:
            directory = out/vector.removesuffix('.json')/arm['id']
            directory.mkdir(parents=True)
            job = dict(copy.deepcopy(stored['job']), mode='ic', config=arm['config'],
                       public_targets=fixture['targets'])
            write(directory/'job.json', job, exclusive=True)
            def run(command, selected, target):
                if args.profile:
                    return execute(command, selected, target, 180, 8*1024**3, cpu)
                require(len(command) == 1, 'native control cannot profile')
                return run_native(Path(command[0]), selected, target, 180)
            result = dict(vector=vector, arm=arm['id'], status='FAILED')
            try:
                inventory = run([str(worker)], dict(job, mode='inventory'), directory/'inventory')
                write(directory/'inventory/process.json', inventory, exclusive=True)
                require(inventory['exit_code'] == 0, 'inventory failed')
                kwargs = dict(job=job, fixture=fixture, report=read(directory/'inventory/stdout.json'),
                    manifest=manifest, metadata=metadata, resources=resources, worker_sha256=digest(worker))
                admitted = make_admission(**kwargs)
                write(directory/'admission.json', admitted, exclusive=True)
                native = run([str(worker)], job, directory/'native')
                write(directory/'native/process.json', native, exclusive=True)
                require(native['exit_code'] == 0, 'native control failed')
                report = read(directory/'native/stdout.json')
                profile = costs = None
                if args.profile:
                    command = ['valgrind', '--tool=callgrind', '--cache-sim=no', '--branch-sim=no',
                        '--separate-threads=no', '--collect-atstart=yes', '--instr-atstart=yes',
                        '--callgrind-out-file='+str(directory/'profile/callgrind.out'), str(worker)]
                    process = run(command, job, directory/'profile')
                    write(directory/'profile/process.json', process, exclusive=True)
                    require(process['exit_code'] == 0, 'profile control failed')
                    profile = read(directory/'profile/stdout.json')
                    costs = parse_profiles(directory/'profile', phase_schema=3)
                record = run_record(admitted, number=0, host_id=host_id, status='complete',
                    native=report, process_wall_ns=native['process_wall_ns'], profile=profile, costs=costs)
                write(directory/'run.json', record, exclusive=True)
                # A declared different row policy must fail before gaining a
                # second candidate identity from the same executed inventory.
                changed = copy.deepcopy(kwargs)
                changed['job']['config']['row_kernel'] = 'word' if arm['config']['row_kernel'] != 'word' else 'full'
                try:
                    make_admission(**changed)
                except InvalidEvidence:
                    pass
                else:
                    raise InvalidEvidence('ignored row policy was admitted')
                result.update(status='VERIFIED', candidate_id=record['candidate_id'])
            except Exception as exc:
                result['reason'] = f'{type(exc).__name__}: {exc}'
            write(directory/'receipt.json', result, exclusive=True)
            outcomes.append(result)
            print(json.dumps(result), flush=True)
    summary = dict(scope='fixed archived toy vectors; wiring and correctness only',
        promotion_eligible=False, scheduled=len(outcomes),
        verified=sum(row['status'] == 'VERIFIED' for row in outcomes), outcomes=outcomes)
    write(out/'summary.json', summary, exclusive=True)
    require(summary['verified'] == summary['scheduled'] == 30, 'retained candidate control failures')


if __name__ == '__main__':
    main()
