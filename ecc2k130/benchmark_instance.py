#!/usr/bin/env python3
"""Measure concurrent ECC2K-130 clients over one setup-inclusive host window.

This does not publish results, launch cloud instances, or change CUDA defaults.
Run IDs and physical GPU UUIDs are distinct within every requested experiment.
"""
import argparse
import csv
from dataclasses import asdict, dataclass
import hashlib
import json
import math
import os
from pathlib import Path
import re
import signal
import statistics
import subprocess
import time

from benchmark_hardware import platform_identity

COMPILERS = {'nvcc', 'ptxas', 'cicc', 'cc1plus', 'nvvm', 'fatbinary', 'cudafe++'}


def digest(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def save(path, value):
    with Path(path).open('x') as stream:
        json.dump(value, stream, indent=2)
        stream.write('\n')


def inventory():
    fields = 'index,uuid,name,driver_version,memory.total,power.limit,compute_cap'
    result = subprocess.run(['nvidia-smi', '--query-gpu=' + fields,
                             '--format=csv,noheader,nounits'], capture_output=True,
                            text=True, check=True, timeout=20)
    devices = []
    for values in csv.reader(result.stdout.splitlines()):
        if len(values) != 7:
            raise ValueError('Unexpected GPU inventory row')
        index, uuid, name, driver, memory, power, capability = [v.strip() for v in values]
        if not re.fullmatch(r'GPU-[0-9a-fA-F-]+', uuid):
            raise ValueError('This benchmark requires physical GPU UUIDs; MIG is not supported')
        devices.append(dict(index=int(index), uuid=uuid, name=name, driver=driver,
                            memory_mib=int(memory), power_limit_watts=float(power),
                            compute_capability=capability))
    if not devices:
        raise ValueError('No physical GPUs found')
    return devices


def select_devices(devices, selection):
    selectors = [str(d['index']) for d in devices] if selection == 'all' else selection.split(',')
    selected = []
    for selector in selectors:
        matches = [d for d in devices if selector.strip() in {str(d['index']), d['uuid']}]
        if len(matches) != 1:
            raise ValueError('Unknown or ambiguous physical GPU: ' + selector)
        if any(d['uuid'] == matches[0]['uuid'] for d in selected):
            raise ValueError('Each physical GPU may appear only once')
        selected.append(matches[0])
    if not selected:
        raise ValueError('Select at least one GPU')
    return selected


@dataclass(frozen=True)
class Work:
    workers: int = 524288
    batch: int = 16
    steps: int = 1024
    launches: int = 32
    collection: bool = False
    verify: int = 0
    dp_weight: int = 34

    def check(self):
        for value in (self.workers, self.batch, self.steps, self.launches):
            if not 0 < value < 2**31:
                raise ValueError('Worker, batch, step and launch counts must be positive int32 values')
        # eccSeedFor stores only the low 32 bits of each walk index.
        if self.workers * self.batch > 2**32:
            raise ValueError('Walk population exceeds the distinct 32-bit seed namespace')
        if self.iterations >= 2**63 or self.verify < 0 or not 0 <= self.dp_weight <= 131:
            raise ValueError('Invalid work or verification budget')

    @property
    def iterations(self):
        return self.workers * self.batch * self.steps * self.launches


def run_ids(base, gpu_count, repetitions):
    if gpu_count < 1 or repetitions < 1 or base < 0:
        raise ValueError('Invalid run-ID allocation')
    # First row is the independent CPU-replay admission check.
    count = gpu_count * (repetitions + 1)
    if base + count > 65536:
        raise ValueError('Distinct 16-bit run IDs would overflow')
    return [list(range(base + r * gpu_count, base + (r + 1) * gpu_count))
            for r in range(repetitions + 1)]


def parse_client(text, work):
    geometry = re.findall(r'backend cuda-packed131: (\d+) threads x (\d+) slots x (\d+) lanes = (\d+) walks, dp weight (\d+), (\d+) steps per launch', text)
    wanted_weight = work.dp_weight if work.collection else 0
    wanted = (work.workers, work.batch, 1, work.workers * work.batch, wanted_weight, work.steps)
    if len(geometry) != 1 or tuple(map(int, geometry[0])) != wanted:
        raise ValueError('Packed backend geometry or collection mode does not match the requested work')
    progress = re.findall(r'\b(\d+) iterations\b', text)
    finished = re.findall(r'finished: ([0-9.]+) M it/s, (\d+) distinguished points \((\d+) verified against the reference, (\d+) dropped\)', text)
    if not progress or int(progress[-1]) != work.iterations or len(finished) != 1:
        raise ValueError('Client did not report its full requested iteration budget')
    rate, points, verified, dropped = finished[0]
    rate = float(rate)
    if 'MISMATCH' in text or int(dropped) or int(verified) != min(work.verify, int(points)):
        raise ValueError('Client verification or dropped-record check failed')
    if not math.isfinite(rate) or rate <= 0 or (work.verify and int(points) < work.verify):
        raise ValueError('Missing rate or insufficient independent CPU replays')
    return dict(iterations=work.iterations, client_million_iterations_per_second=rate,
                distinguished_points=int(points), cpu_replays=int(verified), dropped=int(dropped))


def gpu_processes(_owned=None):
    result = subprocess.run(['nvidia-smi', '--query-compute-apps=gpu_uuid,pid,process_name',
                             '--format=csv,noheader'], capture_output=True, text=True,
                            check=True, timeout=15)
    rows = []
    for values in csv.reader(result.stdout.splitlines()):
        if len(values) != 3:
            raise ValueError('Unexpected GPU process row')
        uuid, pid, name = [v.strip() for v in values]
        rows.append(dict(uuid=uuid, pid=int(pid), name=name))
    return rows


def exited_unreaped(process):
    return os.waitid(os.P_PID, process.pid, os.WEXITED | os.WNOHANG | os.WNOWAIT) is not None


def stop_owned(processes):
    for process in processes:
        if process.returncode is None and not exited_unreaped(process):
            try:
                os.killpg(process.pid, signal.SIGTERM)
            except ProcessLookupError:
                pass
    deadline = time.monotonic() + 2
    while time.monotonic() < deadline and any(p.returncode is None and not exited_unreaped(p) for p in processes):
        time.sleep(0.02)
    for process in processes:
        if process.returncode is None and not exited_unreaped(process):
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
        process.wait()


def run_cohort(binary, devices, ids, work, output, timeout=900, omp_threads=1,
               observer=gpu_processes, poll_seconds=0.2):
    """Own all child PIDs until observations finish, including exited zombies."""
    work.check()
    if len(devices) != len(ids) or len(set(ids)) != len(ids) or len({d['uuid'] for d in devices}) != len(devices):
        raise ValueError('Device and run-ID assignments must be one-to-one')
    if not devices or any(not 0 <= value <= 65535 for value in ids):
        raise ValueError('Invalid device or run-ID assignment')
    if not math.isfinite(timeout) or timeout <= 0 or omp_threads < 1 or poll_seconds <= 0:
        raise ValueError('A finite positive timeout and polling interval are required')
    binary = Path(binary).resolve()
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    binary_sha = digest(binary)
    selected = {d['uuid'] for d in devices}
    save(output / 'plan.json', dict(work=asdict(work), devices=devices, run_ids=ids,
                                   binary=str(binary), binary_sha256=binary_sha,
                                   omp_threads_per_worker=omp_threads, timeout_seconds=timeout))
    try:
        preflight_processes = observer({})
    except Exception as exc:
        save(output / 'preflight.json', dict(passed=False, observation_error=repr(exc)))
        raise RuntimeError('GPU preflight observation failed; no client was launched') from exc
    busy = [p for p in preflight_processes if p['uuid'] in selected]
    save(output / 'preflight.json', dict(passed=not busy, gpu_processes=preflight_processes,
                                        external_gpu_processes=busy))
    if busy:
        raise RuntimeError('A selected GPU is already in use; no client was launched')
    children, streams, metadata, errors, observations, seen = [], [], [], [], [], set()
    start_ns = time.monotonic_ns()
    timed_out = False
    try:
        for index, (device, run_id) in enumerate(zip(devices, ids)):
            log = output / f'gpu{index}.log'
            command = [str(binary), '--packed', '--device', '0', '--threads', str(work.workers),
                       '--steps', str(work.steps), '--launches', str(work.launches),
                       '--run-id', str(run_id), '--verify', str(work.verify)]
            if work.collection:
                command += ['--dp-weight', str(work.dp_weight), '--dp-file', str(output.resolve() / f'gpu{index}.dp')]
            else:
                command += ['--bench']
            env = dict(os.environ, CUDA_VISIBLE_DEVICES=device['uuid'],
                       OMP_NUM_THREADS=str(omp_threads), CUDA_DISABLE_PTX_JIT='1')
            stream = log.open('x')
            streams.append(stream)
            process = subprocess.Popen(command, stdout=stream, stderr=subprocess.STDOUT,
                                       env=env, start_new_session=True)
            children.append(process)
            metadata.append(dict(pid=process.pid, device=device, run_id=run_id,
                                 command=command, log=log.name))
        save(output / 'owned-processes.json', metadata)
        owned = {r['pid']: r['device']['uuid'] for r in metadata}
        while True:
            row = dict(monotonic_ns=time.monotonic_ns(), owned_pids=list(owned))
            try:
                processes = [p for p in observer(owned) if p['uuid'] in selected]
                row['gpu_processes'] = processes
                outside = [p for p in processes if owned.get(p['pid']) != p['uuid'] or p['name'] not in {str(binary), '[No data]'}]
                if outside:
                    errors.append(dict(external_gpu_processes=outside))
                seen.update(p['pid'] for p in processes if owned.get(p['pid']) == p['uuid'])
            except Exception as exc:
                row['observation_error'] = repr(exc)
                errors.append(dict(observation_error=repr(exc)))
            observations.append(row)
            if all(exited_unreaped(p) for p in children):
                break
            if (time.monotonic_ns() - start_ns) / 1e9 >= timeout:
                timed_out = True
                break
            time.sleep(poll_seconds)
    except BaseException as exc:
        errors.append(dict(launch_or_control_error=repr(exc)))
    finally:
        stop_owned(children)
        end_ns = time.monotonic_ns()
        for stream in streams:
            stream.close()
    save(output / 'gpu-observations.json', observations)
    results = []
    for process, row in zip(children, metadata):
        row.update(returncode=process.returncode, log_sha256=digest(output / row['log']))
        try:
            if process.returncode != 0:
                raise ValueError('Client exited unsuccessfully')
            row.update(parse_client((output / row['log']).read_text(), work))
            if work.collection:
                dp = output / f'gpu{len(results)}.dp'
                if dp.stat().st_size != row['distinguished_points'] * 32:
                    raise ValueError('DP file size does not match completed records')
                row.update(dp_file=dp.name, dp_sha256=digest(dp))
            row['passed'] = True
        except (ValueError, OSError) as exc:
            row.update(passed=False, error=str(exc))
        results.append(row)
    passed = (len(results) == len(devices) and all(r['passed'] for r in results)
              and not errors and not timed_out and seen == {p.pid for p in children}
              and digest(binary) == binary_sha)
    total = sum(r['iterations'] for r in results) if passed else None
    elapsed = (end_ns - start_ns) / 1e9
    result = dict(passed=passed, measured_gpu_count=len(devices), elapsed_seconds=elapsed,
                  monotonic_start_ns=start_ns, monotonic_end_ns=end_ns,
                  completed_scalar_iterations=total,
                  cohort_iterations_per_second=total / elapsed if passed else None,
                  timing_scope='First process launch through every process exit and observation/cleanup; includes all client setup. No per-device rates are summed.',
                  observation_scope='Sampled GPU processes; all child PIDs remain owned and unreaped until observations finish. Activity between queries is not excluded.',
                  timed_out=timed_out, observation_errors=errors,
                  observed_worker_pids=sorted(seen), binary_sha256=binary_sha, workers=results)
    save(output / 'result.json', result)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', type=Path, default=Path(__file__).parent / 'build/ecc2k130-local-packed')
    parser.add_argument('--devices', required=True, help='Physical nvidia-smi indices/UUIDs separated by commas, or all')
    parser.add_argument('--out', type=Path, required=True, help='New evidence directory; existing paths are rejected')
    parser.add_argument('--workers', type=int, default=524288)
    parser.add_argument('--batch', type=int, default=16, help='Must match the compiled packed backend')
    parser.add_argument('--steps', type=int, default=1024)
    parser.add_argument('--launches', type=int, default=32)
    parser.add_argument('--repetitions', type=int, default=3)
    parser.add_argument('--run-id-base', type=int, default=52000)
    parser.add_argument('--collect-dp34', action='store_true')
    parser.add_argument('--timeout', type=float, default=900)
    args = parser.parse_args()
    devices = select_devices(inventory(), args.devices)
    ids = run_ids(args.run_id_base, len(devices), args.repetitions)
    work = Work(args.workers, args.batch, args.steps, args.launches, args.collect_dp34)
    work.check()
    if not math.isfinite(args.timeout) or args.timeout <= 0:
        parser.error('timeout must be finite and positive')
    compilers = set(subprocess.check_output(['ps', '-eo', 'comm='], text=True).split()) & COMPILERS
    if compilers:
        raise RuntimeError('Finish active compilers before benchmarking: ' + ', '.join(sorted(compilers)))
    args.out.mkdir(parents=True, exist_ok=False)
    affinity = len(os.sched_getaffinity(0))
    omp = max(1, min(8, affinity // len(devices)))
    identity = dict(platform_identity(), devices=devices, cpu_affinity_count=affinity,
                    measurement_scope='concurrent_process_cohort', measured_gpu_count=len(devices))
    save(args.out / 'experiment.json', dict(hardware=identity, work=asdict(work),
                                          source_sha256=digest(Path(__file__)), binary_sha256=digest(args.binary),
                                          run_ids=ids, repetitions=args.repetitions,
                                          target_iterations_per_second=40e9, full_dlp_s=None,
                                          generic_work_ratio=1))
    rows, admission = [], None
    try:
        admission = run_cohort(args.binary, devices, ids[0], Work(257, args.batch, 64, 2, True, 16, 52),
                               args.out / 'admission', args.timeout, omp)
        if not admission['passed']:
            raise RuntimeError('Independent CPU-replay admission failed; see saved evidence')
        for rep in range(args.repetitions):
            result = run_cohort(args.binary, devices, ids[rep + 1], work,
                                args.out / f'repetition{rep}', args.timeout, omp)
            rows.append(result)
            if not result['passed']:
                raise RuntimeError('Cohort failed; no samples will be removed or replaced')
            print(f"repetition {rep}: {result['cohort_iterations_per_second']/1e9:.6f} B complete iterations/s across {len(devices)} GPU(s), setup included", flush=True)
    except BaseException as exc:
        save(args.out / 'summary.json', dict(passed=False, target_met=False, hardware=identity,
                                            admission=admission, rows=rows, error=repr(exc),
                                            full_dlp_s=None, generic_work_ratio=1))
        raise
    rates = [r['cohort_iterations_per_second'] for r in rows]
    save(args.out / 'summary.json', dict(passed=True, hardware=identity, rows=rows, admission=admission,
                                        median_cohort_iterations_per_second=statistics.median(rates),
                                        minimum_cohort_iterations_per_second=min(rates),
                                        target_met=len(rows) >= 3 and min(rates) >= 40e9,
                                        target_rule='At least three complete cohorts, every cohort at least 40 B/s. This is an absolute throughput observation, not a baseline/candidate speedup claim.',
                                        full_dlp_s=None, generic_work_ratio=1))


if __name__ == '__main__':
    main()
