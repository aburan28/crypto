"""Run the immutable full-corpus WDSat regression contract. Linux/Python 3.10+."""
import argparse
from datetime import datetime, timezone
import hashlib
import itertools
import json
import math
import os
from pathlib import Path
import platform
import random
import re
import resource
import shutil
import statistics
import subprocess
import sys
import time
from urllib.request import urlopen

HERE = Path(__file__).resolve().parent
PILOT = HERE.parent / 'pilot'
sys.path.insert(0, str(PILOT))
from run_pilot import GF, anf_check, check_sat, little


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fetch_inputs(cache, offline):
    for item in json.loads((HERE / 'corpus_manifest.json').read_text()):
        path = cache / item['path']
        if not path.exists() and not offline:
            with urlopen(item['url'], timeout=30) as response:
                data = response.read()
            if hashlib.sha256(data).hexdigest() != item['sha256']:
                raise ValueError(f'download hash mismatch: {path}')
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(data)
        if digest(path) != item['sha256']:
            raise ValueError(f'input hash mismatch: {path}')


def child_limits(memory):
    resource.setrlimit(resource.RLIMIT_AS, (memory, memory))
    resource.setrlimit(resource.RLIMIT_CORE, (0, 0))


def percentile90(values):
    return sorted(values)[math.ceil(.9 * len(values)) - 1]


def aggregate(records, contract):
    cells = []
    for variant in contract['variants']:
        for case in contract['cases']:
            rs = [r for r in records if r['variant'] == variant['id'] and r['case'] == case]
            valid = len(rs) == contract['repetitions'] and all(r['verified'] for r in rs)
            counters = [r['conflicts'] for r in rs if r['conflicts'] is not None]
            stable = valid and len({(r['status'], r['conflicts']) for r in rs}) == 1
            cells.append({'variant': variant['id'], 'case': case, 'completed': len(rs),
                          'verified': valid, 'repeatable': stable,
                          'status': rs[0]['status'] if stable else 'INCOMPLETE_OR_UNSTABLE',
                          'median_conflicts': statistics.median(counters) if len(counters) == len(rs) and rs else None,
                          'median_process_wall_s': statistics.median(r['process_wall_s'] for r in rs) if rs else None,
                          'p90_process_wall_s': percentile90([r['process_wall_s'] for r in rs]) if rs else None})
    groups = []
    for variant in contract['variants']:
        for n, l in [(15, 5), (17, 6), (19, 6)]:
            cs = [c for c in cells if c['variant'] == variant['id'] and c['case'].startswith(f'n{n}l{l}-')]
            complete = all(c['verified'] and c['repeatable'] for c in cs)
            groups.append({'variant': variant['id'], 'n': n, 'l': l, 'cases': len(cs),
                           'verified_cases': sum(c['verified'] for c in cs),
                           'sat_cases': sum(c['status'] == 'SAT' for c in cs),
                           'unsat_cases': sum(c['status'] == 'UNSAT' for c in cs),
                           'algebraic_only_cases': sum(c['status'] == 'SAT_ALGEBRAIC_ONLY' for c in cs),
                           'sum_median_conflicts': sum(c['median_conflicts'] for c in cs) if complete else None,
                           'sum_median_process_wall_s': sum(c['median_process_wall_s'] for c in cs) if complete else None,
                           'full_dlp_S': None, 'cost_over_rho': None, 'cost_over_floor': None})
    return cells, groups


def run(args):
    contract = json.loads((HERE / 'contract.json').read_text())
    if args.output.exists():
        raise ValueError('output must not already exist')
    fetch_inputs(args.cache_dir, args.offline)
    solver = args.solver.resolve()
    solver_hash = digest(solver)
    dest = args.output.resolve()
    dest.mkdir(parents=True, exist_ok=False)
    (dest / 'normalized').mkdir()
    cpu = next((s.split(':', 1)[1].strip() for s in Path('/proc/cpuinfo').read_text().splitlines()
                if s.startswith('model name')), 'unknown')
    metadata = {'suite_id': contract['id'], 'classification': 'accounting',
                'started_at': datetime.now(timezone.utc).isoformat(),
                'contract_sha256': digest(HERE / 'contract.json'),
                'corpus_sha256': digest(HERE / 'corpus_manifest.json'),
                'runner_sha256': digest(Path(__file__)),
                'checker_sha256': digest(PILOT / 'run_pilot.py'),
                'solver_sha256': solver_hash, 'source_revision': args.source_revision,
                'counter_unit': contract['solver_counter_unit'],
                'cpu': cpu, 'platform': platform.platform(), 'python': platform.python_version(),
                'cpu_affinity': sorted(os.sched_getaffinity(0)),
                'load_average_at_start': os.getloadavg(),
                'command': sys.argv, 'peak_rss_bytes': None,
                'memory_measurement': 'Address-space cap enforced per child; peak RSS not measured.',
                'full_dlp_S': None, 'cost_over_rho': None, 'cost_over_floor': None}
    for name in ['build.json', 'build.stdout', 'build.stderr']:
        src = solver.parent / name
        if src.exists():
            shutil.copyfile(src, dest / name)
    if (dest / 'build.json').exists():
        build = json.loads((dest / 'build.json').read_text())
        if build['binary_sha256'] != solver_hash or build['exit_code'] != 0:
            raise ValueError('build metadata does not describe this binary')
    (dest / 'metadata.json').write_text(json.dumps(metadata, indent=2) + '\n')
    shutil.copyfile(HERE / 'contract.json', dest / 'contract.json')
    shutil.copyfile(HERE / 'corpus_manifest.json', dest / 'corpus_manifest.json')
    data = {}
    for case in contract['cases']:
        info = (args.cache_dir / 'instances' / ('INFO' + case + '.dimacs')).read_text().splitlines()
        n, l = map(int, info[0].split())
        raw = (args.cache_dir / 'instances' / ('X' + case + '.anf')).read_text()
        normalized = '\n'.join(s for s in raw.splitlines() if s.strip()) + '\n'
        path = dest / 'normalized' / ('X' + case + '.anf')
        path.write_text(normalized)
        data[case] = (n, l, little(info[1]), little(info[2]), normalized, path)
    records, unsat_proofs = [], {}
    rng = random.Random(contract['shuffle_seed'])
    started = time.perf_counter()
    with (dest / 'raw.jsonl').open('x') as stream:
        for repetition in range(contract['repetitions']):
            cases = list(contract['cases'])
            rng.shuffle(cases)
            for case in cases:
                variants = list(contract['variants'])
                rng.shuffle(variants)
                n, l, modulus, target, raw, inp = data[case]
                for variant in variants:
                    field = GF(n, modulus)
                    cmd = [str(solver), '-i', str(inp), '-n', str(n), '-l', str(l), '-m', str(contract['m']), *variant['args']]
                    row = {'case': case, 'variant': variant['id'], 'repetition': repetition,
                           'command': cmd, 'normalized_sha256': digest(inp),
                           'status': 'ERROR', 'verified': False, 'conflicts': None,
                           'point_relation_verified': False,
                           'field_modulus': modulus, 'target_x': target,
                           'sat_verification_s': 0.0, 'unsat_oracle_s': 0.0}
                    before = resource.getrusage(resource.RUSAGE_CHILDREN)
                    start = time.perf_counter()
                    try:
                        process = subprocess.run(cmd, capture_output=True, text=True,
                                                 timeout=contract['timeout_s'],
                                                 preexec_fn=lambda: child_limits(contract['address_space_limit_bytes']))
                        row.update(returncode=process.returncode, stdout=process.stdout, stderr=process.stderr)
                    except subprocess.TimeoutExpired as exc:
                        row.update(status='TIMEOUT', returncode=None,
                                   stdout=(exc.stdout or b'').decode(errors='replace'),
                                   stderr=(exc.stderr or b'').decode(errors='replace'))
                    row['process_wall_s'] = time.perf_counter() - start
                    after = resource.getrusage(resource.RUSAGE_CHILDREN)
                    row['child_cpu_s'] = after.ru_utime + after.ru_stime - before.ru_utime - before.ru_stime
                    if row['returncode'] == 0:
                        try:
                            nv = int(raw.splitlines()[0].split()[2])
                            bits = next((s for s in row['stdout'].splitlines() if len(s) == nv and set(s) <= {'0', '1'}), None)
                            row['conflicts'] = int(row['stdout'].strip().splitlines()[-1])
                            if row['conflicts'] < 0:
                                raise ValueError('negative conflict count')
                            if bits is not None:
                                t = time.perf_counter()
                                if case in contract['algebraic_only_controls']:
                                    anf_check(raw, bits)
                                    xs = [little(bits[i*l:(i+1)*l]) for i in range(3)]
                                    if field.s4(*xs, target) != 0 or field.lift(target):
                                        raise ValueError('algebraic-only control failed certification')
                                    row.update(status='SAT_ALGEBRAIC_ONLY', verified=True, assignment=bits,
                                               x_witness=xs, certificate={'reason': 'target_x_does_not_lift',
                                               'target_lift_count': 0, 'summand_lift_counts': [len(field.lift(x)) for x in xs]})
                                else:
                                    xs, points = check_sat(field, raw, bits, l, target)
                                    row.update(status='SAT', verified=True, point_relation_verified=True,
                                               assignment=bits, x_witness=xs, point_witness=points)
                                row['sat_verification_s'] = time.perf_counter()-t
                            elif re.search(r'^UNSAT$', row['stdout'], re.M):
                                if case not in unsat_proofs:
                                    t = time.perf_counter(); count = 0
                                    for xs in itertools.combinations_with_replacement(range(1 << l), 3):
                                        count += 1
                                        if field.s4(*xs, target) == 0:
                                            raise ValueError(f'false UNSAT; S4 root {xs}')
                                    elapsed = time.perf_counter() - t
                                    unsat_proofs[case] = {'checks': count, 'wall_s': elapsed,
                                                         'field_modulus': modulus, 'target_x': target}
                                    row['unsat_oracle_s'] = elapsed
                                row.update(status='UNSAT', verified=True, unsat_proof=case)
                            else:
                                raise ValueError('unrecognized solver output')
                        except (AssertionError, ValueError, IndexError, ZeroDivisionError) as exc:
                            row.update(status='VERIFICATION_ERROR', verified=False, error=repr(exc))
                    records.append(row)
                    stream.write(json.dumps(row) + '\n'); stream.flush()
                    if len(records) % 20 == 0 or not row['verified']:
                        print(json.dumps({'finished': len(records), 'case': case, 'variant': variant['id'],
                                          'status': row['status'], 'elapsed_s': round(time.perf_counter()-started, 2)}), flush=True)
    cells, groups = aggregate(records, contract)
    complete = all(c['verified'] and c['repeatable'] for c in cells)
    complete = complete and all(len({r['status'] for r in records if r['case'] == case}) == 1
                                for case in contract['cases'])
    summary = {'status': 'complete' if complete else 'incomplete_or_unstable', 'metadata': metadata,
               'processes': len(records), 'verified_processes': sum(r['verified'] for r in records),
               'cases': len(contract['cases']), 'variants': len(contract['variants']),
               'repetitions': contract['repetitions'], 'cells': cells, 'groups': groups,
               'unsat_proofs': unsat_proofs, 'exhaustive_s4_checks': sum(p['checks'] for p in unsat_proofs.values()),
               'raw_sha256': digest(dest / 'raw.jsonl'),
               'total_loop_wall_s': time.perf_counter()-started,
               'sum_process_wall_s': sum(r['process_wall_s'] for r in records),
               'sum_sat_verification_s': sum(r['sat_verification_s'] for r in records),
               'sum_unsat_oracle_s': sum(r['unsat_oracle_s'] for r in records),
               'finished_at': datetime.now(timezone.utc).isoformat(),
               'limitations': contract['scope']}
    (dest / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    # Normalized third-party inputs are reproducible in the external cache; do not commit them.
    print(json.dumps({'status': summary['status'], 'processes': len(records),
                      'verified': summary['verified_processes'], 'exhaustive_checks': summary['exhaustive_s4_checks']}), flush=True)
    return 0 if complete else 1


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cache-dir', type=Path, required=True)
    parser.add_argument('--solver', type=Path, required=True)
    parser.add_argument('--source-revision', required=True, help='Actual solver commit/revision; retained in results.')
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--offline', action='store_true')
    args = parser.parse_args()
    try:
        raise SystemExit(run(args))
    except (ValueError, OSError) as exc:
        parser.error(str(exc))


if __name__ == '__main__':
    main()
