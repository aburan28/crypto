#!/usr/bin/env python3
"""Exclusive full-pipeline Ir counts; archive raw profiles and replay checks."""
import argparse
import lzma
import hashlib
import itertools
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import time
from run import make_workloads

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--binary', type=Path, required=True)
    ap.add_argument('--valgrind', type=Path, required=True)
    ap.add_argument('--valgrind-lib', type=Path, required=True)
    ap.add_argument('--revision', required=True)
    ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--confirm', action='store_true')
    args = ap.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    source = HERE.parent/'weil_factor_composition_20260914/contract.json'
    contract = json.loads(source.read_text())
    workloads = make_workloads(contract['stage_cases'], args.confirm)
    cpu = min(os.sched_getaffinity(0))
    os.sched_setaffinity(0, {cpu})
    env = {k:v for k,v in os.environ.items() if not k.startswith('IC_')}
    env.update(VALGRIND_LIB=str(args.valgrind_lib.resolve()), RAYON_NUM_THREADS='1',
               IC_PREPROCESS_CACHE='off', IC_REDUCTION_CACHE='off')
    metadata = dict(revision=args.revision, binary_sha256=hashlib.sha256(args.binary.read_bytes()).hexdigest(),
                    harness_sha256=hashlib.sha256((ROOT/'examples/rho_parity_cost.rs').read_bytes()).hexdigest(),
                    profiler=subprocess.check_output([str(args.valgrind), '--version'],env=env,text=True).strip(),
                    cpu_affinity=[cpu], workloads=workloads, confirmation=args.confirm, instruction_unit='Callgrind Ir',
                    timeout_seconds=60, repetitions=1,
                    limitation='user-mode instruction count; kernel execution and memory latency unpriced')
    (args.output/'metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
    (args.output/'protocol.md').write_bytes((HERE/'COST_PROTOCOL.md').read_bytes())
    summary, count = [], 0
    # Publish the archive only after its compressor has closed. A partial
    # profile must never look like an accepted, complete evidence file.
    with tempfile.TemporaryDirectory() as tmp, lzma.open(Path(tmp)/'records.jsonl.xz','wt') as archive:
        def run(mode, case, variant, seed, log=5):
            nonlocal count
            count += 1
            profile = Path(tmp)/'profile'
            for p in Path(tmp).glob('profile*'):
                p.unlink()
            command = [str(args.valgrind.resolve()), '--tool=callgrind', '--instr-atstart=no',
                       '--fn-skip=*', '--separate-threads=no', '--cache-sim=no', '--branch-sim=no',
                       '--callgrind-out-file='+str(profile), str(args.binary.resolve()), mode,
                       str(case), variant, str(seed), str(log)]
            started = time.monotonic()
            try:
                r = subprocess.run(command, env=env, capture_output=True, timeout=60)
                status = 'finished' if r.returncode == 0 else 'error'
                stdout, stderr = r.stdout.decode(), r.stderr.decode()
            except subprocess.TimeoutExpired as e:
                status = 'timeout'
                stdout, stderr = (e.stdout or b'').decode(), (e.stderr or b'').decode()
            profiles = {p.name:p.read_text() for p in sorted(Path(tmp).glob('profile*'))}
            scopes = {}
            for name, text in profiles.items():
                totals = re.findall(r'^totals: (\d+)',text,re.M)
                trigger = re.findall(r'^desc: Trigger: Client Request: (.+)',text,re.M)
                if trigger:
                    assert len(totals) == 1 and trigger[0] not in scopes
                    scopes[trigger[0]] = int(totals[0])
                elif totals:
                    assert int(totals[0]) == 0, 'unscoped instructions'
            rows = [json.loads(x) for x in stdout.splitlines() if x.startswith('{')]
            terminal = next((x for x in rows if x['phase'] == mode),None)
            checks = [a for x in rows for a in x.get('attempts',[])]
            for a in checks:
                assert a['correct']
            total = sum(scopes.values())
            collected = re.search(r'Collected\s*:\s*(\d+)',stderr)
            if status == 'finished':
                assert collected and int(collected[1]) == total
                expected = {'curve_setup','target_generation','rho_and_final_verification'} if mode == 'rho' else (
                    {'curve_setup','calibration_loop'} if mode.startswith('calibrate') else
                    {'curve_setup','target_generation','factor_base_plan_index','driver_and_final_verification'})
                assert set(scopes) == expected, scopes
                assert total > 0
            result = dict(mode=mode, case=case, variant=variant, seed=seed, known_log=log,
                          status=status, verified=terminal.get('verified',False) if terminal else False,
                          scopes=scopes, total_instructions=total,
                          oracle_checks=len(checks), wall_seconds=time.monotonic()-started)
            summary.append(result)
            archive.write(json.dumps(dict(**result, command=command, stdout=stdout, stderr=stderr, profiles=profiles),separators=(',',':'))+'\n')
            archive.flush()
            if count % 20 == 0 or status != 'finished':
                print(count, mode, case, variant, status, flush=True)

        for ci, wi, variant in itertools.product(range(10), range(8), ['charts_linear','pair_table']):
            seed, log = workloads[ci][wi]
            run('dlp', ci, variant, seed, log)
        seen = set()
        for ci, c in enumerate(contract['stage_cases']):
            curve = tuple(c[k] for k in ['n','k','a','b'])
            if curve in seen:
                continue
            seen.add(curve)
            for seed, log in workloads[ci]:
                run('rho', ci, 'rho', seed, log)
            for mode in ['calibrate','calibrate_empty']:
                run(mode, ci, 'calibration', 0)
        archive.close()
        with lzma.open(Path(tmp)/'records.jsonl.xz', 'rt') as check:
            assert sum(1 for line in check if json.loads(line)) == count
        shutil.copyfile(Path(tmp)/'records.jsonl.xz', args.output/'profiles.jsonl.xz')
    (args.output/'counts.json').write_text(json.dumps(summary,indent=2)+'\n')
    (args.output/'complete.json').write_text(json.dumps(dict(processes=count,complete=True))+'\n')


if __name__ == '__main__':
    main()
