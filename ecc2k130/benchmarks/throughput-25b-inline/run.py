#!/usr/bin/env python3
"""Run from ecc2k130/: python3 benchmarks/throughput-25b-inline/run.py."""
import argparse
import hashlib
import json
import pathlib
import subprocess
import sys

sys.path.insert(0, '.')
from codegen.benchreport import benchResult, summarizeSamples, reportsVerified
parser = argparse.ArgumentParser()
parser.add_argument('--output', required=True, type=pathlib.Path)
args = parser.parse_args()
out = args.output.resolve()
if (out / 'results.json').exists():
    raise SystemExit('Refusing to overwrite existing measurements; use a new build directory')
source = out / 'source'
result = {'targetB': 25, 'class': 'engineering', 'unit': 'B complete scalar updates/s',
          'verification': {}, 'runs': {}, 'sourceHashes': {}}
for name in ['include/packed131.h', 'include/packedkernels.cuh', 'Makefile']:
    result['sourceHashes'][name] = hashlib.sha256((source / name).read_bytes()).hexdigest()

def save():
    (out/'results.json').write_text(json.dumps(result, indent=2))

def run(label, args):
    busy = subprocess.check_output(['nvidia-smi', '--query-compute-apps=pid,process_name', '--format=csv,noheader'], text=True).strip()
    if busy:
        raise RuntimeError('GPU occupied: '+busy)
    print(label, flush=True)
    with (out/(label+'.log')).open('w') as log:
        proc = subprocess.Popen(args, stdout=log, stderr=subprocess.STDOUT, text=True)
        contaminated = False
        try:
            for _ in range(300):
                try:
                    proc.wait(timeout=1)
                    break
                except subprocess.TimeoutExpired:
                    active = subprocess.check_output(['nvidia-smi', '--query-compute-apps=pid', '--format=csv,noheader'], text=True).splitlines()
                    if any(pid.strip() != str(proc.pid) for pid in active if pid.strip()):
                        contaminated = True
            else:
                proc.kill()
                proc.wait()
                raise RuntimeError('Run timed out')
        finally:
            if proc.poll() is None:
                proc.kill()
                proc.wait()
    p = subprocess.CompletedProcess(args, proc.returncode, (out/(label+'.log')).read_text())
    if contaminated:
        result.setdefault('contention', []).append(label)
        save()
        if label.startswith('bench-'):
            raise RuntimeError('Another GPU process appeared; timing sample rejected')
    result['runs'][label] = {'args': args, 'returncode': p.returncode}
    save()
    return p

result['gpu'] = subprocess.check_output(['nvidia-smi', '--query-gpu=name,uuid,clocks.sm,power.limit', '--format=csv'], text=True)
samples = {str(i): [] for i in range(4)}
for mode in range(4):
    binary = str(out/f'client-{mode}')
    result.setdefault('binaryHashes', {})[str(mode)] = hashlib.sha256(pathlib.Path(binary).read_bytes()).hexdigest()
    for seed in [0, 1937]:
        p = run(f'verify-{mode}-{seed}', [binary, '--curve', '131', '--packed', '--dp-weight', '50', '--steps', '16', '--launches', '6', '--dp-cap', '262144', '--threads', '1024', '--verify', '300', '--run-id', str(seed)])
        ok = reportsVerified(p.returncode, p.stdout, required=300)
        result['verification'][f'{mode}-{seed}'] = ok
        save()
        if not ok:
            raise RuntimeError('Report verification failed')
for rep in range(3):
    for mode in (range(4) if rep % 2 == 0 else reversed(range(4))):
        args = [str(out/f'client-{mode}'), '--curve', '131', '--packed', '--bench', '--steps', '1024', '--launches', '32', '--verify', '0', '--run-id', str(1937 if rep == 2 else 0)]
        p = run(f'bench-{mode}-{rep}', args)
        samples[str(mode)].append(benchResult(args, p.returncode, p.stdout))
        result['samples'] = samples
        save()
result['summary'] = {m: summarizeSamples(s) for m,s in samples.items()}
base = result['summary']['0']['rate']
result['table'] = {m: {'B/s': s['rate']/1000, 'ratioToBaseline': s['rate']/base if base else None, 'ratioToTarget': s['rate']/25000, 'valid': s['valid']} for m,s in result['summary'].items()}
save()
print(json.dumps(result['table'], indent=2), flush=True)
