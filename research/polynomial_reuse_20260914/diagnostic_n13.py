#!/usr/bin/env python3
"""Diagnostic outside the frozen contract: can the preprocessing cache pay at n = 13?

run-002 found every cache mode slower than `off` at the frozen degrees (7, 9),
and the cache's own counters showed why: 2.5-3.5 lookups per process cannot
repay one cold miss.  n = 13 is the one larger degree the frozen invocation
reaches (n = 11 has no usable prime-order subgroup in the constructor; n = 15
exhausts 500 trials with no relation), and there a process makes 14-23 lookups.

Same conventions as run.py: IC_* stripped, RAYON_NUM_THREADS=1, 30 s and 2 GiB
per process, raw stdout/stderr kept, every status recorded.  Only `off` and
`preprocess-local` are compared -- the reduction cache is not reached by the
default engine, and Redis only adds cost.  Five repetitions, the mode order
alternating between repetitions so that drift favours neither.  A diagnostic,
not a frozen-suite result and not a speedup claim.
"""
import argparse, json, os, pathlib, resource, subprocess, time
ROOT = pathlib.Path(__file__).resolve().parents[2]

def limit():
    resource.setrlimit(resource.RLIMIT_AS, (2 * 1024**3, 2 * 1024**3))

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--output', type=pathlib.Path, required=True)
    a = p.parse_args(); a.output.mkdir(parents=True, exist_ok=False)
    env = {k: v for k, v in os.environ.items() if not k.startswith('IC_')}
    env['RAYON_NUM_THREADS'] = '1'
    ic = ROOT / 'target/release/ic'
    modes = {'off': {}, 'preprocess-local': {'IC_PREPROCESS_CACHE': 'local'}}
    records = []
    for rep in range(5):
        order = ['off', 'preprocess-local'] if rep % 2 == 0 else ['preprocess-local', 'off']
        for seed in (101, 202, 307, 409):
            for solver in ('groebner', 'sat'):
                for mode in order:
                    name = f'n13-{solver}-s{seed}-r{rep}-{mode}'
                    cmd = [str(ic), 'run', '--degree', '13', '--solver', solver, '--seed', str(seed),
                           '--known-log', '5', '--batch', '1', '--max-trials', '500', '--json']
                    start = time.monotonic(); status = 'finished'; code = None
                    try:
                        r = subprocess.run(cmd, env=env | modes[mode], capture_output=True, timeout=30, preexec_fn=limit)
                        out, err, code = r.stdout, r.stderr, r.returncode
                    except subprocess.TimeoutExpired as e:
                        status, out, err = 'timeout', e.stdout or b'', e.stderr or b''
                    (a.output / (name + '.stdout')).write_bytes(out)
                    (a.output / (name + '.stderr')).write_bytes(err)
                    records.append({'name': name, 'command': cmd, 'mode': mode, 'status': status,
                                    'returncode': code, 'wall_seconds': time.monotonic() - start})
                    (a.output / 'processes.json').write_text(json.dumps(records, indent=2) + '\n')
    print('Saved', len(records), 'process records', flush=True)

if __name__ == '__main__':
    main()
