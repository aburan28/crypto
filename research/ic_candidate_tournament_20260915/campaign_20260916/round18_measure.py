#!/usr/bin/env python3
"""Development measurement for round 0018: verify every arm against the oracle
on every round-0017 confirmation fixture, then measure Ir and native wall per
cell for each arm against the incumbent and against rho.

    python3 campaign_20260916/round18_measure.py ARM_DIR [reps]

ARM_DIR holds the three arm workers built from the trees
`round18_candidates.py` writes -- either as `ARM_DIR/<arm>` or in the cargo
layout `ARM_DIR/<arm>/target/.../examples/ic_tournament_worker`.

Arms, all built from `runs/round-0017/source_candidates/orbits/source`:
  block   -- round18-block.patch: the scan block's clamp ceiling, 64 -> 16.
  column  -- round18-column.patch: a column is the orbit representative itself
             rather than `[h]` of it, with the matching factor of `h` gone from
             the row and the descent, and the report declaring
             `column_convention: representative`.
  both    -- the two patches together.

Every arm must recover exactly the incumbent's logarithms; the arms differ in
what they cost, never in what they answer, and this script refuses to report a
measurement if that is not so.

Raw rows are written under the system temp directory and the path is printed;
the tables below are the record.
"""
import collections
import json
import random
import re
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
import tournament as T  # noqa: E402
from oracle import verify  # noqa: E402

ARM_NAMES = ('block', 'column', 'both')
INCUMBENT = ROOT / 'runs/round-0017/source_candidates/orbits/worker'
RHO = ROOT / 'runs/round-0017/worker'
FIXTURES = ROOT / 'runs/round-0017/fixtures.json'
CANDIDATES = ROOT / 'runs/round-0017/candidates.json'
OUT = Path(tempfile.gettempdir()) / 'round18_measure.json'
HOLDOUT = {'n19a1', 'n29a1'}


def ir(binary, job, env):
    r = subprocess.run(
        ['valgrind', '--tool=callgrind', '--callgrind-out-file=/dev/null', str(binary)],
        input=json.dumps(job), text=True, capture_output=True, env=env, timeout=900)
    m = re.search(r'Collected\s*:\s*(\d+)', r.stderr)
    if not m:
        raise RuntimeError(f'no Ir collected: {r.stderr[-400:]}')
    return int(m.group(1))


def run(binary, job, env):
    r = subprocess.run([str(binary)], input=json.dumps(job), text=True,
                       capture_output=True, env=env, timeout=900)
    if r.returncode != 0:
        raise RuntimeError(f'worker exit {r.returncode}: {r.stderr[-400:]}')
    return json.loads(r.stdout)


def main():
    if len(sys.argv) < 2:
        raise SystemExit('usage: round18_measure.py ARM_DIR [reps]')
    arm_dir = Path(sys.argv[1]).resolve()
    # One directory holding the three arm workers, however they were built:
    # either <dir>/<arm> or the cargo layout <dir>/arm18-<arm>/target/.../examples/.
    arms = {}
    for a in ARM_NAMES:
        for cand in (arm_dir / a,
                     arm_dir / f'arm18-{a}',
                     arm_dir / a / 'target/x86_64-unknown-linux-musl/release/examples/ic_tournament_worker',
                     arm_dir / f'arm18-{a}/target/x86_64-unknown-linux-musl/release/examples/ic_tournament_worker'):
            if cand.is_file():
                arms[a] = cand
                break
        else:
            raise SystemExit(f'no worker for arm {a!r} under {arm_dir}')
    ARMS = arms
    reps = int(sys.argv[2]) if len(sys.argv) > 2 else 5
    fx = json.loads(FIXTURES.read_text())['confirmation']
    cfg = json.loads(CANDIDATES.read_text())[0]['config']
    env = T.child_env()
    cells = sorted({c['cell'] for c in fx}, key=lambda s: (int(s[1:s.index('a')]), s))

    # 1. Oracle-verify every arm on all 96 confirmation fixtures and require
    #    the incumbent's own solutions back.
    bad = collections.Counter()
    conv = collections.Counter()
    for case in fx:
        job = dict(case['job'], mode='ic', config=cfg)
        want = verify(run(INCUMBENT, job, env), case['fixture'], expected_mode='ic', summands=3)
        for name, binary in ARMS.items():
            rep = run(binary, job, env)
            got = verify(rep, case['fixture'], expected_mode='ic', summands=3)
            conv[(name, rep.get('column_convention', 'cofactor'))] += 1
            if got['solutions'] != want['solutions'] or got['factor_base_sha256'] != want['factor_base_sha256']:
                bad[name] += 1
    print(f'verified {len(fx)} fixtures x {len(ARMS)} arms against the oracle')
    for (name, c), n in sorted(conv.items()):
        print(f"  {name:7} column_convention={c:15} on {n} fixtures")
    if bad:
        raise SystemExit(f'arms disagree with the incumbent: {dict(bad)}')
    print('  every arm returned the incumbent\'s logarithms and base hash on every fixture\n')

    # 2. Ir and native wall, four fixtures a cell, arms interleaved.
    rows = collections.defaultdict(list)
    for cell in cells:
        for case in [x for x in fx if x['cell'] == cell][:4]:
            jobs = {'incumbent': (INCUMBENT, dict(case['job'], mode='ic', config=cfg)),
                    'rho': (RHO, dict(case['job'], mode='rho', config=cfg))}
            for name, binary in ARMS.items():
                jobs[name] = (binary, dict(case['job'], mode='ic', config=cfg))
            irs = {k: ir(b, j, env) for k, (b, j) in jobs.items()}
            walls = {k: [] for k in jobs}
            for _ in range(reps):
                order = list(jobs)
                random.shuffle(order)
                for k in order:
                    b, j = jobs[k]
                    payload = json.dumps(j).encode()
                    t0 = time.perf_counter()
                    subprocess.run([str(b)], input=payload, capture_output=True, env=env)
                    walls[k].append(time.perf_counter() - t0)
            rows[cell].append((irs, {k: statistics.median(v) for k, v in walls.items()}))

    gm = statistics.geometric_mean
    for metric, idx in (('instructions', 0), ('native wall', 1)):
        print(f'{metric}: each arm over the incumbent, and over rho  (H = holdout)')
        names = list(ARMS)
        head = '  '.join(f'{n + "/inc":>12}' for n in names) + '  ' + '  '.join(f'{n + "/rho":>12}' for n in names)
        print(f"{'cell':8} {'inc/rho':>9}  {head}")
        for cell in cells:
            R = rows[cell]
            f = lambda a, b: gm([r[idx][a] / r[idx][b] for r in R])
            tag = 'H' if cell in HOLDOUT else ' '
            vals = '  '.join(f'{f(n, "incumbent"):12.4f}' for n in names)
            vals2 = '  '.join(f'{f(n, "rho"):12.4f}' for n in names)
            print(f'{cell:7}{tag} {f("incumbent", "rho"):9.4f}  {vals}  {vals2}')
        dev = [c for c in cells if c not in HOLDOUT]
        for label, group in (('development', dev), ('all eight', cells)):
            line = []
            for n in names:
                line.append(gm([gm([r[idx][n] / r[idx]['incumbent'] for r in rows[c]]) for c in group]))
            print(f'  {label:12} gm over cells, arm/incumbent: ' +
                  '  '.join(f'{n}={v:.4f}' for n, v in zip(names, line)))
        print()

    OUT.write_text(json.dumps({c: [(i, w) for i, w in rows[c]] for c in cells}, indent=1))
    print('raw rows written to', OUT)


if __name__ == '__main__':
    main()
