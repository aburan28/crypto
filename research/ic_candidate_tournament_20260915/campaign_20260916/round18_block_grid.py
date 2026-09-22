#!/usr/bin/env python3
"""Measure a grid of fixed scan-block sizes at every cell of the round-0017
panel, to decide what the block rule should be.

    python3 campaign_20260916/round18_block_grid.py PROBE_WORKER [fixtures_per_cell]

`round18_block_sweep.py` showed the shipped rule oversizes the block: it aims
at one expected witness per block, and at n23a1 the value it asks for (102)
costs 1.04x while 16 costs 0.971x.  That sweep moved the block by lowering the
clamp ceiling, which cannot push a block *below* the rule's own value, so at
the three cells where the rule already asks for the floor of 8 it could only
ever measure 8.  This script sets the block outright and so covers every cell
with the same grid, which is what choosing a constant requires.

Instructions are deterministic for a fixed input, so the spread across the
fixtures of a cell is real per-fixture variation, not measurement noise; the
geometric mean of the per-fixture ratio against the shipped block is reported
alongside the median so a single odd fixture is visible.

DEVELOPMENT CELLS DECIDE.  The panel holds out n19a1 and n29a1, and any rule
fitted here is fitted on the six development cells only; the holdout columns
are printed so the pre-registration can state a prediction for them that the
tournament then tests.

Raw rows are written under the system temp directory and the path is printed;
the tables below are the record.
"""
import collections
import json
import re
import statistics
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
import tournament as T  # noqa: E402

ORBITS = ROOT / 'runs/round-0017/source_candidates/orbits/worker'
FIXTURES = ROOT / 'runs/round-0017/fixtures.json'
CANDIDATES = ROOT / 'runs/round-0017/candidates.json'
OUT = Path(tempfile.gettempdir()) / 'round18_block_grid.json'
BLOCKS = [4, 8, 12, 16, 24, 32, 48, 64]
HOLDOUT = {'n19a1', 'n29a1'}


def ir(binary, job, env):
    r = subprocess.run(
        ['valgrind', '--tool=callgrind', '--callgrind-out-file=/dev/null', str(binary)],
        input=json.dumps(job), text=True, capture_output=True, env=env, timeout=900)
    m = re.search(r'Collected\s*:\s*(\d+)', r.stderr)
    if not m:
        raise RuntimeError(f'no Ir collected: {r.stderr[-400:]}')
    return int(m.group(1))


def env_for(env, block=None):
    e = dict(env)
    for k in ('IC_SCAN_CEIL', 'IC_SCAN_BLOCK', 'IC_TABLE_ROWS'):
        e.pop(k, None)
    if block is not None:
        e['IC_SCAN_BLOCK'] = str(block)
    return e


def run(binary, job, e):
    r = subprocess.run([str(binary)], input=json.dumps(job), text=True,
                       capture_output=True, env=e, timeout=900)
    if r.returncode != 0:
        raise RuntimeError(f'worker exit {r.returncode}: {r.stderr[-400:]}')
    return json.loads(r.stdout)


def main():
    probe = Path(sys.argv[1]).resolve()
    per_cell = int(sys.argv[2]) if len(sys.argv) > 2 else 6
    fx = json.loads(FIXTURES.read_text())['confirmation']
    cfg = json.loads(CANDIDATES.read_text())[0]['config']
    env = T.child_env()
    cells = sorted({c['cell'] for c in fx}, key=lambda s: (int(s[1:s.index('a')]), s))

    drop = lambda d: {k: v for k, v in d.items() if k != 'elapsed_seconds'}
    for cell in cells:
        job = dict([x for x in fx if x['cell'] == cell][0]['job'], mode='ic', config=cfg)
        if drop(run(probe, job, env_for(env))) != drop(run(ORBITS, job, env_for(env))):
            raise SystemExit(f'probe differs from the frozen orbits worker at {cell}; refusing to measure')
    print(f'probe verified byte-identical to the frozen orbits worker on {len(cells)} cells '
          f'(every field but elapsed_seconds)\n')

    per = collections.defaultdict(dict)   # cell -> block -> [ir per fixture]
    meta = {}
    raw = []
    for cell in cells:
        cases = [x for x in fx if x['cell'] == cell][:per_cell]
        base = run(probe, dict(cases[0]['job'], mode='ic', config=cfg), env_for(env))
        n = int(cell[1:cell.index('a')])
        size = len(base['factor_base_orbits']) * 2 * n
        shipped = base['pair_table']['scan_block']
        meta[cell] = {'size': size, 'shipped': shipped, 'rows': base['pair_table']['rows']}
        for block in sorted({*BLOCKS, shipped}):
            if block > size:
                continue
            vals = []
            for case in cases:
                job = dict(case['job'], mode='ic', config=cfg)
                e = env_for(env, block)
                rep = run(probe, job, e)
                assert rep['status'] == 'complete', (cell, block, rep['status'])
                assert rep['pair_table']['scan_block'] == min(block, size), (cell, block, rep['pair_table'])
                v = ir(probe, job, e)
                vals.append(v)
                raw.append({'cell': cell, 'block': block, 'fixture': case['id'],
                            'ir': v, 'trials': rep['trials']})
            per[cell][block] = vals

    gm = statistics.geometric_mean
    print(f'per-cell Ir against the shipped block, {per_cell} fixtures each '
          f'(gm = geometric mean of the per-fixture ratio; * = shipped, H = holdout)')
    header = '  '.join(f'{b:>7}' for b in BLOCKS)
    print(f"\n{'cell':8} {'ship':>4}  {header}")
    for cell in cells:
        ship = meta[cell]['shipped']
        base_vals = per[cell][ship]
        row = []
        for b in BLOCKS:
            if b not in per[cell]:
                row.append(f'{"-":>7}')
                continue
            ratio = gm([a / c for a, c in zip(per[cell][b], base_vals)])
            row.append(f'{ratio:7.4f}' + ('*' if b == ship else ' '))
        tag = 'H' if cell in HOLDOUT else ' '
        print(f'{cell:7}{tag} {ship:>4}  ' + ' '.join(row))

    print('\ndevelopment cells only (the six the rule may be fitted on):')
    dev = [c for c in cells if c not in HOLDOUT]
    print(f"{'block':>6} {'gm over dev cells':>18} {'worst dev cell':>16} {'worst value':>12}")
    best_block, best_score = None, None
    for b in BLOCKS:
        ratios = {}
        for cell in dev:
            if b not in per[cell]:
                continue
            ratios[cell] = gm([x / y for x, y in zip(per[cell][b], per[cell][meta[cell]['shipped']])])
        if len(ratios) != len(dev):
            continue
        score = gm(list(ratios.values()))
        worst = max(ratios, key=ratios.get)
        print(f'{b:6} {score:18.4f} {worst:>16} {ratios[worst]:12.4f}')
        if best_score is None or score < best_score:
            best_block, best_score = b, score
    print(f'\nbest fixed block on the development cells: {best_block} '
          f'(gm {best_score:.4f} against the shipped rule)')

    print('\nwhat that block would do at the two holdout cells:')
    for cell in sorted(HOLDOUT):
        if best_block in per[cell]:
            ratio = gm([x / y for x, y in zip(per[cell][best_block], per[cell][meta[cell]['shipped']])])
            print(f'  {cell}: {ratio:.4f}  (shipped block {meta[cell]["shipped"]})')

    OUT.write_text(json.dumps(raw, indent=1))
    print('\nraw rows written to', OUT)


if __name__ == '__main__':
    main()
