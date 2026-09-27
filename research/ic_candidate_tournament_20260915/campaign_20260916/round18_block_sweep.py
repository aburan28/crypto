#!/usr/bin/env python3
"""Measure the cost of the scan-block clamp ceiling at every cell of the
round-0017 panel.

    python3 campaign_20260916/round18_block_sweep.py PROBE_WORKER [fixtures_per_cell]

`decompose` scans the base in blocks; each block pays one batch inversion and
is discarded on an early exit, so the shipped rule sizes a block at about one
expected witness, `chunk = (1/hit + 1).clamp(8, 64).min(|F|)`, with
`hit = lambda * cov(t)`.

The ceiling of 64 in that clamp binds at exactly two cells of the panel, the
two with the largest subgroups: n23a0 wants 71 and n23a1 wants 101, and every
other cell's own value is already below 64.  n23a1 is the cell the campaign
keeps losing, so the question this sweep answers is whether the ceiling is
costing anything, or whether 64 happens to be at least as good as the value
the rule asks for.

`PROBE_WORKER` is a build of `runs/round-0017/source_candidates/orbits/source`
with `campaign_20260916/round18-rows-probe.patch` applied; with none of its
environment variables set it is byte-identical to the frozen `orbits` worker
except for `elapsed_seconds`, which this script checks before measuring.

Raw rows are written under the system temp directory and the path is printed;
the table below is the record.
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
OUT = Path(tempfile.gettempdir()) / 'round18_block_sweep.json'
# The shipped ceiling, then the values the rule actually asks for at the two
# clamped cells, then past them to find the true shape.
CEILINGS = [8, 16, 32, 64, 96, 128, 192, 256, 384]


def ir(binary, job, env):
    r = subprocess.run(
        ['valgrind', '--tool=callgrind', '--callgrind-out-file=/dev/null', str(binary)],
        input=json.dumps(job), text=True, capture_output=True, env=env, timeout=900)
    m = re.search(r'Collected\s*:\s*(\d+)', r.stderr)
    if not m:
        raise RuntimeError(f'no Ir collected: {r.stderr[-400:]}')
    return int(m.group(1))


def run(binary, job, env, ceil=None):
    e = dict(env)
    e.pop('IC_SCAN_CEIL', None)
    e.pop('IC_SCAN_BLOCK', None)
    e.pop('IC_TABLE_ROWS', None)
    if ceil is not None:
        e['IC_SCAN_CEIL'] = str(ceil)
    r = subprocess.run([str(binary)], input=json.dumps(job), text=True,
                       capture_output=True, env=e, timeout=900)
    if r.returncode != 0:
        raise RuntimeError(f'worker exit {r.returncode}: {r.stderr[-400:]}')
    return json.loads(r.stdout), e


def main():
    probe = Path(sys.argv[1]).resolve()
    per_cell = int(sys.argv[2]) if len(sys.argv) > 2 else 3
    fx = json.loads(FIXTURES.read_text())['confirmation']
    cfg = json.loads(CANDIDATES.read_text())[0]['config']
    env = T.child_env()
    cells = sorted({c['cell'] for c in fx}, key=lambda s: (int(s[1:s.index('a')]), s))

    drop = lambda d: {k: v for k, v in d.items() if k != 'elapsed_seconds'}
    for cell in cells:
        case = [x for x in fx if x['cell'] == cell][0]
        job = dict(case['job'], mode='ic', config=cfg)
        if drop(run(probe, job, env)[0]) != drop(run(ORBITS, job, env)[0]):
            raise SystemExit(f'probe differs from the frozen orbits worker at {cell}; refusing to measure')
    print(f'probe verified byte-identical to the frozen orbits worker on {len(cells)} cells '
          f'(every field but elapsed_seconds)\n')

    out = collections.defaultdict(dict)
    raw = []
    for cell in cells:
        cases = [x for x in fx if x['cell'] == cell][:per_cell]
        base, _ = run(probe, dict(cases[0]['job'], mode='ic', config=cfg), env)
        shipped_block = base['pair_table']['scan_block']
        n = int(cell[1:cell.index('a')])
        size = len(base['factor_base_orbits']) * 2 * n
        seen = set()
        for ceil in CEILINGS:
            irs, blocks, trials = [], [], []
            for case in cases:
                job = dict(case['job'], mode='ic', config=cfg)
                rep, e = run(probe, job, env, ceil=ceil)
                assert rep['status'] == 'complete', (cell, ceil, rep['status'])
                blocks.append(rep['pair_table']['scan_block'])
                trials.append(rep['trials'])
                irs.append(ir(probe, job, e))
                raw.append({'cell': cell, 'ceiling': ceil, 'block': blocks[-1],
                            'fixture': case['id'], 'ir': irs[-1], 'trials': trials[-1]})
            blk = blocks[0]
            # Once the ceiling stops binding the block stops changing; record
            # each distinct block once so the table shows the real sweep.
            if blk in seen:
                continue
            seen.add(blk)
            out[cell][blk] = {'ir': statistics.median(irs), 'trials': statistics.median(trials),
                              'ceiling': ceil}
        out[cell]['meta'] = {'size': size, 'shipped_block': shipped_block,
                             'rows': base['pair_table']['rows']}

    print(f"{'cell':7} {'|F|':>4} {'t':>2} {'shipped':>7} {'best':>5} {'best/shipped':>13} {'clamped':>8}")
    for cell in cells:
        meta = out[cell]['meta']
        ts = {b: v['ir'] for b, v in out[cell].items() if isinstance(b, int)}
        best = min(ts, key=ts.get)
        clamped = 'YES' if meta['shipped_block'] == 64 else 'no'
        print(f"{cell:7} {meta['size']:4} {meta['rows']:2} {meta['shipped_block']:7} {best:5} "
              f"{ts[best] / ts[meta['shipped_block']]:13.4f} {clamped:>8}")

    print(f"\nper-cell detail (median Ir over {per_cell} fixtures; * marks the shipped block)")
    for cell in cells:
        meta = out[cell]['meta']
        print(f"\n  {cell}  |F|={meta['size']}  rows={meta['rows']}")
        print(f"    {'block':>5} {'Ir':>10} {'Ir/shipped':>11} {'trials':>7}")
        for blk in sorted(b for b in out[cell] if isinstance(b, int)):
            v = out[cell][blk]
            mark = '*' if blk == meta['shipped_block'] else ' '
            print(f"   {mark}{blk:>5} {v['ir']:10.0f} "
                  f"{v['ir'] / out[cell][meta['shipped_block']]['ir']:11.4f} {v['trials']:7.0f}")

    OUT.write_text(json.dumps(raw, indent=1))
    print('\nraw rows written to', OUT)


if __name__ == '__main__':
    main()
