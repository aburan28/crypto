#!/usr/bin/env python3
"""Measure the true cost of `t` pair-table rows at every cell of the round-0017
panel, against the analytic rule the shipped worker uses.

    python3 campaign_20260916/round18_row_sweep.py PROBE_WORKER [fixtures_per_cell]

`PROBE_WORKER` is a build of `runs/round-0017/source_candidates/orbits/source`
with `campaign_20260916/round18-rows-probe.patch` applied, which lets
`IC_TABLE_ROWS` override the row count (and, through it, the scan block).  With
the variable unset the probe is byte-identical to the frozen `orbits` worker on
every confirmation fixture except `elapsed_seconds`; this script checks that
first and refuses to measure if it does not hold.

Round 0017 falsified the analytic rule's scan side at n23a1: it predicted two
rows would cost 0.973x three, and two rows actually cost 1.272x three.  The
rule's build term `t*|F|` is exact, so the error is in the scan term
`(K + TABLE_PROBES)/(lambda*cov(t))`.  This sweep measures the real minimum so
a corrected rule can be fitted to it rather than to an unchecked model.

Raw rows are written under the system temp directory and the path is printed;
the table below is the record.
"""
import collections
import json
import os
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
OUT = Path(tempfile.gettempdir()) / 'round18_row_sweep.json'
TABLE_PROBES = 3.0  # the constant the shipped rule uses


def ir(binary, job, env):
    """Valgrind Ir of one complete cold job."""
    r = subprocess.run(
        ['valgrind', '--tool=callgrind', '--callgrind-out-file=/dev/null', str(binary)],
        input=json.dumps(job), text=True, capture_output=True, env=env, timeout=600)
    m = re.search(r'Collected\s*:\s*(\d+)', r.stderr)
    if not m:
        raise RuntimeError(f'no Ir collected: {r.stderr[-400:]}')
    return int(m.group(1))


def run(binary, job, env, rows=None):
    e = dict(env)
    if rows is None:
        e.pop('IC_TABLE_ROWS', None)
    else:
        e['IC_TABLE_ROWS'] = str(rows)
    r = subprocess.run([str(binary)], input=json.dumps(job), text=True,
                       capture_output=True, env=e, timeout=600)
    if r.returncode != 0:
        raise RuntimeError(f'worker exit {r.returncode}: {r.stderr[-400:]}')
    return json.loads(r.stdout)


def analytic(size, orbits, r, rows):
    """The shipped rule's expected work, for comparison with the measurement."""
    lam = size * (size + 1.0) / 2.0 / r
    rest = (orbits - rows) / orbits
    cov = 1.0 - rest * rest
    return rows * size + (orbits + TABLE_PROBES) / max(lam * cov, 1e-12)


def main():
    probe = Path(sys.argv[1]).resolve()
    per_cell = int(sys.argv[2]) if len(sys.argv) > 2 else 3
    fx = json.loads(FIXTURES.read_text())['confirmation']
    cfg = json.loads(CANDIDATES.read_text())[0]['config']
    env = T.child_env()
    cells = sorted({c['cell'] for c in fx}, key=lambda s: (int(s[1:s.index('a')]), s))

    # 1. The probe must be the frozen arm when the override is absent.
    drop = lambda d: {k: v for k, v in d.items() if k != 'elapsed_seconds'}
    for cell in cells:
        case = [x for x in fx if x['cell'] == cell][0]
        job = dict(case['job'], mode='ic', config=cfg)
        if drop(run(probe, job, env)) != drop(run(ORBITS, job, env)):
            raise SystemExit(f'probe differs from the frozen orbits worker at {cell}; refusing to measure')
    print(f'probe verified byte-identical to the frozen orbits worker on {len(cells)} cells '
          f'(every field but elapsed_seconds)\n')

    rows_out = collections.defaultdict(dict)
    raw = []
    for cell in cells:
        cases = [x for x in fx if x['cell'] == cell][:per_cell]
        shipped = run(probe, dict(cases[0]['job'], mode='ic', config=cfg), env)
        K = shipped['pair_table']['rows_possible']
        n = int(cell[1:cell.index('a')])
        size = len(shipped['factor_base_orbits']) * 2 * n
        order = cases[0]['job']['curve']['subgroup_order'] if 'curve' in cases[0]['job'] else None
        for t in range(1, K + 1):
            irs, trials, blocks = [], [], []
            for case in cases:
                job = dict(case['job'], mode='ic', config=cfg)
                e = dict(env, IC_TABLE_ROWS=str(t))
                rep = run(probe, job, env, rows=t)
                assert rep['pair_table']['rows'] == t, (cell, t, rep['pair_table'])
                assert rep['status'] == 'complete', (cell, t, rep['status'])
                irs.append(ir(probe, job, e))
                trials.append(rep['trials'])
                blocks.append(rep['pair_table']['scan_block'])
                raw.append({'cell': cell, 'rows': t, 'fixture': case['id'],
                            'ir': irs[-1], 'trials': trials[-1], 'scan_block': blocks[-1]})
            rows_out[cell][t] = {'ir': statistics.median(irs), 'trials': statistics.median(trials),
                                 'scan_block': blocks[0]}
        rows_out[cell]['meta'] = {'K': K, 'size': size, 'shipped': shipped['pair_table']['rows'],
                                  'order': order}

    print(f"{'cell':7} {'K':>2} {'|F|':>4} {'shipped':>7} {'best':>4} {'cost(best)/cost(shipped)':>24}")
    for cell in cells:
        meta = rows_out[cell]['meta']
        ts = {t: v['ir'] for t, v in rows_out[cell].items() if isinstance(t, int)}
        best = min(ts, key=ts.get)
        print(f"{cell:7} {meta['K']:2} {meta['size']:4} {meta['shipped']:7} {best:4} "
              f"{ts[best] / ts[meta['shipped']]:24.4f}")

    print(f"\nper-cell detail (median Ir over {per_cell} fixtures; * marks the shipped choice)")
    for cell in cells:
        meta = rows_out[cell]['meta']
        print(f"\n  {cell}  K={meta['K']}  |F|={meta['size']}")
        print(f"    {'t':>2} {'Ir':>10} {'Ir/shipped':>11} {'trials':>7} {'block':>6}")
        for t in sorted(k for k in rows_out[cell] if isinstance(k, int)):
            v = rows_out[cell][t]
            mark = '*' if t == meta['shipped'] else ' '
            print(f"   {mark}{t:>2} {v['ir']:10.0f} "
                  f"{v['ir'] / rows_out[cell][meta['shipped']]['ir']:11.4f} "
                  f"{v['trials']:7.0f} {v['scan_block']:6}")

    OUT.write_text(json.dumps(raw, indent=1))
    print('\nraw rows written to', OUT)


if __name__ == '__main__':
    main()
