#!/usr/bin/env python3
"""Measure the factor-base size curve, to check the derivation round 0019's
headline rests on.

    python3 campaign_20260916/round19_base_sweep.py [cases]

`round19_model.py` derives, from the solver's own published cost model, that
the base the panel builds is already at the model's optimum at every cell and
that enlarging it is a loss.  That derivation is what says round 0019 has no
factor-base lever to take, so it is measured here rather than taken on faith.

**The parameter is the orbit count, not the requested point count.**  The first
version of this sweep varied `factor_base.points` from `2*degree` to `12*degree`
and measured a flat line: instruction counts identical to six figures at every
value.  The reason is in `build_subgroup_orbit_factor_base`, which samples
orbits in batches of eight and stops at the first rebuild that reaches the
target, so every request at or below one batch returns that whole batch.  The
contract's `6*degree` is below one batch at every cell on the panel.  **The
panel's factor base is seven or eight orbits everywhere, and the contract's
`points` has been inert since round 0002.**  The base can be enlarged in whole
batches; it cannot be shrunk below one.  This sweep therefore requests point
counts that cross batch boundaries and records the orbit count each one
actually built.

Three controls:

  * rho never reads the factor base, so its Ir must not move with it.  The
    residual permitted here is 0.1%: the request is part of the job document,
    so a longer `points` field costs a few dozen instructions to parse, and
    that -- tens of instructions in millions -- is the whole of the observed
    drift.  Anything larger means the sweep is measuring the harness.
  * every run must return the incumbent's own logarithm, verified against
    `oracle.py`.  A cheaper base that answers differently is not a base.
  * the base the panel actually builds is in the sweep, so the curve is
    anchored at the value every round from 0002 on has used.

Changing the base size changes the instance IC is given, not the discrete
logarithm being solved, so this is a diagnostic and not a tournament stage.
Raw rows go under the system temp directory; the printed tables are the record.
"""
import collections
import json
import re
import statistics as st
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
import tournament as T  # noqa: E402
from oracle import verify  # noqa: E402

WORKER = ROOT / 'runs/round-0018b/worker'
FIXTURES = ROOT / 'runs/round-0018b/fixtures.json'
CANDIDATES = ROOT / 'runs/round-0018b/candidates.json'
OUT = Path(tempfile.gettempdir()) / 'round19_base_sweep.json'
MODEL = json.loads(Path(__file__).with_name('round-0019-base-size-model.json').read_text())
BATCH = 8
BATCHES = (1, 2, 3, 4, 5)   # requested multiples of one batch of orbits
RHO_TOLERANCE = 1e-3


def _ir(job, env, timeout=3600):
    r = subprocess.run(
        ['valgrind', '--tool=callgrind', '--callgrind-out-file=/dev/null', str(WORKER)],
        input=json.dumps(job), text=True, capture_output=True, env=env, timeout=timeout)
    m = re.search(r'Collected\s*:\s*(\d+)', r.stderr)
    if not m:
        raise RuntimeError(f'no Ir collected: {r.stderr[-400:]}')
    return int(m.group(1))


def _run(job, env, timeout=3600):
    r = subprocess.run([str(WORKER)], input=json.dumps(job), text=True,
                       capture_output=True, env=env, timeout=timeout)
    if r.returncode != 0:
        raise RuntimeError(f'worker exit {r.returncode}: {r.stderr[-300:]}')
    return json.loads(r.stdout)


def main():
    cases_per_cell = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    fx = json.loads(FIXTURES.read_text())['confirmation']
    cfg = json.loads(CANDIDATES.read_text())[0]['config']
    env = T.child_env()
    cells = sorted({c['cell'] for c in fx}, key=lambda s: (int(s[1:s.index('a')]), s))

    rows = collections.defaultdict(lambda: collections.defaultdict(list))
    rho = collections.defaultdict(list)
    bad = []
    for cell in cells:
        degree = int(cell[1:cell.index('a')])
        for case in [x for x in fx if x['cell'] == cell][:cases_per_cell]:
            want = None
            for batches in BATCHES:
                # One orbit is `degree` abscissae and at most 2*degree signed
                # points; ask for just past the batch boundary we want.
                points = 1 if batches == 1 else 2 * degree * BATCH * (batches - 1) + 1
                job = dict(case['job'], mode='ic', config=cfg)
                job['factor_base'] = dict(job['factor_base'], points=points)
                try:
                    report = _run(job, env)
                except Exception as exc:                        # noqa: BLE001
                    bad.append((cell, case['id'], points, str(exc)[:100]))
                    continue
                built = len(report.get('factor_base_orbits') or report.get('factor_base') or [])
                got = verify(report, case['fixture'], expected_mode='ic', summands=3)
                if want is None:
                    want = got['solutions']
                elif got['solutions'] != want:
                    bad.append((cell, case['id'], points, 'different logarithm'))
                rows[cell][built].append(_ir(job, env))
            for points in (1, 2 * degree * BATCH * 2 + 1):
                job = dict(case['job'], mode='rho', config=cfg)
                job['factor_base'] = dict(job['factor_base'], points=points)
                rho[case['id']].append(_ir(job, env))

    drift = {k: v for k, v in rho.items() if max(v) / min(v) - 1 > RHO_TOLERANCE}
    worst = max((max(v) / min(v) - 1 for v in rho.values()), default=0.0)
    print(f'rho control: {len(rho)} fixtures, worst spread {worst:.2e} '
          f'(tolerance {RHO_TOLERANCE:.0e}) -> {"invariant" if not drift else f"DRIFTED {drift}"}')
    print(f'oracle: {"same logarithm at every base size" if not bad else bad}\n')

    gm = st.geometric_mean
    built_at = {c: min(rows[c], key=lambda k: k) for c in cells if rows[c]}
    print('measured Ir over the orbit count, relative to the base the panel builds')
    ks = sorted({k for c in cells for k in rows[c]})
    print(f'{"cell":8} {"built":>5} ' + ' '.join(f'{("K=" + str(k)):>8}' for k in ks))
    for cell in cells:
        b = gm(rows[cell][built_at[cell]])
        print(f'{cell:8} {built_at[cell]:5d} ' +
              ' '.join(f'{gm(rows[cell][k]) / b:8.4f}' if rows[cell].get(k) else f'{"-":>8}'
                       for k in ks))

    print('\nthe model over the same orbit counts (round19_model.py)')
    print(f'{"cell":8} {"built":>5} ' + ' '.join(f'{("K=" + str(k)):>8}' for k in ks))
    for cell in cells:
        curve = MODEL[cell]['curve']
        print(f'{cell:8} {built_at[cell]:5d} ' +
              ' '.join(f'{curve[str(k)]:8.4f}' if str(k) in curve else f'{"-":>8}' for k in ks))

    print('\nmeasured / model (1.000 = the derivation got it right)')
    print(f'{"cell":8} {"built":>5} ' + ' '.join(f'{("K=" + str(k)):>8}' for k in ks))
    for cell in cells:
        b = gm(rows[cell][built_at[cell]])
        curve = MODEL[cell]['curve']
        print(f'{cell:8} {built_at[cell]:5d} ' +
              ' '.join(f'{(gm(rows[cell][k]) / b) / curve[str(k)]:8.4f}'
                       if rows[cell].get(k) and str(k) in curve else f'{"-":>8}' for k in ks))

    OUT.write_text(json.dumps({c: {str(k): v for k, v in rows[c].items()} for c in cells}, indent=1))
    print('\nraw rows written to', OUT)
    if bad or drift:
        raise SystemExit('controls failed; the sweep above is not evidence')


if __name__ == '__main__':
    main()
