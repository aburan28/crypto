#!/usr/bin/env python3
"""Lean rho against matched rho on round 0024's fixtures.

    python3 campaign_20260916/round25_rho_lean_check.py \
        MATCHED_WORKER LEAN_WORKER [--per-cell 6] [--cells n13a0,...] \
        [--ic-check 1] [--json OUT]

MATCHED_WORKER is `round24-sources/matched` built as it is (round 0024's rho:
Frobenius classes named by the IC arm's normal-basis rotation). LEAN_WORKER is
the same tree with `round25-rho-lean.patch` applied. Both are the tree's
`examples/ic_tournament_worker`, run in `rho` mode with the round's config.

Fixtures are round 0024's development and selection cases for every cell that
has them, and the first confirmation cases for the cells that have none there
(n19a1, n29a1, n59a0, n61a1), `--per-cell` of each (default 6).

Per fixture and arm, one callgrind run (`--callgrind-out-file` in a temporary
directory, so the worker's phase dumps are kept): its report goes through
`oracle.verify(report, fixture, expected_mode='rho')` -- which accepts only the
unique d in [0, r) with [d]G = Q, i.e. the fixture's logarithm -- and the two
arms must report the same logarithm. Recorded: total instructions (callgrind
`Collected`), the `rho_solve` phase's instructions (the worker's own
client-request dump around the rho call), iterations and walk group additions.

Per cell: the lean/matched instruction ratio (geometric mean over fixtures,
95% band from the t distribution of the log ratios), the same for the
rho_solve phase, median iterations of each arm and the geometric mean of the
per-fixture iteration ratio with its band, and instructions per walk addition
(rho_solve instructions over walk additions, pooled over the cell's fixtures).

`--ic-check K` also runs the first K fixtures of each cell in `ic` mode
(pair_table, 3 summands) through both workers and reports whether the reports
(without `elapsed_seconds`) are identical -- the patch must not change what
the IC arm computes -- and the instruction difference by worker phase (callgrind
counts of one binary vary by a few instructions from run to run, and a changed
build can move code generation of shared functions slightly).

Exit status 0 only if every report verified and the arms agreed on every
logarithm.
"""
import argparse
import json
import math
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

ROUND = ROOT / 'runs/round-0024'
CELLS = ('n13a0', 'n17a1', 'n19a0', 'n19a1', 'n23a0', 'n23a1', 'n29a1', 'n31a0',
         'n37a0', 'n43a1', 'n59a0', 'n61a1')
RHO = {'batch_trials': 1, 'linear_algebra': 'sparse', 'max_trials': 65536,
       'solver': 'pair_table', 'summands': 3}
IC = dict(RHO)
# Two-sided 97.5% t quantiles by degrees of freedom.
T975 = {1: 12.706, 2: 4.303, 3: 3.182, 4: 2.776, 5: 2.571, 6: 2.447, 7: 2.365,
        8: 2.306, 9: 2.262, 10: 2.228, 11: 2.201, 12: 2.179, 15: 2.131, 20: 2.086,
        30: 2.042, 60: 2.000}


def t975(df):
    keys = sorted(k for k in T975 if k <= df)
    return T975[keys[-1]] if keys else float('nan')


def callgrind(worker, job, env, timeout=3600):
    """(report, total Ir, {phase label: Ir}) for one run under callgrind."""
    with tempfile.TemporaryDirectory() as tmp:
        done = subprocess.run(
            ['valgrind', '--tool=callgrind', f'--callgrind-out-file={tmp}/cg.%p', worker],
            input=json.dumps(job), text=True, capture_output=True, env=env, timeout=timeout)
        found = re.search(r'Collected\s*:\s*(\d+)', done.stderr)
        if not found:
            raise RuntimeError(f'no Ir collected: {done.stderr[-300:]}')
        phases = {}
        for dump in Path(tmp).iterdir():
            text = dump.read_text(errors='replace')
            label = re.search(r'^desc: Trigger: Client Request: (\S+)', text, re.M)
            total = re.search(r'^totals:\s*(\d+)', text, re.M) or re.search(r'^summary:\s*(\d+)', text, re.M)
            if total:
                # A label the worker dumps more than once (the IC arm's
                # collection rounds) is summed; the tail after the last dump
                # is `(rest)`.
                key = label.group(1) if label else '(rest)'
                phases[key] = phases.get(key, 0) + int(total.group(1))
    report = json.loads(done.stdout) if done.stdout.strip() else {}
    return report, int(found.group(1)), phases


def band(logs):
    """Geometric mean ratio and its 95% band from the log ratios."""
    ratio = math.exp(st.mean(logs))
    if len(logs) < 2:
        return ratio, float('nan'), float('nan')
    half = t975(len(logs) - 1) * st.stdev(logs) / math.sqrt(len(logs))
    return ratio, ratio * math.exp(-half), ratio * math.exp(half)


def fmt_band(values):
    r, lo, hi = values
    return f'{r:6.3f} [{lo:6.3f}, {hi:6.3f}]'


def cases_for(fixtures, cell, per_cell):
    early = [dict(c, label=f'{s}/{c["id"]}') for s in ('development', 'selection')
             for c in fixtures[s] if c['cell'] == cell]
    if not early:
        early = [dict(c, label=f'confirmation/{c["id"]}')
                 for c in fixtures['confirmation'] if c['cell'] == cell]
    return early[:per_cell]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('matched')
    ap.add_argument('lean')
    ap.add_argument('--per-cell', type=int, default=6)
    ap.add_argument('--cells', default=','.join(CELLS))
    ap.add_argument('--ic-check', type=int, default=0)
    ap.add_argument('--json')
    args = ap.parse_args()
    workers = {'matched': str(Path(args.matched).resolve()), 'lean': str(Path(args.lean).resolve())}
    env = T.child_env()
    fixtures = json.loads((ROUND / 'fixtures.json').read_text())
    cells = args.cells.split(',')
    print(f'lean/matched rho on {ROUND.name}, up to {args.per_cell} fixtures a cell')
    print(f'  matched = {workers["matched"]}\n  lean    = {workers["lean"]}\n')
    rows, ok, summary = [], True, []
    for cell in cells:
        for case in cases_for(fixtures, cell, args.per_cell):
            job = dict(case['job'], mode='rho', config=RHO)
            row = dict(cell=cell, case=case['label'])
            answers = {}
            for arm, worker in workers.items():
                report, total, phases = callgrind(worker, job, env)
                if report.get('status') != 'complete':
                    print(f'  {case["label"]} {arm}: status {report.get("status")}')
                    ok = False
                    break
                try:
                    proof = verify(report, case['fixture'], expected_mode='rho')
                except Exception as err:  # the oracle's own rejection
                    print(f'  {case["label"]} {arm}: ORACLE REJECTED: {err}')
                    ok = False
                    break
                answers[arm] = proof['solutions']
                sol = report['solutions']
                row[arm] = dict(ir=total, rho_solve_ir=phases.get('rho_solve'),
                                iterations=sum(s['iterations'] for s in sol),
                                walk_group_additions=sum(s['walk_group_additions'] for s in sol),
                                restarts=sum(s['restarts'] for s in sol),
                                recovered=proof['solutions'])
            else:
                if answers['matched'] != answers['lean']:
                    print(f'  {case["label"]}: ARMS DISAGREE ON THE LOGARITHM')
                    ok = False
                    continue
                rows.append(row)
                m, l = row['matched'], row['lean']
                print(f'  {case["label"]:28} it {m["iterations"]:>6} / {l["iterations"]:>6}  '
                      f'Ir {m["ir"]:>11,} / {l["ir"]:>11,}  rho_solve {m["rho_solve_ir"]:>11,} / '
                      f'{l["rho_solve_ir"]:>11,}  d={l["recovered"][0]}', flush=True)
    print()
    head = (f'{"cell":6} {"n":>2} {"Ir lean/matched":>26} {"rho_solve lean/matched":>26} '
            f'{"med it m":>8} {"med it l":>8} {"iter ratio (geomean)":>26} '
            f'{"Ir/add m":>8} {"Ir/add l":>8}')
    print(head)
    all_logs = []
    for cell in cells:
        rs = [r for r in rows if r['cell'] == cell]
        if not rs:
            print(f'{cell:6}  0  --')
            continue
        logs = [math.log(r['lean']['ir'] / r['matched']['ir']) for r in rs]
        solve_logs = [math.log(r['lean']['rho_solve_ir'] / r['matched']['rho_solve_ir']) for r in rs]
        it_logs = [math.log(r['lean']['iterations'] / r['matched']['iterations']) for r in rs]
        all_logs += logs
        med_m = st.median(r['matched']['iterations'] for r in rs)
        med_l = st.median(r['lean']['iterations'] for r in rs)
        per_add = {arm: sum(r[arm]['rho_solve_ir'] for r in rs) / sum(r[arm]['walk_group_additions'] for r in rs)
                   for arm in ('matched', 'lean')}
        cell_row = dict(cell=cell, fixtures=len(rs), ir_ratio=band(logs), rho_solve_ratio=band(solve_logs),
                        median_iterations_matched=med_m, median_iterations_lean=med_l,
                        iteration_ratio=band(it_logs), ir_per_addition=per_add,
                        identical_walks=all(r['lean']['iterations'] == r['matched']['iterations']
                                            and r['lean']['walk_group_additions'] == r['matched']['walk_group_additions']
                                            and r['lean']['restarts'] == r['matched']['restarts'] for r in rs))
        summary.append(cell_row)
        print(f'{cell:6} {len(rs):>2} {fmt_band(cell_row["ir_ratio"]):>26} {fmt_band(cell_row["rho_solve_ratio"]):>26} '
              f'{med_m:>8g} {med_l:>8g} {fmt_band(cell_row["iteration_ratio"]):>26} '
              f'{per_add["matched"]:>8.0f} {per_add["lean"]:>8.0f}', flush=True)
    if all_logs:
        print(f'\nall {len(all_logs)} fixtures: Ir lean/matched {fmt_band(band(all_logs))}')
    same = all(c['identical_walks'] for c in summary)
    print('walks: iterations, walk additions and restarts identical on every fixture' if same else
          'walks: NOT identical on every fixture (see the iteration ratio column)')

    ic_rows = []
    if args.ic_check:
        print(f'\nIC arm unchanged? first {args.ic_check} fixture(s) a cell, mode ic, pair_table, 3 summands')
        for cell in cells:
            for case in cases_for(fixtures, cell, args.per_cell)[:args.ic_check]:
                job = dict(case['job'], mode='ic', config=IC)
                seen = {}
                for arm, worker in workers.items():
                    report, total, phases = callgrind(worker, job, env)
                    try:
                        verify(report, case['fixture'], expected_mode='ic', summands=3)
                        good = True
                    except Exception as err:
                        print(f'  {case["label"]} {arm}: ORACLE REJECTED: {err}')
                        good = False
                        ok = False
                    report.pop('elapsed_seconds', None)
                    seen[arm] = (json.dumps(report, sort_keys=True), total, good, phases)
                identical = seen['matched'][0] == seen['lean'][0]
                pm, pl = seen['matched'][3], seen['lean'][3]
                deltas = {k: pl.get(k, 0) - pm.get(k, 0) for k in sorted(set(pm) | set(pl))
                          if pl.get(k, 0) != pm.get(k, 0)}
                ic_rows.append(dict(cell=cell, case=case['label'], report_identical=identical,
                                    ir_matched=seen['matched'][1], ir_lean=seen['lean'][1],
                                    phase_delta_lean_minus_matched=deltas))
                ok = ok and identical
                print(f'  {case["label"]:28} report identical: {identical}  Ir {seen["matched"][1]:>12,} / '
                      f'{seen["lean"][1]:>12,}  (lean/matched {seen["lean"][1] / seen["matched"][1]:.6f})', flush=True)
                print(f'  {"":28} phase Ir, lean - matched: '
                      + (', '.join(f'{k} {v:+,}' for k, v in deltas.items()) or 'none'), flush=True)

    print('\nEVERY REPORT VERIFIED AND THE ARMS AGREED ON EVERY LOGARITHM' if ok else '\nNOT CLEAN -- see above')
    out = dict(rows=rows, cells=summary, ic_check=ic_rows, clean=ok)
    if args.json:
        Path(args.json).write_text(json.dumps(out, indent=1) + '\n')
    print(json.dumps({'rows': rows}))
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
