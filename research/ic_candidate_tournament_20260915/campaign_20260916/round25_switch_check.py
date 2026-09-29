#!/usr/bin/env python3
"""The switch collector against the two collectors it chooses between.

    python3 campaign_20260916/round25_switch_check.py SWITCH_WORKER \
        [PAIR_WORKER [TRIPLE_WORKER [per_cell [cells]]]]

SWITCH_WORKER is `round24-sources/counted` + `round25-switch.patch`, run as
`solver: pair_or_triple`, `summands: 4`. PAIR_WORKER defaults to round 0024's
frozen incumbent (`runs/round-0024/worker`, the `matched` tree: pair table,
three summands). TRIPLE_WORKER defaults to round 0024's frozen `counted`
worker (`solver: triple_counted`, four summands).

Per fixture:

1. The switch's report names the collector it ran (`collector`). With
   `elapsed_seconds`, `collector` and `collector_model` removed, it must equal
   that collector's report from its own worker with `elapsed_seconds` removed.
2. The switch's report verifies under the amended checker
   (`evaluator-r25/oracle.py`) at the arm's configured `summands: 4`, and
   recovers the same logarithm as both collectors.
3. Instructions, callgrind `Collected:`, for all three workers, and where
   the switch ran the pair table, for the triple worker's own tree run as
   `pair_table` at three summands (the same `Collector::Pair` code the
   switch runs, without the switch). The workers
   are copied to one directory under names of equal length first, so the
   start-up's path handling costs the same in each.

Fixtures are round 0024's closed ones: development then selection, then
confirmation (which alone holds the holdout cells), `per_cell` (default 4) a
cell. No round-0024 measurement is used by the switch: its rule is fixed in
the source (`koblitz_tiny_ic::choose_collector`) and reads only `(n, r,
points)`.
"""
import hashlib
import json
import math
import re
import shutil
import statistics as st
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE / 'evaluator-r25'))
import oracle as NEW  # noqa: E402  the amended checker

ROUND = ROOT / 'runs/round-0024'
CELLS = ('n13a0', 'n17a1', 'n19a0', 'n19a1', 'n23a0', 'n23a1', 'n29a1', 'n31a0',
         'n37a0', 'n43a1', 'n59a0', 'n61a1')
BASE = {'batch_trials': 1, 'linear_algebra': 'sparse', 'max_trials': 65536}
PAIR = dict(BASE, solver='pair_table', summands=3)
TRIPLE = dict(BASE, solver='triple_counted', summands=4)
SWITCH = dict(BASE, solver='pair_or_triple', summands=4)
# Round 0024's confirmation, triple/pair in instructions (RESULTS.md). Used
# only to score the choice after the fact, never by the switch.
R24_TRIPLE_OVER_PAIR = dict(n13a0=1.967, n17a1=2.205, n19a0=2.327, n19a1=2.115, n23a0=1.542,
                            n23a1=1.323, n29a1=3.835, n31a0=3.260, n37a0=0.524, n43a1=0.291,
                            n59a0=0.648, n61a1=0.314)
SWITCH_ONLY = ('elapsed_seconds', 'collector', 'collector_model')
Z = 1.959963985


def child_env():
    # As the frozen tournament's `child_env`.
    import os
    keep = ('PATH', 'LANG', 'LC_ALL', 'TZ', 'CARGO_HOME', 'RUSTUP_HOME', 'HOME')
    env = {k: os.environ[k] for k in keep if k in os.environ}
    env.update(RAYON_NUM_THREADS='1', IC_ARTIFACT_CACHE='off', IC_F2_BACKEND='cpu', PYTHONHASHSEED='0')
    return env


def run(worker, job, env, timeout=1800):
    done = subprocess.run([worker], input=json.dumps(job), text=True,
                          capture_output=True, env=env, timeout=timeout)
    return json.loads(done.stdout) if done.stdout else {}


def instructions(worker, job, env, timeout=7200):
    done = subprocess.run(
        ['valgrind', '--tool=callgrind', '--callgrind-out-file=/dev/null', worker],
        input=json.dumps(job), text=True, capture_output=True, env=env, timeout=timeout)
    found = re.search(r'Collected\s*:\s*(\d+)', done.stderr)
    if not found:
        raise RuntimeError(f'no Ir collected: {done.stderr[-300:]}')
    return int(found.group(1))


def band(logs):
    ratio = math.exp(st.mean(logs))
    err = st.stdev(logs) / math.sqrt(len(logs)) if len(logs) > 1 else float('nan')
    return ratio, ratio * math.exp(-Z * err), ratio * math.exp(Z * err)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def check(case, workers, env):
    """One fixture: equivalence, verification, and three instruction counts."""
    job = lambda config: dict(case['job'], mode='ic', config=config)
    switch = run(workers['switch'], job(SWITCH), env)
    out = {'case': case['label'], 'cell': case['cell'], 'problems': []}
    if switch.get('status') != 'complete':
        out['problems'].append(f'switch: {switch.get("status")} {switch.get("reason", "")}')
        return out
    chosen = switch['collector']
    out['collector'] = chosen
    out['model'] = switch['collector_model']
    reference = run(workers['pair' if chosen == 'pair_table' else 'triple'],
                    job(PAIR if chosen == 'pair_table' else TRIPLE), env)
    stripped = {k: v for k, v in switch.items() if k not in SWITCH_ONLY}
    ref = {k: v for k, v in reference.items() if k != 'elapsed_seconds'}
    out['identical'] = stripped == ref
    if not out['identical']:
        diff = sorted(k for k in set(stripped) | set(ref) if stripped.get(k) != ref.get(k))
        out['problems'].append(f'report differs from {chosen} in {diff}')
    try:
        proof = NEW.verify(switch, case['fixture'], expected_mode='ic', summands=4)
        out['summands'] = switch['summands']
        out['relation_summands'] = proof.get('relation_summands', 4)
    except NEW.InvalidEvidence as exc:
        out['problems'].append(f'switch report refused: {exc}')
        return out
    other = run(workers['triple' if chosen == 'pair_table' else 'pair'],
                job(TRIPLE if chosen == 'pair_table' else PAIR), env)
    answers = {json.dumps(proof['solutions'])}
    for name, report in ((chosen, reference), ('the other collector', other)):
        try:  # each collector at its own configured count
            answers.add(json.dumps(NEW.verify(report, case['fixture'], expected_mode='ic',
                                              summands=report.get('summands'))['solutions']))
        except NEW.InvalidEvidence as exc:
            out['problems'].append(f'{name} refused: {exc}')
    out['same_logarithm'] = len(answers) == 1
    if not out['same_logarithm']:
        out['problems'].append('collectors disagree on the logarithm')
    for arm, config in (('switch', SWITCH), ('pair', PAIR), ('triple', TRIPLE)):
        out['ir_' + arm] = instructions(workers[arm], job(config), env)
    if chosen == 'pair_table':
        # The switch's pair path is the counted tree's `Collector::Pair`, not the
        # matched tree's binary: the same tree's own pair run separates what the
        # switch adds from what the two trees differ by.
        out['ir_counted_tree_pair'] = instructions(workers['triple'], job(PAIR), env)
    return out


def main():
    if len(sys.argv) < 2:
        raise SystemExit(__doc__)
    given = [Path(p).resolve() if p else None for p in sys.argv[1:4]] + [None, None]  # '' = default
    sources = {'switch': given[0],
               'pair': given[1] or ROUND / 'worker',
               'triple': given[2] or ROUND / 'source_candidates/counted/worker'}
    per_cell = int(sys.argv[4]) if len(sys.argv) > 4 else 4
    cells = sys.argv[5].split(',') if len(sys.argv) > 5 else CELLS
    env = child_env()
    fixtures = json.loads((ROUND / 'fixtures.json').read_text())
    cases = [dict(c, label=f'{stage}/{c["id"]}')
             for stage in ('development', 'selection', 'confirmation') for c in fixtures[stage]]
    print('workers (sha256):')
    for arm, path in sources.items():
        print(f'  {arm:6} {sha(path)}  {path}')
    print(f'{per_cell} fixtures a cell from {ROUND.name} (development, selection, then confirmation)\n')
    rows, ok = [], True
    with tempfile.TemporaryDirectory() as tmp:
        workers = {}
        for arm, path in sources.items():
            workers[arm] = str(Path(tmp) / f'w-{arm[:2]}')  # equal-length names
            shutil.copy2(path, workers[arm])
        chosen = [c for cell in cells for c in [x for x in cases if x['cell'] == cell][:per_cell]]
        with ThreadPoolExecutor(max_workers=3) as pool:
            results = list(pool.map(lambda c: check(c, workers, env), chosen))
    print(f'{"cell":7} {"n":>2} {"collector":>15} {"W3/W2":>6} {"identical":>9} {"verified@4":>10} '
          f'{"switch/pair, instr":>26} {"switch/triple, instr":>26} {"switch/chosen":>13} {"switch/same-tree":>16} {"r24 better":>10} {"match":>5}')
    summary = {}
    for cell in cells:
        mine = [r for r in results if r['cell'] == cell]
        for r in mine:
            for p in r['problems']:
                print(f'  {r["case"]}: {p}')
                ok = False
        good = [r for r in mine if not r['problems']]
        picks = {r.get('collector') for r in mine}
        pick = picks.pop() if len(picks) == 1 else 'MIXED'
        if pick == 'MIXED':
            ok = False
        model = good[0]['model'] if good else {}
        q = model['triple_work'] / model['pair_work'] if model else float('nan')
        cols = []
        for other in ('pair', 'triple'):
            logs = [math.log(r['ir_switch'] / r['ir_' + other]) for r in good]
            if len(logs) >= 2:
                ratio, lo, hi = band(logs)
                cols.append(f'{ratio:7.4f} [{lo:6.4f}, {hi:6.4f}]')
            else:
                cols.append(f'{"--":>26}')
        own = [r['ir_switch'] / r['ir_' + ('pair' if r['collector'] == 'pair_table' else 'triple')] for r in good]
        tree = [r['ir_switch'] / r.get('ir_counted_tree_pair', r['ir_triple']) for r in good]
        better = 'triple' if R24_TRIPLE_OVER_PAIR[cell] < 1 else 'pair'
        agrees = (pick == 'triple_counted') == (better == 'triple')
        ok &= agrees or pick == 'MIXED'
        print(f'{cell:7} {len(good):>2} {pick:>15} {q:6.3f} {sum(r["identical"] for r in good):>5}/{len(mine):<3} '
              f'{sum(1 for r in good):>6}/{len(mine):<3} {cols[0]:>26} {cols[1]:>26} '
              f'{min(own) if own else float("nan"):.4f}-{max(own) if own else float("nan"):.4f} '
              f'{min(tree) if tree else float("nan"):>7.4f}-{max(tree) if tree else float("nan"):.4f} '
              f'{better:>10} {"yes" if agrees else "NO":>5}')
        summary[cell] = {'collector': pick, 'model': model, 'W3_over_W2': q, 'r24_better': better,
                         'matches_r24_better_of_two': agrees}
        rows.extend(mine)
    total = len(rows)
    identical = sum(1 for r in rows if r.get('identical'))
    verified = sum(1 for r in rows if 'relation_summands' in r)
    print(f'\n{identical}/{total} switch reports equal the chosen collector\'s report '
          f'(apart from {", ".join(SWITCH_ONLY)}); {verified}/{total} verify under the amended checker at summands 4')
    print('EVERY FIXTURE EQUIVALENT, VERIFIED AND AGREED; EVERY CHOICE MATCHES ROUND 0024\'S BETTER OF TWO'
          if ok else 'NOT CLEAN -- see above')
    print(json.dumps({'summary': summary, 'rows': rows}))
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
