#!/usr/bin/env python3
"""Before round 0024 as amended: the pair arm, the triple arm and two rhos,
on round 0023's closed fixtures.

    python3 campaign_20260916/round24_amended_probe.py \
        MATCHED_WORKER COUNTED_WORKER OLD_WORKER [per_cell [cells]]

MATCHED_WORKER is the round's incumbent build (`matched`), whose rho mode is
the round's rho: Frobenius classes named by the IC arm's normal-basis
rotation. OLD_WORKER is round 0023's frozen `scaled` worker, whose rho mode
is the rho every earlier round was scored against. COUNTED_WORKER is the
triple-sum collector with counted sizing.

Per fixture: the pair arm (`matched`, pair_table, 3 summands), the triple arm
(`counted`, triple_counted, 4 summands), matched rho and old rho. Every
report goes through `oracle.py`; all four must recover the same logarithm,
and a fixture any of them fails to finish is reported rather than dropped
silently. Instructions are callgrind `Collected:`.

Fixtures are round 0023's development and selection cases, then confirmation
(which alone holds the holdout cells): nothing here draws from round 0024's
seed.
"""
import json
import math
import re
import statistics as st
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
import tournament as T  # noqa: E402
from oracle import verify  # noqa: E402

ROUND = ROOT / 'runs/round-0023'
CELLS = ('n13a0', 'n17a1', 'n19a0', 'n19a1', 'n23a0', 'n23a1', 'n29a1', 'n31a0',
         'n37a0', 'n43a1', 'n59a0', 'n61a1')
PAIR = {'batch_trials': 1, 'linear_algebra': 'sparse', 'max_trials': 65536,
        'solver': 'pair_table', 'summands': 3}
COUNTED = dict(PAIR, solver='triple_counted', summands=4)
Z = 1.959963985
RATIOS = (('counted/pair', 'counted', 'pair'), ('pair/rho', 'pair', 'rho'),
          ('counted/rho', 'counted', 'rho'), ('rho/old_rho', 'rho', 'old_rho'))


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


def main():
    if len(sys.argv) < 4:
        raise SystemExit(__doc__)
    matched, counted, old = (str(Path(p).resolve()) for p in sys.argv[1:4])
    per_cell = int(sys.argv[4]) if len(sys.argv) > 4 else 4
    cells = sys.argv[5].split(',') if len(sys.argv) > 5 else CELLS
    env = T.child_env()
    fixtures = json.loads((ROUND / 'fixtures.json').read_text())
    cases = [dict(c, label=f'{stage}/{c["id"]}')
             for stage in ('development', 'selection', 'confirmation') for c in fixtures[stage]]
    arms = {'pair': (matched, PAIR, 'ic', 3), 'counted': (counted, COUNTED, 'ic', 4),
            'rho': (matched, PAIR, 'rho', None), 'old_rho': (old, PAIR, 'rho', None)}
    print(f'{per_cell} fixtures a cell from {ROUND.name}\n')
    print(f'{"cell":7} {"n":>2} ' + ' '.join(f'{name:>26}' for name, *_ in RATIOS))
    rows, ok = [], True
    for cell in cells:
        logs = {name: [] for name, *_ in RATIOS}
        chosen = [c for c in cases if c['cell'] == cell][:per_cell]
        for case in chosen:
            ir, answers = {}, set()
            for arm, (worker, config, mode, summands) in arms.items():
                job = dict(case['job'], config=config, mode=mode)
                report = run(worker, job, env)
                if report.get('status') != 'complete':
                    print(f'  {case["label"]} {arm}: {report.get("status")}')
                    ok = False
                    break
                proof = (verify(report, case['fixture'], expected_mode='rho') if mode == 'rho' else
                         verify(report, case['fixture'], expected_mode='ic', summands=summands))
                answers.add(json.dumps(proof['solutions'], sort_keys=True))
                ir[arm] = instructions(worker, job, env)
            else:
                if len(answers) != 1:
                    print(f'  {case["label"]}: ARMS DISAGREE ON THE LOGARITHM')
                    ok = False
                    continue
                rows.append(dict(case=case['label'], cell=cell, **ir))
                for name, a, b in RATIOS:
                    logs[name].append(math.log(ir[a] / ir[b]))
        cols = []
        for name, *_ in RATIOS:
            if len(logs[name]) >= 2:
                r, lo, hi = band(logs[name])
                cols.append(f'{r:7.3f} [{lo:6.3f}, {hi:6.3f}]')
            else:
                cols.append(f'{"--":>26}')
        print(f'{cell:7} {len(logs["pair/rho"]):>2} ' + ' '.join(f'{c:>26}' for c in cols), flush=True)
    print('\nEVERY ARM FINISHED AND AGREED ON EVERY FIXTURE' if ok else '\nNOT CLEAN -- see above')
    print(json.dumps({'rows': rows}))
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
