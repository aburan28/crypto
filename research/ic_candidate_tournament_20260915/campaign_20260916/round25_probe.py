#!/usr/bin/env python3
"""Before round 0025: the pair arm, the switch and two rhos, on round 0024's
closed fixtures.

    python3 campaign_20260916/round25_probe.py LEAN_WORKER SWITCH_WORKER MATCHED_WORKER [per_cell [cells]]

LEAN_WORKER is the round's incumbent build (`round25-sources/lean`); its rho
mode is the round's rho. MATCHED_WORKER is round 0024's frozen incumbent,
whose rho mode is the matched rho round 0024 was scored against. Every report
is checked by the amended checker the round will run (`evaluator-r25`), at
each arm's own configured summand count. All four must recover the same
logarithm; a fixture any of them fails to finish is reported, not dropped.
The lean and matched rho must also report identical walks: iterations,
walk additions, restarts and the recovered logarithm, per target.

Fixtures are round 0024's development and selection cases, then confirmation
(which alone holds the holdout cells): nothing here draws from round 0025's
seed.
"""
import importlib.util
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

spec = importlib.util.spec_from_file_location('oracle_r25', Path(__file__).resolve().parent / 'evaluator-r25/oracle.py')
oracle = importlib.util.module_from_spec(spec)
spec.loader.exec_module(oracle)

ROUND = ROOT / 'runs/round-0024'
CELLS = ('n13a0', 'n17a1', 'n19a0', 'n19a1', 'n23a0', 'n23a1', 'n29a1', 'n31a0',
         'n37a0', 'n43a1', 'n59a0', 'n61a1')
PAIR = {'batch_trials': 1, 'linear_algebra': 'sparse', 'max_trials': 65536,
        'solver': 'pair_table', 'summands': 3}
SWITCH = dict(PAIR, solver='pair_or_triple', summands=4)
Z = 1.959963985
RATIOS = (('pair/rho', 'pair', 'rho'), ('switch/rho', 'switch', 'rho'),
          ('switch/pair', 'switch', 'pair'), ('rho/matched', 'rho', 'matched'))


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
    lean, switch, matched = (str(Path(p).resolve()) for p in sys.argv[1:4])
    per_cell = int(sys.argv[4]) if len(sys.argv) > 4 else 4
    cells = sys.argv[5].split(',') if len(sys.argv) > 5 else CELLS
    env = T.child_env()
    fixtures = json.loads((ROUND / 'fixtures.json').read_text())
    cases = [dict(c, label=f'{stage}/{c["id"]}')
             for stage in ('development', 'selection', 'confirmation') for c in fixtures[stage]]
    arms = {'pair': (lean, PAIR, 'ic'), 'switch': (switch, SWITCH, 'ic'),
            'rho': (lean, PAIR, 'rho'), 'matched': (matched, PAIR, 'rho')}
    print(f'{per_cell} fixtures a cell from {ROUND.name}; checker evaluator-r25\n')
    print(f'{"cell":7} {"n":>2} {"picks":>15} ' + ' '.join(f'{name:>26}' for name, *_ in RATIOS))
    rows, ok = [], True
    for cell in cells:
        logs = {name: [] for name, *_ in RATIOS}
        picks = set()
        for case in [c for c in cases if c['cell'] == cell][:per_cell]:
            ir, answers, reports = {}, set(), {}
            for arm, (worker, config, mode) in arms.items():
                job = dict(case['job'], config=config, mode=mode)
                report = run(worker, job, env)
                if report.get('status') != 'complete':
                    print(f'  {case["label"]} {arm}: {report.get("status")}')
                    ok = False
                    break
                proof = oracle.verify(report, case['fixture'], expected_mode=mode,
                                      **({} if mode == 'rho' else {'summands': config['summands']}))
                answers.add(json.dumps(proof['solutions'], sort_keys=True))
                reports[arm] = report
                ir[arm] = instructions(worker, job, env)
            else:
                picks.add(reports['switch'].get('collector'))
                if len(answers) != 1:
                    print(f'  {case["label"]}: ARMS DISAGREE ON THE LOGARITHM')
                    ok = False
                    continue
                if json.dumps(reports['rho']['solutions'], sort_keys=True) != json.dumps(reports['matched']['solutions'], sort_keys=True):
                    print(f'  {case["label"]}: lean and matched rho walked differently')
                    ok = False
                rows.append(dict(case=case['label'], cell=cell, collector=reports['switch'].get('collector'), **ir))
                for name, a, b in RATIOS:
                    logs[name].append(math.log(ir[a] / ir[b]))
        cols = []
        for name, *_ in RATIOS:
            if len(logs[name]) >= 2:
                r, lo, hi = band(logs[name])
                cols.append(f'{r:7.3f} [{lo:6.3f}, {hi:6.3f}]')
            else:
                cols.append(f'{"--":>26}')
        print(f'{cell:7} {len(logs["pair/rho"]):>2} {",".join(sorted(p or "?" for p in picks)):>15} '
              + ' '.join(f'{c:>26}' for c in cols), flush=True)
    print('\nEVERY ARM FINISHED AND AGREED ON EVERY FIXTURE' if ok else '\nNOT CLEAN -- see above')
    print(json.dumps({'rows': rows}))
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
