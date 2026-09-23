#!/usr/bin/env python3
"""Is `pairinv` the same computation as `scaled`, and what does it save?

    python3 campaign_20260916/round24_pairinv_check.py SCALED_WORKER PAIRINV_WORKER [per_cell [cells]]

`cells` is a comma list; the default is the four cells the design read first.

`round24-pairinv.patch` changes only how the decomposition scan batches its
field inversions: +P and -P share an abscissa, hence the denominator
`target.x + P.x`, so it is inverted once and read twice.  Nothing about which
rests are computed, in which order, or what is done with them changes -- so the
two workers must return the SAME report on every job.  This script checks that
first, on the whole parsed report except its wall clock (`elapsed_seconds`,
the one field that differs between two runs of the SAME worker), and only then
counts instructions.

Fixtures are round 0023's cases, development and selection first and then
confirmation (which alone holds the holdout cells), from a closed round:
nothing here draws from round 0024's seed.  Both workers run the arm config
round 0023 froze for `scaled`.  Every `pairinv` report is also checked by
`oracle.py` against its fixture.
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
CELLS = ('n13a0', 'n23a1', 'n37a0', 'n43a1')
CLOCK = 'elapsed_seconds'
Z = 1.959963985


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


def main():
    if len(sys.argv) < 3:
        raise SystemExit('usage: round24_pairinv_check.py SCALED_WORKER PAIRINV_WORKER [per_cell]')
    scaled, pairinv = (str(Path(p).resolve()) for p in sys.argv[1:3])
    per_cell = int(sys.argv[3]) if len(sys.argv) > 3 else 6
    cells = sys.argv[4].split(',') if len(sys.argv) > 4 else CELLS
    env = T.child_env()
    arms = json.loads((ROUND / 'candidates.json').read_text())
    config = next(a for a in arms if a['id'] == 'scaled')['config']
    fixtures = json.loads((ROUND / 'fixtures.json').read_text())
    cases = [dict(c, label=f'{stage}/{c["id"]}')
             for stage in ('development', 'selection', 'confirmation') for c in fixtures[stage]]
    print(f'config {json.dumps(config, sort_keys=True)}; {per_cell} fixtures a cell '
          f'from {ROUND.name}\n')
    print(f'{"cell":7} {"n":>3} {"same":>5} {"pairinv/scaled":>15} {"95% band":>18} '
          f'{"min":>7} {"max":>7}')
    ok = True
    for cell in cells:
        chosen = [c for c in cases if c['cell'] == cell][:per_cell]
        same, logs = 0, []
        for case in chosen:
            job = dict(case['job'], config=config, mode='ic')
            a, b = run(scaled, job, env), run(pairinv, job, env)
            if a.get('status') != 'complete' or b.get('status') != 'complete':
                print(f'  {case["label"]}: incomplete ({a.get("status")}, {b.get("status")})')
                ok = False
                continue
            a.pop(CLOCK, None), b.pop(CLOCK, None)
            if json.dumps(a, sort_keys=True) != json.dumps(b, sort_keys=True):
                keys = sorted(k for k in set(a) | set(b) if a.get(k) != b.get(k))
                print(f'  {case["label"]}: REPORTS DIFFER in {keys}')
                ok = False
                continue
            verify(b, case['fixture'], expected_mode='ic', summands=config['summands'])
            same += 1
            logs.append(math.log(instructions(pairinv, job, env) / instructions(scaled, job, env)))
        if len(logs) < 2:
            print(f'{cell:7} {len(chosen):>3} {same:>5}   NO RATIO')
            continue
        ratio, err = math.exp(st.mean(logs)), st.stdev(logs) / math.sqrt(len(logs))
        lo, hi = ratio * math.exp(-Z * err), ratio * math.exp(Z * err)
        print(f'{cell:7} {len(chosen):>3} {same:>5} {ratio:>15.4f} '
              f'{f"[{lo:.4f}, {hi:.4f}]":>18} {math.exp(min(logs)):>7.4f} {math.exp(max(logs)):>7.4f}')
    print('\nIDENTICAL REPORTS ON EVERY JOB' if ok else '\nNOT EQUIVALENT -- see above')
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
