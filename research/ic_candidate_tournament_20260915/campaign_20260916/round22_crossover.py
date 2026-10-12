#!/usr/bin/env python3
"""The crossover again, with rho actually completing.

    python3 campaign_20260916/round22_crossover.py WORKER [fixtures]

**This supersedes round 0021's `n37a0` number.**  `round21_crossover.py`
checked that the IC arm completed and then measured both arms:

    report = run(worker, dict(job, mode='ic'), env)
    if report.get('status') != 'complete':
        continue
    ...
    rho_ir.append(instructions(worker, dict(job, mode='rho'), env))

`instructions()` parses callgrind's `Collected:` line and never sees the
worker's JSON status, so a rho run cut off at `max_trials` contributed its
TRUNCATED instruction count as the denominator.  On round 0021's own draw,
**4 of the 64 rho runs at `n37a0` did not complete** (64 of 64 did at
`n23a1`).  An under-charged denominator inflates IC/rho, so the published
1.533 [1.238, 1.899] is an upper-biased reading and the correction runs
downward.  `n23a1`'s 0.831 is unaffected.

That is the campaign's own rule -- an arm that does not complete has no cost,
and a budget exhaustion is never evidence about an algorithm.  Round 0021
enforced it at `n41a0`, by hand, and missed it at `n37a0`.

**Dropping the four fixtures is not the fix.**  They are exactly the fixtures
where rho is expensive, so excluding them biases the ratio down, in IC's
favour -- the same error with the opposite sign.  The fix is to let rho finish.

## The budget, and why this is a re-scoped protocol rather than an amendment

`max_trials` was 4096 in the frozen config.  It is not a ceiling: the worker
itself accepts up to 65,536 (`ic_tournament_worker.rs`, "invalid collection
limits").  rho needs about 2,200 steps at `n37a0` and about 9,200 at `n43a1`,
which is why `n37a0` sat just under 4096 and looked like the boundary.

It would be convenient if raising it were INERT wherever both arms already
finished -- a cap that is never reached cannot be observed, so every published
measurement would read unchanged and this would be an additive amendment.
**That was predicted and it is false.**  `round22_budget_effect.py` measures
identical answers and identical factor bases at every cell, and a DIFFERENT
instruction count at every cell on both arms, because

    max_iterations_per_restart: cfg.max_trials,

makes `max_trials` rho's iterations-per-RESTART -- how often it abandons a walk
and starts over, which is a parameter of the algorithm being measured -- and
the candidate takes it at `TinyIc::new`.  It is protocol, not a safety valve.

So this ladder is measured at ONE cap throughout and is not interchangeable
with the panel measured at 4096.  `n23a1` is carried at BOTH caps as the
anchor: it is the cell where the published panel, round 0021 and this ladder
all meet, and the two readings of it are what any cross-cap statement has to
rest on.  Note also what the change cannot do -- a longer restart interval can
only help rho -- so it cannot flatter the candidate.

**Round 0021's four truncated fixtures bias the ratio DOWNWARD, which is the
opposite of what this script's first version asserted.**  The assertion was
that a truncated rho is charged too little, raising the ratio.  It is charged
too MUCH: an `incomplete` rho has exhausted its restarts, so it did a great
deal of work and returned nothing, where a completed run often finds its
collision early.  1.533 is therefore a lower bound on `n37a0`, and this script
measures 1.763 [1.547, 2.009] with every rho run required to finish.

## What this measures

Three cells, all certified by `oracle.py`, spanning three decades of `r`:

    n23a1   r = 4,196,903          the panel's largest cell
    n37a0   r = 230,603,167        round 0021's second point, corrected
    n43a1   r = 4,644,189,029      new -- a third point, so the rate gets
                                   curvature instead of being a two-point line

A fixture counts only if BOTH arms complete and the IC report verifies.  A cell
that loses any fixture prints no ratio at all: a partial panel is the bias this
script exists to remove, not a result to report with a caveat.
"""
import json
import math
import random
import re
import statistics as st
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
import tournament as T  # noqa: E402
from oracle import verify  # noqa: E402

FROZEN = json.loads((ROOT / 'runs/round-0020/candidates.json').read_text())[0]['config']
# The one difference from the frozen config, and it is a protocol difference.
CONFIG = dict(FROZEN, max_trials=65536)
# (label, degree, curve_a, r, base points, configs to measure at).  `n23a1` is
# carried at both caps: it is the only cell where the published panel and this
# ladder can be compared, so it is the anchor rather than just the low point.
CELLS = [('n23a1', 23, 1, 4_196_903, 138, ('frozen', 'raised')),
         ('n37a0', 37, 0, 230_603_167, 222, ('raised',)),
         ('n43a1', 43, 1, 4_644_189_029, 222, ('raised',))]
SEEDS = (20260922, 4242)
Z = 1.959963985


def run(worker, job, env, timeout=1800):
    done = subprocess.run([worker], input=json.dumps(job), text=True,
                          capture_output=True, env=env, timeout=timeout)
    return json.loads(done.stdout) if done.stdout else {}


def instructions(worker, job, env, timeout=3600):
    done = subprocess.run(
        ['valgrind', '--tool=callgrind', '--callgrind-out-file=/dev/null', worker],
        input=json.dumps(job), text=True, capture_output=True, env=env, timeout=timeout)
    found = re.search(r'Collected\s*:\s*(\d+)', done.stderr)
    if not found:
        raise RuntimeError(f'no Ir collected: {done.stderr[-300:]}')
    return int(found.group(1))


def main():
    if len(sys.argv) < 2:
        raise SystemExit('usage: round22_crossover.py WORKER [fixtures]')
    worker = str(Path(sys.argv[1]).resolve())
    fixtures = int(sys.argv[2]) if len(sys.argv) > 2 else 64
    env = T.child_env()

    probe = run(worker, {'mode': 'fixture', 'degree': 43, 'curve_a': 1, 'target_seeds': [1],
                         'algorithm_seed': 1, 'config': CONFIG,
                         'factor_base': {'kind': 'subgroup_orbits', 'seed': 43, 'points': 222}}, env)
    if probe.get('status') != 'fixture':
        raise SystemExit(f'this worker refuses degree 43 ({probe.get("reason")!r}); '
                         'apply round21-wide-pair-table.patch and rebuild')

    caps = {'frozen': FROZEN, 'raised': CONFIG}
    print(f'max_trials: frozen {FROZEN["max_trials"]}, raised {CONFIG["max_trials"]}; '
          f'{fixtures} fixtures a cell over seeds {SEEDS}')
    print('a fixture counts only if BOTH arms complete and both reports verify\n')
    print(f'{"cell":8} {"r":>16} {"cap":>7} {"IC/rho":>8} {"sd(log)":>8} '
          f'{"95% band":>18} {"complete":>10}')
    table = {}
    for label, degree, curve_a, r, points, tags in CELLS:
        for tag in tags:
            cfg = caps[tag]
            ic_ir, rho_ir, dropped = [], [], []
            for stream in SEEDS:
                rng = random.Random(stream)
                for _ in range(fixtures // len(SEEDS)):
                    job = {'degree': degree, 'curve_a': curve_a,
                           'target_seeds': [rng.getrandbits(64)],
                           'algorithm_seed': rng.getrandbits(64), 'config': cfg,
                           'factor_base': {'kind': 'subgroup_orbits', 'seed': 43,
                                           'points': points}}
                    fixture = run(worker, dict(job, mode='fixture'), env)
                    if fixture.get('status') != 'fixture':
                        dropped.append('no fixture')
                        continue
                    report = run(worker, dict(job, mode='ic'), env)
                    rho = run(worker, dict(job, mode='rho'), env)
                    # Both, before either is charged. This is the check round
                    # 0021 made for IC and not for rho.
                    if report.get('status') != 'complete':
                        dropped.append('IC incomplete')
                        continue
                    if rho.get('status') != 'complete':
                        dropped.append('rho incomplete')
                        continue
                    verify(report, fixture['fixture'], expected_mode='ic', summands=3)
                    verify(rho, fixture['fixture'], expected_mode='rho')
                    ic_ir.append(instructions(worker, dict(job, mode='ic'), env))
                    rho_ir.append(instructions(worker, dict(job, mode='rho'), env))
            n = len(ic_ir)
            if dropped or n < 2:
                # Not a caveat to print under a number. The fixtures that drop
                # are the ones where an arm is expensive, so a partial panel is
                # biased toward whichever arm survived -- which is the defect
                # being corrected here, reappearing with a different sign.
                reasons = {x: dropped.count(x) for x in sorted(set(dropped))}
                print(f'{label:8} {r:>16,} {tag:>7}   NO RATIO -- {reasons} '
                      f'({n} usable)')
                continue
            logs = [math.log(a / b) for a, b in zip(ic_ir, rho_ir)]
            ratio, spread = math.exp(st.mean(logs)), st.stdev(logs)
            err = spread / math.sqrt(n)
            lo, hi = ratio * math.exp(-Z * err), ratio * math.exp(Z * err)
            table[(label, tag)] = (r, ratio, lo, hi)
            print(f'{label:8} {r:>16,} {tag:>7} {ratio:>8.3f} {spread:>8.3f} '
                  f'{f"[{lo:.3f}, {hi:.3f}]":>18} {f"{n}/{fixtures}":>10}')

    anchor = [table.get(('n23a1', tag)) for tag in ('frozen', 'raised')]
    if all(anchor):
        (_, frozen_v, *_), (_, raised_v, *_) = anchor
        print(f'\nanchor: n23a1 reads {frozen_v:.3f} at the frozen cap and '
              f'{raised_v:.3f} at the raised one')
        print(f'  the cap is worth {raised_v / frozen_v:.3f}x in the RATIO at this cell.')
        print('  That number, not an argument about caps, is what relates the ladder')
        print('  below to the panel published at the frozen cap.')

    ladder = [(cell[0], table[(cell[0], 'raised')])
              for cell in CELLS if (cell[0], 'raised') in table]
    if len(ladder) >= 2:
        print('\nrate between consecutive cells, all at the raised cap')
        for (c0, (r0, v0, *_)), (c1, (r1, v1, *_)) in zip(ladder, ladder[1:]):
            print(f'  {c0} -> {c1}: r^{math.log(v1 / v0) / math.log(r1 / r0):.3f} '
                  f'({r1 / r0:.0f}x in r, {v1 / v0:.2f}x in the ratio)')
        print('round19_model.py derives r^(1/6) = r^0.167 for the balanced optimum.')
        if len(ladder) >= 3:
            print('Three cells give the first CURVATURE this campaign has measured:')
            print('two consecutive rates that can disagree with each other.')
    print('\nEvery ratio above is between two arms that both completed. A fixture')
    print('where either did not is not a cost and not evidence; it drops the cell.')


if __name__ == '__main__':
    main()
