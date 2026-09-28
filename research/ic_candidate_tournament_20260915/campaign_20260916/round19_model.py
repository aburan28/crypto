#!/usr/bin/env python3
"""The solver's own cost model over factor-base size, on the base it really builds.

    python3 campaign_20260916/round19_model.py

`koblitz_tiny_ic.rs` publishes its cost model and minimises it at run time:

    W(t) = t*size + (K + TABLE_PROBES) / (lambda * cov(t))
    lambda = size*(size+1) / (2r)        cov(t) = 1 - ((K-t)/K)^2

`size` is the number of signed base points and `K` the column count after the
signed-Frobenius orbits are folded, so `size = 2*n*K` for a base of `K` orbits
on a degree-`n` curve.  The first term builds the pair table; the second scans
for relations, and it carries `r` in the numerator through `lambda`.

**What `K` actually is.**  The contract asks for `6*degree` points, and every
round from 0002 on has been read as running with that base.  It does not:
`build_subgroup_orbit_factor_base` samples orbits in batches of eight and stops
at the first rebuild that reaches the target, so any request at or below one
batch returns that whole batch.  Measured on the panel (`round19_base_sweep.py`),
`points=1` and `points=368` both return eight orbits at degree 23, and the
contract's 138 returns eight as well.  **The panel's factor base is seven or
eight orbits at every cell, `size = 16*degree` signed points, and the contract's
`points` has been inert since round 0002** -- not a chosen size, a batch
boundary.  The tables below are therefore computed over `K`, the orbit count,
which is the parameter that exists, and the achievable values of `K` are
multiples of the batch: a base can be enlarged, never shrunk below one batch.

**The boundary.**  Treating `K` as free, the scan term is
`2r(K+3) / (4 n^2 K^2 cov(t))`.  At fixed `t` and large `K`, `cov(t) -> 2t/K`
and the scan term tends to `r/(4 n^2 t)`: growing the base buys nothing,
because coverage thins exactly as fast as the pair density grows, while the
table cost `2 n t K` grows linearly.  Progress needs `t` to grow with `K`; at
`t = aK` the table costs `2 n a K^2` and the scan `~ c r / K`, so

    K* ~ r^(1/3) / n,   size* ~ r^(1/3),   W* ~ r^(2/3)

against rho's `Theta((r/2n)^(1/2))`.  **At its own optimal factor-base size
this method is `r^(2/3)` where rho is `r^(1/2)`: the ratio grows as `r^(1/6)`
and must eventually exceed one.**  That is a boundary in the AGENTS.md sec. 1
sense -- derived, not measured, and not movable by tuning constants.  It is a
statement about this pair-table collector, not about index calculus.

No measurement here and no runs; the units are the solver's own and only
ratios within a cell are meaningful.  `round19_base_sweep.py` measures it.
"""
import json
from pathlib import Path

TABLE_PROBES = 3.0  # koblitz_tiny_ic.rs
BATCH = 8           # build_subgroup_orbit_factor_base
CELLS = [('n13a0', 13, 2003), ('n29a1', 29, 42457), ('n17a1', 17, 65587),
         ('n19a0', 19, 130873), ('n19a1', 19, 262543), ('n31a0', 31, 1439393),
         ('n23a0', 23, 2095853), ('n23a1', 23, 4196903)]
# What the panel really runs: one batch, less any orbit the sampler rejected.
# Measured per cell by round19_base_sweep.py; 8 everywhere but n13a0, where the
# closure at degree 13 gives 7.
BUILT = {'n13a0': 7, 'n17a1': 8, 'n19a0': 8, 'n19a1': 8,
         'n23a0': 8, 'n23a1': 8, 'n29a1': 8, 'n31a0': 8}
KS = [8, 16, 24, 32, 40, 56, 72, 104, 152, 216]


def work(columns, degree, r):
    """(W, t) for a base of `columns` orbits, with the row count the solver picks."""
    size = 2 * degree * columns
    lam = size * (size + 1) / 2.0 / r
    best = None
    for t in range(1, columns + 1):
        rest = (columns - t) / columns
        w = t * size + (columns + TABLE_PROBES) / max(lam * (1 - rest * rest), 1e-12)
        if best is None or w < best[0]:
            best = (w, t)
    return best


def main():
    print('## What the panel runs, against the continuous optimum\n')
    print('| cell | r | orbits built | size | r^(1/3) | size/optimal | model K* | W(K*)/W(built) |')
    print('|:--|--:|--:|--:|--:|--:|--:|--:|')
    out = {}
    for cell, degree, r in CELLS:
        built = BUILT[cell]
        base = work(built, degree, r)[0]
        grid = sorted({built} | set(KS))
        ws = {k: work(k, degree, r) for k in grid}
        star = min(ws, key=lambda k: ws[k][0])
        size = 2 * degree * built
        print(f'| `{cell}` | {r:,} | {built} | {size} | {r ** (1 / 3):.1f} '
              f'| {size / r ** (1 / 3):.2f} | {star} | {ws[star][0] / base:.3f} |')
        out[cell] = {'orbits_built': built, 'size': size, 'K_star': star,
                     'W_ratio_at_K_star': round(ws[star][0] / base, 4),
                     't_at_built': ws[built][1], 't_at_K_star': ws[star][1],
                     'curve': {k: round(ws[k][0] / base, 4) for k in grid}}

    print('\n## W(K)/W(built), in the solver\'s own units\n')
    head = ' '.join(f'{("K=" + str(k)):>7}' for k in KS)
    print(f"{'cell':7} {'r':>9} {'built':>5} {head}  {'K*':>4} {'t@built':>7} {'t@K*':>5}")
    for cell, degree, r in CELLS:
        c = out[cell]
        print(f'{cell:7} {r:9d} {c["orbits_built"]:5d} ' +
              ' '.join(f'{c["curve"][k]:7.3f}' for k in KS) +
              f'  {c["K_star"]:4d} {c["t_at_built"]:7d} {c["t_at_K_star"]:5d}')

    print('\nThe base can only be enlarged -- one batch is the floor -- and at every cell on')
    print('this panel enlarging it is a loss or a wash.  There is no factor-base lever here,')
    print('which is why round 0019 carries no new arm.  The `size/optimal` column is the')
    print('crossover mechanism: it falls monotonically with r, because `size` grows linearly')
    print('in the degree while `r` grows exponentially in it.')
    Path(__file__).with_name('round-0019-base-size-model.json').write_text(
        json.dumps(out, indent=1) + '\n')


if __name__ == '__main__':
    main()
