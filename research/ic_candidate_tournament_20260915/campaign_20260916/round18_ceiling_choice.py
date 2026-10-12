#!/usr/bin/env python3
"""Choose the scan-block clamp ceiling from the grid measurement, on the
development cells only.

    python3 campaign_20260916/round18_ceiling_choice.py [raw_rows.json]

Reads the raw rows `round18_block_grid.py` wrote (default: its output path
under the system temp directory) and, for each candidate ceiling, works out
what block each cell would then use -- `min(rule_value, ceiling)`, where the
rule's own value is the block the shipped worker reports -- and what that block
costs against the shipped one.

A ceiling never raises a cell's block above what the rule asks for, so a cell
whose rule value is already at or below the ceiling is untouched and scores
exactly 1. The choice is made on the six development cells; the two holdout
cells are printed but take no part in it.
"""
import collections
import json
import statistics
import sys
import tempfile
from pathlib import Path

DEFAULT = Path(tempfile.gettempdir()) / 'round18_block_grid.json'
HOLDOUT = {'n19a1', 'n29a1'}
# The block each cell's own rule asks for, as the shipped worker reports it.
# n23a0 and n23a1 are the two the shipped ceiling of 64 already caps (the rule
# asks 71 and 102 there); the rest are the rule's value outright.
SHIPPED = {'n13a0': 8, 'n17a1': 8, 'n19a0': 13, 'n19a1': 25,
           'n23a0': 64, 'n23a1': 64, 'n29a1': 8, 'n31a0': 50}
UNCAPPED = dict(SHIPPED, n23a0=71, n23a1=102)
CEILINGS = [8, 12, 16, 24, 32, 48, 64]


def main():
    raw = json.loads(Path(sys.argv[1] if len(sys.argv) > 1 else DEFAULT).read_text())
    by = collections.defaultdict(dict)          # cell -> block -> {fixture: ir}
    for row in raw:
        by[row['cell']].setdefault(row['block'], {})[row['fixture']] = row['ir']
    cells = sorted(by, key=lambda s: (int(s[1:s.index('a')]), s))
    dev = [c for c in cells if c not in HOLDOUT]
    gm = statistics.geometric_mean

    def ratio(cell, block):
        """Cost of `block` against the shipped block at `cell`, paired by fixture."""
        ship = SHIPPED[cell]
        if block == ship:
            return 1.0
        a, b = by[cell].get(block), by[cell].get(ship)
        if not a or not b:
            return None
        shared = sorted(set(a) & set(b))
        return gm([a[f] / b[f] for f in shared]) if shared else None

    print('what each ceiling would cost, per cell (1.0000 = the ceiling does not '
          'bind there, so nothing changes)\n')
    head = '  '.join(f'{c:>8}' for c in CEILINGS)
    print(f"{'cell':9} {'rule':>5} {'ship':>5}  {head}")
    table = {}
    for cell in cells:
        row = []
        for ceil in CEILINGS:
            blk = min(UNCAPPED[cell], ceil)
            r = ratio(cell, blk) if blk != SHIPPED[cell] else 1.0
            table[(cell, ceil)] = r
            row.append(f'{r:8.4f}' if r is not None else f'{"n/a":>8}')
        tag = 'H' if cell in HOLDOUT else ' '
        print(f'{cell:8}{tag} {UNCAPPED[cell]:5} {SHIPPED[cell]:5}  ' + '  '.join(row))

    print('\ngeometric mean over cells:')
    print(f"{'ceiling':>8} {'development (decides)':>22} {'all eight':>12} {'worst dev cell':>16}")
    best = None
    for ceil in CEILINGS:
        d = [table[(c, ceil)] for c in dev]
        a = [table[(c, ceil)] for c in cells]
        if any(x is None for x in a):
            continue
        gd, ga = gm(d), gm(a)
        worst = max(dev, key=lambda c: table[(c, ceil)])
        print(f'{ceil:8} {gd:22.4f} {ga:12.4f} {worst:>10} {table[(worst, ceil)]:5.4f}')
        if best is None or gd < best[1]:
            best = (ceil, gd)
    print(f'\nchosen on the development cells: ceiling {best[0]} (gm {best[1]:.4f})')
    print('at the holdout cells that ceiling gives: ' +
          ', '.join(f'{c} {table[(c, best[0])]:.4f}' for c in sorted(HOLDOUT)))
    print('\nno cell regresses under it: ' +
          str(all(table[(c, best[0])] <= 1.0 + 1e-9 for c in cells)))


if __name__ == '__main__':
    main()
