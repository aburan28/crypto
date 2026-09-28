#!/usr/bin/env python3
"""Recreate the arm trees of round 0024 as amended (ROUND24-pair-inverse.md §7).

    python3 campaign_20260916/round24_amended_candidates.py

Requires the restored round-0020 evidence (`evidence/restore.py --archive
round-0020`): candidate sources are archived, not committed. Trees land in
`campaign_20260916/round24-sources/` (ignored by git).

Every tree carries `rho-normal-basis.patch`, so the round's rho -- the
incumbent build in rho mode -- names Frobenius classes with the IC arm's own
normal-basis rotation instead of walking the orbit by squaring
(research/ic_triple_counted_20260923/RESULTS.md §3).

  matched = round 0023's scaled arm (round-0020 `both`
            + round21-wide-pair-table.patch + round23-scaled-base.patch)
            + rho-normal-basis.patch; the round's incumbent and its rho
  counted = scaled + triple-table.patch + counted-sizing.patch
            + rho-normal-basis.patch; the triple-sum collector with counted
            sizing, run as `solver: triple_table`, `summands: 4`
"""
import json, shutil, subprocess
from pathlib import Path

WORK = Path(__file__).resolve().parent
ROOT = WORK.parent
REPO = ROOT.parent.parent
BASE = ROOT / 'runs/round-0020/source_candidates/both/source'
SCALED = [WORK / 'round21-wide-pair-table.patch', WORK / 'round23-scaled-base.patch']
TRIPLE = REPO / 'research/ic_triple_table_20260923/triple-table.patch'
COUNTED = REPO / 'research/ic_triple_counted_20260923/counted-sizing.patch'
RHO = REPO / 'research/ic_triple_counted_20260923/rho-normal-basis.patch'
PATCHES = {'matched': SCALED + [RHO], 'counted': SCALED + [TRIPLE, COUNTED, RHO]}


def main():
    assert BASE.is_dir(), 'restore round-0020 first'
    out = WORK / 'round24-sources'
    for name, patches in PATCHES.items():
        tree = out / name
        if tree.exists():
            shutil.rmtree(tree)
        tree.mkdir(parents=True)
        for item in ('Cargo.toml', 'Cargo.lock', 'build.rs', '.cargo/config.toml'):
            if (BASE / item).exists():
                (tree / item).parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(BASE / item, tree / item)
        for folder in ('src', 'examples', 'gpu'):
            if (BASE / folder).is_dir():
                shutil.copytree(BASE / folder, tree / folder)
        for patch in patches:
            subprocess.run(['patch', '-p1', '-s', '-i', str(patch)], cwd=tree, check=True)
        print(json.dumps({'candidate': name, 'source_root': str(tree),
                          'patches': [str(p.relative_to(REPO)) for p in patches]}))


if __name__ == '__main__':
    main()
