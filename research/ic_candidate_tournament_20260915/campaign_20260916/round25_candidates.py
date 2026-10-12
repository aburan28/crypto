#!/usr/bin/env python3
"""Recreate the arm trees of round 0025 (ROUND25-lean-rho-switch.md).

    python3 campaign_20260916/round25_candidates.py

Requires the restored round-0020 evidence (`evidence/restore.py --archive
round-0020`): candidate sources are archived, not committed. Trees land in
`campaign_20260916/round25-sources/` (ignored by git).

Every tree carries the matched rho (rho-normal-basis.patch) and then the
lean rho (round25-rho-lean.patch), so the round's rho -- the incumbent build
in rho mode -- takes the matched rho's exact walk at the lower per-step cost
of the IC arm's own single-word arithmetic.

  lean   = round 0024's `matched` (round-0020 `both` + round21-wide-pair-table
           + round23-scaled-base + rho-normal-basis) + round25-rho-lean.patch;
           the incumbent, pair_table at three summands, and the round's rho
  switch = round 0024's `counted` (the same + triple-table + counted-sizing)
           + round25-switch.patch + round25-rho-lean.patch; pair_or_triple at
           four summands
"""
import json, shutil, subprocess
from pathlib import Path

WORK = Path(__file__).resolve().parent
ROOT = WORK.parent
REPO = ROOT.parent.parent
BASE = ROOT / 'runs/round-0020/source_candidates/both/source'
SCALED = [WORK / 'round21-wide-pair-table.patch', WORK / 'round23-scaled-base.patch']
TRIPLE = [REPO / 'research/ic_triple_table_20260923/triple-table.patch',
          REPO / 'research/ic_triple_counted_20260923/counted-sizing.patch']
MATCHED = REPO / 'research/ic_triple_counted_20260923/rho-normal-basis.patch'
LEAN = WORK / 'round25-rho-lean.patch'
SWITCH = WORK / 'round25-switch.patch'
PATCHES = {'lean': SCALED + [MATCHED, LEAN],
           'switch': SCALED + TRIPLE + [MATCHED, SWITCH, LEAN]}


def main():
    assert BASE.is_dir(), 'restore round-0020 first'
    out = WORK / 'round25-sources'
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
