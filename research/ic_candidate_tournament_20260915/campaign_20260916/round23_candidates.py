#!/usr/bin/env python3
"""Recreate round 0023's two arm trees from the frozen round-0020 winner and
committed patches, then check the registry against them.

    python3 campaign_20260916/round23_candidates.py

Requires the restored round-0020 evidence (`evidence/restore.py --archive
round-0020`): candidate sources are archived, not committed. Trees land in
`campaign_20260916/round23-sources/` (ignored by git).

  incumbent = round-0020 `both` + round21-wide-pair-table.patch
              (pair table widened to u64, degree bounds 31 -> 61; round 0021
              measured it at 1.0007x the promoted worker with identical
              logarithms and factor bases at all eight panel cells)
  scaled    = incumbent + round23-scaled-base.patch
              (orbit count set from the subgroup order by the rule frozen from
              round 0022's measured optima; one batch of eight on every panel
              cell, 16 at n37a0, 24 at n43a1, n59a0 and n61a1)
"""
import json, shutil, subprocess
from pathlib import Path

WORK = Path(__file__).resolve().parent
ROOT = WORK.parent
BASE = ROOT / 'runs/round-0020/source_candidates/both/source'
PATCHES = {'incumbent': ['round21-wide-pair-table.patch'],
           'scaled': ['round21-wide-pair-table.patch', 'round23-scaled-base.patch']}


def main():
    assert BASE.is_dir(), 'restore round-0020 first'
    out = WORK / 'round23-sources'
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
            subprocess.run(['patch', '-p1', '-s', '-i', str(WORK / patch)], cwd=tree, check=True)
        print(json.dumps({'candidate': name, 'source_root': str(tree), 'patches': patches}))
    registry = json.loads((WORK / 'round-0023-single-candidates.json').read_text())
    for arm in registry:
        if 'source_root' in arm:
            assert Path(arm['source_root']).is_dir(), arm['source_root']
    print(json.dumps({'registry': str(WORK / 'round-0023-single-candidates.json'), 'arms': len(registry)}))


if __name__ == '__main__':
    main()
