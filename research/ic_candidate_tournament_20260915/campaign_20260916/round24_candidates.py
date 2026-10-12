#!/usr/bin/env python3
"""Recreate round 0024's arm trees from the frozen round-0020 winner and
committed patches.

    python3 campaign_20260916/round24_candidates.py

Requires the restored round-0020 evidence (`evidence/restore.py --archive
round-0020`): candidate sources are archived, not committed. Trees land in
`campaign_20260916/round24-sources/` (ignored by git).

  incumbent = round 0023's incumbent (round-0020 `both`
              + round21-wide-pair-table.patch), retained by round 0023
  scaled  = round 0023's scaled arm (round-0020 `both`
            + round21-wide-pair-table.patch + round23-scaled-base.patch)
  pairinv = scaled + round24-pairinv.patch
            (+P and -P share an abscissa, so the decomposition scan inverts
            their common denominator once instead of twice; the scan order,
            every rest, witness and relation are unchanged)
"""
import json, shutil, subprocess
from pathlib import Path

WORK = Path(__file__).resolve().parent
ROOT = WORK.parent
BASE = ROOT / 'runs/round-0020/source_candidates/both/source'
SCALED = ['round21-wide-pair-table.patch', 'round23-scaled-base.patch']
PATCHES = {'incumbent': SCALED[:1], 'scaled': SCALED,
           'pairinv': SCALED + ['round24-pairinv.patch']}


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
            subprocess.run(['patch', '-p1', '-s', '-i', str(WORK / patch)], cwd=tree, check=True)
        print(json.dumps({'candidate': name, 'source_root': str(tree), 'patches': patches}))
    registry = json.loads((WORK / 'round-0024-single-candidates.json').read_text())
    for arm in registry:
        if 'source_root' in arm:
            assert Path(arm['source_root']).is_dir(), arm['source_root']
    print(json.dumps({'registry': str(WORK / 'round-0024-single-candidates.json'), 'arms': len(registry)}))


if __name__ == '__main__':
    main()
