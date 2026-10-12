#!/usr/bin/env python3
"""Recreate the round-0018 challenger trees from the frozen round-0017 winner
source and the committed patches, then print the registry to freeze.

    python3 campaign_20260916/round18_candidates.py

Requires the restored round-0017 evidence (`evidence/restore.py --archive
round-0017`), because candidate sources are archived, not committed. The trees
land in `campaign_20260916/round18-sources/<arm>` (ignored by git);
`round-0018-single-candidates.json` names them by absolute path.

The baseline is the tree round 0017 sealed as its winner, `orbits`, unchanged.

  block   = orbits + round18-block.patch
  column  = orbits + round18-column.patch
  both    = orbits + round18-block.patch + round18-column.patch

`round18-oracle-convention.patch` is already applied to `oracle.py`; it is the
checker amendment `column` and `both` need and is recorded beside these.
"""
import json
import shutil
import subprocess
from pathlib import Path

WORK = Path(__file__).resolve().parent
ROOT = WORK.parent
BASE = ROOT / 'runs/round-0017/source_candidates/orbits/source'
PATCHES = {
    'block': ['round18-block.patch'],
    'column': ['round18-column.patch'],
    'both': ['round18-block.patch', 'round18-column.patch'],
}


def main():
    assert BASE.is_dir(), 'restore round-0017 first (evidence/restore.py --archive round-0017)'
    out = WORK / 'round18-sources'
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
    registry = json.loads((WORK / 'round-0018-single-candidates.json').read_text())
    for arm in registry:
        if 'source_root' in arm:
            assert Path(arm['source_root']).is_dir(), arm['source_root']
    print(json.dumps({'registry': str(WORK / 'round-0018-single-candidates.json'),
                      'arms': len(registry)}))


if __name__ == '__main__':
    main()
