#!/usr/bin/env python3
"""Recreate the round-0009 baseline source tree from the frozen round-0008
incumbent source and the committed patches, then print the registry to freeze.

    python3 campaign_20260916/round9_candidates.py

Requires the restored evidence (`evidence/restore.py --archive round-0008`), because the baseline sources are archived, not
committed. Trees land in `campaign_20260916/round9-sources/<id>` (ignored by
git); `round-0009-single-candidates.json` names them by absolute path.
"""
import json
import shutil
import subprocess
from pathlib import Path

WORK = Path(__file__).resolve().parent
ROOT = WORK.parent
BASE = ROOT / 'runs/round-0008/source'
PATCHES = {'fastcurve': 'round9-fastcurve.patch'}


def main():
    assert BASE.is_dir(), 'restore round-0006 first'
    out = WORK / 'round9-sources'
    for name, patch in PATCHES.items():
        tree = out / name
        if tree.exists():
            shutil.rmtree(tree)
        tree.mkdir(parents=True)
        for item in ('Cargo.toml', 'Cargo.lock', 'build.rs'):
            if (BASE / item).exists():
                shutil.copy2(BASE / item, tree / item)
        for folder in ('src', 'examples', 'gpu'):
            shutil.copytree(BASE / folder, tree / folder)
        subprocess.run(['patch', '-p1', '-s', '-i', str(WORK / patch)], cwd=tree, check=True)
        print(json.dumps({'candidate': name, 'source_root': str(tree)}))
    registry = json.loads((WORK / 'round-0009-single-candidates.json').read_text())
    for arm in registry:
        if 'source_root' in arm:
            assert Path(arm['source_root']).is_dir(), arm['source_root']
    print(json.dumps({'registry': str(WORK / 'round-0009-single-candidates.json'), 'arms': len(registry)}))


if __name__ == '__main__':
    main()
