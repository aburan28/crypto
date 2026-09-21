#!/usr/bin/env python3
"""Recreate the round-0017 challenger tree from the frozen round-0016 scan_io
source and the committed patch, then print the registry to freeze.

    python3 campaign_20260916/round17_candidates.py

Requires the restored round-0016 evidence (`evidence/restore.py --archive
round-0016`), because candidate sources are archived, not committed. The tree
lands in `campaign_20260916/round17-sources/orbits` (ignored by git);
`round-0017-single-candidates.json` names it by absolute path.

The baseline and the control are the trees round 0016 sealed, unchanged.
`orbits` is `scan_io` plus `round17-orbits.patch`: the certificate names the
factor base by its orbit representatives and the checker
(`round17-oracle-orbits.patch`, already applied to `oracle.py`) expands them.
"""
import json, shutil, subprocess
from pathlib import Path

WORK = Path(__file__).resolve().parent
ROOT = WORK.parent
BASE = ROOT / 'runs/round-0016/source_candidates/scan_io/source'
PATCHES = {'orbits': ['round17-orbits.patch'],
           'orbits_rows': ['round17-orbits.patch', 'round17-rows.patch']}


def main():
    assert BASE.is_dir(), 'restore round-0016 first'
    out = WORK / 'round17-sources'
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
    registry = json.loads((WORK / 'round-0017-single-candidates.json').read_text())
    for arm in registry:
        if 'source_root' in arm:
            assert Path(arm['source_root']).is_dir(), arm['source_root']
    print(json.dumps({'registry': str(WORK / 'round-0017-single-candidates.json'), 'arms': len(registry)}))


if __name__ == '__main__':
    main()
