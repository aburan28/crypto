#!/usr/bin/env python3
"""Recreate the round-0015 source trees from the frozen round-0014 source and
the committed patches, then print the registry to freeze.

    python3 campaign_20260916/round15_candidates.py

Requires the restored evidence (`evidence/restore.py --archive round-0014`),
because the baseline sources are archived, not committed. Trees land in
`campaign_20260916/round15-sources/<id>` (ignored by git);
`round-0015-single-candidates.json` names them by absolute path.

The round changes nothing shared: the baseline is the round-0014 winner source
exactly as it was sealed, so any difference measured here is the challengers'
own. `scan` is the IC-only change to the tiny pipeline (`round15-scan.patch`):
the orbit name of an abscissa found from the longest circular run of zeros,
and the block scan of a decomposition keeping only each rest's abscissa.
`scan_io` adds the report serialisation change (`round15-fastio.patch`).
Every tree carries `.cargo/config.toml`, which the evaluator snapshots and seals.
"""
import json
import shutil
import subprocess
from pathlib import Path

WORK = Path(__file__).resolve().parent
ROOT = WORK.parent
BASE = ROOT / 'runs/round-0014/source'
PATCHES = {
    'scan': ['round15-scan.patch'],
    'scan_io': ['round15-scan.patch', 'round15-fastio.patch'],
}


def main():
    assert BASE.is_dir(), 'restore round-0014 first'
    out = WORK / 'round15-sources'
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
            shutil.copytree(BASE / folder, tree / folder)
        for patch in patches:
            subprocess.run(['patch', '-p1', '-s', '-i', str(WORK / patch)], cwd=tree, check=True)
        print(json.dumps({'candidate': name, 'source_root': str(tree), 'patches': patches}))
    registry = json.loads((WORK / 'round-0015-single-candidates.json').read_text())
    for arm in registry:
        if 'source_root' in arm:
            assert Path(arm['source_root']).is_dir(), arm['source_root']
    print(json.dumps({'registry': str(WORK / 'round-0015-single-candidates.json'), 'arms': len(registry)}))


if __name__ == '__main__':
    main()
