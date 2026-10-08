#!/usr/bin/env python3
"""Recreate the round-0013 source trees from the frozen round-0012 source and
the committed patches, then print the registry to freeze.

    python3 campaign_20260916/round13_candidates.py

Requires the restored evidence (`evidence/restore.py --archive round-0012`),
because the baseline sources are archived, not committed. Trees land in
`campaign_20260916/round13-sources/<id>` (ignored by git);
`round-0013-single-candidates.json` names them by absolute path.

`ld` is the round's baseline (`--source-root`): the general and the
single-word scalar multiplications in López–Dahab coordinates
(`round13-ld.patch`), shared by every arm including rho. `ld_canon` adds
the IC-only orbit naming from the longest zero run (`round13-canon.patch`);
`ld_ic` adds to that the IC-only serialisation change of round 0011
(`round13-fastio.patch`, rebased). Every tree carries `.cargo/config.toml`,
which the evaluator snapshots and seals.
"""
import json
import shutil
import subprocess
from pathlib import Path

WORK = Path(__file__).resolve().parent
ROOT = WORK.parent
BASE = ROOT / 'runs/round-0012/source'
PATCHES = {
    'ld': ['round13-ld.patch'],
    'ld_canon': ['round13-ld.patch', 'round13-canon.patch'],
    'ld_ic': ['round13-ld.patch', 'round13-canon.patch', 'round13-fastio.patch'],
}


def main():
    assert BASE.is_dir(), 'restore round-0012 first'
    out = WORK / 'round13-sources'
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
    registry = json.loads((WORK / 'round-0013-single-candidates.json').read_text())
    for arm in registry:
        if 'source_root' in arm:
            assert Path(arm['source_root']).is_dir(), arm['source_root']
    print(json.dumps({'registry': str(WORK / 'round-0013-single-candidates.json'), 'arms': len(registry)}))


if __name__ == '__main__':
    main()
