#!/usr/bin/env python3
"""Recreate the round-0010 source trees from the frozen round-0009 incumbent
source and the committed patches, then print the registry to freeze.

    python3 campaign_20260916/round10_candidates.py

Requires the restored evidence (`evidence/restore.py --archive round-0009`),
because the baseline sources are archived, not committed. Trees land in
`campaign_20260916/round10-sources/<id>` (ignored by git);
`round-0010-single-candidates.json` names them by absolute path. `lean` is the
round's baseline (`--source-root`); `lean_stdprobe` is its ablation control.
Both carry `.cargo/config.toml`, which the evaluator snapshots and seals.
"""
import json
import shutil
import subprocess
from pathlib import Path

WORK = Path(__file__).resolve().parent
ROOT = WORK.parent
BASE = ROOT / 'runs/round-0009/source'
PATCHES = {'lean': 'round10-lean.patch', 'lean_stdprobe': 'round10-lean_stdprobe.patch'}


def main():
    assert BASE.is_dir(), 'restore round-0009 first'
    out = WORK / 'round10-sources'
    for name, patch in PATCHES.items():
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
        subprocess.run(['patch', '-p1', '-s', '-i', str(WORK / patch)], cwd=tree, check=True)
        print(json.dumps({'candidate': name, 'source_root': str(tree)}))
    registry = json.loads((WORK / 'round-0010-single-candidates.json').read_text())
    for arm in registry:
        if 'source_root' in arm:
            assert Path(arm['source_root']).is_dir(), arm['source_root']
    print(json.dumps({'registry': str(WORK / 'round-0010-single-candidates.json'), 'arms': len(registry)}))


if __name__ == '__main__':
    main()
