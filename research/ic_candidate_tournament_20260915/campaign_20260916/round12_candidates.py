#!/usr/bin/env python3
"""Recreate the round-0012 source trees from the frozen round-0011 source and
the committed patches, then print the registry to freeze.

    python3 campaign_20260916/round12_candidates.py

Requires the restored evidence (`evidence/restore.py --archive round-0011`),
because the baseline sources are archived, not committed. Trees land in
`campaign_20260916/round12-sources/<id>` (ignored by git);
`round-0012-single-candidates.json` names them by absolute path.

`arena` is the round's baseline (`--source-root`): the musl static executable
without relocations (`round12-musl.patch`, the cargo config only) and the
worker's bump-pointer arena allocator (`round12-arena.patch`). `arena_fastio`
adds the IC-only serialisation change of round 0011 (`round12-fastio.patch`, the
same change rebased on the arena worker).
`arena_glibc` and `musl_sysalloc` are the two ablation controls: each carries
one of the two baseline changes without the other. The library source is
identical in all four trees. Every tree carries `.cargo/config.toml`, which
the evaluator snapshots and seals.
"""
import json
import shutil
import subprocess
from pathlib import Path

WORK = Path(__file__).resolve().parent
ROOT = WORK.parent
BASE = ROOT / 'runs/round-0011/source'
PATCHES = {
    'arena': ['round12-musl.patch', 'round12-arena.patch'],
    'arena_fastio': ['round12-musl.patch', 'round12-arena.patch', 'round12-fastio.patch'],
    'arena_glibc': ['round12-arena.patch'],
    'musl_sysalloc': ['round12-musl.patch'],
}


def main():
    assert BASE.is_dir(), 'restore round-0011 first'
    out = WORK / 'round12-sources'
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
    registry = json.loads((WORK / 'round-0012-single-candidates.json').read_text())
    for arm in registry:
        if 'source_root' in arm:
            assert Path(arm['source_root']).is_dir(), arm['source_root']
    print(json.dumps({'registry': str(WORK / 'round-0012-single-candidates.json'), 'arms': len(registry)}))


if __name__ == '__main__':
    main()
