#!/usr/bin/env python3
"""Recreate round 0020's single challenger tree from the frozen baseline and
the committed round-0018 patches, then check the registry.

    python3 campaign_20260916/round19_candidates.py

The tree lands in `campaign_20260916/round20-sources/both` (ignored by git);
`round-0020-single-candidates.json` names it by absolute path.

The baseline is `runs/round-0019/source` -- the tree round 0018b retained as
its winner, byte-identical to the tree round 0017 sealed as `orbits`. Requires
the restored evidence (`evidence/restore.py --archive round-0018b`), because
candidate sources are archived, not committed.

  both = baseline + round18-block.patch + round18-column.patch

Round 0018b measured that arm at 0.9652 [0.9448, 0.9828] against the incumbent
in instructions and 0.9613 [0.9288, 0.9895] natively. It passed the
no-regression gate on both final stages and failed only the strict rho gate, at
`n23a1` alone, where ROUND19-single-target.md sec. 2 shows twelve fixtures
cannot resolve the cell either way. Round 0019 carries no new arm: sec. 3's
base-size model says there is nothing left to take at that cell.

`round18-oracle-convention.patch` is already applied to `oracle.py`; it is the
checker amendment `both` needs, and it reads either column convention.
"""
import json
import shutil
import subprocess
from pathlib import Path

WORK = Path(__file__).resolve().parent
ROOT = WORK.parent
BASE = ROOT / 'runs/round-0019/source'
# The same tree under its other name. Both must be present and identical; the
# round-0018 baseline error was exactly a source root that looked right.
TWIN = ROOT / 'runs/round-0017/source_candidates/orbits/source'
PATCHES = {'both': ['round18-block.patch', 'round18-column.patch']}
COPY_FILES = ('Cargo.toml', 'Cargo.lock', 'build.rs', '.cargo/config.toml')
COPY_DIRS = ('src', 'examples', 'gpu')


def digest(path):
    return subprocess.run(['sha256sum', str(path)], capture_output=True, text=True,
                          check=True).stdout.split()[0]


def main():
    assert BASE.is_dir(), f'restore round-0018b first: {BASE}'
    # The baseline must be the round-0017 winner, not the round-0017 incumbent.
    # Round 0018 was run against the wrong one of those two and had to be
    # discarded; this is the check that would have caught it.
    worker = ROOT / 'runs/round-0019/worker'
    twin_worker = ROOT / 'runs/round-0017/source_candidates/orbits/worker'
    if twin_worker.exists():
        assert digest(worker) == digest(twin_worker), \
            f'baseline worker {digest(worker)} is not the round-0017 winner {digest(twin_worker)}'
    marker = (BASE / 'examples/ic_tournament_worker.rs').read_text()
    assert 'factor_base_orbits' in marker, \
        'baseline does not carry the round-0017 orbit certificate; wrong source root'
    if TWIN.is_dir():
        assert (TWIN / 'examples/ic_tournament_worker.rs').read_text() == marker, \
            'the two names for the baseline tree disagree'

    out = WORK / 'round20-sources'
    for name, patches in PATCHES.items():
        tree = out / name
        if tree.exists():
            shutil.rmtree(tree)
        tree.mkdir(parents=True)
        for item in COPY_FILES:
            if (BASE / item).exists():
                (tree / item).parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(BASE / item, tree / item)
        for folder in COPY_DIRS:
            if (BASE / folder).is_dir():
                shutil.copytree(BASE / folder, tree / folder)
        for patch in patches:
            subprocess.run(['patch', '-p1', '-s', '-i', str(WORK / patch)], cwd=tree, check=True)
        assert 'representative' in (tree / 'examples/ic_tournament_worker.rs').read_text(), \
            f'{name}: the column patch did not reach the worker'
        print(json.dumps({'candidate': name, 'source_root': str(tree), 'patches': patches}))

    registry = json.loads((WORK / 'round-0020-single-candidates.json').read_text())
    for arm in registry:
        if 'source_root' in arm:
            assert Path(arm['source_root']).is_dir(), arm['source_root']
    print(json.dumps({'registry': str(WORK / 'round-0020-single-candidates.json'),
                      'arms': [a['id'] for a in registry],
                      'baseline_worker_sha256': digest(worker)}))


if __name__ == '__main__':
    main()
