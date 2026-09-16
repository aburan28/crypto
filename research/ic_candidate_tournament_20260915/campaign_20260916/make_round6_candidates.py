"""Create isolated, reviewable round-0006 candidates from the round-0005 winner.

Each candidate is a copy of the frozen `combined_descent` source with exactly
one mechanism applied (plus one declared combination), a unified diff, and a
registry entry. No measurement happens here; `tournament.py prepare` freezes
and builds the trees it is given."""
import copy
import difflib
import json
from pathlib import Path
import shutil
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from round6_mechanisms import CANDIDATES, HYPOTHESES  # noqa: E402

WORK = Path(__file__).resolve().parent
ROOT = WORK.parent
PARENT = ROOT / 'runs/round-0005-batch16'
SOURCES = WORK / 'round6-sources'
FALSIFICATION = ('Any incorrect/missing target rejects; promotion requires >=20% lower '
                 'instructions and native wall, upper paired 95% limits <1 and every cell '
                 '<=1.10, on confirmation and replay; rho parity is reported separately.')


def main():
    decision = json.loads((PARENT / 'decision.json').read_text())
    assert decision['status'] == 'promoted' and decision['winner'] == 'combined_descent'
    arm = next(a for a in json.loads((PARENT / 'candidates.json').read_text()) if a['id'] == 'combined_descent')
    source = PARENT / arm['source_directory']
    config = copy.deepcopy(arm['config'])
    registry = [{'id': 'incumbent', 'config': config, 'parent': 'combined_descent',
                 'parent_round': str(PARENT),
                 'hypothesis': 'Round-0005 winner, freshly measured on the same complete cold 16-target jobs.'}]
    for name, transforms in CANDIDATES.items():
        destination = SOURCES / name
        if destination.exists():
            raise SystemExit(f'refusing to replace an existing candidate: {destination}')
        shutil.copytree(source, destination)
        diff = []
        for relative, steps in transforms.items():
            original = (source / relative).read_text()
            changed = original
            for step in steps:
                changed = step(changed)
            (destination / relative).write_text(changed)
            diff.extend(difflib.unified_diff(original.splitlines(True), changed.splitlines(True),
                                             fromfile='a/' + relative, tofile='b/' + relative))
        (WORK / f'round6-{name}.patch').write_text(''.join(diff))
        registry.append({'id': name, 'source_root': str(destination), 'config': copy.deepcopy(config),
                         'parent': 'incumbent', 'hypothesis': HYPOTHESES[name],
                         'falsification': FALSIFICATION})
    (WORK / 'round-0006-candidates.json').write_text(json.dumps(registry, indent=2) + '\n')
    print(json.dumps({'baseline_source_root': str(source), 'candidates': len(registry) - 1,
                      'registry': str(WORK / 'round-0006-candidates.json')}))


if __name__ == '__main__':
    main()
