"""Prepare the next fixed round from selection; run only after parent promotion."""
import copy
import difflib
import fcntl
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

from next_mechanisms import fast_lift, folded

WORK = Path(__file__).resolve().parent
ROOT = WORK.parent
PARENT = ROOT / 'runs/round-0003b'
ROUND = ROOT / 'runs/round-0004'
REL = Path('src/cryptanalysis/koblitz_index_calculus.rs')


def read(path):
    return json.loads(path.read_text())


def main():
    deadline = time.monotonic() + 3 * 3600
    while not (PARENT / 'summaries/selection.json').exists():
        if time.monotonic() > deadline:
            raise TimeoutError('parent selection did not complete')
        time.sleep(5)
    selected = read(PARENT / 'summaries/selection.json')['provisional_challenger']
    parent_arm = next(a for a in read(PARENT / 'candidates.json') if a['id'] == selected)
    source = PARENT / parent_arm['source_directory']
    config = parent_arm['config']
    original = (source / REL).read_text()
    registry = [{'id': 'incumbent', 'config': config,
                 'parent': selected, 'parent_round': str(PARENT),
                 'hypothesis': 'Previous selected implementation, admitted only after parent promotion.'}]
    for name, transform in [('fast_lift', fast_lift), ('folded', folded),
                            ('folded_lift', lambda s: folded(fast_lift(s)))]:
        destination = WORK / 'round4-sources' / name
        shutil.copytree(source, destination)
        changed = transform(original)
        (destination / REL).write_text(changed)
        (WORK / ('round4-' + name + '.patch')).write_text(''.join(difflib.unified_diff(
            original.splitlines(True), changed.splitlines(True),
            fromfile='a/' + str(REL), tofile='b/' + str(REL))))
        registry.append({'id': name, 'config': copy.deepcopy(config),
                         'source_root': str(destination), 'parent': 'incumbent',
                         'hypothesis': {'fast_lift': 'Preserve exact factor-base points while avoiding allocating half-trace arithmetic.',
                                        'folded': 'Build the same pair-sum coverage from signed Frobenius orbit representatives.',
                                        'folded_lift': 'Measure both independently validated mechanisms together.'}[name],
                         'falsification': 'Any bad or unmatched certificate; or failure of both-metric confirmation/replay gate.'})
    for name, delta in [('folded_lift_dense', {'linear_algebra': 'dense'}),
                        ('folded_lift_batch4', {'batch_trials': 4})]:
        candidate = copy.deepcopy(registry[-1] if len(registry) == 4 else registry[3])
        candidate.update(id=name, parent='folded_lift',
                         hypothesis='Test one residual scheduling/linear-algebra cost on the combined source.')
        candidate['config'].update(delta)
        registry.append(candidate)
    registry_path = WORK / 'round-0004-candidates.json'
    registry_path.write_text(json.dumps(registry, indent=2) + '\n')
    (WORK / 'ROUND4.md').write_text(f'''# Round 0004 pre-registration

Parent: round-0003b selection {selected}; source `{source}`. No parent holdout
data was used to select these mechanisms. Execution requires that same candidate
to pass the parent's full confirmation and replay.

Fresh seed 2026091604. Five challengers: fast factor-base lifting, folded pair
table, both, both with dense scalar algebra, both with batch4. Base support,
decomposition size, target count (one), all verification and full costs remain
fixed. Mathematical coverage and weak K-instruction floor stay unchanged.
Class: engineering. Same >=20% instruction and native-time promotion gates,
paired 95% limits, per-cell limits, confirmation and replay; same explicit
candidate/rho parity criterion. Budget 1800 paired jobs.

Builds are restricted to CPUs 0–1 while the parent measures on CPU 7. Round4
measurements start after the parent finishes, with every arm pinned to CPU 1.
The parent's full replay supplies a timing check after these builds complete.
Each round's runtime claims use its own fresh matched reference/candidate runs.

Preflight evidence: `preflight-next-tests.log` compares fast lifting against
general arithmetic (including point order) and folded decomposition against an
explicit full table on all five curve cells. These are correctness tests, not
performance claims. Frozen E2E oracle checks still cover every measured input.
''')
    os.sched_setaffinity(0, {0, 1})
    with (WORK / 'prepare-round-0004.log').open('w') as log:
        subprocess.run([sys.executable, str(ROOT / 'tournament.py'), 'prepare',
                        '--out', str(ROUND), '--source-root', str(source),
                        '--candidates', str(registry_path), '--seed', '2026091604',
                        '--cpu', '1', '--require-native-progress'],
                       stdout=log, stderr=subprocess.STDOUT, check=True)
    print('Round4 prepared; waiting for the parent final decision.', flush=True)
    with (PARENT / 'operation.lock').open('a+') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        decision = read(PARENT / 'decision.json')
        if decision['status'] != 'promoted' or decision['winner'] != selected:
            raise RuntimeError('Parent was not promoted: preserve this prepared proposal without running it.')
    with (ROUND / 'operation.jsonl').open('w') as log:
        subprocess.run([sys.executable, str(ROUND / 'evaluator/tournament.py'),
                        'run', '--round', str(ROUND)], stdout=log,
                       stderr=subprocess.STDOUT, check=True)
    print('Round4 completed.', flush=True)


if __name__ == '__main__':
    main()
