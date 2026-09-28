#!/usr/bin/env python3
"""Run the registered final toy improvement panel through the frozen tournament.

Rounds one and two stay sealed. This wrapper binds the accepted cold/online
references, restores both completed rounds, excludes their exposed points and
supplemental fixtures, and freezes the third panel. It does not retune on prior
confirmation or replay outcomes.
"""
import argparse
import json
from pathlib import Path
import shutil
import subprocess
import sys

import campaign_rules_v2 as rules
from oracle import require
from tournament import digest, read, write

HERE = Path(__file__).resolve().parent
PANEL = HERE / 'goal_20260924/improvement-v2/round3.json'
PANEL_SHA256 = '30b04a370e77e2a2bb7c01b7e96fd1ad8b1785104240113ac71f85eac28beca8'
HISTORY = HERE / 'goal_20260924/improvement/target-history.json'
GOAL = HERE / 'goal_20260924'
EXPOSED = (
    GOAL / 'generic-reference-readiness/fixtures.json',
    GOAL / 'generic-reference-qualification/fixtures.json',
    GOAL / 'generic-adapter-control/fixtures.json',
)


def registry(panel, source):
    require(panel['round'] == 3 and panel['candidate_panel'] == 'round3-v1',
            'unknown candidate panel')
    rows = panel['candidates']
    require(len(rows) == 11 and rows[0]['id'] == 'incumbent' and
            rows[0]['source'] == 'qualified-pairinv' and rows[0]['config'] == rules.CONFIG,
            'changed registered incumbent or arm count')
    require(all(row['source'] == 'round3-v1' for row in rows[1:]), 'unknown candidate source')
    return [dict(id=row['id'], config=row['config'],
                 **({'source_root': str(source)} if index else {}))
            for index, row in enumerate(rows)]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    require(digest(PANEL) == PANEL_SHA256, 'changed premeasurement candidate panel')
    panel = read(PANEL)
    write(out / 'registered-panel.json', panel, exclusive=True)
    commands = []

    def run(command, *, capture_json=False):
        commands.append([str(value) for value in command])
        write(out / 'commands.json', commands)
        print(json.dumps(dict(command=len(commands) - 1, argv=commands[-1])), flush=True)
        log_path = out / f'command-{len(commands) - 1}.log'
        with log_path.open('x') as log:
            completed = subprocess.run(commands[-1], stdout=log,
                                       stderr=subprocess.STDOUT, check=False, text=True)
        write(out / f'command-{len(commands) - 1}.status.json',
              dict(exit_code=completed.returncode), exclusive=True)
        require(completed.returncode == 0,
                f'command {len(commands) - 1} failed; retained output in {log_path}')
        if not capture_json:
            return None
        for line in reversed(log_path.read_text().splitlines()):
            line = line.strip()
            if line.startswith('{'):
                return json.loads(line)
        raise AssertionError('command omitted JSON receipt')

    require(digest(HISTORY) == rules.HISTORY_SHA256, 'changed target history')
    for path in EXPOSED:
        require(path.is_file(), f'missing exposed fixture corpus: {path.name}')

    run([sys.executable, HERE / 'target_history.py', '--history', HISTORY,
         '--repository', HERE.parents[1]])

    references = out / 'reference-evidence'
    run([sys.executable, HERE / 'evidence/restore.py', '--archive',
         'ic-generic-reference-qualification-20260926', '--out', references])
    bundle = references / 'ic-generic-reference-qualification'
    declared = run([sys.executable, HERE / 'reference_registry_v2.py', '--bundle', bundle,
                    '--out', out / 'reference-registry.json'], capture_json=True)
    require(declared['status'] == 'DECLARED' and declared['workers_executed'] == 0
            and declared['targets_generated'] == 0, 'reference declaration was not dry')
    write(out / 'reference-declaration.json', declared, exclusive=True)
    require([arm['id'] for arm in read(out / 'reference-registry.json')] ==
            list(rules.REFERENCE_ROLES), 'reference registry omitted a required role')

    prior_root = out / 'prior-round1'
    run([sys.executable, HERE / 'evidence/restore.py', '--archive',
         'ic-improvement-round1-20260925', '--out', prior_root])
    prior = prior_root / 'ic-improvement-round1/round/tournament'
    require((prior / 'contract.json').is_file() and (prior / 'decision.json').is_file(),
            'restored round-one tournament is incomplete')
    prior_root2 = out / 'prior-round2'
    run([sys.executable, HERE / 'evidence/restore.py', '--archive',
         'ic-improvement-round2-20260928', '--out', prior_root2])
    prior2 = prior_root2 / 'ic-improvement-round2/round/tournament'
    require((prior2 / 'contract.json').is_file() and (prior2 / 'decision.json').is_file(),
            'restored round-two tournament is incomplete')

    incumbent = Path(declared['incumbent_source'])
    require(read(incumbent.parent / 'preparation.json')['source_manifest_sha256'] == rules.IC_SOURCE,
            'restored qualified incumbent differs')
    require((references / 'ic-generic-qualification-build/worker').is_file(),
            'missing sealed generic reference worker')

    candidate = out / 'candidate-source'
    run([sys.executable, HERE / 'producer/prepare.py', '--reference', 'pairinv',
         '--candidate-panel', panel['candidate_panel'], '--out', candidate])
    require(read(candidate / 'preparation.json')['source_manifest_sha256'] ==
            panel['candidate_source_sha256'], 'candidate differs from premeasurement registration')

    for source in (incumbent, candidate / 'source'):
        run(['cargo', 'fetch', '--locked', '--manifest-path', source / 'Cargo.toml'])

    write(out / 'candidates.json', registry(panel, candidate / 'source'), exclusive=True)
    shutil.copy2(HERE / 'goal_20260924/improvement-v2/PROTOCOL.md', out / 'PROTOCOL.md')
    campaign = out / 'tournament'
    command = [sys.executable, HERE / 'tournament.py', 'prepare',
               '--campaign-version', '2', '--attempt-number', '3',
               '--source-root', incumbent,
               '--reference-registry', out / 'reference-registry.json',
               '--qualified-report', declared['qualified_report'],
               '--qualified-observer', declared['qualified_observer'],
               '--target-history', HISTORY,
               '--prior-round', prior, '--prior-round', prior2,
               '--out', campaign,
               '--candidates', out / 'candidates.json',
               '--cells', '17a1,19a0,23a0,23a1,31a0',
               '--holdout-cells', '29a1',
               '--profile', 'pilot', '--seed', '2026092553',
               '--timeout', '180', '--max-processes', '3500',
               '--selection-width', '6', '--exploration-slots', '1',
               '--comparison-kind', 'factor-base-policy',
               '--require-native-progress']
    for path in EXPOSED:
        command.extend(['--exposed-fixtures', path])
    run(command)
    run([sys.executable, campaign / 'evaluator/tournament.py', 'run', '--round', campaign])
    run([sys.executable, campaign / 'evaluator/tournament.py', 'verify', '--round', campaign])
    report = read(campaign / 'decision.json')
    require(report['attempt_number'] == 3 and report['familywise_rule'] == rules.RULE,
            'decision omitted frozen improvement rules')
    write(out / 'summary.json', report, exclusive=True)
    print(json.dumps(report), flush=True)


if __name__ == '__main__':
    main()
