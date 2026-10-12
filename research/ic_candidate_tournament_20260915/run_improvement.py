#!/usr/bin/env python3
"""Run the versioned first toy improvement panel through the existing tournament."""
import argparse
import json
from pathlib import Path
import shutil
import subprocess
import sys

import campaign_rules as rules
from oracle import require
from tournament import digest, read, write

HERE = Path(__file__).resolve().parent
PANEL = HERE/'goal_20260924/improvement/round1.json'


def registry(panel, source):
    require(panel['round'] == 1 and panel['candidate_panel'] == 'round1-v1', 'unknown candidate panel')
    rows = panel['candidates']
    require(len(rows) == 16 and rows[0]['id'] == 'incumbent' and
            rows[0]['source'] == 'qualified-pairinv' and rows[0]['config'] == rules.CONFIG,
            'changed registered incumbent or arm count')
    require(all(row['source'] == 'round1-v1' for row in rows[1:]), 'unknown candidate source')
    return [dict(id=row['id'], config=row['config'],
                 **({'source_root':str(source)} if index else {})) for index,row in enumerate(rows)]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    out = args.out.resolve(); out.mkdir(parents=True, exist_ok=False)
    panel = read(PANEL)
    write(out/'registered-panel.json', panel, exclusive=True)
    commands = []

    def run(command):
        commands.append([str(value) for value in command])
        write(out/'commands.json', commands)
        print(json.dumps(dict(command=len(commands)-1, argv=commands[-1])), flush=True)
        with (out/f'command-{len(commands)-1}.log').open('x') as log:
            subprocess.run(commands[-1], stdout=log, stderr=subprocess.STDOUT, check=True)

    require(digest(HERE/'goal_20260924/improvement/target-history.json') == rules.HISTORY_SHA256,
            'changed target history')
    run([sys.executable, HERE/'target_history.py', '--history',
         HERE/'goal_20260924/improvement/target-history.json', '--repository', HERE.parents[1]])
    restored = out/'reference-evidence'
    run([sys.executable, HERE/'evidence/restore.py', '--archive',
         'ic-reference-qualification-20260925', '--out', restored])
    previous = restored/'ic-reference-qualification/full/ic-reference-qualification-36127866931-1/tournament'
    incumbent, cold_rho = previous/'source_candidates/pairinv/source', previous/'source'
    candidate = out/'candidate-source'
    run([sys.executable, HERE/'producer/prepare.py', '--reference', 'pairinv',
         '--candidate-panel', panel['candidate_panel'], '--out', candidate])
    require(read(candidate/'preparation.json')['source_manifest_sha256'] == panel['candidate_source_sha256'],
            'candidate differs from premeasurement registration')
    for source, expected in ((incumbent, rules.IC_SOURCE), (cold_rho, rules.COLD_RHO_SOURCE)):
        require(read(source.parent/'preparation.json')['source_manifest_sha256'] == expected,
                'restored qualified source differs')
    for source in (incumbent, cold_rho, candidate/'source'):
        run(['cargo', 'fetch', '--locked', '--manifest-path', source/'Cargo.toml'])
    write(out/'candidates.json', registry(panel, candidate/'source'), exclusive=True)
    write(out/'rho-config.json', dict(rules.CONFIG, rho_parallel_walks=4), exclusive=True)
    shutil.copy2(HERE/'goal_20260924/improvement/PROTOCOL.md', out/'PROTOCOL.md')
    campaign = out/'tournament'
    run([sys.executable, HERE/'tournament.py', 'prepare', '--attempt-number', '1',
         '--source-root', incumbent, '--rho-source-root', cold_rho, '--rho-config', out/'rho-config.json',
         '--qualified-report', previous/'qualification.json', '--target-history',
         HERE/'goal_20260924/improvement/target-history.json', '--out', campaign,
         '--candidates', out/'candidates.json', '--cells', '17a1,19a0,23a0,23a1,31a0',
         '--holdout-cells', '29a1', '--profile', 'pilot', '--seed', '2026092551',
         '--timeout', '180', '--max-processes', '3500', '--selection-width', '6',
         '--exploration-slots', '1', '--comparison-kind', 'factor-base-policy', '--require-native-progress'])
    run([sys.executable, campaign/'evaluator/tournament.py', 'run', '--round', campaign])
    run([sys.executable, campaign/'evaluator/tournament.py', 'verify', '--round', campaign])
    report = read(campaign/'decision.json')
    require(report['attempt_number'] == 1 and report['familywise_rule'] == rules.RULE,
            'decision omitted frozen improvement rules')
    write(out/'summary.json', report, exclusive=True)
    print(json.dumps(report), flush=True)


if __name__ == '__main__':
    main()
