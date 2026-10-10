"""Replay the retained comparison and observer records without running workers."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent
EVIDENCE = HERE.parents[1] / 'evidence'
ARCHIVE = 'ic-generic-reference-qualification-20260926.tar.zst'


def read(path):
    return json.loads(path.read_text())


def require(condition, message):
    if not condition:
        raise ValueError(message)


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def restore(destination):
    entry = next(row for row in read(EVIDENCE / 'manifest.json')['archives'] if row['file'] == ARCHIVE)
    module('generic_qualification_restore', EVIDENCE / 'restore.py').restore(entry, destination.resolve())
    root = destination / 'ic-generic-reference-qualification'
    require((root / 'tournament/fixtures.json').read_bytes() == (HERE / 'fixtures.json').read_bytes(),
            'changed exposed-point export')
    return root


def checked_json(command):
    print('Checking ' + ' '.join(map(str, command)), file=sys.stderr, flush=True)
    result = subprocess.run(list(map(str, command)), capture_output=True, text=True, timeout=1200)
    require(result.returncode == 0, f'frozen audit failed: {command!r}\n' + result.stdout + result.stderr)
    return json.loads(result.stdout)


def replay_controls(root):
    controls = root / 'controls'
    reports = []
    for reference in ('both', 'scaled', 'pairinv'):
        directory = controls / f'ic-producer-{reference}-36290704597-1'
        report = checked_json([sys.executable, '-I',
            directory / 'ic-producer-evidence/evaluator/producer/audit.py', directory])
        expected = (9, 3) if reference == 'both' else (15, 5)
        require(report['status'] == 'VERIFIED' and report['reference'] == reference
                and (report['verified_pairs'], report['rho_verified_pairs']) == expected,
                'changed prerequisite producer census')
        reports.append(dict(kind='producer', reference=reference, audit=report))
    driver = controls / 'ic-driver-36290704597-1'
    for name, tool, expected in (
            ('native', 'autolab.py', 12), ('native-failure', 'autolab.py', 6),
            ('tournament', 'tournament.py', 15), ('qualification-control', 'tournament.py', 66)):
        directory = driver / 'ic-driver-evidence' / name
        report = checked_json([sys.executable, directory / 'evaluator' / tool, 'verify', '--round', directory])
        require(report['status'] == 'VERIFIED' and report['trial_receipts'] == expected,
                'changed prerequisite driver census')
        reports.append(dict(kind=name, audit=report))
    directory = driver / 'ic-generic-driver-control/tournament'
    report = checked_json([sys.executable, directory / 'evaluator/tournament.py', 'verify', '--round', directory])
    require(report['status'] == 'VERIFIED' and report['trial_receipts'] == 102,
            'changed mixed-adapter control census')
    reports.append(dict(kind='mixed-adapter', audit=report))
    directory = driver / 'ic-generic-driver-observer'
    report = checked_json([sys.executable, directory / 'evaluator/generic_observer.py', 'verify', '--out', directory])
    require(report == dict(status='VERIFIED', observer_pairs=45, workers_executed=0),
            'changed observer control census')
    summary = read(directory / 'summary.json')
    require((summary['complete_pairs'], summary['incomplete_pairs'], summary['failed_pairs']) == (36, 9, 0),
            'expected preparation failures changed')
    reports.append(dict(kind='observer-control', audit=report))
    require(len(reports) == 9, 'partial prerequisite replay')
    return reports


def replay(root):
    root = root.resolve()
    campaign, observer = root / 'tournament', root / 'observer'
    compared = checked_json([sys.executable, campaign / 'evaluator/tournament.py', 'verify', '--round', campaign])
    require(compared['status'] == 'VERIFIED' and compared['trial_receipts'] == 1350,
            'comparison replay is incomplete')
    observed = checked_json([sys.executable, observer / 'evaluator/generic_observer.py', 'verify', '--out', observer])
    require(observed == dict(status='VERIFIED', observer_pairs=360, workers_executed=0),
            'observer replay is incomplete')
    program = ("import json,runpy,sys; from pathlib import Path; root=Path(sys.argv[1]); "
               "sys.path.insert(0,str(root/'tournament/evaluator')); "
               "runner=runpy.run_path(str(root/'registered-runner.py')); "
               "print(json.dumps(runner['panel_result'](root/'tournament')))")
    result = checked_json([sys.executable, '-c', program, root])
    require(result == read(root / 'comparison-summary.json'), 'changed frozen panel summary')
    summary = read(observer / 'summary.json')
    result.update(observer_summary_sha256=hashlib.sha256((observer / 'summary.json').read_bytes()).hexdigest(),
                  observer_pairs=360, observer_complete_pairs=summary['complete_pairs'],
                  observer_incomplete_pairs=summary['incomplete_pairs'], accepted_reference_binding_changed=False)
    require(result == read(root / 'summary.json'), 'changed final registered summary')
    parent_ids, observer_ids = set(), set()
    for path in campaign.glob('runs/**/receipt.json'):
        run = read(path)['measurement']
        parent_ids.add(run['run_id'])
        if run.get('adapter') == 'generic-v1':
            parent_ids.add(run['profile_execution']['run_id'])
    for path in observer.glob('runs/*/receipt.json'):
        observer_ids.update(read(path)['run_ids'].values())
    require(len(parent_ids) == 1830 and len(observer_ids) == 720
            and not parent_ids.intersection(observer_ids),
            'observer executions collide with the parent campaign')
    controls = replay_controls(root)
    return dict(status='VERIFIED', comparison=compared, observer=observed, control_audits=controls,
                frozen_audits=len(controls) + 2,
                distinct_record_ids=len(parent_ids | observer_ids), workers_executed=0,
                promotion_eligible=False, accepted_reference_binding_changed=False)


def exported(root):
    return module('generic_qualification_export', HERE / 'export_report.py').export(root)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bundle', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    require(sys.version_info[:2] == (3, 12), 'use Python 3.12 for frozen summaries')
    result = replay(args.bundle)
    require(exported(args.bundle) == read(HERE / 'RESULTS.json'), 'changed result export')
    with args.out.open('x') as stream:
        json.dump(result, stream, indent=2)
        stream.write('\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
