"""Replay the retained driver controls without executing new measured jobs."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

ARCHIVE_SHA256 = '972f1b6e16cc9a09960fccc80ba5f64331c3d12a5a312581521c53cea53d334c'
RUNS = {
    'initial': '36118318243',
    'intermediate': '36119164868',
    'clockfix': '36120078732',
    'final': '36121443713',
}


def checked_json(command):
    completed = subprocess.run(command, capture_output=True, text=True)
    if completed.returncode:
        raise RuntimeError(f'Audit failed: {command!r}\n{completed.stdout}\n{completed.stderr}')
    report = json.loads(completed.stdout)
    if report.get('status') != 'VERIFIED':
        raise ValueError(f'Audit did not verify: {report!r}')
    return report


def replay(bundle, also_python311):
    support = bundle / 'initial-replay-support'
    provenance = json.loads((support / 'provenance.json').read_text())
    dependency = (support / 'driver_admission.py').read_bytes()
    if hashlib.sha256(dependency).hexdigest() != provenance['dependency_sha256']:
        raise ValueError('Initial replay dependency differs from its retained provenance')
    original = bundle / 'initial' / provenance['source_artifact']
    if dependency != original.read_bytes():
        raise ValueError('Initial replay dependency differs from the same CI run sealed copy')

    results = []
    for phase, run in RUNS.items():
        for reference in ('both', 'scaled', 'pairinv'):
            root = bundle / phase / f'ic-producer-{reference}-{run}-1'
            script = root / 'ic-producer-evidence/evaluator/producer/audit.py'
            if phase == 'initial':
                # Only this original bundle lacks the transitive import. Supply
                # its same-CI sealed copy without changing any archived file.
                code = (f'import sys,runpy;sys.path.insert(0,{str(support)!r});'
                        f'sys.argv=[{str(script)!r},{str(root)!r}];'
                        f'runpy.run_path({str(script)!r},run_name="__main__")')
                command = [sys.executable, '-I', '-c', code]
            else:
                command = [sys.executable, '-I', str(script), str(root)]
            report = checked_json(command)
            expected = (9, 3) if reference == 'both' else (15, 5)
            if (report['reference'] != reference or
                    (report['verified_pairs'], report['rho_verified_pairs']) != expected):
                raise ValueError(f'Unexpected producer census: {report!r}')
            record = {key: report[key] for key in (
                'status', 'reference', 'source_manifest_sha256', 'worker_sha256',
                'verified_pairs', 'rho_verified_pairs')}
            record.update(phase=phase, kind='producer', supplemented_dependency=phase == 'initial')
            results.append(record)

        driver = bundle / phase / f'ic-driver-{run}-1'
        for name, tool in (('native', 'autolab.py'), ('native-failure', 'autolab.py'),
                           ('tournament', 'tournament.py')):
            if phase == 'initial' and name == 'native-failure':
                continue
            root = driver / name
            versions = [(sys.executable, '3.12')]
            if phase == 'final' and also_python311:
                versions.append(('python3.11', '3.11'))
            for executable, version in versions:
                report = checked_json([executable, str(root / 'evaluator' / tool),
                                       'verify', '--round', str(root)])
                results.append(dict(phase=phase, kind=name, python=version, audit=report))

    for name in ('native', 'failure'):
        root = bundle / 'final-native' / name
        report = checked_json([sys.executable, str(root / 'evaluator/autolab.py'),
                               'verify', '--round', str(root)])
        results.append(dict(phase='final-native', kind=name, python='3.12', audit=report))
    if len(results) != (28 if also_python311 else 25):
        raise ValueError('Incomplete audit list')
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bundle', type=Path, required=True, help='Freshly restored bundle root')
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--also-python311', action='store_true')
    args = parser.parse_args()
    if sys.version_info[:2] != (3, 12):
        parser.error('Use Python 3.12 for the earlier sealed floating-point summaries')
    archive = (Path(__file__).resolve().parents[2] / 'evidence' /
               'ic-driver-linux-controls-20260925.tar.zst')
    with archive.open('rb') as stream:
        actual = hashlib.file_digest(stream, 'sha256').hexdigest()
    if actual != ARCHIVE_SHA256:
        raise ValueError('Retained archive hash differs from the reviewed bundle')
    results = replay(args.bundle.resolve(), args.also_python311)
    args.out.write_text(json.dumps(dict(
        archive_sha256=actual, files=83113, expanded_file_bytes=372870950,
        audits=results), indent=2) + '\n')
    print(json.dumps(dict(status='VERIFIED', fresh_restore_audits=len(results),
                          final_python_versions=['3.11', '3.12'] if args.also_python311 else ['3.12'])))


if __name__ == '__main__':
    main()
