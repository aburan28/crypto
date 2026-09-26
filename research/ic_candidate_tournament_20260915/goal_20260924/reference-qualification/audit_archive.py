"""Strictly replay qualification evidence; never execute measured workers."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys


ARCHIVE = 'ic-reference-qualification-20260925.tar.zst'
RUNS = {'initial': '36126765057', 'final': '36127830907', 'full': '36127866931'}


def checked_json(command):
    result = subprocess.run(command, capture_output=True, text=True)
    if result.returncode:
        raise RuntimeError(f'Audit failed: {command!r}\n{result.stdout}\n{result.stderr}')
    report = json.loads(result.stdout)
    if report.get('status') != 'VERIFIED':
        raise ValueError(f'Audit did not verify: {report!r}')
    return report


def replay(bundle):
    results = []
    for phase, run in RUNS.items():
        for reference in ('both', 'scaled', 'pairinv'):
            root = bundle / phase / f'ic-producer-{reference}-{run}-1'
            report = checked_json([sys.executable, '-I',
                str(root / 'ic-producer-evidence/evaluator/producer/audit.py'), str(root)])
            expected = (9, 3) if reference == 'both' else (15, 5)
            if (report['reference'] != reference or
                    (report['verified_pairs'], report['rho_verified_pairs']) != expected):
                raise ValueError('Producer census changed')
            results.append(dict(phase=phase, kind='producer', reference=reference, audit=report))
        for name, tool, expected in (
                ('native', 'autolab.py', 12), ('native-failure', 'autolab.py', 6),
                ('tournament', 'tournament.py', 15),
                ('qualification-control', 'tournament.py', 66)):
            root = bundle / phase / f'ic-driver-{run}-1' / name
            report = checked_json([sys.executable, str(root / 'evaluator' / tool),
                                   'verify', '--round', str(root)])
            if report['trial_receipts'] != expected:
                raise ValueError('Driver census changed')
            results.append(dict(phase=phase, kind=name, audit=report))
        print(json.dumps(dict(phase=phase, audits=len(results))), flush=True)

    root = bundle / 'full/ic-reference-qualification-36127866931-1/tournament'
    report = checked_json([sys.executable, str(root / 'evaluator/tournament.py'),
                           'verify', '--round', str(root)])
    if report['trial_receipts'] != 1290:
        raise ValueError('Full qualification census changed')
    qualification = json.loads((root / 'qualification.json').read_text())
    if (qualification['status'] != 'DEVELOPMENT_REFERENCES_SELECTED' or
            qualification['promotion_eligible'] or qualification['improvement_rounds_used'] or
            qualification['heldout_data_used'] or len(qualification['table']) != 21 or
            any(not row['qualified'] or row['verified_runs'] != 45
                for row in qualification['table'])):
        raise ValueError('Qualification completion/scope changed')
    results.append(dict(phase='full', kind='reference-qualification', audit=report))
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bundle', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if sys.version_info[:2] != (3, 12):
        parser.error('Use Python 3.12; strict floating summaries also require the Linux runtime')
    evidence = Path(__file__).resolve().parents[2] / 'evidence'
    entry = next(row for row in json.loads((evidence / 'manifest.json').read_text())['archives']
                 if row['file'] == ARCHIVE)
    with (evidence / ARCHIVE).open('rb') as stream:
        digest = hashlib.file_digest(stream, 'sha256').hexdigest()
    if digest != entry['sha256']:
        raise ValueError('Retained archive hash differs from its manifest')
    results = replay(args.bundle.resolve())
    if len(results) != 22:
        raise ValueError('Incomplete audit list')
    args.out.write_text(json.dumps(dict(status='VERIFIED', archive=entry,
        platform=sys.platform, python=sys.version, audits=results), indent=2) + '\n')
    print(json.dumps(dict(status='VERIFIED', fresh_restore_audits=len(results),
        full_qualification_pairs=1290, promotion_eligible=False)), flush=True)


if __name__ == '__main__':
    main()
