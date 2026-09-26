"""Replay a transported bounded round and its controls without executing workers."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys


def read(path):
    return json.loads(path.read_text())


def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def control_audit(bundle, campaign):
    # The main audit has already checked every byte in this sealed evaluator.
    # Use its code for the earlier fixed-vector controls, too.
    sys.path.insert(0, str(campaign/'evaluator'))
    from driver_admission import make_admission, run_record
    from identity import sha256
    from oracle import require
    from tournament import parse_profiles
    controls = bundle/'controls/ic-candidate-controls'
    prepared = bundle/'controls/ic-candidates'
    preparation = read(prepared/'preparation.json')
    manifest = read(prepared/'source-manifest.json')
    panel = read(bundle/'round/registered-panel.json')
    source = sha256(manifest)
    require(source == preparation['source_manifest_sha256'] == panel['candidate_source_sha256'],
            'control source differs from registered candidate')
    for name, expected in manifest.items():
        require(digest(prepared/'source'/name) == expected, 'changed control source '+name)
    expected_inventory = {str(p.relative_to(prepared/'source')) for p in (prepared/'source').rglob('*') if p.is_file()}
    require(expected_inventory == set(manifest), 'control source inventory differs')
    outcomes = []
    for vector in ('n13-public', 'n23'):
        for arm in panel['candidates'][1:]:
            directory = controls/vector/arm['id']
            job, admitted = read(directory/'job.json'), read(directory/'admission.json')
            require(job['config'] == arm['config'], 'changed registered control policy')
            require(admitted['build']['compiler'].startswith('rustc 1.94.1 ') and
                    admitted['build']['target'] == 'x86_64-unknown-linux-musl', 'changed control build')
            resources = admitted['workload']['record']['resource_envelope']
            metadata = dict(preparation=preparation, build=admitted['build'])
            expected = make_admission(job=job, fixture=admitted['fixture'],
                report=read(directory/'inventory/stdout.json'), manifest=manifest, metadata=metadata,
                resources=resources, worker_sha256=admitted['worker_sha256'])
            require(expected == admitted, 'control admission differs')
            for mode in ('inventory', 'native', 'profile'):
                process = read(directory/mode/'process.json')
                require(process['exit_code'] == 0 and process['process_status'] == 'EXITED',
                        'failed control process')
                require(process['cpu'] == resources['cpu'] and
                        process['memory_cap_bytes'] == resources['memory_bytes'], 'changed control resources')
            native, profile = [read(directory/mode/'stdout.json') for mode in ('native', 'profile')]
            costs = parse_profiles(directory/'profile',
                compressed=bool(list((directory/'profile').glob('callgrind.out*.gz'))), phase_schema=3)
            stored = read(directory/'run.json')
            actual = run_record(admitted, number=0, host_id=stored['provenance']['host_id'], status='complete',
                native=native, process_wall_ns=read(directory/'native/process.json')['process_wall_ns'],
                profile=profile, costs=costs)
            require(actual == stored, 'control record differs from independent replay')
            outcome = dict(vector=vector+'.json', arm=arm['id'], status='VERIFIED',
                           candidate_id=stored['candidate_id'])
            require(read(directory/'receipt.json') == outcome, 'changed control receipt')
            outcomes.append(outcome)
    summary = read(controls/'summary.json')
    require(summary['scheduled'] == summary['verified'] == len(outcomes) == 30 and
            summary['outcomes'] == outcomes and summary['promotion_eligible'] is False,
            'changed control summary')
    return dict(status='VERIFIED', pairs=len(outcomes),
                scope='fixed-vector mathematical certificates, policy admission and phase closure; no performance selection')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bundle', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if sys.platform != 'linux' or sys.version_info[:2] != (3, 12):
        parser.error('Exact archived floating-summary replay requires Linux/Python 3.12')
    bundle = args.bundle.resolve()
    campaign = bundle/'round/tournament'
    result = subprocess.run([sys.executable, str(campaign/'evaluator/tournament.py'),
                             'verify', '--round', str(campaign)], capture_output=True, text=True)
    if result.returncode:
        raise RuntimeError(f'Frozen tournament audit failed:\n{result.stdout}\n{result.stderr}')
    audit = json.loads(result.stdout)
    contract = read(campaign/'contract.json')
    decision = read(campaign/'decision.json')
    expected_stages = ['aa', 'smoke', 'development', 'selection', 'confirmation', 'replay']
    if (audit['status'] != 'VERIFIED' or contract['stages'] != expected_stages or
            contract['attempt_number'] != decision['attempt_number'] or
            not all((campaign/'summaries'/f'{stage}.json').is_file() for stage in expected_stages)):
        raise ValueError('Round did not complete every declared stage')
    controls = control_audit(bundle, campaign)
    report = dict(status='VERIFIED', round=contract['attempt_number'],
        contract_sha256=digest(campaign/'contract.json'), decision_sha256=digest(campaign/'decision.json'),
        tournament=audit, controls=controls, decision_status=decision['status'],
        promotion_eligible=decision['promotion_eligible'], python=sys.version, platform=sys.platform)
    with args.out.open('x') as stream:
        json.dump(report, stream, indent=2); stream.write('\n')
    print(json.dumps(report), flush=True)


if __name__ == '__main__':
    main()
