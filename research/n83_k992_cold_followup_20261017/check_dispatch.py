#!/usr/bin/env python3
"""Check a cloud dispatch against the sealed K992 follow-up protocol."""

import argparse
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent

ARMS = (
    'IC1N83Ckb1a1fb164672PDP4rootRCguidedLAgaussTDdirectISO0h829ed476f094',
    'IC1N83Ckb1a1fb169984PDP4rootRCguidedLAgaussTDdirectISO0h3eb055926289',
    'IC1N83Ckb1a1fb183762PDP4rootRCguidedLAgaussTDdirectISO0h82329aa062fa',
)
SOURCE_SHA256 = '8a2fe55dd90dc33892c197d27260182bf783fc8212fb485d8297daf2bbba8e12'
TARGETS_SHA256 = '89c7208dd4fc7a5621cf727b3fac9a6714f04f66a312e516b418c31764e7547e'
ANSWERS_SHA256 = '363adc7f9bded3f4a25a80c7b27da5a2305afb30ee6741932cbe76f455c2865b'


def sha256(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def expected():
    files = ['run.py', 'verify.py', 'select.py', 'aggregate.py',
             'test_protocol.py', 'generate.py', 'materialize_candidates.py',
             'fixture_manifest.json',
             'targets.jsonl', 'answers.jsonl']
    files += [f'candidate_{candidate_id}.json' for candidate_id in ARMS]
    return {
        'schema': 'n83.k992-cold-followup-dispatch/v1',
        'source_sha256': SOURCE_SHA256,
        'targets_sha256': TARGETS_SHA256,
        'answers_sha256': ANSWERS_SHA256,
        'candidate_ids': list(ARMS),
        'runner_os': 'ubuntu-24.04',
        'host_count': 2,
        'points_per_host': 5,
        'ic_rho_pairs_per_host': 15,
        'calibrations_per_host': 30,
        'ic_cap_seconds': 480,
        'rho_cap_seconds': 120,
        'producer_cap_seconds_per_host': 7200,
        'workflow_sha256': sha256(
            HERE.parents[1] / '.github/workflows/n83-k992-cold-followup.yml'),
        'files_sha256': {name: sha256(HERE / name) for name in files},
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('dispatch', type=Path, nargs='?')
    parser.add_argument('--emit', action='store_true')
    args = parser.parse_args()
    frozen = expected()
    if args.emit:
        if args.dispatch:
            parser.error('--emit takes no dispatch path')
        print(json.dumps(frozen, sort_keys=True, indent=2))
        return
    if not args.dispatch:
        parser.error('dispatch path required')
    observed = json.loads(args.dispatch.read_text())
    if observed != frozen:
        missing = sorted(set(frozen) - set(observed))
        extra = sorted(set(observed) - set(frozen))
        differing = sorted(key for key in set(frozen) & set(observed)
                           if frozen[key] != observed[key])
        raise ValueError(f'dispatch differs: missing={missing}, extra={extra}, '
                         f'differing={differing}')
    print(json.dumps({'status': 'PASS_FROZEN_DISPATCH',
                      'dispatch_sha256': sha256(args.dispatch),
                      'candidate_ids': list(ARMS)}, sort_keys=True))


if __name__ == '__main__':
    main()
