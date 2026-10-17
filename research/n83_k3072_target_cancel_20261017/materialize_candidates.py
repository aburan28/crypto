#!/usr/bin/env python3
"""Materialize exact K3072 serial/cancellation IC1 identities."""

import copy
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
TEMPLATE = ROOT / ('research/n83_k3072_current_source_20261017/'
                   'candidate_IC1N83Ckb1a1fb509952PDP4rootRCguidedLAgaussTDdirectISO0hec7ae07d8b9a.json')
SOURCE = ROOT / 'examples/koblitz_orbit_dlp_fast.rs'
SOURCE_SHA256 = '8109fb33374ce2aa66c6a5a0b6e144d3b0143168ac756f9e0bd70af60a911829'
SOURCE_COMMIT = '99b7d7c8fcb201f11dc910d15e031ac4d137c5a5'
PREFIX = 'IC1N83Ckb1a1fb509952PDP4rootRCguidedLAgaussTDdirectISO0h'
POLICIES = {'serial': (1, False), 'parallel_cancel': (2, True)}


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'),
                      ensure_ascii=False).encode()


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def materialize():
    if digest(SOURCE) != SOURCE_SHA256:
        raise ValueError('executed target-cancellation source differs')
    template = json.loads(TEMPLATE.read_text())
    old_id = template.pop('candidate_id')
    if old_id != PREFIX + hashlib.sha256(canonical(template)).hexdigest()[:12]:
        raise ValueError('K3072 template identity differs')
    result = {}
    for policy, (threads, cancel) in POLICIES.items():
        value = copy.deepcopy(template)
        implementation = value['implementation']
        implementation['source_commit'] = SOURCE_COMMIT
        implementation['source_sha256']['ic'] = SOURCE_SHA256
        implementation['algorithm_flags']['target_threads'] = threads
        implementation['algorithm_flags']['target_cancel'] = cancel
        value['point_decomposition']['solver_source_sha256'] = SOURCE_SHA256
        value['relation_collection']['source_digest'] = SOURCE_SHA256
        value['relation_linear_algebra']['source_digest'] = SOURCE_SHA256
        value['target_descent']['source_digest'] = SOURCE_SHA256
        value['target_descent']['worker_policy'] = policy
        candidate_id = PREFIX + hashlib.sha256(canonical(value)).hexdigest()[:12]
        value['candidate_id'] = candidate_id
        path = HERE / f'candidate_{policy}.json'
        blob = json.dumps(value, sort_keys=True, indent=2) + '\n'
        if path.exists() and path.read_text() != blob:
            raise ValueError(f'candidate manifest differs: {path}')
        path.write_text(blob)
        result[policy] = {'candidate_id': candidate_id,
                          'manifest_sha256': digest(path)}
    output = HERE / 'candidate_ids.json'
    blob = json.dumps(result, sort_keys=True, indent=2) + '\n'
    if output.exists() and output.read_text() != blob:
        raise ValueError('candidate ID receipt differs')
    output.write_text(blob)
    return result


if __name__ == '__main__':
    print(json.dumps(materialize(), sort_keys=True))
