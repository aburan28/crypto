#!/usr/bin/env python3
"""Recover registered SAT Python sources from their measured Git commits."""
import argparse
import gzip
import hashlib
import io
import json
from pathlib import Path
import subprocess
import tarfile

from identity import sha256, write_immutable
from oracle import require

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
REGISTRATIONS = {
    'v1': ('static-sat-full', 'fa1e9f8fcc68cbc63ce75861f3989bbe16b5862d'),
    'v2': ('static-sat-paired-v2', '6958067b53f110be39a7cbcd9627aab43908f906'),
}
SUPPLEMENTAL_ROLES = (
    'research/ic_candidate_tournament_20260915/producer/evidence.py',
    'research/ic_candidate_tournament_20260915/producer/timing.py',
)


def freeze(out):
    out = Path(out)
    require(not out.exists(), 'frozen SAT runtime output already exists')
    files, registrations = {}, {}
    for version, (directory, commit) in REGISTRATIONS.items():
        source = json.loads((HERE/'goal_20260924'/directory/'source-manifest.json').read_text())
        roles = {item['role']: item['sha256'] for item in source['components']
                 if item['role'].endswith('.py')}
        require(roles, 'registered SAT runtime has no Python sources')
        for role, expected in roles.items():
            require(role.startswith(('research/', 'scripts/'))
                    and '..' not in Path(role).parts, 'unsafe registered source role')
            data = subprocess.check_output(['git', 'show', commit+':'+role], cwd=ROOT)
            require(hashlib.sha256(data).hexdigest() == expected,
                    'measured commit differs from registered source: '+role)
            require(role not in files or files[role] == data,
                    'shared registered SAT source differs between versions')
            files[role] = data
        supplemental = {}
        for role in SUPPLEMENTAL_ROLES:
            require(role not in roles, 'supplemental source is already registered')
            data = subprocess.check_output(['git', 'show', commit+':'+role], cwd=ROOT)
            require(role not in files or files[role] == data,
                    'shared supplemental SAT source differs between versions')
            files[role] = data
            supplemental[role] = hashlib.sha256(data).hexdigest()
        registrations[version] = dict(
            directory=directory, source_commit=commit,
            source_manifest_sha256=sha256(source), files=roles,
            supplemental_files=supplemental,
            complete_preexecution_python_manifest=False,
            supplemental_provenance='post-execution recovery from recorded Git commit')
    out.mkdir(parents=True)
    archive = out/'runtime.tar.gz'
    with archive.open('xb') as raw, gzip.GzipFile(fileobj=raw, mode='wb', mtime=0,
                                                filename='') as compressed:
        with tarfile.open(fileobj=compressed, mode='w|') as tar:
            for role, data in sorted(files.items()):
                item = tarfile.TarInfo(role)
                item.size, item.mode, item.mtime = len(data), 0o644, 0
                tar.addfile(item, io.BytesIO(data))
    receipt = dict(schema_version=1,
                   scope='post-execution source recovery; registered hashes unchanged',
                   archive_sha256=hashlib.sha256(archive.read_bytes()).hexdigest(),
                   files={role: hashlib.sha256(data).hexdigest()
                          for role, data in sorted(files.items())},
                   registrations=registrations)
    write_immutable(out/'receipt.json', receipt)
    return receipt


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    receipt = freeze(args.out)
    print(json.dumps({'archive_sha256': receipt['archive_sha256'],
                      'source_files': len(receipt['files'])}, sort_keys=True))
