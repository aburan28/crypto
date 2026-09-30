#!/usr/bin/env python3
"""Retain the closed local diagnostic panel; never rerun a measured arm."""
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
AUDITOR_COMMIT = '6958067b53f110be39a7cbcd9627aab43908f906'
AUDITORS = ('audit_static_sat_full.py', 'audit_static_sat_full_v2.py',
            'audit_static_cms_s4_natural.py')
RUNS = {
    'sat-v1': 'ic-static-sat-full-20260929-run0',
    'sat-v2': 'ic-static-sat-paired-v2-20260929-run0',
    'f5': 'ic-fresh-paired-f5-20260929-run0',
    'incumbent': 'ic-fresh-paired-pairinv-20260929-run0',
    'rho-selected': 'ic-fresh-paired-pairinv-rho-20260929-run0',
    'rho-generic': 'ic-fresh-paired-rho-20260929-run0',
}
EXTRAS = {
    'sat-v1-audit.json': 'ic-static-sat-full-20260929-audit.json',
    'sat-v2-audit.json': 'ic-static-sat-paired-v2-20260929-audit.json',
    'sat-v2-admissibility-review.json':
        'ic-static-sat-paired-v2-20260929-admissibility-review.json',
    'sat-v2-load.jsonl': 'ic-static-sat-paired-v2-20260929-host-load-observations.jsonl',
    'f5-load.jsonl': 'ic-fresh-paired-f5-20260929-host-load-observations.jsonl',
}


def pack(raw_root, prepared, output):
    raw_root, prepared, output = map(Path, (raw_root, prepared, output))
    require(not output.exists(), 'paired evidence output already exists')
    files = {}
    for alias, name in RUNS.items():
        directory = raw_root/name
        terminal = directory/('summary.json' if alias.startswith('sat-')
                              else 'result.json')
        result = json.loads(terminal.read_text())
        require(result['status'] in {'COMPLETE', 'AUDITED_COMPLETE',
                                     'AUDITED_BOUNDED_INCOMPLETE'},
                'paired run is not terminal: '+alias)
        for path in sorted(directory.rglob('*')):
            require(not path.is_symlink(), 'symlinked raw paired evidence')
            if path.is_file():
                data = path.read_bytes()
                require(path.read_bytes() == data, 'raw paired evidence changed')
                files[alias+'/'+path.relative_to(directory).as_posix()] = data
    for name, original in EXTRAS.items():
        files['reviews/'+name] = (raw_root/original).read_bytes()
    build = json.loads(files['incumbent/build-record.json'])
    for name, digest in build['source_manifest'].items():
        path = prepared/'source'/name
        require(path.is_file() and not path.is_symlink(),
                'prepared pairinv source missing or symlinked')
        data = path.read_bytes()
        require(hashlib.sha256(data).hexdigest() == digest,
                'prepared pairinv source differs from build receipt')
        files['pairinv-source/'+name] = data
    files['pairinv-preparation.json'] = (prepared/'preparation.json').read_bytes()
    require(sha256(json.loads(files['pairinv-preparation.json']))
            == sha256(build['prepared_source_receipt']),
            'pairinv preparation differs from measured build')
    output.mkdir(parents=True)
    archive = output/'evidence.tar.gz'
    with archive.open('xb') as raw, gzip.GzipFile(
            fileobj=raw, mode='wb', mtime=0, filename='') as compressed:
        with tarfile.open(fileobj=compressed, mode='w|') as tar:
            for name, data in sorted(files.items()):
                item = tarfile.TarInfo(name)
                item.size, item.mode, item.mtime = len(data), 0o644, 0
                tar.addfile(item, io.BytesIO(data))
    inventory = {name: dict(bytes=len(data),
                           sha256=hashlib.sha256(data).hexdigest())
                 for name, data in sorted(files.items())}
    receipt = dict(schema_version=1, archive_file='evidence.tar.gz',
                   archive_bytes=archive.stat().st_size,
                   archive_sha256=hashlib.sha256(archive.read_bytes()).hexdigest(),
                   inventory=inventory, retention='all closed raw run files',
                   claim_scope='correctness and paired natural-query diagnostics',
                   promotion_eligible=False, online_speedup=None)
    write_immutable(output/'receipt.json', receipt)
    audit_dir = output/'audit-source'
    audit_dir.mkdir()
    hashes = {}
    for name in AUDITORS:
        role = 'research/ic_candidate_tournament_20260915/'+name
        data = subprocess.check_output(['git', 'show', AUDITOR_COMMIT+':'+role],
                                       cwd=ROOT)
        (audit_dir/name).write_bytes(data)
        hashes[name] = hashlib.sha256(data).hexdigest()
    write_immutable(output/'auditor-receipt.json', dict(
        schema_version=1, source_commit=AUDITOR_COMMIT, files=hashes,
        stage='post-execution independent auditor source recovery',
        original_candidate_manifests_changed=False))
    return receipt


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--raw-root', type=Path, required=True)
    parser.add_argument('--prepared', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    receipt = pack(args.raw_root, args.prepared, args.out)
    print(json.dumps({key: receipt[key] for key in
                      ('archive_bytes', 'archive_sha256')}, sort_keys=True))
