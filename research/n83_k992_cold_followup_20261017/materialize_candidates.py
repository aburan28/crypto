#!/usr/bin/env python3
"""Materialize current-source IC1 manifests from replayed lower-K controls."""

import hashlib
import importlib.util
import json
import tarfile
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SCREEN = HERE / 'local_screen'
ARCHIVE_SHA256 = '9269fde6cd09fffd8866396be119c05356e75e3975f660d4d2c88a54cf8f0bc0'
SOURCE_COMMIT = '416489e09c3d8885b6c12f262935e6ba29373460'
SOURCE_SHA256 = '8a2fe55dd90dc33892c197d27260182bf783fc8212fb485d8297daf2bbba8e12'
EXPECTED = {
    992: 'IC1N83Ckb1a1fb164672PDP4rootRCguidedLAgaussTDdirectISO0h829ed476f094',
    1024: 'IC1N83Ckb1a1fb169984PDP4rootRCguidedLAgaussTDdirectISO0h3eb055926289',
    1107: 'IC1N83Ckb1a1fb183762PDP4rootRCguidedLAgaussTDdirectISO0h82329aa062fa',
}


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    value = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(value)
    return value


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def write_once(path, record):
    blob = json.dumps(record, sort_keys=True, indent=2) + '\n'
    if path.exists() and path.read_text() != blob:
        raise ValueError(f'candidate manifest changed: {path.name}')
    path.write_text(blob)


def main():
    archive = SCREEN / 'raw.tar.gz'
    if (digest(archive) != ARCHIVE_SHA256 or
            digest(ROOT / 'examples/koblitz_orbit_dlp_fast.rs') != SOURCE_SHA256):
        raise ValueError('archive or executed source differs')
    replay = json.loads((SCREEN / 'replay.json').read_text())
    if replay['status'] != 'PASS_INDEPENDENT_ARITHMETIC_REPLAY':
        raise ValueError('local arithmetic replay missing')
    rows = {row['k']: row for row in replay['rows']}
    finegrid = module('k992_finegrid_manifest',
                      ROOT / 'research/n83_e2e_finegrid_20261016/materialize.py')
    dynamic = finegrid.dynamic
    dynamic.SOURCE_COMMIT = SOURCE_COMMIT
    dynamic.panel.SOURCE_COMMIT = SOURCE_COMMIT
    dynamic.panel.SOURCE_SHA256 = SOURCE_SHA256
    field, curve = dynamic.panel.curve_record()
    with tempfile.TemporaryDirectory(prefix='n83-k992-manifests-') as temp:
        with tarfile.open(archive, 'r:gz') as saved:
            for k, expected_id in EXPECTED.items():
                folder = f'k{k}root11'
                summary_blob = saved.extractfile(f'{folder}/stdout.jsonl').read()
                base_blob = saved.extractfile(f'{folder}/base_dump.jsonl').read()
                if (hashlib.sha256(summary_blob).hexdigest() != rows[k]['stdout_sha256'] or
                        hashlib.sha256(base_blob).hexdigest() != rows[k]['base_dump_sha256']):
                    raise ValueError(f'K{k}: archived inputs differ from replay')
                summary = json.loads(summary_blob)
                base_path = Path(temp) / f'k{k}-base.jsonl'
                base_path.write_bytes(base_blob)
                record = dynamic.candidate(k, 'root11_dynamic', summary,
                                           base_path, field, curve)
                if (record['candidate_id'] != expected_id or
                        record['factor_base']['actual_usable_points'] != 166 * k or
                        record['factor_base']['source_base_hash_blake3'] != rows[k]['base_hash']):
                    raise ValueError(f'K{k}: IC1 candidate identity differs')
                write_once(HERE / f'candidate_{expected_id}.json', record)
    print(json.dumps({'status': 'PASS', 'candidate_ids': list(EXPECTED.values())},
                     sort_keys=True))


if __name__ == '__main__':
    main()
