#!/usr/bin/env python3
"""Bind the immutable first-run ZIP and replay it with a versioned path-order correction."""
from __future__ import annotations

import hashlib
import importlib.util
import json
import shutil
import stat
import subprocess
import tempfile
import zipfile
from pathlib import Path

import corrected_replay as successor

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent / 'symbolic_dag_dimacs_linux_v2_20260928'
REPO = HERE.parents[3]
BUNDLE = HERE / 'evidence/first_outcome_36534548089.zip'
ROOT_NAME = 'dag-linux-v2-36534548089'
ORDER_ERROR = 'manifest omits or invents an archive file'


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as stream:
        while chunk := stream.read(1 << 20):
            h.update(chunk)
    return h.hexdigest()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def freeze_gate() -> dict:
    frozen = json.loads((HERE / 'FROZEN.json').read_text())
    require(frozen['schema'] == 'dag-linux-v2-first-outcome-replay-v1', 'successor schema')
    require(frozen['original_pr'] == 921 and
            frozen['original_head'] == 'd70396a649cf4a12ece39f01bf60f121c6b2a9ae',
            'original exact head')
    check = subprocess.run(['git', 'merge-base', '--is-ancestor',
                            frozen['original_head'], 'HEAD'], cwd=REPO,
                           capture_output=True)
    require(check.returncode == 0, 'successor does not descend from original head')
    paths = {
        'original_freeze_sha256': SOURCE / 'FROZEN.json',
        'original_verifier_sha256': SOURCE / 'ci_replay.py',
        'successor_verifier_sha256': HERE / 'corrected_replay.py',
        'gate_sha256': HERE / 'gate.py',
        'protocol_sha256': HERE / 'PROTOCOL.md',
        'workflow_sha256': REPO / '.github/workflows/ecc2k130-dag-linux-v2-outcome-replay.yml',
        'artifact_zip_sha256': BUNDLE,
    }
    for key, path in paths.items():
        require(sha(path) == frozen[key], f'freeze drift: {key}')
    require(frozen['artifact']['run_id'] == 36534548089 and
            frozen['artifact']['artifact_id'] == 11017624684 and
            frozen['artifact']['name'] == ROOT_NAME and
            frozen['artifact']['digest'] == 'sha256:' + frozen['artifact_zip_sha256'] and
            frozen['artifact']['size_bytes'] == BUNDLE.stat().st_size,
            'source artifact metadata drift')
    original = (SOURCE / 'ci_replay.py').read_text()
    old_here = 'HERE = Path(__file__).resolve().parent\n'
    new_here = ("HERE = Path(__file__).resolve().parent.parent / "
                "'symbolic_dag_dimacs_linux_v2_20260928'\n")
    old_order = "require([row['path'] for row in manifest] == sorted(actual),"
    new_order = ("require([row['path'] for row in manifest] == "
                 "sorted(actual, key=lambda rel: Path(rel).parts),")
    require(original.count(old_here) == original.count(old_order) == 1,
            'original patch anchor changed')
    expected = original.replace(old_here, new_here).replace(old_order, new_order)
    require((HERE / 'corrected_replay.py').read_text() == expected,
            'successor verifier changes more than source location and manifest order')
    return frozen


def extract_original(root: Path, frozen: dict) -> Path:
    with zipfile.ZipFile(BUNDLE) as archive:
        members = archive.infolist()
        require(len(members) == 115 and len({m.filename for m in members}) == 115,
                'ZIP member count/uniqueness')
        require(archive.testzip() is None, 'ZIP CRC failed')
        for member in members:
            parts = Path(member.filename).parts
            mode = member.external_attr >> 16
            require(len(parts) >= 2 and parts[0] == ROOT_NAME and
                    not Path(member.filename).is_absolute() and
                    '..' not in parts and not member.is_dir() and
                    not stat.S_ISLNK(mode), 'unsafe ZIP member')
            target = root.joinpath(*parts)
            target.parent.mkdir(parents=True, exist_ok=True)
            with archive.open(member) as source, target.open('wb') as dest:
                shutil.copyfileobj(source, dest)
    outcome = root / ROOT_NAME
    receipt = outcome / 'receipt.json'
    manifest = outcome / 'MANIFEST.json'
    require(sha(receipt) == frozen['receipt_sha256'] and
            sha(manifest) == frozen['manifest_sha256'],
            'original receipt/manifest byte drift')
    rows = json.loads(manifest.read_text())
    require(len(rows) == frozen['payload_count'] == 113 and
            sum(row['bytes'] for row in rows) == frozen['payload_bytes'] == 7444986,
            'original payload inventory drift')
    return receipt


def original_module():
    spec = importlib.util.spec_from_file_location('original_frozen_replay',
                                                  SOURCE / 'ci_replay.py')
    require(spec is not None and spec.loader is not None, 'original verifier import')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def expect_order_rejection(verifier, receipt: Path, original_frozen: dict) -> None:
    try:
        verifier.check_archive(receipt, original_frozen)
    except AssertionError as exc:
        require(str(exc) == ORDER_ERROR, f'unexpected archive rejection: {exc}')
    else:
        raise AssertionError('altered/omitted manifest path accepted')


def path_mutation_controls(receipt: Path, original_frozen: dict) -> list[str]:
    manifest_path = receipt.parent / 'MANIFEST.json'
    manifest_raw, receipt_raw = manifest_path.read_bytes(), receipt.read_bytes()
    rows = json.loads(manifest_raw)
    labels = []
    try:
        for label, altered in (
            ('omitted_path', rows[:-1]),
            ('invented_path', [dict(rows[0], path='invented/file')] + rows[1:]),
        ):
            manifest_path.write_text(json.dumps(altered, sort_keys=True, indent=2) + '\n')
            r = json.loads(receipt_raw)
            r['manifest_sha256'] = sha(manifest_path)
            receipt.write_text(json.dumps(r, sort_keys=True, indent=2) + '\n')
            expect_order_rejection(successor, receipt, original_frozen)
            labels.append(label)
    finally:
        manifest_path.write_bytes(manifest_raw)
        receipt.write_bytes(receipt_raw)
    return labels


def main() -> None:
    frozen = freeze_gate()
    original = original_module()
    original_frozen = successor.check_static()
    with tempfile.TemporaryDirectory(prefix='dag-linux-v2-replay-') as tmp:
        receipt = extract_original(Path(tmp), frozen)
        expect_order_rejection(original, receipt, original_frozen)
        controls = path_mutation_controls(receipt, original_frozen)
        require(sha(receipt) == frozen['receipt_sha256'] and
                sha(receipt.parent / 'MANIFEST.json') == frozen['manifest_sha256'],
                'control did not restore original archive')
        result = successor.check_archive(receipt, original_frozen)
    require(result['decision'] == 'PASS_INDEPENDENT_REPLAY',
            'successor semantic replay did not pass')
    print(json.dumps({'decision': 'PASS_SUCCESSOR_SEMANTIC_REPLAY',
                      'original_verifier': 'FAIL_MANIFEST_ORDER_ONLY',
                      'successor': result, 'path_mutation_controls': controls,
                      'artifact_sha256': frozen['artifact_zip_sha256'],
                      'receipt_sha256': frozen['receipt_sha256'],
                      'manifest_sha256': frozen['manifest_sha256']}, sort_keys=True))


if __name__ == '__main__':
    main()
