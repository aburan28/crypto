#!/usr/bin/env python3
"""Read-only replay of the sealed, harmless Linux preparation artifact."""
from __future__ import annotations

import hashlib
import json
import stat
import subprocess
from pathlib import Path

from check_parent import main as check_parent

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
ARCHIVE = HERE / 'preparation_run_1'
PREPARATION_HEAD = 'ba8aeb651edba2f7a4ff20f336f499d0c57adf44'
RUN_ID = 36220351050


def require(condition: bool, message: str) -> None:
    if not condition:
        raise SystemExit('NOT_ADMITTED: ' + message)


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(path: Path) -> dict:
    return json.loads(path.read_text())


def main() -> None:
    check_parent()
    freeze = load(HERE / 'LINUX_PREPARATION_FREEZE.json')
    prepare = load(HERE / 'PREPARE.json')
    require(freeze['schema'] == 'k0-dag-dimacs-linux-preparation-freeze-v1', 'freeze schema')
    require(freeze['status'] == 'PASS_HARMLESS_PREPARATION_ONLY', 'freeze status')
    require(freeze['measured_attempt2_admitted'] is False, 'measured release flag')
    require(freeze['v2_release_head'] is None and freeze['measured_archive'] is None,
            'measured release fields')
    require(prepare['measured_attempt2_admitted'] is False and
            prepare['cadical']['linux_binary_sha256'] is None and
            prepare['cadical']['linux_build_receipt_sha256'] is None,
            'design snapshot was modified')
    require(sha(HERE / 'PREPARE.json') == freeze['prepare_sha256'], 'PREPARE.json hash')
    for name in ('attempt1_record_head', 'attempt1_freeze_sha256',
                 'attempt1_receipt_sha256', 'attempt1_manifest_sha256'):
        require(freeze[name] == prepare[name], f'{name} mismatch')
    require(freeze['preparation_checkout_head'] == PREPARATION_HEAD,
            'preparation event head')
    require(subprocess.run(['git', 'merge-base', '--is-ancestor', PREPARATION_HEAD, 'HEAD'],
                           cwd=ROOT, capture_output=True).returncode == 0,
            'preparation event head is not an ancestor')

    archive_meta = freeze['archive']
    require(archive_meta['path'] == ARCHIVE.name and
            archive_meta['github_actions_run_id'] == RUN_ID and
            archive_meta['github_actions_job_id'] == 108344308874 and
            archive_meta['github_artifact_id'] == 10899366900 and
            archive_meta['github_artifact_name'] == 'cadical-linux-prep-36220351050',
            'preparation run metadata')
    require(archive_meta['github_reported_artifact_digest'] ==
            'sha256:cad1826b51f58d98d2042449d3afd85d19f35e2c5da192807cd26bffe0b8f6a2',
            'GitHub-reported artifact digest')
    manifest_path = ARCHIVE / 'MANIFEST.json'
    require(sha(manifest_path) == archive_meta['manifest_sha256'], 'manifest hash')
    manifest = load(manifest_path)
    require(manifest['schema'] == 'k0-dag-dimacs-linux-preparation-manifest-v1',
            'manifest schema')
    entries = manifest['files']
    require(len(entries) == archive_meta['raw_file_count'] == 18,
            'raw artifact count')
    names = [entry['path'] for entry in entries]
    require(names == sorted(set(names)), 'manifest paths not unique and sorted')
    actual = {p.relative_to(ARCHIVE).as_posix() for p in ARCHIVE.rglob('*') if p.is_file()}
    require(actual == set(names) | {'MANIFEST.json'}, 'missing or extra archive file')
    for entry in entries:
        name = entry['path']
        rel = Path(name)
        require(not rel.is_absolute() and '..' not in rel.parts and name != 'MANIFEST.json',
                'unsafe archive path')
        path = ARCHIVE / rel
        require(not path.is_symlink(), f'symlink in archive: {name}')
        require(path.stat().st_size == entry['bytes'] and sha(path) == entry['sha256'] and
                format(stat.S_IMODE(path.stat().st_mode), '04o') == entry['mode'],
                f'archive byte/mode drift: {name}')

    probe_path = ARCHIVE / 'rlimit_probe.json'
    build_path = ARCHIVE / 'cadical_prepare/PREPARE_RECEIPT.json'
    binary_path = ARCHIVE / 'cadical_prepare/cadical'
    require(sha(probe_path) == freeze['cap_probe_receipt_sha256'], 'cap receipt hash')
    require(sha(build_path) == freeze['build_receipt_sha256'], 'build receipt hash')
    require(freeze['linux_binary_path'] == binary_path.relative_to(HERE).as_posix() and
            sha(binary_path) == freeze['linux_binary_sha256'] and
            binary_path.stat().st_size == freeze['linux_binary_bytes'],
            'committed Linux binary hash/size')
    binary = binary_path.read_bytes()
    require(binary[:6] == b'\x7fELF\x02\x01' and
            int.from_bytes(binary[18:20], 'little') == 62,
            'binary is not ELF64 x86-64')

    probe = load(probe_path)
    build = load(build_path)
    require(probe['schema'] == 'k0-dag-dimacs-linux-rlimit-probe-v1' and
            probe['decision'] == 'PASS_HARMLESS_CAPS_ONLY', 'cap decision')
    require(build['schema'] == 'k0-dag-dimacs-linux-cadical-build-v1' and
            build['decision'] == 'PASS_BUILD_ONLY_NO_MEASURED_CHILD', 'build decision')
    for receipt in (probe, build):
        require(receipt['checkout_head'] == PREPARATION_HEAD and
                receipt['machine'] == freeze['machine'] == 'x86_64' and
                receipt['runner_image_os'] == freeze['runner_image_os'] == 'ubuntu24' and
                receipt['runner_image_version'] == freeze['runner_image_version'] and
                receipt['python_executable'].endswith('/python3'),
                'runner or checkout provenance')
    require(probe['python'].startswith('3.12.') and
            build['python_version'] == probe['python'], 'Python version')
    require(build['prepare_sha256'] == freeze['prepare_sha256'] and
            build['probe_receipt_sha256'] == freeze['cap_probe_receipt_sha256'] and
            build['attempt1_receipt_sha256'] == freeze['attempt1_receipt_sha256'],
            'build input provenance')
    caps = prepare['target']['hard_rlimit_as_caps_bytes']
    require(caps == freeze['caps_bytes'] == [536870912, 1073741824, 2147483648],
            'cap schedule changed')
    require(len(probe['caps']) == len(caps), 'cap observation count')
    for expected, observed in zip(caps, probe['caps']):
        require(observed['cap_bytes'] == expected and
                observed['true_exit_code'] == observed['python_exit_code'] == 0 and
                observed['observation']['limit'] == [expected, expected] and
                observed['observation']['over_cap_allocation'] == 'REJECTED' and
                observed['observation']['exception'] == 'OSError',
                f'hard cap not observed at {expected}')

    cadical = prepare['cadical']
    require((freeze['source_commit'], freeze['source_tree'],
             freeze['source_manifest_sha256']) ==
            (cadical['commit'], cadical['tree'], cadical['source_manifest_sha256']),
            'source identity')
    require(len(build['builds']) == 2 and
            [b['label'] for b in build['builds']] == ['a', 'b'],
            'two independent builds not recorded')
    for name, item in zip(('build_a_binary_sha256', 'build_b_binary_sha256'),
                          build['builds']):
        source = item['source']
        require(source['commit'] == cadical['commit'] and
                source['tree'] == cadical['tree'] and
                source['tracked_files'] == cadical['tracked_files'] and
                source['source_manifest_sha256'] == cadical['source_manifest_sha256'] and
                source['version'] == item['version_output'] == '3.0.1',
                f'build {item["label"]} source/version')
        require(item['binary_sha256'] == freeze[name] == freeze['linux_binary_sha256'] and
                item['binary_bytes'] == freeze['linux_binary_bytes'],
                f'build {item["label"]} binary')
    require(build['linux_binary_sha256'] == freeze['linux_binary_sha256'] and
            build['linux_binary_bytes'] == freeze['linux_binary_bytes'],
            'selected binary receipt')
    for name in ('gcc', 'gxx', 'ar', 'make', 'ldd', 'os_release', 'uname'):
        require(isinstance(build['toolchain'].get(name), str) and
                bool(build['toolchain'][name]), f'missing toolchain {name}')
    require('Ubuntu 24.04.5 LTS' in build['toolchain']['os_release'] and
            '13.3.0' in build['toolchain']['gcc'] and
            'GLIBC 2.39' in build['toolchain']['ldd'],
            'unexpected toolchain')
    print(json.dumps({'decision': 'PASS_HARMLESS_PREPARATION_ONLY',
                      'run_id': RUN_ID,
                      'raw_files': len(entries),
                      'manifest_sha256': sha(manifest_path),
                      'binary_sha256': sha(binary_path),
                      'measured_attempt2_admitted': False}, sort_keys=True))


if __name__ == '__main__':
    main()
