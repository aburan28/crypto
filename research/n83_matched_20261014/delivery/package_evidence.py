#!/usr/bin/env python3
"""Package existing evidence bytes. No experiments or statistics are computed."""
import hashlib
import json
import shutil
import zipfile
from pathlib import Path


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def main():
    repo = Path(__file__).resolve().parents[3]
    study = repo / 'research/n83_matched_20261014'
    run = study / 'reports/run-20261014T1940Z'
    output = Path('/workspace/report_output')
    archive = output / 'n83_matched_ic_rho_evidence_bundle.zip'
    committed_archive = run / archive.name
    manifest = run / 'HASH_MANIFEST.json'
    receipt = study / 'delivery/EVIDENCE_ARCHIVE_RECEIPT.json'
    excluded = {manifest, receipt, committed_archive}
    if archive.exists() or committed_archive.exists():
        raise SystemExit('Refusing to overwrite an existing evidence archive')
    entries = {}

    def add(path, name):
        if path in excluded or '__pycache__' in path.parts or path.suffix == '.pyc':
            return
        if (path.parent.name in {'preview', 'combined_preview'}
                and path.name.startswith('page_') and path.suffix == '.png'):
            return  # Regenerable page rasters; retain PDFs, figure PNGs and QA receipts.
        if name in entries:
            raise ValueError(f'Duplicate archive name: {name}')
        entries[name] = path

    for path in sorted(study.rglob('*')):
        if path.is_file():
            add(path, str(path.relative_to(repo)))
    benchmark = repo / 'benchmarks/20261014-n83-matched-ic-rho'
    for path in sorted(benchmark.rglob('*')):
        if path.is_file():
            add(path, str(path.relative_to(repo)))
    for relative in ['AGENTS.md', 'benchmarks/README.md',
                     '.github/workflows/n83-scalar-arithmetic-lean.yml',
                     'scripts/check_curve_names.py',
                     'scripts/test_curve_names_archived_labels.py',
                     'docs/ic/progress-timeline.json',
                     'docs/index-calculus-scoreboard.html',
                     'scripts/build_ic_leaderboard.py',
                     'docs/ic-leaderboard.html', 'docs/ic/LEADERBOARD.md',
                     'docs/ic/leaderboard.json', 'docs/browser/data.json',
                     'docs/ic-current-state.html', 'scripts/site/panel_index.py']:
        add(repo / relative, relative)
    binaries = Path('/workspace/n83_matched/frozen_binaries')
    for name in ['runner', 'analysis', 'ic', 'rho_matched']:
        add(binaries / name, 'runtime_binaries/' + name)
    for version in ['v2', 'v3']:
        path = Path('/workspace/n83_matched/verification/verify_matched_campaign_' + version)
        add(path, 'runtime_binaries/' + path.name)
    add(Path('/workspace/n83_evaluation/build/release/n83-matched-supplement'),
        'runtime_binaries/n83-matched-supplement')
    add(Path('/workspace/n83_matched/collector_only_build/release/n83-matched-collect-only'),
        'runtime_binaries/n83-matched-collect-only')
    add(Path('/workspace/n83_evaluation/frozen_binaries/isolated_bench'),
        'runtime_binaries/isolated_bench')
    historical = output / 'n83_faster_rho_evidence_bundle.zip'
    with zipfile.ZipFile(historical) as previous:
        if previous.testzip() is not None:
            raise ValueError('Historical evidence CRC failure')
        historical_count = len(previous.infolist())
    add(historical, 'historical/' + historical.name)
    old_pdf = output / 'unified_curve_basis_report_faster_rho.pdf'
    if sha(old_pdf) != '856edd55dd4172374f0f60b23e8605230ce7e63ee41d025c2bc7e9e6d0a8fb0c':
        raise ValueError('Historical 70-page PDF identity changed')
    add(old_pdf, 'historical/' + old_pdf.name)
    payload = {
        'kind': 'recorded_evidence_hash_manifest',
        'research_results_recomputed': False,
        'scope': 'Existing files only; manifest and external archive receipt exclude themselves',
        'historical_inner_archive_crc_checked': True,
        'historical_inner_archive_file_count': historical_count,
        'files': [{'archive_path': name, 'bytes': path.stat().st_size,
                   'sha256': sha(path)} for name, path in sorted(entries.items())],
    }
    manifest.write_text(json.dumps(payload, indent=2) + '\n')
    # add() excludes the manifest during file enumeration; include it once here.
    entries[str(manifest.relative_to(repo))] = manifest
    with zipfile.ZipFile(archive, 'x', compression=zipfile.ZIP_DEFLATED,
                         compresslevel=6) as bundle:
        for name, path in sorted(entries.items()):
            bundle.write(path, name)
    with zipfile.ZipFile(archive) as bundle:
        bad = bundle.testzip()
        if bad is not None:
            raise ValueError(f'Evidence archive CRC failure: {bad}')
        if len(bundle.infolist()) != len(entries):
            raise ValueError('Archive file count mismatch')
        # Bind every archived payload to its already recorded manifest identity.
        for item in payload['files']:
            if hashlib.sha256(bundle.read(item['archive_path'])).hexdigest() != item['sha256']:
                raise ValueError('Archived file differs from manifest: ' + item['archive_path'])
    shutil.copyfile(archive, committed_archive)
    result = {
        'kind': 'recorded_evidence_packaging_receipt', 'status': 'PASS',
        'research_results_recomputed': False, 'crc_check': 'PASS',
        'manifest_entry_hashes_checked': len(payload['files']),
        'archive_files': len(entries), 'archive_bytes': archive.stat().st_size,
        'archive_sha256': sha(archive), 'manifest_sha256': sha(manifest),
        'historical_inner_archive_crc_check': 'PASS',
        'committed_copy_sha256': sha(committed_archive),
        'scope': 'File packaging and byte identity only; scientific results are native recorded JSON',
    }
    receipt.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
