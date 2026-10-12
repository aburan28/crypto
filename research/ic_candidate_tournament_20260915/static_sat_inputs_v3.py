"""Admit the retained n17 static-SAT native build, without live-path fallbacks.

This adapter recognizes the existing physical macOS ARM64 research build.
A different platform or build needs its own reviewed source/build adapter.
"""
import hashlib
import io
import json
from pathlib import PurePosixPath
import tarfile

from generic_bases import construct
from identity import curve_record, factor_base_inventory, sha256
from oracle import Curve, require

EXPORTER_SHA256 = '4da7f5781da1aba2f76d6b9ce8f6e4ea43845f2465cf8edf0612615931d0629a'
EXPORTER_SOURCE_SHA256 = '67acaff696500679a6239d47e74bba52b41c0c64fd0c7413ae7fbd7e00505b11'
RUST_MANIFEST_SHA256 = 'c64e4b3102bface63a2305efbff4bd85810cc112cb43546da2992ba48e9e85b7'
RUST_COMMIT = '765c3c5f19032bd852163805f257c56babef2040'
CMS_SHA256 = '6c509f09622f103d8a3ad90afc151e1c4031275c052d7f8af481465b4f27f2af'
CMS_RECEIPT_SHA256 = 'c2a9be07fb510e378b935433c8c091dd4f515e2f05e24154f8ef00e83fc191f4'
ROLES = {'fixture.json', 'base.json', 'bin/exporter', 'bin/cms',
         'exporter/source.rs', 'exporter/build-record.json', 'exporter/build.log',
         'exporter/build-exit.json', 'rust/source-manifest.json',
         'rust/root-source.tar.gz', 'cms/receipt.json', 'cms/source.tar',
         'cms/cadical.tar', 'cms/cadiback.tar'}


def digest(data):
    return hashlib.sha256(data).hexdigest()


def native_admission(files):
    require(set(files) == ROLES, 'SAT native input inventory differs from adapter')
    read = lambda role: json.loads(files[role])
    fixture = read('fixture.json')
    require(fixture['targets'] == [] and fixture['target_seeds'] == []
            and fixture['target_scalar_constructed'] is False,
            'target data leaked into reusable SAT inputs')
    curve = Curve(fixture)
    require(curve_record(fixture)['curve']['curve_id'] == 'EC1N17Ckb1hbbe2b5b6b1e6',
            'static SAT adapter requires the admitted exact n17a1 curve')
    base, _ = construct(curve, dict(kind='standard_subspace', dimension=6))
    report = read('base.json')
    require(set(report) == {'factor_base', 'columns', 'column_convention'}
            and report['factor_base'] == [list(point) for point in base]
            and report['columns'] == 29 and report['column_convention'] == 'cofactor',
            'retained SAT base differs from independent construction')
    inventory = factor_base_inventory(report, fixture)
    require(inventory['geometric_point_count'] == 63
            and inventory['usable_point_count'] == 62,
            'SAT geometric/usable factor-base count differs')
    build = read('exporter/build-record.json')
    source = read('rust/source-manifest.json')
    require(digest(files['bin/exporter']) == EXPORTER_SHA256
            and digest(files['exporter/source.rs']) == EXPORTER_SOURCE_SHA256
            and sha256(source) == RUST_MANIFEST_SHA256
            and build['source_commit'] == RUST_COMMIT
            and build['parent_source_manifest_sha256'] == RUST_MANIFEST_SHA256
            and build['exporter_source_sha256'] == EXPORTER_SOURCE_SHA256
            and build['exporter_executable_sha256'] == EXPORTER_SHA256
            and build['build_log_sha256'] == digest(files['exporter/build.log'])
            and read('exporter/build-exit.json')['exit_code'] == 0
            and build['build_command'] == ['cargo', 'build', '--locked', '--offline',
                                          '--release', '--no-default-features',
                                          '--example', 'koblitz_pdp_export']
            and build['platform'] == {'system': 'Darwin', 'machine': 'arm64'},
            'SAT exporter source/build receipt differs from accepted build')
    with tarfile.open(fileobj=io.BytesIO(files['rust/root-source.tar.gz']), mode='r:gz') as tar:
        members = tar.getmembers()
        require(len(members) == len(source['root_files'])
                and {m.name for m in members} == set(source['root_files'])
                and all(m.isfile() and not PurePosixPath(m.name).is_absolute()
                        and '..' not in PurePosixPath(m.name).parts for m in members),
                'unsafe or incomplete retained Rust source archive')
        for item in members:
            require(digest(tar.extractfile(item).read()) == source['root_files'][item.name],
                    'retained Rust source bytes differ from accepted manifest')
    cms = read('cms/receipt.json')
    require(digest(files['bin/cms']) == CMS_SHA256
            and digest(files['cms/receipt.json']) == CMS_RECEIPT_SHA256
            and cms['binaries']['cryptominisat']['sha256'] == CMS_SHA256
            and cms['status'] == 'completed' and cms['source_clean'] is True
            and cms['platform']['system'] == 'Darwin'
            and cms['platform']['machine'] == 'arm64',
            'SAT solver binary/build receipt differs from accepted static build')
    for name in ('source', 'cadical', 'cadiback'):
        receipt = cms['source_archives'][name]
        data = files['cms/'+name+'.tar']
        require(len(data) == receipt['bytes'] and digest(data) == receipt['sha256'],
                'retained CMS source archive differs from build receipt')
    admission = dict(schema_version=3, platform=dict(system='Darwin', machine='arm64'),
                     curve_id=curve_record(fixture)['curve']['curve_id'],
                     inventory=inventory, exporter_binary_sha256=EXPORTER_SHA256,
                     exporter_source_sha256=EXPORTER_SOURCE_SHA256,
                     rust_source_manifest_sha256=RUST_MANIFEST_SHA256,
                     rust_source_commit=RUST_COMMIT, cms_binary_sha256=CMS_SHA256,
                     cms_build_receipt_sha256=CMS_RECEIPT_SHA256,
                     scope='accepted local source/build receipts; not remote attestation')
    return fixture, report, curve, base, admission
