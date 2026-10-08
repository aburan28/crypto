"""Reusable, source/build-bound inputs for one n17 F5 development pipeline.

Only the retained physical macOS ARM64 worker is admitted here. No rebuild,
worker execution, target generation or new query feasibility oracle runs.
"""
import argparse
import copy
import hashlib
import io
import json
from pathlib import Path
import platform
import tarfile

from generic_bases import verify_base
from generic_build import verify_build_record
from identity import sha256
from oracle import Curve, require
from register_paired_generic import inputs as retained_worker_inputs
from replay_paired_n17_evidence import retained_files
from static_sat_assets_v3 import freeze_assets


SOURCE_MANIFEST = 'c64e4b3102bface63a2305efbff4bd85810cc112cb43546da2992ba48e9e85b7'
WORKER_SHA256 = '94caf3d67e57dde09488763ec19e791aca54bbb35f1968c73e42d37a210b436a'
DEPS_BUNDLE_SHA256 = '79f60c5b969a374c5595735498dddbb5677436085ce643a36a43d9886e3384ef'
HERE = Path(__file__).resolve().parent
ROLES = {'bin/worker', 'build/build-record.json', 'build/build-policy.json',
         'build/build-exit.json', 'build/build.log', 'rust/source-manifest.json',
         'rust/root-source.tar.gz', 'rust/dependency-source.tar.gz',
         'fixture.json', 'inventory.json', 'retention.json'}


def digest(data):
    return hashlib.sha256(data).hexdigest()


def source_archive(data, expected):
    observed = {}
    with tarfile.open(fileobj=io.BytesIO(data), mode='r:gz') as archive:
        for item in archive:
            require(item.isfile() and item.name not in observed,
                    'non-file or duplicate F5 source member')
            observed[item.name] = digest(archive.extractfile(item).read())
    require(observed == expected, 'F5 retained source bytes differ from build manifest')


def native_admission(files, *, check_host=True):
    require(set(files) == ROLES, 'F5 native input roles missing or extra')
    if check_host:
        require(platform.system() == 'Darwin' and platform.machine() == 'arm64',
                'F5 v1 admits the retained physical macOS ARM64 build only')
    record = json.loads(files['build/build-record.json'])
    source = json.loads(files['rust/source-manifest.json'])
    identity = verify_build_record(record, source)
    require(sha256(source) == SOURCE_MANIFEST
            and digest(files['bin/worker']) == WORKER_SHA256 == record['worker_sha256']
            and record['build']['target_arch'] == 'aarch64'
            and record['build']['target_os'] == 'macos'
            and json.loads(files['build/build-policy.json']) == record['build']
            and json.loads(files['build/build-exit.json'])['exit_code'] == 0,
            'F5 accepted worker or build receipt changed')
    source_archive(files['rust/root-source.tar.gz'], source['root_files'])
    dependencies = {item['package']+'-'+item['version']+'/'+name: expected
                    for item in source['dependencies'] for name, expected in item['files'].items()}
    source_archive(files['rust/dependency-source.tar.gz'], dependencies)
    fixture = json.loads(files['fixture.json'])
    inventory = json.loads(files['inventory.json'])
    require(fixture['degree'] == 17 and fixture['curve_a'] == 1
            and int(fixture['subgroup_order']) == 65587
            and int(fixture['group_order']) == 131174 and int(fixture['cofactor']) == 2
            and fixture['generator'] == ['43693', '23339']
            and fixture['targets'] == [] and fixture['target_seeds'] == []
            and fixture['target_scalar_constructed'] is False
            and inventory['fixture'] == fixture,
            'F5 reusable fixture contains a target or differs from the admitted curve')
    job = dict(degree=17, curve_a=1, factor_base=dict(kind='standard_subspace', dimension=6), config={})
    base = verify_base(inventory, fixture, job)
    require(len(inventory['factor_base']) == 63 and inventory['columns'] == 29
            and inventory['collector_dispatch'] == dict(strategy='Groebner', field_kernel='pmull',
                    pair_table=False, query_rule='trial-keyed-sample', collection_window=None),
            'F5 geometry or declared worker inventory differs')
    receipt = json.loads(files['retention.json'])
    require(receipt == dict(schema_version=1, source_commit='765c3c5f19032bd852163805f257c56babef2040',
            parent_archive_sha256='e0ce19fc28c58e2dc1cae9649a16af74099016fff8183e5ce69f65b39c804f02',
            dependency_retention_archive_sha256=DEPS_BUNDLE_SHA256,
            source_manifest_sha256=SOURCE_MANIFEST, worker_sha256=WORKER_SHA256,
            scope='accepted local source/build receipts; no remote attestation or new hardware claim'),
            'F5 source/build retention provenance differs')
    native = dict(build_identity=identity, source_manifest_sha256=SOURCE_MANIFEST,
                  worker_sha256=WORKER_SHA256, build_record_sha256=sha256(record),
                  source_commit=receipt['source_commit'], host_class='physical-macos-arm64')
    return fixture, inventory, Curve(fixture), base, record, source, native


def retained_assets():
    """Read and validate reusable retained bytes without native execution."""
    old_files, record, source, old, _ = retained_worker_inputs()
    bundle = HERE/'goal_20260924/f5-boolean-system-control-20260930/native-v2'
    require(json.loads((bundle/'receipt.json').read_text())['archive_sha256'] == DEPS_BUNDLE_SHA256,
            'accepted F5 dependency retention archive differs')
    deps = retained_files(bundle)['dependency-source.tar.gz']
    fixture = dict(old['fixture'], targets=[], target_seeds=[], target_scalar_constructed=False)
    # Keep only the independently checked static inventory, never old query
    # witnesses, scalar answers, rates, costs, target seeds or measured ranks.
    inventory = {key: copy.deepcopy(old[key]) for key in (
        'factor_base', 'columns', 'effective_config', 'effective_factor_base', 'collector_dispatch')}
    inventory['fixture'] = fixture
    files = {'bin/worker':old_files['build/worker'],
             'rust/source-manifest.json':old_files['build/source-manifest.json'],
             'rust/root-source.tar.gz':old_files['build/root-source.tar.gz'],
             'rust/dependency-source.tar.gz':deps,
             'fixture.json':json.dumps(fixture, sort_keys=True, separators=(',', ':')).encode(),
             'inventory.json':json.dumps(inventory, sort_keys=True, separators=(',', ':')).encode()}
    for name in ('build-record.json', 'build-policy.json', 'build-exit.json', 'build.log'):
        files['build/'+name] = old_files['build/'+name]
    files['retention.json'] = json.dumps(dict(schema_version=1,
        source_commit='765c3c5f19032bd852163805f257c56babef2040',
        parent_archive_sha256='e0ce19fc28c58e2dc1cae9649a16af74099016fff8183e5ce69f65b39c804f02',
        dependency_retention_archive_sha256=DEPS_BUNDLE_SHA256,
        source_manifest_sha256=SOURCE_MANIFEST, worker_sha256=WORKER_SHA256,
        scope='accepted local source/build receipts; no remote attestation or new hardware claim'),
        sort_keys=True, separators=(',', ':')).encode()
    native_admission(files, check_host=False)
    return files


def stage(output):
    files = retained_assets()
    native_admission(files)
    manifest, seal = freeze_assets(files, {'bin/worker'}, output)
    return dict(manifest_sha256=sha256(manifest), seal=seal)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(stage(args.out), sort_keys=True))
