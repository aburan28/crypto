"""Synthetic build receipts and a retained report; no native worker execution."""
import copy
import gzip
import hashlib
import io
import json
from pathlib import Path
import tarfile
import tempfile
import unittest
from unittest.mock import patch

from generic_stages import DEFAULTS
from identity import sha256
from oracle import InvalidEvidence
from prepared_f5_runtime_v1 import (CONSUMED_SOURCE_MANIFEST, CONSUMED_WORKER_SHA256,
                                   WORKER_ROLE, admit_native, audit,
                                   mathematical_registration, package_assets, run)
from static_sat_assets_v3 import verified_assets
from prepared_ic_state_v1 import accepted_files
from prepared_target_v1 import CERTIFICATE_SEALS, native_job
from sat_runtime_execution_v3 import binding

HERE = Path(__file__).resolve().parent


def gzip_tar(members):
    raw = io.BytesIO()
    with gzip.GzipFile(fileobj=raw, mode='wb', mtime=0) as compressed:
        with tarfile.open(fileobj=compressed, mode='w') as archive:
            for name, data in sorted(members.items()):
                item = tarfile.TarInfo(name)
                item.size, item.mtime = len(data), 0
                archive.addfile(item, io.BytesIO(data))
    return raw.getvalue()


def digest(data):
    return hashlib.sha256(data).hexdigest()


class PreparedF5RuntimeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.doc = json.loads((HERE/'goal_20260924/prepared-ic-state-v1/f5-preparation.json').read_text())
        cls.accepted = accepted_files('f5')
        report = cls.accepted['execution/entry-output/pipeline.stdout']
        fixture = copy.deepcopy(report['fixture'])
        fixture.update(targets=[], target_seeds=[], target_scalar_constructed=False)
        inventory = {key: copy.deepcopy(report[key]) for key in (
            'factor_base', 'columns', 'effective_factor_base', 'collector_dispatch')}
        inventory['fixture'] = fixture
        cls.fixture = fixture
        cls.inventory = inventory
        cls.worker_source = b'fn run_prepared_target() {}\n'
        cls.worker = b'prepared-f5-synthetic-worker\n'
        cls.files = cls.asset_files(cls.worker_source, cls.worker, 'linux', 'x86_64')
        cls.native = admit_native(cls.files)[-1]
        cls.spec = dict(entrypoint=dict(module='prepared_f5_runtime_v1', callable='run'),
            runtime_watchdog_seconds=600,
            runtime_seal=dict(manifest_sha256='a'*64, archive_sha256='b'*64),
            asset_manifest=dict(schema_version=3, components=[]),
            asset_seal=dict(schema_version=3, manifest_sha256='c'*64, archive_sha256='d'*64,
                            registration_stage='before-execution'),
            interpreter={'scope': 'registration control only'})
        cls.spec['binding'] = binding(cls.spec)
        cls.panel = dict(question='prepared-development-source-control', run_number=0,
            algorithm_seed=2026093032, max_descent_queries=8,
            target_input=dict(point=[52411, 72106], seed=None,
                input_law='one-disclosed-public-point; fixture-construction-excluded',
                point_was_previously_supplied=True, known_scalar_supplied=False),
            resources=dict(host_class='linux-x86_64', cpu_workers=1, target_count=1,
                memory_limit_bytes=None, total_wall_limit_seconds=600))

    @classmethod
    def asset_files(cls, worker_source, worker, target_os, target_arch):
        cargo, lock = b'[package]\nname="crypto"\n', b'# lock\n'
        root = {WORKER_ROLE: worker_source, 'Cargo.toml': cargo, 'Cargo.lock': lock}
        source = dict(schema_version=1, root_files={name: digest(data) for name, data in root.items()},
                      dependencies=[])
        policy = dict(schema_version=1, source_manifest_sha256=sha256(source),
            rustc='synthetic', cargo='synthetic', native_tools={},
            arguments=['cargo', 'build', '--example', 'ic_tournament_worker'],
            target_os=target_os, target_arch=target_arch,
            flags=dict(rustflags='', incremental=False, features=[], profile='release'),
            scope='synthetic registration control; not a compiled worker')
        identity = dict(schema_version=1, source_manifest_sha256=sha256(source),
                        build_sha256=sha256(policy), target_arch=target_arch, target_os=target_os)
        record = dict(schema_version=1, source_manifest_sha256=sha256(source),
                      build_sha256=sha256(policy), build=policy, identity=identity,
                      worker_sha256=digest(worker), builder_sha256='11'*32)
        files = {
            'bin/worker': worker,
            'build/build-record.json': json.dumps(record).encode(),
            'build/build-policy.json': json.dumps(policy).encode(),
            'build/build-exit.json': b'{"exit_code":0}',
            'build/build.log': b'synthetic build log\n',
            'rust/source-manifest.json': json.dumps(source).encode(),
            'rust/root-source.tar.gz': gzip_tar(root),
            'rust/dependency-source.tar.gz': gzip_tar({}),
            'fixture.json': json.dumps(cls.fixture).encode(),
            'inventory.json': json.dumps(cls.inventory).encode(),
        }
        return files

    def registration(self, panel=None, spec=None, files=None, doc=None):
        return mathematical_registration(panel or self.panel, spec or self.spec,
            files or self.files, doc or self.doc, CERTIFICATE_SEALS['f5'])

    def prepared_report(self):
        report = copy.deepcopy(self.accepted['execution/entry-output/pipeline.stdout'])
        job = native_job(self.doc, CERTIFICATE_SEALS['f5'], point=[52411, 72106],
                         algorithm_seed=2026093032)
        report.update(preparation_mode='imported-certified-log-table-v1',
            preparation_mathematical_state_sha256=job['prepared']['mathematical_state_sha256'],
            reusable_symbolic_template_prepared=True,
            trials=0, relations=[], collection_reports=[], solve_attempts=0,
            effective_config=dict(copy.deepcopy(DEFAULTS), **job['config']),
            scalar_verified=True, generic_build=self.native['build_identity'],
            factor_base=job['prepared']['factor_base'],
            column_logs=job['prepared']['columns'])
        return report

    def test_code_binding_excludes_asset_seals_seeds_and_certificate_history(self):
        first = self.registration()
        self.assertTrue(first['candidate']['candidate_id'].startswith('IC1N17Ckb1fb62PDP3f5'))
        self.assertEqual(first['candidate']['record']['factor_base']['inventory']['usable_point_count'], 62)
        record = json.dumps(first['candidate']['record'])
        self.assertNotIn(CERTIFICATE_SEALS['f5'], record)
        self.assertNotIn(self.doc['provenance']['archive_sha256'], record)
        self.assertNotIn(self.spec['binding']['asset_archive_sha256'], record)
        self.assertNotIn(self.spec['binding']['asset_manifest_sha256'], record)
        self.assertIn(self.native['worker_sha256'], record)
        self.assertIn(self.native['worker_source_sha256'], record)
        for key in ('algorithm_seed', 'run_number'):
            panel = copy.deepcopy(self.panel)
            panel[key] += 1
            changed = self.registration(panel)
            self.assertEqual(first['candidate'], changed['candidate'])
            self.assertNotEqual(first['seal']['run_id'], changed['seal']['run_id'])
        spec = copy.deepcopy(self.spec)
        spec['asset_seal'] = dict(spec['asset_seal'], archive_sha256='e'*64)
        spec['binding'] = binding(spec)
        self.assertEqual(first['candidate'], self.registration(spec=spec)['candidate'])
        spec = copy.deepcopy(self.spec)
        spec['runtime_seal'] = dict(spec['runtime_seal'], manifest_sha256='f'*64)
        spec['binding'] = binding(spec)
        self.assertNotEqual(first['candidate'], self.registration(spec=spec)['candidate'])

    def test_consumed_worker_fresh_targets_and_missing_marker_are_rejected(self):
        with patch('prepared_f5_runtime_v1.CONSUMED_WORKER_SHA256', digest(self.worker)), \
             self.assertRaises(InvalidEvidence):
            admit_native(self.files)
        with patch('prepared_f5_runtime_v1.CONSUMED_SOURCE_MANIFEST',
                   self.native['source_manifest_sha256']), \
             self.assertRaises(InvalidEvidence):
            admit_native(self.files)
        self.assertNotEqual(self.native['worker_sha256'], CONSUMED_WORKER_SHA256)
        self.assertNotEqual(self.native['source_manifest_sha256'], CONSUMED_SOURCE_MANIFEST)
        for change in ('fresh', 'point', 'scalar', 'cap', 'host', 'watchdog', 'marker'):
            panel, spec, files = copy.deepcopy(self.panel), copy.deepcopy(self.spec), self.files
            if change == 'fresh':
                panel['question'] = 'fresh-paired-qualification'
            elif change == 'point':
                panel['target_input']['point'] = [471, 57570]
            elif change == 'scalar':
                panel['target_input']['known_scalar_supplied'] = True
            elif change == 'cap':
                panel['max_descent_queries'] = 9
            elif change == 'host':
                panel['resources']['host_class'] = 'macos-aarch64'
            elif change == 'watchdog':
                spec['runtime_watchdog_seconds'] += 1
            else:
                files = self.asset_files(b'fn ordinary_collection() {}\n', self.worker, 'linux', 'x86_64')
            with self.assertRaises(InvalidEvidence):
                self.registration(panel, spec, files)
        swapped = json.loads((HERE/'goal_20260924/prepared-ic-state-v1/sat-preparation.json').read_text())
        with self.assertRaises(InvalidEvidence):
            mathematical_registration(self.panel, self.spec, self.files, swapped, CERTIFICATE_SEALS['sat'])

    def test_bound_entry_audits_the_retained_relation_and_rejects_tampering(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            out = root/'entry-output'
            out.mkdir()
            args = self.registration()
            spec = dict(copy.deepcopy(self.spec), arguments=args)
            (root/'execution.json').write_text(json.dumps(spec))
            report = self.prepared_report()
            preflight = dict(returncode=0, timed_out=False, scope='mocked build identity')
            pipeline = dict(returncode=0, timed_out=False, scope='mocked prepared worker')

            def fake_meter(_execution, _role, _arguments, directory, name, _seconds, stdin_argument=None):
                if name == 'build_identity':
                    self.assertIsNone(stdin_argument)
                    (directory/'build_identity.stdout').write_text(
                        json.dumps(self.native['build_identity']))
                    return preflight
                self.assertEqual(stdin_argument, 'job')
                (directory/'pipeline.stdout').write_text(json.dumps(report))
                return pipeline

            with patch('prepared_f5_runtime_v1.check_extracted_assets', return_value=self.files), \
                 patch('prepared_f5_runtime_v1.platform.system', return_value='Linux'), \
                 patch('prepared_f5_runtime_v1.platform.machine', return_value='x86_64'), \
                 patch('prepared_f5_runtime_v1.meter', side_effect=fake_meter):
                result = run(args, out)
                self.assertEqual(result['status'], 'COMPLETE')
                self.assertEqual(result['run_id'], args['seal']['run_id'])
                summary = json.loads((out/'summary.json').read_text())
                self.assertEqual(summary['recovered_scalar'], 24886)
                self.assertEqual(summary['ordinary_queries_executed'], 0)
                self.assertFalse(summary['source_bound_execution_admitted'])
                self.assertIsNone(summary['online_speedup'])
                with self.assertRaises(InvalidEvidence):
                    run(args, out)
            with patch('prepared_f5_runtime_v1.audit_execution',
                       return_value={'entrypoint_succeeded': True, 'scope': 'mocked-source-control'}), \
                 patch('prepared_f5_runtime_v1.check_extracted_assets', return_value=self.files), \
                 patch('prepared_f5_runtime_v1.audit_meter', side_effect=[preflight, pipeline]):
                admitted = audit(root, spec)
                self.assertEqual(admitted['status'], 'ADMITTED_COMPLETE_PREPARED_F5_CONTROL')
                self.assertEqual(admitted['recovered_scalar'], 24886)
                self.assertTrue(admitted['source_bound_execution_admitted'])
                self.assertFalse(admitted['headline_online_admissible'] or admitted['promotion_eligible'])
                self.assertIsNone(admitted['online_speedup'])
            summary['online_wall_ns'] += 1
            (out/'summary.json').write_text(json.dumps(summary))
            with patch('prepared_f5_runtime_v1.audit_execution',
                       return_value={'entrypoint_succeeded': True}), \
                 patch('prepared_f5_runtime_v1.check_extracted_assets', return_value=self.files), \
                 patch('prepared_f5_runtime_v1.audit_meter', side_effect=[preflight, pipeline]):
                with self.assertRaises(InvalidEvidence):
                    audit(root, spec)

    def test_packaged_assets_roundtrip_to_the_same_registration(self):
        with tempfile.TemporaryDirectory() as temporary:
            snapshot = Path(temporary)/'assets'
            package_assets(self.files, snapshot)
            files = verified_assets(snapshot, json.loads((snapshot/'manifest.json').read_text()),
                                    json.loads((snapshot/'seal.json').read_text()))
            self.assertEqual(self.registration(), self.registration(files=files))

    def test_wrong_platform_and_ordinary_collection_do_not_pass_the_entrypoint(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            out = root/'entry-output'
            out.mkdir()
            macos = self.asset_files(self.worker_source, self.worker, 'macos', 'aarch64')
            panel = copy.deepcopy(self.panel)
            panel['resources']['host_class'] = 'macos-aarch64'
            args = self.registration(panel, files=macos)
            spec = dict(copy.deepcopy(self.spec), arguments=args)
            (root/'execution.json').write_text(json.dumps(spec))
            with patch('prepared_f5_runtime_v1.check_extracted_assets', return_value=macos), \
                 patch('prepared_f5_runtime_v1.meter', side_effect=AssertionError('worker ran')):
                with self.assertRaises(InvalidEvidence):
                    run(args, out)
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            out = root/'entry-output'
            out.mkdir()
            args = self.registration()
            spec = dict(copy.deepcopy(self.spec), arguments=args)
            (root/'execution.json').write_text(json.dumps(spec))
            report = self.prepared_report()
            report['trials'] = 1

            def fake_meter(_execution, _role, _arguments, directory, name, _seconds, stdin_argument=None):
                if name == 'build_identity':
                    (directory/'build_identity.stdout').write_text(
                        json.dumps(self.native['build_identity']))
                    return dict(returncode=0, timed_out=False)
                (directory/'pipeline.stdout').write_text(json.dumps(report))
                return dict(returncode=0, timed_out=False)

            with patch('prepared_f5_runtime_v1.check_extracted_assets', return_value=self.files), \
                 patch('prepared_f5_runtime_v1.platform.system', return_value='Linux'), \
                 patch('prepared_f5_runtime_v1.platform.machine', return_value='x86_64'), \
                 patch('prepared_f5_runtime_v1.meter', side_effect=fake_meter):
                with self.assertRaises(InvalidEvidence):
                    run(args, out)


if __name__ == '__main__':
    unittest.main()
