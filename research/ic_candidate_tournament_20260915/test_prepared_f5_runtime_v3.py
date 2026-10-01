"""Real rebuilt assets and retained mathematics; mocked native/source calls.

These tests are registration/rejection controls, never new F5 execution,
hardware validation, natural-yield observations or competitive measurements.
"""
import copy
import io
import json
from pathlib import Path
import tempfile
import tarfile
import unittest
from unittest.mock import patch

import generic_build
from identity import sha256
from oracle import InvalidEvidence
from prepared_f5_inputs_v2 import (MATHEMATICS, MATHEMATICS_SHA256, WORKER,
    WORKER_SOURCE_SHA256, archive_members, digest, native_admission, source_archive)
from prepared_f5_runtime_v3 import (audit, mathematical_registration, register,
    run, transport)
from prepared_target_v1 import CERTIFICATE_SEALS
from prepared_runtime_transport_v2 import HELPER_ROLE, validate_admission
from sat_runtime_execution_v3 import audit_execution, binding, execute, register as register_runtime
from static_sat_assets_v3 import verified_assets
import test_prepared_target_v1 as controls

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


class PreparedF5RuntimeV3Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.docs = {family:json.loads((HERE/'goal_20260924/prepared-ic-state-v1'/
            (family+'-preparation.json')).read_text()) for family in ('f5','sat')}
        cls.assets = HERE/'goal_20260924/prepared-f5-runtime-v3/native-inputs-macos-arm64'
        manifest,seal = [json.loads((cls.assets/name).read_text()) for name in ('manifest.json','seal.json')]
        cls.files = verified_assets(cls.assets,manifest,seal)
        cls.build,cls.source,cls.native = native_admission(cls.files)
        cls.spec = dict(entrypoint=dict(module='prepared_f5_runtime_v3',callable='run'),
            runtime_watchdog_seconds=600,runtime_seal=dict(manifest_sha256='a'*64,archive_sha256='b'*64),
            asset_manifest=manifest,asset_seal=seal,interpreter={'scope':'registration control only'})
        cls.spec['binding'] = binding(cls.spec)
        cls.panel = dict(question='prepared-development-source-control',algorithm_seed=2026093032,
            max_attempts=8,run_number=0,target_input=dict(point=[52411,72106],seed=None,
                input_law='one-disclosed-public-point; fixture-construction-excluded',
                point_was_previously_supplied=True,known_scalar_supplied=False),
            resources=dict(host_class='physical-macos-arm64',cpu_workers=1,target_count=1,
                memory_limit_bytes=None,total_wall_limit_seconds=600))
        controls.PreparedTargetTests.setUpClass()

    def registration(self, panel=None, spec=None, document=None):
        return mathematical_registration(panel or self.panel,spec or self.spec,self.files,
            document or self.docs['f5'],CERTIFICATE_SEALS['f5'])

    def test_mathematical_fixture_is_shared_and_contains_no_preparation_history(self):
        fixture = json.loads((ROOT/MATHEMATICS).read_text())
        self.assertEqual(digest((ROOT/MATHEMATICS).read_bytes()),MATHEMATICS_SHA256)
        self.assertEqual(set(fixture),{'mathematical_state_sha256','factor_base','columns'})
        for document in self.docs.values():
            self.assertEqual(fixture,dict(mathematical_state_sha256=document['record_sha256'],
                factor_base=[[str(v) for v in point] for point in document['record']['factor_base']['points']],
                columns=[dict(point=[str(v) for v in item['point']],log=str(item['log']))
                         for item in document['record']['column_logs']]))
        # This pin identifies the retained v1 worker, not today's mutable tree.
        with tarfile.open(fileobj=io.BytesIO(self.files['rust/root-source.tar.gz']), mode='r:gz') as archive:
            retained_worker = archive.extractfile(WORKER).read()
        self.assertEqual(digest(retained_worker),WORKER_SOURCE_SHA256)

    def test_self_consistent_changed_source_still_cannot_enter_the_old_admission_gate(self):
        files = copy.deepcopy(self.files)
        source = copy.deepcopy(self.source)
        with tarfile.open(fileobj=io.BytesIO(files['rust/root-source.tar.gz']), mode='r:gz') as archive:
            root_files = {item.name:archive.extractfile(item).read() for item in archive}
        root_files[WORKER] += b'\n// prospective source change; not a built worker\n'
        source['root_files'][WORKER] = digest(root_files[WORKER])
        policy = copy.deepcopy(self.build['build'])
        policy['source_manifest_sha256'] = sha256(source)
        record = copy.deepcopy(self.build)
        record.update(source_manifest_sha256=sha256(source),build=policy,build_sha256=sha256(policy))
        record['identity'].update(source_manifest_sha256=sha256(source),build_sha256=sha256(policy))
        files['rust/root-source.tar.gz'] = source_archive(root_files)
        for role, value in [('rust/source-manifest.json',source),
                            ('build/build-policy.json',policy),('build/build-record.json',record)]:
            files[role] = json.dumps(value).encode()
        generic_build.verify_build_record(record,source)
        with self.assertRaisesRegex(InvalidEvidence, '^prepared worker source lacks the admitted mode'):
            native_admission(files)

    def test_conservative_rust_manifest_retains_only_mathematics_not_certificate(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root/'src').mkdir()
            for name,data in [('Cargo.toml',b''),('Cargo.lock',b''),
                              (WORKER,(ROOT/WORKER).read_bytes()),(MATHEMATICS,(ROOT/MATHEMATICS).read_bytes())]:
                path = root/name; path.parent.mkdir(parents=True,exist_ok=True); path.write_bytes(data)
            source = generic_build.source_manifest(root,{'packages':[]})
        self.assertIn(MATHEMATICS,source['root_files'])
        self.assertFalse(any(name.endswith('-preparation.json') for name in source['root_files']))
        self.assertFalse(any(name.endswith('-preparation.json') for name in self.source['root_files']))

    def test_corrected_assets_cannot_enter_the_immutable_old_input_gate(self):
        from prepared_f5_inputs_v1 import native_admission as old_admission
        with self.assertRaisesRegex(InvalidEvidence, '^prepared worker source lacks the admitted mode'):
            old_admission(self.files)
        old_assets = HERE/'goal_20260924/prepared-f5-runtime-v1/native-inputs-macos-arm64'
        old_manifest,old_seal = [json.loads((old_assets/name).read_text())
                                for name in ('manifest.json','seal.json')]
        old_files = verified_assets(old_assets,old_manifest,old_seal)
        old_admission(old_files)
        with self.assertRaisesRegex(InvalidEvidence, '^prepared worker source lacks the admitted mode'):
            native_admission(old_files)

    def test_binary_recipe_and_complete_dependency_retention_are_required(self):
        for change in ('binary','recipe','dependency','old_source'):
            files = copy.deepcopy(self.files)
            if change == 'binary': files['bin/worker'] += b'changed'
            elif change == 'recipe': files['build/build-policy.json'] = b'{}'
            elif change == 'dependency': files['rust/dependency-source.tar.gz'] = source_archive({})
            else:
                source = copy.deepcopy(self.source)
                source['root_files'][WORKER] = '0'*64
                files['rust/source-manifest.json'] = json.dumps(source).encode()
            with self.assertRaises(InvalidEvidence): native_admission(files)
        with self.assertRaises(InvalidEvidence):
            archive_members(source_archive({'../escape':b'code'}),{'../escape':digest(b'code')})

    def test_candidate_binds_code_math_not_seeds_resources_or_transport_receipts(self):
        first = self.registration()
        for key in ('algorithm_seed','run_number','wall'):
            panel = copy.deepcopy(self.panel)
            if key == 'wall': panel['resources']['total_wall_limit_seconds'] += 1
            else: panel[key] += 1
            spec = copy.deepcopy(self.spec)
            spec['runtime_watchdog_seconds'] = panel['resources']['total_wall_limit_seconds']
            changed = self.registration(panel,spec)
            self.assertEqual(first['candidate'],changed['candidate'])
            self.assertNotEqual(first['seal']['run_id'],changed['seal']['run_id'])
        for value in (CERTIFICATE_SEALS['f5'],sha256(self.build),self.build['build_sha256'],
                      self.spec['asset_seal']['archive_sha256'],self.spec['asset_seal']['manifest_sha256']):
            self.assertNotIn(value,json.dumps(first['candidate']['record']))
        spec = copy.deepcopy(self.spec)
        spec['runtime_seal']['manifest_sha256'] = 'c'*64; spec['binding'] = binding(spec)
        self.assertNotEqual(first['candidate'],self.registration(spec=spec)['candidate'])
        self.assertTrue(first['candidate']['candidate_id'].startswith('IC1N17Ckb1fb62PDP3f5'))

    def test_control_caps_fresh_point_bool_and_changed_certificate_fail_closed(self):
        for change in ('fresh','point','bool','scalar','extra','cap','wall','watchdog','certificate'):
            panel,spec,doc = copy.deepcopy(self.panel),copy.deepcopy(self.spec),copy.deepcopy(self.docs['f5'])
            if change == 'fresh': panel['question'] = 'fresh-paired-qualification'
            elif change == 'point': panel['target_input']['point'] = [471,57570]
            elif change == 'bool': panel['resources']['cpu_workers'] = True
            elif change == 'scalar': panel['target_input']['known_scalar_supplied'] = 0
            elif change == 'extra': panel['ordinary_queries'] = 1
            elif change == 'cap': panel['max_attempts'] = 9
            elif change == 'wall': panel['resources']['total_wall_limit_seconds'] = 30
            elif change == 'watchdog': spec['runtime_watchdog_seconds'] += 1
            else: doc['provenance']['archive_sha256'] = '0'*64
            with self.assertRaises(InvalidEvidence): self.registration(panel,spec,doc)

    def control_execution(self, root, *, timeout=False):
        """Construct an explicitly mocked runtime around real retained witnesses."""
        report,_ = controls.PreparedTargetTests().native_math_control()
        arguments = self.registration()
        spec = dict(copy.deepcopy(self.spec),arguments=arguments)
        report['fixture'] = arguments['fixture']
        report['generic_build'] = self.native['build_identity']
        report['effective_config'] = dict(copy.deepcopy(report['effective_config']),max_trials=8)
        root.mkdir()
        out = root/'entry-output'; out.mkdir()
        binary = root/'asset-files/bin/worker'; binary.parent.mkdir(parents=True)
        binary.write_bytes(self.files['bin/worker'])
        (root/'execution.json').write_text(json.dumps(spec))
        process = dict(returncode=-9 if timeout else 0,timed_out=timeout,
            native_wall_ns=report['online_wall_ns']+1000000,
            metrics=dict(wall_seconds=1.0,user_seconds=0.9,system_seconds=0.1,
                total_core_seconds=1.0,single_core_seconds=1.0,peak_rss_bytes=1234,meter='mock unit control'))

        def native_meter(execution,role,argv,directory,name,seconds,**kwargs):
            self.assertEqual(role,'bin/worker')
            if name == 'build_identity':
                self.assertEqual(argv,['--build-identity'])
                (directory/'build_identity.stdout').write_text(json.dumps(self.native['build_identity']))
                return dict(returncode=0,timed_out=False)
            self.assertEqual((name,argv,seconds,kwargs),('pipeline',[],570,{'stdin_argument':'job'}))
            if not timeout: (directory/'pipeline.stdout').write_text(json.dumps(report))
            (directory/'pipeline.metrics.json').write_text(json.dumps(process))
            return process

        with patch('prepared_f5_runtime_v3.check_extracted_assets',return_value=self.files),\
             patch('prepared_f5_inputs_v2.platform.system',return_value='Darwin'),\
             patch('prepared_f5_inputs_v2.platform.machine',return_value='arm64'),\
             patch('prepared_f5_runtime_v3.meter',side_effect=native_meter) as calls:
            result = run(arguments,out)
            self.assertEqual(calls.call_count,2)  # One preflight, exactly one native attempt; no retry.
        return arguments,spec,report,process,result

    def audit_control(self,root,spec,process):
        with patch('prepared_f5_runtime_v3.audit_execution',return_value={'entrypoint_succeeded':True,
                        'scope':'mock source gates; no new execution admitted by this test'}),\
             patch('prepared_f5_runtime_v3.check_extracted_assets',return_value=self.files),\
             patch('prepared_f5_runtime_v3.audit_meter',side_effect=lambda *a,**k:
                   dict(returncode=0,timed_out=False) if a[2]=='build_identity' else process):
            return audit(root,spec)

    def test_entrypoint_and_auditor_replay_retained_negatives_witness_and_clock(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)/'execution'
            arguments,spec,report,process,result = self.control_execution(root)
            self.assertEqual(result['status'],'COMPLETE')
            audited = self.audit_control(root,spec,process)
            validate_admission(audited, dict(spec, runtime_manifest={
                'components': [{'role': HELPER_ROLE}]}))
            self.assertEqual(audited['recovered_scalar'],24886)
            self.assertEqual(audited['target_attempt_count'],3)
            self.assertEqual(audited['ordinary_queries_executed'],0)
            self.assertFalse(audited['headline_online_admissible'] or audited['promotion_eligible'])
            for change in ('engine','clock','ordinary','summary'):
                changed = copy.deepcopy(report)
                if change == 'engine': changed['solutions'][0]['attempts'][0]['pdp']['stats']['engine']={'MatrixF4':{'max_degree':3}}
                elif change == 'clock': changed['online_wall_ns'] += 1
                elif change == 'ordinary': changed['trials'] = 1
                else: (root/'entry-output/summary.json').write_text('{}')
                (root/'entry-output/pipeline.stdout').write_text(json.dumps(changed))
                with self.assertRaises(InvalidEvidence): self.audit_control(root,spec,process)

    def test_native_timeout_keeps_unknown_target_cost_without_retry(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)/'execution'
            _,spec,_,process,result = self.control_execution(root,timeout=True)
            self.assertEqual(result['status'],'NATIVE_TIMEOUT')
            audited = self.audit_control(root,spec,process)
            validate_admission(audited, dict(spec, runtime_manifest={
                'components': [{'role': HELPER_ROLE}]}))
            self.assertFalse(audited['scalar_verified'])
            self.assertIsNone(audited['online_wall_ns'])
            self.assertIsNone(audited['target_attempt_count'])

    def test_real_source_registration_is_one_use_and_unexecuted_transport_rejects(self):
        with tempfile.TemporaryDirectory() as temporary,\
             patch('prepared_f5_inputs_v2.platform.system',return_value='Darwin'),\
             patch('prepared_f5_inputs_v2.platform.machine',return_value='arm64'):
            root = Path(temporary)
            spec = register(ROOT,self.assets,self.panel,self.docs['f5'],CERTIFICATE_SEALS['f5'],root/'registered')
            self.assertEqual(spec['arguments']['seal']['registration_stage'],'before-execution')
            with self.assertRaises(InvalidEvidence):
                register(ROOT,self.assets,self.panel,self.docs['f5'],CERTIFICATE_SEALS['f5'],root/'registered')
            with patch('prepared_runtime_transport_v2.subprocess.run',side_effect=AssertionError('no audit child may run')):
                with self.assertRaises((InvalidEvidence,FileNotFoundError)):
                    transport(root/'registered',sha256(spec),root/'transport')
            self.assertFalse((root/'transport').exists())

    def test_new_adapters_import_from_the_real_fresh_isolated_source_snapshot(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            modules = ['prepared_f5_runtime_v3','prepared_runtime_transport_v2']
            spec = register_runtime(ROOT,root/'registration',module='sat_runtime_execution_v3',
                action='import_probe',arguments=modules,timeout_seconds=60)
            process = execute(root/'registration',root/'execution',expected_spec=spec,timeout_seconds=60)
            self.assertEqual(process['exit_code'],0,(root/'execution/stderr.txt').read_text())
            checked = audit_execution(root/'execution',spec)
            self.assertTrue(checked['complete_source_gates'] and checked['entrypoint_succeeded'])
            terminal = json.loads((root/'execution/after.json').read_text())
            for module in modules:
                self.assertIn(module,terminal['loaded_modules'])
            self.assertEqual(terminal['result'],{'imported':modules,'measured_solver_executed':False})
            self.assertFalse(checked['promotion_eligible'])


if __name__ == '__main__':
    unittest.main()
