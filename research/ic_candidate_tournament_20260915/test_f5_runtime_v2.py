"""Source/identity controls and independent old-report audits; no new F5 job."""
import copy
import json
from pathlib import Path
import tempfile
import unittest

from audit_f5_runtime_v2 import natural_queries
from f5_runtime_inputs_v1 import native_admission, retained_assets
from f5_runtime_registration_v2 import mathematical_registration, validate_panel
from identity import sha256
from oracle import InvalidEvidence
from register_paired_generic import inputs
from sat_runtime_execution_v3 import binding, register as register_runtime
from static_sat_assets_v3 import freeze_assets

ROOT = Path(__file__).resolve().parents[2]


class F5RuntimeV2Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.files = retained_assets()
        cls.fixture, cls.inventory, cls.curve, *_ = native_admission(cls.files, check_host=False)
        cls.old = inputs()[3]
        # This is only a registration control. These settings have no measured
        # invocation and do not authorize executing a historical protocol.
        cls.panel = dict(question='development-source-control', run_number=0,
            algorithm_seed=2026093032, config=dict(solver='f5', linear_algebra='dense', summands=3,
                groebner_degree=3, node_budget=8192, conflict_budget=100000,
                batch_trials=8, max_trials=512),
            resources=dict(host_class='physical-macos-arm64', cpu_workers=1, target_count=1,
                memory_limit_bytes=None, total_wall_limit_seconds=7200),
            target_input=dict(point=[52411,72106], seed=None,
                input_law='one-supplied-public-point; seed-is-provenance',
                point_was_previously_supplied=True, known_scalar_supplied=False))
        cls.spec = dict(entrypoint=dict(module='f5_runtime_pipeline_v2', callable='run'),
            runtime_watchdog_seconds=7200, runtime_seal=dict(manifest_sha256='a'*64,
                archive_sha256='b'*64), interpreter={'registration_control': True})
        cls.spec['binding'] = binding(cls.spec)

    def registration(self, panel=None, spec=None):
        return mathematical_registration(panel or self.panel, spec or self.spec,
                                         self.files, check_host=False)

    def test_target_and_seed_change_workload_without_changing_candidate(self):
        first = self.registration()
        self.assertEqual(first['job']['target_seeds'], [])
        self.assertEqual(first['fixture']['target_seeds'], [None])
        changed = copy.deepcopy(self.panel)
        changed['algorithm_seed'] += 1
        second = self.registration(changed)
        self.assertEqual(first['candidate'], second['candidate'])
        self.assertNotEqual(first['workload'], second['workload'])
        changed['target_input']['point'] = [853, 39791]
        third = self.registration(changed)
        self.assertEqual(first['candidate'], third['candidate'])
        self.assertNotEqual(second['workload'], third['workload'])
        self.assertEqual(first['candidate']['record']['factor_base']['inventory']['usable_point_count'], 62)
        self.assertEqual(first['candidate']['record']['factor_base']['inventory']['effective_columns'], 29)
        self.assertIn('PDP3f5RCsampleLAgaussTDpdpISO0h', first['candidate']['candidate_id'])

    def test_algorithm_and_controller_source_change_candidate(self):
        first = self.registration()
        changed = copy.deepcopy(self.panel)
        changed['config']['solver'] = 'f4'
        second = self.registration(changed)
        self.assertNotEqual(first['candidate'], second['candidate'])
        self.assertEqual(first['workload'], second['workload'])
        changed = copy.deepcopy(self.panel)
        changed['config']['node_budget'] += 1
        self.assertNotEqual(first['candidate'], self.registration(changed)['candidate'])
        spec = copy.deepcopy(self.spec)
        spec['runtime_seal']['manifest_sha256'] = 'c'*64
        spec['binding'] = binding(spec)
        self.assertNotEqual(first['candidate'], self.registration(spec=spec)['candidate'])
        self.assertEqual(first['workload'], self.registration(spec=spec)['workload'])

    def test_new_registration_seals_controller_packages_and_assets_without_worker(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            assets = root / 'assets'
            freeze_assets(self.files, {'bin/worker'}, assets)
            spec = register_runtime(ROOT, root / 'registration', module='f5_runtime_pipeline_v2',
                action='run', arguments=None, timeout_seconds=7200, asset_snapshot=assets,
                arguments_factory=lambda value: mathematical_registration(
                    self.panel, value, self.files, check_host=False))
            self.assertEqual(spec['arguments'], self.registration(spec=spec))
            roles = {item['role'] for item in spec['runtime_manifest']['components']}
            for name in ('f5_runtime_pipeline_v2.py', 'static_sat_native_v3.py',
                         'producer/evidence.py', 'producer/timing.py'):
                self.assertIn('research/ic_candidate_tournament_20260915/' + name, roles)
            self.assertFalse((root / 'registration/process.json').exists())
            self.assertEqual(spec['arguments']['seal']['registration_stage'], 'before-execution')
            self.assertIn('asset_archive_sha256', spec['binding'])
            self.assertEqual(spec['arguments']['source']['execution_binding'], spec['binding'])
            self.assertEqual(sha256(spec['arguments']['candidate']),
                             spec['arguments']['seal']['candidate_sha256'])

    def test_known_scalar_freshness_wrong_subgroup_and_boolean_resources_fail_closed(self):
        for key, value in [('known_scalar_supplied', True), ('point_was_previously_supplied', False),
                           ('point', [0, 1])]:
            changed = copy.deepcopy(self.panel)
            changed['target_input'][key] = value
            with self.assertRaises(InvalidEvidence):
                validate_panel(changed, self.curve)
        changed = copy.deepcopy(self.panel)
        changed['question'] = 'fresh-paired-qualification'
        with self.assertRaises(InvalidEvidence):
            validate_panel(changed, self.curve)
        for key in ('cpu_workers', 'target_count'):
            changed = copy.deepcopy(self.panel)
            changed['resources'][key] = True
            with self.assertRaises(InvalidEvidence):
                validate_panel(changed, self.curve)
        for key in ('node_budget', 'conflict_budget', 'batch_trials', 'max_trials'):
            changed = copy.deepcopy(self.panel)
            changed['config'][key] = 2**64
            with self.assertRaisesRegex(InvalidEvidence, 'native 64-bit'):
                validate_panel(changed, self.curve)

    def test_changed_native_worker_source_or_target_leak_fails_admission(self):
        changed = dict(self.files)
        changed['bin/worker'] += b'changed'
        with self.assertRaisesRegex(InvalidEvidence, 'worker or build'):
            native_admission(changed, check_host=False)
        changed = dict(self.files)
        source = json.loads(changed['rust/source-manifest.json'])
        source['root_files'][next(iter(source['root_files']))] = '0'*64
        changed['rust/source-manifest.json'] = json.dumps(source).encode()
        with self.assertRaises(InvalidEvidence):
            native_admission(changed, check_host=False)
        changed = dict(self.files)
        fixture = json.loads(changed['fixture.json'])
        fixture['targets'] = [['52411', '72106']]
        changed['fixture.json'] = json.dumps(fixture).encode()
        with self.assertRaisesRegex(InvalidEvidence, 'reusable fixture'):
            native_admission(changed, check_host=False)

    def test_retained_natural_queries_keep_failed_feasible_cases_and_reject_bad_witness(self):
        result = natural_queries(self.old, self.curve)
        self.assertEqual(result['attempts'], 104)
        self.assertEqual(result['status_mix'], {'incomplete': 75, 'witness': 29})
        self.assertEqual(result['exact_group_feasible'], 33)
        self.assertEqual(result['feasible_but_no_witness'], 4)
        self.assertEqual(result['verified_native_witnesses'], 29)
        changed = copy.deepcopy(self.old)
        for batch in changed['collection_reports']:
            witness = next((a for a in batch['attempts'] if a['pdp']['outcome'] == 'witness'), None)
            if witness is not None:
                witness['pdp']['points'][0] = -1
                break
        with self.assertRaisesRegex(InvalidEvidence, 'index leaves base'):
            natural_queries(changed, self.curve)
        changed = copy.deepcopy(self.old)
        changed['collection_reports'][0]['attempts'][0]['b'] = 1
        with self.assertRaisesRegex(InvalidEvidence, 'depends on supplied target'):
            natural_queries(changed, self.curve)


if __name__ == '__main__':
    unittest.main()
