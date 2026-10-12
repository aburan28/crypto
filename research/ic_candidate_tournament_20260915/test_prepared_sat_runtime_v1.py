"""Real source assets and preparation, mocked native calls; no SAT dispatch."""
import copy
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from identity import sha256
from oracle import InvalidEvidence
from prepared_ic_state_v1 import accepted_files
from prepared_sat_runtime_v1 import audit, mathematical_registration, run
from prepared_runtime_transport_v1 import HELPER_ROLE, validate_admission
from prepared_target_v1 import CERTIFICATE_SEALS
from sat_runtime_execution_v3 import binding
from static_sat_assets_v3 import verified_assets

HERE = Path(__file__).resolve().parent


class PreparedSatRuntimeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.doc = json.loads((HERE/'goal_20260924/prepared-ic-state-v1/sat-preparation.json').read_text())
        assets = HERE/'goal_20260924/static-sat-runtime-v3/native-inputs-macos-arm64'
        manifest, seal = [json.loads((assets/name).read_text()) for name in ('manifest.json','seal.json')]
        cls.files = verified_assets(assets, manifest, seal)
        cls.spec = dict(entrypoint=dict(module='prepared_sat_runtime_v1',callable='run'),
            runtime_watchdog_seconds=600, runtime_seal=dict(manifest_sha256='a'*64,archive_sha256='b'*64),
            asset_manifest=manifest, asset_seal=seal, interpreter={'scope':'registration control only'})
        cls.spec['binding'] = binding(cls.spec)
        cls.panel = dict(question='prepared-development-source-control', run_number=0,
            descent_query_seed=2026093012, export_nonce=2026093013, max_descent_queries=8,
            cms_conflict_budget=100000, cms_timeout_seconds=60, export_timeout_seconds=30,
            target_input=dict(point=[52411,72106],seed=None,
                input_law='one-disclosed-public-point; fixture-construction-excluded',
                point_was_previously_supplied=True,known_scalar_supplied=False),
            resources=dict(host_class='physical-macos-arm64',cpu_workers=1,target_count=1,
                memory_limit_bytes=None,total_wall_limit_seconds=600))
        cls.source_row = accepted_files('sat')['execution/entry-output/summary.json']['target_attempts'][0]

    def registration(self, panel=None, spec=None, doc=None):
        return mathematical_registration(panel or self.panel, spec or self.spec, self.files,
                                          doc or self.doc, CERTIFICATE_SEALS['sat'])

    def test_source_and_math_enter_candidate_while_seeds_and_certificate_history_stay_in_run(self):
        first = self.registration()
        for key in ('descent_query_seed','export_nonce','run_number'):
            panel = copy.deepcopy(self.panel)
            panel[key] += 1
            changed = self.registration(panel)
            self.assertEqual(first['candidate'], changed['candidate'])
            self.assertNotEqual(first['seal']['run_id'], changed['seal']['run_id'])
        for seal in CERTIFICATE_SEALS.values():
            self.assertNotIn(seal, json.dumps(first['candidate']['record']))
        for key in ('asset_archive_sha256','asset_manifest_sha256'):
            self.assertNotIn(self.spec['binding'][key], json.dumps(first['candidate']['record']))
        self.assertEqual(first['seal']['preparation_certificate_sha256'], CERTIFICATE_SEALS['sat'])
        self.assertEqual(first['candidate']['record']['factor_base']['inventory']['usable_point_count'],62)
        self.assertTrue(first['candidate']['candidate_id'].startswith('IC1N17Ckb1fb62PDP3sat'))
        spec = copy.deepcopy(self.spec)
        spec['runtime_seal']['manifest_sha256'] = 'c'*64
        spec['binding'] = binding(spec)
        self.assertNotEqual(first['candidate'], self.registration(spec=spec)['candidate'])

    def test_freshness_flags_caps_and_unknown_fields_do_not_authorize_dispatch(self):
        for change in ('fresh','point','scalar','count','bool','extra','watchdog'):
            panel, spec = copy.deepcopy(self.panel), copy.deepcopy(self.spec)
            if change == 'fresh': panel['question'] = 'fresh-paired-qualification'
            elif change == 'point': panel['target_input']['point'] = [471,57570]
            elif change == 'scalar': panel['target_input']['known_scalar_supplied'] = 0
            elif change == 'count': panel['resources']['target_count'] = True
            elif change == 'bool': panel['max_descent_queries'] = True
            elif change == 'extra': panel['max_relation_queries'] = 1
            else: spec['runtime_watchdog_seconds'] += 1
            with self.assertRaises(InvalidEvidence):
                self.registration(panel,spec)

    def test_altered_or_swapped_preparation_is_not_silently_accepted(self):
        doc = copy.deepcopy(self.doc)
        doc['certificate']['inputs']['attempts'][0]['scalar'] += 1
        with self.assertRaises(InvalidEvidence):
            self.registration(doc=doc)
        f5 = json.loads((HERE/'goal_20260924/prepared-ic-state-v1/f5-preparation.json').read_text())
        with self.assertRaises(InvalidEvidence):
            mathematical_registration(self.panel,self.spec,self.files,f5,CERTIFICATE_SEALS['f5'])

    def test_bound_entry_and_auditor_math_controls_preserve_ids_and_reject_summary_tampering(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            out = root/'entry-output'
            out.mkdir()
            args = self.registration()
            spec = dict(copy.deepcopy(self.spec),arguments=args)
            (root/'execution.json').write_text(json.dumps(spec))
            preflight = dict(returncode=0,timed_out=False,scope='mocked native preflight')

            def fake_meter(*unused):
                (out/'cms_preflight.stdout').write_text('CryptoMiniSat version 5.14.7\n')
                return preflight

            def replay_query(panel,item,execution,curve,base,directory):
                self.assertEqual(item['point'],self.source_row['public_point'])
                row = copy.deepcopy(self.source_row)
                row.update(trial=item['trial'],verification_wall_ns=0)
                return row

            with patch('prepared_sat_runtime_v1.check_extracted_assets',return_value=self.files), \
                 patch('prepared_sat_runtime_v1.platform.system',return_value='Darwin'), \
                 patch('prepared_sat_runtime_v1.platform.machine',return_value='arm64'), \
                 patch('prepared_sat_runtime_v1.meter',side_effect=fake_meter), \
                 patch('prepared_sat_runtime_v1.one_query',side_effect=replay_query):
                result = run(args,out)
                self.assertEqual(result['status'],'COMPLETE')
                self.assertEqual(result['run_id'],args['seal']['run_id'])
                summary = json.loads((out/'summary.json').read_text())
                self.assertEqual(summary['recovered_scalar'],24886)
                self.assertFalse(summary['source_bound_execution_admitted'])
                with self.assertRaises(InvalidEvidence):
                    run(args,out)
            witness = self.source_row['point_witness']['point_indices']
            with patch('prepared_sat_runtime_v1.audit_execution',return_value={'entrypoint_succeeded':True,'scope':'mocked-source-control'}), \
                 patch('prepared_sat_runtime_v1.check_extracted_assets',return_value=self.files), \
                 patch('prepared_sat_runtime_v1.verify_query',return_value=(True,witness,'VALID_POINT_WITNESS')), \
                 patch('prepared_sat_runtime_v1.audit_meter',return_value=preflight):
                admitted = audit(root,spec)
                validate_admission(admitted, dict(spec, runtime_manifest={
                    'components': [{'role': HELPER_ROLE}]}))
                self.assertEqual(admitted['status'],'ADMITTED_COMPLETE_PREPARED_SAT_CONTROL')
                self.assertEqual(admitted['recovered_scalar'],24886)
                self.assertFalse(admitted['headline_online_admissible'] or admitted['promotion_eligible'])
                self.assertIsNone(admitted['online_speedup'])
                for key, value in [('online_wall_ns',summary['online_wall_ns']+1),
                                   ('recovered_scalar',1),('ordinary_queries_executed',1),
                                   ('target_status_mix',{})]:
                    bad = dict(summary, **{key:value})
                    (out/'summary.json').write_text(json.dumps(bad))
                    with self.assertRaises(InvalidEvidence): audit(root,spec)


if __name__ == '__main__':
    unittest.main()
