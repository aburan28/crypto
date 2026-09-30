"""Real admitted native inputs and complete identities; no solver dispatch."""
import copy
import json
from pathlib import Path
import tempfile
import unittest

from oracle import InvalidEvidence
from sat_runtime_execution_v3 import binding,register
from static_sat_assets_v3 import verified_assets
from static_sat_inputs_v3 import native_admission
from static_sat_registration_v3 import mathematical_registration

ROOT=Path(__file__).resolve().parents[2]
ASSETS=Path(__file__).resolve().parent/'goal_20260924/static-sat-runtime-v3/native-inputs-macos-arm64'


class StaticSatRegistrationV3Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        manifest=json.loads((ASSETS/'manifest.json').read_text())
        seal=json.loads((ASSETS/'seal.json').read_text())
        cls.files=verified_assets(ASSETS,manifest,seal)
        cls.spec=dict(entrypoint=dict(module='static_sat_pipeline_v3',callable='run'),
                      runtime_watchdog_seconds=600,runtime_seal=dict(
                          manifest_sha256='a'*64,archive_sha256='b'*64),
                      asset_manifest=manifest,asset_seal=seal,interpreter={'control':True})
        cls.spec['binding']=binding(cls.spec)
        cls.panel=dict(cms_conflict_budget=1_000_000,cms_timeout_seconds=120,
                       export_timeout_seconds=60,max_relation_queries=256,max_descent_queries=64,
                       relation_query_seed=2026093001,descent_query_seed=2026093002,
                       export_nonce=2026093003,run_number=0,question='development-source-control',
                       target_input=dict(point=[52411,72106],seed=None,
                           input_law='one-supplied-public-point-development-control',
                           point_was_previously_supplied=True,known_scalar_supplied=False),
                       resources=dict(host_class='physical-macos-arm64',cpu_workers=1,
                           target_count=1,memory_limit_bytes=None,total_wall_limit_seconds=600))

    def test_seed_changes_workload_and_source_changes_candidate(self):
        first=mathematical_registration(self.panel,self.spec,self.files)
        panel=copy.deepcopy(self.panel)
        panel['descent_query_seed']+=1
        second=mathematical_registration(panel,self.spec,self.files)
        self.assertEqual(first['candidate'],second['candidate'])
        self.assertNotEqual(first['workload'],second['workload'])
        changed=copy.deepcopy(self.spec)
        changed['runtime_seal']['manifest_sha256']='c'*64
        changed['binding']=binding(changed)
        third=mathematical_registration(self.panel,changed,self.files)
        self.assertNotEqual(first['candidate'],third['candidate'])
        self.assertEqual(first['workload'],third['workload'])
        self.assertTrue(first['candidate']['candidate_id'].startswith('IC1N17Ckb1fb62PDP3sat'))
        self.assertEqual(first['candidate']['record']['factor_base']['inventory']['effective_columns'],29)
        self.assertEqual(first['seal']['registration_stage'],'before-execution')

    def test_complete_source_factory_is_sealed_before_any_execution(self):
        with tempfile.TemporaryDirectory() as temporary:
            out=Path(temporary)/'registration'
            spec=register(ROOT,out,module='static_sat_pipeline_v3',action='run',arguments=None,
                          timeout_seconds=600,asset_snapshot=ASSETS,
                          arguments_factory=lambda s:mathematical_registration(self.panel,s,self.files))
            self.assertEqual(spec['arguments'],mathematical_registration(self.panel,spec,self.files))
            self.assertFalse((out/'process.json').exists())
            self.assertTrue((out/'registration-seal.json').is_file())
            self.assertIn('asset_archive_sha256',spec['binding'])

    def test_mutated_native_sources_planted_inputs_and_freshness_flags_fail_closed(self):
        changed=dict(self.files)
        changed['exporter/source.rs']+=b'// changed\n'
        with self.assertRaisesRegex(InvalidEvidence,'source/build'):
            native_admission(changed)
        fixture=json.loads(self.files['fixture.json'])
        fixture['targets']=[[52411,72106]]
        changed=dict(self.files,**{'fixture.json':json.dumps(fixture).encode()})
        with self.assertRaisesRegex(InvalidEvidence,'target data leaked'):
            native_admission(changed)
        panel=copy.deepcopy(self.panel)
        panel['question']='fresh-paired-qualification'
        panel['target_input']['point_was_previously_supplied']=False
        with self.assertRaisesRegex(InvalidEvidence,'campaign adapter'):
            mathematical_registration(panel,self.spec,self.files)


if __name__=='__main__':
    unittest.main()
