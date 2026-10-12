"""Replay the actual one-shot F5 v2 result; native execution is forbidden."""
import copy
import hashlib
import json
from pathlib import Path
import tempfile
import unittest

from audit_f5_runtime_v2 import audit as legacy_audit
from generic_queries_exact_v1 import verify_queries
from oracle import InvalidEvidence
from publish_f5_runtime_v3 import frozen_audit, replay

BUNDLE = Path(__file__).resolve().parent/'goal_20260924/f5-source-bound-runtime-v2/results-20260930'
ARCHIVE_SHA256 = 'd62ff4b0e5848b7e36474d1a6af2be7fef1ee060d6cc05ae1e02c86c2957c1b2'
EXECUTION_SHA256 = '9409cb6c40548816153f5b9978e8fa0637b5aaa2801472575e7cc55a9f61d56d'


class F5V2FullDevelopmentTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        assert hashlib.sha256((BUNDLE/'evidence.tar.gz').read_bytes()).hexdigest() == ARCHIVE_SHA256
        cls.temporary = tempfile.TemporaryDirectory()
        cls.root = Path(cls.temporary.name)/'transport'
        cls.result = replay(BUNDLE, cls.root, EXECUTION_SHA256)
        cls.report = json.loads((cls.root/'execution/entry-output/pipeline.stdout').read_text())

    @classmethod
    def tearDownClass(cls):
        cls.temporary.cleanup()

    def test_actual_complete_pipeline_closes_all_proofs_and_charged_intervals(self):
        d = self.result
        self.assertEqual(d['status'], 'AUDITED_COMPLETE_IC')
        self.assertTrue(d['complete_ic_admitted'])
        self.assertEqual((d['final_rank'],d['verified_target_count']), (29,1))
        self.assertEqual(d['online_wall_ns'],10512454542)
        self.assertEqual(sum(x for x in d['online_phases_ns'].values() if x is not None),d['online_wall_ns'])
        self.assertTrue(d['run']['cold_phase_ledger']['complete'])
        natural = d['natural_query_audit']
        self.assertEqual((natural['attempts'],natural['verified_native_witnesses']), (216,61))
        self.assertEqual(natural['status_mix'], {'proved_unsat':155,'witness':61})
        self.assertEqual(natural['feasible_but_no_witness'],0)
        proofs = d['native_only_admission']['stages']['query_law']['negative_proofs']
        self.assertEqual(proofs['base_size'],63)
        self.assertEqual(proofs['independently_proved_negative_queries'],157)
        self.assertEqual(sum(row['stage']=='target' for row in proofs['rows']),2)
        self.assertEqual(d['native_only_admission']['stages']['certificate']['solutions'],['24886'])
        gate = json.loads((self.root/'frozen-replay-source-gate.json').read_text())
        self.assertFalse(gate['native_execution'])
        self.assertTrue(gate['isolated'] and gate['site_disabled'])
        self.assertIn('generic_exact_negative_v1',gate['after'])
        self.assertFalse(d['headline_online_admissible'] or d['promotion_eligible'])
        self.assertIsNone(d['online_speedup'])

    def test_original_v2_negative_gate_still_rejects_the_unchanged_execution(self):
        spec=json.loads((self.root/'registration/execution.json').read_text())
        with self.assertRaisesRegex(InvalidEvidence,'negative claim needs an independent proof adapter'):
            legacy_audit(self.root/'execution',spec)

    def test_forged_ordinary_or_target_negative_cannot_use_the_new_gate(self):
        for stage in ['ordinary','target']:
            report=copy.deepcopy(self.report)
            attempts=([row for batch in report['collection_reports'] for row in batch['attempts']]
                      if stage=='ordinary' else report['solutions'][0]['attempts'])
            attempt=next(row for row in attempts if row['pdp']['outcome']=='witness')
            attempt['pdp'].update(outcome='proved_unsat',points=None)
            with self.subTest(stage=stage),self.assertRaisesRegex(InvalidEvidence,'false negative'):
                verify_queries(report,report['fixture'],3)

    def test_wrong_external_invocation_hash_cannot_extract_evidence(self):
        out=Path(self.temporary.name)/'wrong-seal'
        with self.assertRaisesRegex(InvalidEvidence,'externally frozen'):
            replay(BUNDLE,out,'0'*64)
        self.assertFalse(out.exists())

    def test_changed_retained_auditor_source_fails_before_mathematical_replay(self):
        path=self.root/'post-auditor-source/research/ic_candidate_tournament_20260915/generic_exact_negative_v1.py'
        original=path.read_bytes()
        try:
            path.chmod(0o644)
            path.write_bytes(original+b'\n# changed retained source\n')
            with self.assertRaisesRegex(InvalidEvidence,'proof-adapter replay failed'):
                frozen_audit(self.root,EXECUTION_SHA256)
            self.assertIn('retained F5 auditor source differs',
                          (self.root/'frozen-replay.stderr').read_text())
        finally:
            path.write_bytes(original)
            path.chmod(0o444)


if __name__=='__main__':
    unittest.main()
