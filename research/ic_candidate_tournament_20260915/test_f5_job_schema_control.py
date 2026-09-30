"""Replay the closed native parser control; never call execute_job."""
import json
import unittest

from f5_job_schema_control import CONTROL
from f5_runtime_registration_v2 import SCHEMA_CONTROL_SHA256
from oracle import InvalidEvidence
from replay_f5_job_schema_control import replay, replay_files, REGISTRATION
from replay_paired_n17_evidence import retained_files


class F5JobSchemaControlTests(unittest.TestCase):
    def test_closed_transport_defaults_and_six_actual_native_rejections(self):
        result = replay(CONTROL/'native', SCHEMA_CONTROL_SHA256)
        self.assertEqual(result['status'], 'PASS_NATIVE_JOB_SCHEMA_ONLY')
        self.assertEqual((result['parsed_cases'], result['rejected_cases']), (2, 6))
        self.assertFalse(result['complete_ic_admitted'])
        self.assertIsNone(result['measured_costs'])
        self.assertFalse(result['promotion_eligible'])
        self.assertIsNone(result['online_speedup'])
        files = retained_files(CONTROL/'native')
        raw = json.loads(files['native.stdout'])
        null = next(row for row in raw['results'] if row['label'] == 'historical-null-seed')
        self.assertIn('invalid type: null, expected u64', null['error'])
        self.assertEqual(raw['results'][0]['observed']['target_seeds'], [])

    def test_source_input_binary_output_and_external_seal_tampering_fail_closed(self):
        files = retained_files(CONTROL/'native')
        for role in ('diagnostic', 'kernel-original.rs', REGISTRATION+'kernel_append.rs',
                     REGISTRATION+'inputs.json', 'native.stdout', 'native.stderr',
                     'analysis-src/research/ic_candidate_tournament_20260915/generic_stages.py'):
            with self.subTest(role=role):
                changed = dict(files)
                changed[role] += b'changed'
                with self.assertRaises((InvalidEvidence, ValueError)):
                    replay_files(changed)
        with self.assertRaisesRegex(InvalidEvidence, 'externally retained archive seal'):
            replay(CONTROL/'native', '0'*64)


if __name__ == '__main__':
    unittest.main()
