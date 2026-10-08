"""Replay the closed complete SAT control; never execute its native producer."""
import json
from pathlib import Path
import tempfile
import unittest

from audit_static_sat_full_v3 import audit
from oracle import InvalidEvidence
from publish_static_sat_v3_control import replay

BUNDLE = (Path(__file__).resolve().parent / 'goal_20260924/static-sat-runtime-v3/'
          'full-development-20260929/results-20260930')
EXPECTED = '4e6547c03cfe1fa54ae9a994c160e5d4afa987659dc72880c651885a0e9e8fa9'


class StaticSatV3FullDevelopmentTests(unittest.TestCase):
    def test_complete_transport_and_target_certificate_faults(self):
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / 'transport'
            result = replay(BUNDLE, output, EXPECTED)
            self.assertEqual(result['status'], 'AUDITED_COMPLETE_SINGLE_ARM')
            self.assertTrue(result['complete_source_bound'])
            self.assertTrue(result['online_endpoint_admissible'])
            self.assertEqual((result['attempts'], result['exact_feasible'],
                              result['verified_relations'], result['final_rank']),
                             (149, 37, 37, 29))
            self.assertEqual(result['statuses'], {'VALID_POINT_WITNESS': 37,
                             'SOURCE_UNSAT': 106, 'CONFLICT_BUDGET_INCONCLUSIVE': 6})
            self.assertEqual(result['recovered_scalar'], '24886')
            self.assertEqual(result['online_wall_ns'], 708643292)
            self.assertEqual(sum(result['online_phases_ns'].values()),
                             result['online_wall_ns'])
            self.assertFalse(result['headline_online_admissible'])
            self.assertFalse(result['fresh_target_qualified'])
            self.assertFalse(result['same_point_rho_audited'])
            self.assertFalse(result['promotion_eligible'])
            self.assertIsNone(result['online_speedup'])

            spec = json.loads((output / 'registration/execution.json').read_text())
            summary = output / 'execution/entry-output/summary.json'
            original = summary.read_text()
            changed = json.loads(original)
            changed['recovered_scalar'] = '24887'
            summary.write_text(json.dumps(changed))
            with self.assertRaises(InvalidEvidence):
                audit(output / 'execution', spec)
            summary.write_text(original)

            # Changing native argv invalidates the source-bound receipt before
            # any target scalar could be used to excuse a wrong computation.
            intent = output / 'execution/entry-output/collection/trial-00/cms.intent.json'
            changed = json.loads(intent.read_text())
            changed['arguments'][changed['arguments'].index('--maxconfl') + 1] = '1000001'
            intent.write_text(json.dumps(changed))
            with self.assertRaisesRegex(InvalidEvidence, 'command or watchdog'):
                audit(output / 'execution', spec)

    def test_closed_invocation_requires_external_seal(self):
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / 'wrong-seal'
            with self.assertRaisesRegex(InvalidEvidence, 'externally frozen'):
                replay(BUNDLE, output, '0' * 64)
            self.assertFalse(output.exists())


if __name__ == '__main__':
    unittest.main()
