"""Replay the complete static-SAT IC identity and previously unseen target."""
import unittest

from frozen_sat_runtime import replay
from static_sat_registration import REGISTRATION
from tournament import read


class FullStaticSatRegistrationTests(unittest.TestCase):
    def test_registered_candidate_workload_and_target_replay(self):
        result = replay('v1')
        self.assertEqual(result['candidate_id'],
                         read(REGISTRATION/'candidate.json')['candidate_id'])
        self.assertEqual(result['workload_id'],
                         read(REGISTRATION/'workload.json')['workload_id'])
        self.assertEqual(result['actual_usable_points'], 62)
        self.assertEqual(result['columns'], 29)
        self.assertEqual(result['target_counter'], 0)
        self.assertEqual(result['target'], [114119, 85674])
        self.assertFalse(result['complete_preexecution_python_manifest'])
        self.assertFalse(result['promotion_eligible'])


if __name__ == '__main__':
    unittest.main()
