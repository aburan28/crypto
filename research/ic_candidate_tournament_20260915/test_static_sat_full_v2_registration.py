"""Replay the corrected SAT candidate's frozen point and exact identity."""
import unittest

from frozen_sat_runtime import replay
from static_sat_registration_v2 import REGISTRATION
from tournament import read


class FullStaticSatV2RegistrationTests(unittest.TestCase):
    def test_identities_and_paired_target_replay(self):
        result=replay('v2')
        candidate=read(REGISTRATION/'candidate.json')
        workload=read(REGISTRATION/'workload.json')
        self.assertEqual(result['candidate_id'],candidate['candidate_id'])
        self.assertEqual(result['workload_id'],workload['workload_id'])
        self.assertEqual(workload['record']['targets'],[[52411,72106]])
        self.assertEqual(result['target_counter'],1)
        self.assertEqual(result['target'],[52411,72106])
        self.assertEqual(result['actual_usable_points'],62)
        self.assertEqual(result['columns'],29)
        self.assertFalse(result['complete_preexecution_python_manifest'])
        self.assertFalse(result['promotion_eligible'])


if __name__=='__main__':
    unittest.main()
