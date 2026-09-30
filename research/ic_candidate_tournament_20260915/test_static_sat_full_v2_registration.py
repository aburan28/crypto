"""Replay the corrected SAT candidate's frozen point and exact identity."""
import unittest

from run_static_sat_full_v2 import PANEL,admit,target_from_seed
from static_sat_registration_v2 import REGISTRATION,identities,source_manifest
from tournament import read


class FullStaticSatV2RegistrationTests(unittest.TestCase):
    def test_identities_and_paired_target_replay(self):
        panel=read(PANEL)
        source,method,candidate,workload=identities(panel)
        self.assertEqual(source,read(REGISTRATION/'source-manifest.json'))
        self.assertEqual(source,source_manifest())
        self.assertEqual(method,read(REGISTRATION/'method.json'))
        self.assertEqual(candidate,read(REGISTRATION/'candidate.json'))
        self.assertEqual(workload,read(REGISTRATION/'workload.json'))
        self.assertEqual(workload['record']['targets'],[[52411,72106]])
        _,curve,_,target,_,matrix,seal=admit(
            panel,require_local_binary=False)
        self.assertEqual(target_from_seed(curve,2026092948,
                                          'ic-paired-target-v1'),
                         (1,(52411,72106)))
        self.assertEqual(target,(52411,72106))
        self.assertEqual(len(matrix.columns),29)
        self.assertEqual(seal['candidate_id'],candidate['candidate_id'])


if __name__=='__main__':
    unittest.main()
