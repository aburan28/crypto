"""Replay the complete static-SAT IC identity and previously unseen target."""
import unittest

from run_static_sat_full import PANEL, admit, target_from_seed
from static_sat_registration import REGISTRATION, identities, source_manifest
from tournament import read


class FullStaticSatRegistrationTests(unittest.TestCase):
    def test_registered_candidate_workload_and_target_replay(self):
        panel = read(PANEL)
        source, method, candidate, workload = identities(panel)
        self.assertEqual(source, read(REGISTRATION/'source-manifest.json'))
        self.assertEqual(source, source_manifest())
        self.assertEqual(method, read(REGISTRATION/'method.json'))
        self.assertEqual(candidate, read(REGISTRATION/'candidate.json'))
        self.assertEqual(workload, read(REGISTRATION/'workload.json'))
        self.assertEqual(candidate['record']['factor_base']['inventory']
                         ['usable_point_count'], 62)
        self.assertEqual(candidate['record']['factor_base']['inventory']
                         ['effective_columns'], 29)
        _, curve, _, target, _, matrix, seal = admit(
            panel, require_local_binary=False)
        self.assertEqual(target_from_seed(curve, 2026092938),
                         (0, (114119, 85674)))
        self.assertEqual(target, (114119, 85674))
        self.assertEqual(len(matrix.columns), 29)
        self.assertEqual(seal['candidate_id'], candidate['candidate_id'])


if __name__ == '__main__':
    unittest.main()
