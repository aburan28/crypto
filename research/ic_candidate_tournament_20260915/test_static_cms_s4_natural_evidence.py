"""Replay all 32 raw SAT rows against the independent exact group oracle."""
import hashlib
import unittest

from audit_static_cms_s4_natural import audit
from run_static_cms_s4_natural import REGISTRATION
from tournament import read


EVIDENCE = REGISTRATION/'evidence.tar.gz'
EVIDENCE_SHA256 = 'e760b906acb441aeefd964f30ff3bb36b0be28a698afab84148350328a56648c'
RESULT_SHA256 = '0fcd6e548d32a8d8c5d56b31c0232b12e4edc6bb3e33d8aa495b95002738b4db'


class StaticCmsNaturalEvidenceTests(unittest.TestCase):
    def test_every_attempt_formula_and_exact_label_replays(self):
        self.assertEqual(hashlib.sha256(EVIDENCE.read_bytes()).hexdigest(),
                         EVIDENCE_SHA256)
        result_path = REGISTRATION/'RESULT.json'
        self.assertEqual(hashlib.sha256(result_path.read_bytes()).hexdigest(),
                         RESULT_SHA256)
        result = read(result_path)
        self.assertEqual(audit(EVIDENCE), result)
        self.assertEqual(result['attempts'], 32)
        self.assertEqual(result['exact_feasible'], 6)
        self.assertEqual(result['verified_witnesses'], 6)
        self.assertEqual(result['feasible_missed'], 0)
        self.assertEqual(result['statuses'],
                         {'SOURCE_UNSAT': 26, 'VALID_POINT_WITNESS': 6})
        self.assertFalse(result['full_sat_ic_admission'])
        self.assertIsNone(result['online_speedup'])


if __name__ == '__main__':
    unittest.main()
