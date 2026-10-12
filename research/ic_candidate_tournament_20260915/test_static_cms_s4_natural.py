"""Before-dispatch checks for the fresh, label-free natural SAT schedule."""
import hashlib
import unittest

from run_static_cms_s4_natural import PANEL, PANEL_SHA256, admit
from tournament import read


class StaticCmsNaturalTests(unittest.TestCase):
    def test_frozen_new_query_law_has_no_outcome_labels(self):
        self.assertEqual(hashlib.sha256(PANEL.read_bytes()).hexdigest(),
                         PANEL_SHA256)
        panel = read(PANEL)
        _, _, base, built = admit(panel, require_local_binary=False)
        self.assertEqual(len(base), 63)
        self.assertEqual(built['status'], 'pass')
        self.assertEqual([row['trial'] for row in panel['schedule']],
                         list(range(32)))
        self.assertTrue(all(set(row) == {'trial', 'probe_scalar', 'point'}
                            for row in panel['schedule']))
        self.assertEqual(len({tuple(row['point'])
                              for row in panel['schedule']}), 32)


if __name__ == '__main__':
    unittest.main()
