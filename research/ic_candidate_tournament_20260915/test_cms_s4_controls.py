"""Pre-execution controls for the disclosed wide-S4 external SAT pilot."""
import unittest

from run_cms_s4_controls import (
    PANEL, PANEL_SHA256, lift_source_assignment, preflight, sha_bytes,
)
from tournament import read


class CmsS4ControlsTests(unittest.TestCase):
    def test_frozen_source_binary_and_disclosed_points(self):
        self.assertEqual(sha_bytes(PANEL.read_bytes()), PANEL_SHA256)
        panel = read(PANEL)
        _, curve, base = preflight(panel, require_local_solver=False)
        self.assertEqual(len(base), 65)
        self.assertEqual([p['trial'] for p in panel['schedule']], [0, 1, 3, 10])
        self.assertEqual([p['exact_relation_exists'] for p in panel['schedule']],
                         [False, False, True, True])
        self.assertEqual(panel['cms_max_models_per_query'], 1)

    def test_coordinate_decoder_requires_a_group_witness(self):
        _, curve, base = preflight(read(PANEL), require_local_solver=False)
        positive = [base[i] for i in (4, 22, 45)]
        xs = [point[0] for point in positive]
        model = [bool(x >> j & 1) for x in xs for j in range(6)]
        manifest = dict(ell=6, factor_base_basis_bitmasks=[str(1 << i) for i in range(6)])
        found = lift_source_assignment(model, manifest, curve, base, curve.mul(curve.g, 81459))
        self.assertTrue(found['group_replay'])
        self.assertEqual(found['x_coordinates'], xs)
        absent = lift_source_assignment(model, manifest, curve, base, curve.mul(curve.g, 38307))
        self.assertFalse(absent['group_replay'])


if __name__ == '__main__':
    unittest.main()
