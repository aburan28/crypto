"""Pre-execution controls for the registered static-source SAT pilot."""
import unittest

from run_cms_s4_controls import lift_source_assignment
from run_static_cms_s4_controls import PANEL, PANEL_SHA256, admit, sha_bytes
from tournament import read


class StaticCmsS4ControlsTests(unittest.TestCase):
    def test_parent_labels_source_receipt_and_panel_replay_portably(self):
        self.assertEqual(sha_bytes(PANEL.read_bytes()), PANEL_SHA256)
        panel = read(PANEL)
        _, curve, base, built = admit(panel, require_local_binary=False)
        self.assertEqual(len(base), 63)
        self.assertEqual(built['status'], 'pass')
        self.assertEqual([item['trial'] for item in panel['schedule']],
                         [0, 3, 4, 67, 71, 75])
        self.assertEqual([item['exact_relation_exists']
                          for item in panel['schedule']],
                         [False, False, True, True, True, True])
        self.assertEqual(curve.mul(curve.g, 63975), (120363, 39325))

    def test_coordinate_decoder_accepts_a_real_group_witness_only(self):
        _, curve, base, _ = admit(read(PANEL), require_local_binary=False)
        triple = [base[i] for i in (44, 58, 9)]
        xs = [point[0] for point in triple]
        model = [bool(x >> j & 1) for x in xs for j in range(6)]
        manifest = dict(ell=6,
                        factor_base_basis_bitmasks=[str(1 << j)
                                                    for j in range(6)])
        found = lift_source_assignment(model, manifest, curve, base,
                                       curve.mul(curve.g, 63975))
        self.assertTrue(found['group_replay'])
        self.assertEqual(found['x_coordinates'], xs)
        absent = lift_source_assignment(model, manifest, curve, base,
                                        curve.mul(curve.g, 19198))
        self.assertFalse(absent['group_replay'])


if __name__ == '__main__':
    unittest.main()
