"""The v3 smoke registration stays locked and undispatched."""
import unittest

from oracle import InvalidEvidence
from run_generic_backend_qualification_v3 import (
    DISPATCH_AUTHORIZED, LOST_V2, PANEL, PANEL_SHA256, exposure_block,
    registration_check, status_report)
from tournament import digest, read


class V3RegistrationTest(unittest.TestCase):
    def test_panel_bytes_and_layout_stay_locked(self):
        self.assertEqual(digest(PANEL), PANEL_SHA256)
        panel = read(PANEL)
        layout = registration_check(panel)
        self.assertEqual(layout['status'], 'PASS_STATIC_LAYOUT_ONLY')
        self.assertEqual(layout['impossible_algebraic_cells'], 0)
        algebraic = [row for row in layout['rows'] if row['solver'] in ('f4', 'f5')]
        self.assertTrue(algebraic)
        self.assertTrue(all(row['boolean_variables'] <= 64 for row in algebraic))

    def test_seed_is_not_dispatched_without_the_v2_census(self):
        self.assertFalse(DISPATCH_AUTHORIZED)
        panel = read(PANEL)
        report = status_report(panel)
        self.assertFalse(report['dispatch_authorized'])
        self.assertEqual(report['measurement'], 'not_run')
        self.assertEqual(report['seed'], 2026093001)
        if not LOST_V2.is_file():
            self.assertEqual(exposure_block(), 'v2 exposure census is not on this checkout')
            self.assertIn('v2 exposure census', report['dispatch_block'])
        else:
            self.assertIsNone(exposure_block())
            self.assertTrue(report['v2_exposure_present'])
            self.assertEqual(report['dispatch_block'], 'dispatch path is not in this checker')

    def test_a_changed_seed_is_rejected(self):
        panel = read(PANEL)
        panel['seed'] = 2026092902
        with self.assertRaises(InvalidEvidence):
            registration_check(panel)


if __name__ == '__main__':
    unittest.main()
