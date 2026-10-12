"""Premeasurement controls for the disclosed F4/F5/SAT recovery panel."""
import unittest

from run_generic_backend_recovery_pilot import (
    CELLS, PANEL, PANEL_SHA256, SOLVERS, declared_job, digest,
    fixtures_from_inventory, resources, static_preflight, validate_panel,
)
from tournament import read


class RecoveryPilotTests(unittest.TestCase):
    def test_registered_source_inputs_and_exact_fixed_schedule(self):
        self.assertEqual(digest(PANEL), PANEL_SHA256)
        panel = read(PANEL)
        validate_panel(panel)
        fixtures = fixtures_from_inventory(panel)
        self.assertEqual(set(fixtures), set(CELLS))
        self.assertEqual(len(CELLS) * len(SOLVERS), 20)
        for cell in CELLS:
            jobs = [declared_job(panel, cell, solver, fixtures[cell])
                    for solver in SOLVERS]
            self.assertEqual(len({str(job['public_targets']) for job in jobs}), 1)
            self.assertEqual({job['algorithm_seed'] for job in jobs},
                             {panel['algorithm_seed']})
            self.assertTrue(all(job['exclusive_phases'] for job in jobs))
            self.assertTrue(all(job['factor_base'] == panel['factor_base'] for job in jobs))
            self.assertTrue(all(job['config']['max_trials'] ==
                                panel['per_cell_limits'][cell]['max_trials'] for job in jobs))
            self.assertEqual(resources(panel, cell)['rayon_threads'], 1)
        self.assertEqual(fixtures['n17a1']['targets'], [['853', '39791']])

    def test_algebraic_layout_preflight_and_sat_separation(self):
        result = static_preflight(read(PANEL))
        self.assertEqual(result['status'], 'PASS_STATIC_LAYOUT_ONLY')
        self.assertEqual(result['impossible_algebraic_cells'], 0)
        algebra = [row for row in result['rows'] if row.get('cell')]
        self.assertEqual(len(algebra), 10)
        self.assertTrue(all(row['status'] == 'WITHIN_LAYOUT_CAP_ONLY'
                            for row in algebra))
        sat = [row for row in result['rows'] if row['candidate'].startswith('sat_')]
        self.assertEqual(len(sat), 2)
        self.assertTrue(all(row['status'] == 'SEPARATE_SAT_PATH_UNCHECKED'
                            for row in sat))


if __name__ == '__main__':
    unittest.main()
