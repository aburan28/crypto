"""Adversarial premeasurement layout and independent base controls."""
import copy
import json
from pathlib import Path
import unittest

from generic_bases import construct, verify_base
from generic_solver_feasibility import (SOURCE_OBJECTS, assess, basis_length,
                                        check_source, check_source_checkout)
from identity import sha256
from oracle import Curve, InvalidEvidence

HERE = Path(__file__).resolve().parent
PANEL = HERE / 'goal_20260924/generic-backend-qualification-v2/panel.json'
STATIC_AUDIT = HERE / 'goal_20260924/generic-backend-qualification-v2/STATIC-FEASIBILITY.json'
EXPOSURES = HERE / 'goal_20260924/generic-backend-qualification-v2/lost-campaign-exposures.json'
INVENTORY = HERE / 'goal_20260924/generic-backend-qualification-v2/standard-subspace-d6-inventory-control.json'


class GenericSolverFeasibilityTests(unittest.TestCase):
    def test_reviewed_source_and_entire_registered_panel_expose_hard_failure(self):
        check_source(HERE.parents[1])
        with self.assertRaisesRegex(InvalidEvidence, 'not the reviewed worker'):
            check_source_checkout(HERE.parents[1])
        result = assess(json.loads(PANEL.read_text()))
        self.assertEqual(result, json.loads(STATIC_AUDIT.read_text()))
        self.assertEqual(result['status'], 'FAIL_STATIC_LAYOUT')
        algebraic = [row for row in result['rows'] if row['status'] == 'ABOVE_LAYOUT_CAP']
        self.assertEqual(len(algebraic), 20)  # four F4/F5 arms × five cells
        self.assertEqual({row['boolean_variables'] for row in algebraic},
                         {68, 76, 92, 124})
        self.assertEqual({row['basis_length'] for row in algebraic}, {17, 19, 23, 31})
        sat = [row for row in result['rows'] if row['solver'] in ('sat_xor', 'sat_cnf')]
        self.assertEqual(len(sat), 2)
        self.assertTrue(all(row['status'] == 'SEPARATE_SAT_PATH_UNCHECKED' for row in sat))

    def test_smaller_base_clears_only_static_cap_and_rejects_bad_layout(self):
        panel = json.loads(PANEL.read_text())
        for candidate in panel['candidates']:
            if candidate['config']['solver'] in ('f4', 'f5', 'inherited_f4'):
                candidate['config']['factor_base'] = dict(recipe='standard_subspace', dimension=6)
        result = assess(panel)
        self.assertEqual(result['status'], 'PASS_STATIC_LAYOUT_ONLY')
        algebraic = [row for row in result['rows'] if row['status'] == 'WITHIN_LAYOUT_CAP_ONLY']
        self.assertEqual(len(algebraic), 20)
        self.assertEqual({row['boolean_variables'] for row in algebraic},
                         {35, 37, 41, 49})
        for bad in (0, True, 21, 31, '6'):
            changed = copy.deepcopy(panel)
            changed['candidates'][4]['config']['factor_base']['dimension'] = bad
            with self.assertRaises(InvalidEvidence):
                assess(changed)
        with self.assertRaises(InvalidEvidence):
            basis_length({'recipe': 'frobenius_union', 'seed_masks': [1]}, 17)
        self.assertEqual(basis_length({'kind': 'standard_subspace', 'dimension': 6}, 17), 6)
        self.assertEqual(basis_length({'kind': 'subgroup_orbits', 'seed': 43, 'points': 68}, 17), 17)
        self.assertEqual(basis_length({'recipe': 'subgroup_orbits', 'seed': 43,
                                       'requested_points': '4n'}, 17), 17)
        with self.assertRaises(InvalidEvidence):
            basis_length({'kind': 'subgroup_orbits', 'recipe': 'standard_subspace',
                          'seed': 43, 'points': 68}, 17)
        changed = copy.deepcopy(panel)
        changed['candidates'] = [c for c in changed['candidates']
                                 if c['config']['solver'] not in ('f4', 'f5', 'inherited_f4')]
        with self.assertRaisesRegex(InvalidEvidence, 'no algebraic arm'):
            assess(changed)
        changed = copy.deepcopy(panel)
        changed['candidates'][4]['id'] = changed['candidates'][5]['id']
        with self.assertRaisesRegex(InvalidEvidence, 'duplicate candidate'):
            assess(changed)

    def test_independent_polynomial_subspace_reconstructs_ordered_points(self):
        fixture = json.loads(EXPOSURES.read_text())['attempts'][0]['fixture']
        curve = Curve(fixture)
        points, detail = construct(curve, dict(kind='standard_subspace', dimension=6))
        self.assertEqual(detail['nominal_dimension'], 6)
        self.assertEqual(detail['basis'], [1, 2, 4, 8, 16, 32])
        self.assertEqual(detail['domain_abscissae'], 64)
        self.assertEqual(points[0], (0, 1))
        self.assertEqual(len(points), len(set(points)))
        self.assertEqual([x for x, _ in points], sorted(x for x, _ in points))
        self.assertTrue(all(x < 64 and curve.decode([str(x), str(y)]) == (x, y)
                            for x, y in points))
        for bad in (0, True, 17, 21, '6'):
            with self.assertRaises(InvalidEvidence):
                construct(curve, dict(kind='standard_subspace', dimension=bad))
        with self.assertRaises(InvalidEvidence):
            construct(curve, dict(kind='standard_subspace', dimension=6, seed=43))

    def test_disclosed_worker_inventories_replay_against_independent_base(self):
        control = json.loads(INVENTORY.read_text())
        self.assertEqual(control['status'], 'BASE_CONSTRUCTION_CONTROL_ONLY')
        self.assertEqual(control['producer_factor_base_source_object'],
                         SOURCE_OBJECTS['src/cryptanalysis/koblitz_index_calculus.rs'])
        self.assertEqual(len(control['rows']), 5)
        exposures = {f"n{item['fixture']['degree']}a{item['fixture']['curve_a']}": item['fixture']
                     for item in json.loads(EXPOSURES.read_text())['attempts'][:5]}
        for row in control['rows']:
            report, job = row['worker_report'], row['job']
            self.assertEqual(sha256(report), row['worker_report_sha256'])
            self.assertEqual(row['cell'], f"n{job['degree']}a{job['curve_a']}")
            original = exposures[row['cell']]
            self.assertEqual(sha256(original), row['fixture_sha256'])
            replay = copy.deepcopy(report['fixture'])
            replay['target_seeds'] = original['target_seeds']
            self.assertEqual(replay, original)
            self.assertEqual(verify_base(report, report['fixture'], job),
                             row['independent_receipt'])
            self.assertGreater(row['independent_receipt']['inventory']['usable_point_count'], 0)
            self.assertGreater(row['independent_receipt']['inventory']['effective_columns'], 0)


if __name__ == '__main__':
    unittest.main()
