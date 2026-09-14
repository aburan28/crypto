"""Numerical boundaries, frozen certificates, and evidence-preserving CLI checks."""
import hashlib
import itertools
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from ec_index_calculus_budget import budget, monomials

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / 'pilot'))
from run_pilot import GF, anf_check, check_sat


class BudgetTests(unittest.TestCase):
    def example(self, **kw):
        args = dict(order=2**60, signed_base=2**20, columns=2**19,
                    m=3, rho_orbit=2, rho_rate=1e6)
        args.update(kw)
        return budget(**args)

    def test_published_numerical_example(self):
        r = self.example()
        self.assertEqual(r['rows_target_not_rank_guarantee'], 576717)
        self.assertAlmostEqual(r['attempt_budget_s'] * 1e6, 253.304, places=3)
        self.assertEqual([monomials(30, d)['column_universe'] for d in (4, 5, 6)],
                         [31931, 174437, 768212])
        frozen = json.loads((ROOT / 'ec_index_calculus_budget_example.json').read_text())
        self.assertEqual({k: frozen[k] for k in r}, r)

    def test_no_success_is_not_a_finite_attack(self):
        r = self.example(success=0, attempt_s=1)
        self.assertIsNone(r['expected_attempts'])
        self.assertIsNone(r['modeled_full_attack_s'])
        self.assertIsNone(r['modeled_rho_over_ic'])

    def test_all_other_costs_are_charged(self):
        r = self.example(other_s=1000, attempt_s=0.1)
        self.assertTrue(r['other_costs_already_exhaust_rho_budget'])
        self.assertEqual(r['attempt_budget_s'], 0)
        self.assertGreater(r['modeled_full_attack_s'], 1000)
        self.assertLess(r['modeled_rho_over_ic'], 1)

    def test_orbit_does_not_compress_columns_implicitly(self):
        plain, quotient = self.example(), self.example(rho_orbit=8)
        self.assertAlmostEqual(plain['rho_expected_steps'] / quotient['rho_expected_steps'], 2)
        self.assertEqual(plain['expected_attempts'], quotient['expected_attempts'])

    def test_count_floor_covers_actual_repeated_sums(self):
        # Compare the exact support in a small cyclic group, including collisions.
        support = [0, 1, 4]
        actual = len({sum(t) % 101 for t in itertools.combinations_with_replacement(support, 3)}) / 101
        r = budget(101, 3, 2, 3, 1, 1, success=actual)
        self.assertLessEqual(actual, r['p_counting_upper_bound_uniform_target'])
        self.assertGreaterEqual(r['expected_attempts'], r['attempt_count_floor_first_relation_mode'])

    def test_invalid_inputs_cannot_make_plausible_forecasts(self):
        for kw in ({'rho_rate': math.inf}, {'other_s': math.nan}, {'success': -1},
                   {'usable': 0}, {'columns': 2**21}, {'m': 2.5}, {'attempt_s': -1}):
            with self.subTest(kw=kw), self.assertRaises(ValueError):
                self.example(**kw)


class EvidenceTests(unittest.TestCase):
    def test_frozen_artifact_integrity(self):
        manifest = json.loads((ROOT / 'artifact_manifest.json').read_text())
        for item in manifest:
            with self.subTest(path=item['path']):
                self.assertEqual(hashlib.sha256((ROOT / item['path']).read_bytes()).hexdigest(), item['sha256'])

    def test_frozen_summary_and_point_certificates(self):
        summary = json.loads((ROOT / 'summary.json').read_text())
        evidence = json.loads((ROOT / 'pilot/results/pilot_results.json').read_text())
        self.assertIsNone(summary['full_dlp_S'])
        self.assertIsNone(summary['cost_over_matched_rho'])
        self.assertEqual(sum(r['exhaustive_s4_checks'] for r in evidence['cases']), 97504)
        for raw, row in zip(evidence['cases'], summary['cases'], strict=True):
            self.assertEqual(raw['case'], row['case'])
            for key in ('status', 'median_wall_s', 'conflicts', 'normalized_sha256'):
                self.assertEqual(raw[key], row[key])
            if raw['status'] == 'SAT':
                field = GF(row['n'], row['field_modulus'])
                pts = [tuple(p) for p in raw['point_witness']['summands']]
                self.assertTrue(all(field.on_curve(p) for p in pts))
                total = field.add(field.add(pts[0], pts[1]), pts[2])
                self.assertEqual(total, tuple(raw['point_witness']['sum']))
                self.assertEqual(total[0], row['target_x'])
                self.assertEqual(field.s4(*(p[0] for p in pts), row['target_x']), 0)

    def test_bad_assignments_and_missing_point_lifts_are_rejected(self):
        anf_check('p anf 1 1\nx 1 0\n', '1')
        for bits in ('0', '2', ''):
            with self.subTest(bits=bits), self.assertRaises(AssertionError):
                anf_check('p anf 1 1\nx 1 0\n', bits)
        class NoLift:
            def s4(self, *args): return 0
            def lift(self, x): return []
        with self.assertRaisesRegex(AssertionError, 'no verified point relation'):
            check_sat(NoLift(), 'p anf 3 1\nx 1 0\n', '100', 1, 0)

    def test_existing_evidence_is_refused_before_fetch_or_build(self):
        with tempfile.TemporaryDirectory() as existing:
            for script, args in (
                ('build_pilot.py', ['--cache-dir', '/unused', '--build-dir', existing]),
                ('run_pilot.py', ['--cache-dir', '/unused', '--solver', '/unused', '--output', existing])):
                result = subprocess.run([sys.executable, str(ROOT / 'pilot' / script), *args],
                                        capture_output=True, text=True)
                self.assertNotEqual(result.returncode, 0)
                self.assertIn('must not already exist', result.stderr)
                self.assertEqual(list(Path(existing).iterdir()), [])

    def test_optimized_python_cannot_disable_certificate_checks(self):
        result = subprocess.run([sys.executable, '-O', str(ROOT / 'pilot/run_pilot.py'), '--help'],
                                capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('certificate checks must remain enabled', result.stderr)


if __name__ == '__main__':
    unittest.main()
