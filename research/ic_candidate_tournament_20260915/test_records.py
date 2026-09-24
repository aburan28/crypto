"""Identity/accounting attacks against a real independent IC certificate.

The method digest and integer costs below are synthetic test data. They are
never exported as a measured candidate or a performance comparison.
"""
import copy
from pathlib import Path
import tempfile
import unittest

from identity import (STAGE_FIELDS, candidate_manifest, canonical, curve_record,
                      factor_base_inventory, run_id, sha256, workload_manifest,
                      write_immutable)
from measurement import (PHASES, PDP_OUTCOMES, exclusive_ledger, legacy_ledger,
                         measured_run, query_diagnostics, report_sha256)
from oracle import Curve, InvalidEvidence
from test_certificate import FIXTURE, REPORT


def method_fixture():
    method = {stage: {key: 'test-only' for key in names.split()}
              for stage, names in STAGE_FIELDS.items()}
    for stage in STAGE_FIELDS:
        method[stage]['source_sha256'] = '1'*64
    method['point_decomposition'].update(summands=3, solver='pair')
    method['relation_collection']['collector'] = 'sample'
    method['relation_linear_algebra'].update(solver='gauss', modulus=2003)
    method['target_descent']['method'] = 'pdp'
    method.update(isogeny='none', endomorphism={'order_conductor': None,
                  'frobenius_order_conductor': None, 'volcano_levels': []},
                  factor_base={'construction': {'test_only': 'retained-certificate'}, 'nominal_bound': 64},
                  implementation={'source_manifest_sha256': '1'*64,
                      'components': [{'role': 'test-only-worker', 'sha256': '2'*64}], 'flags': {}})
    return method


class IdentityTests(unittest.TestCase):
    def setUp(self):
        self.method = method_fixture()

    def test_target_seeds_do_not_change_curve_or_candidate(self):
        other = copy.deepcopy(FIXTURE)
        other['target_seeds'] = [4, 5]
        other['targets'].reverse()
        self.assertEqual(curve_record(FIXTURE), curve_record(other))
        self.assertEqual(candidate_manifest(FIXTURE, REPORT, self.method),
                         candidate_manifest(other, REPORT, self.method))

    def test_full_identity_binds_base_order_code_and_method(self):
        original = candidate_manifest(FIXTURE, REPORT, self.method)
        self.assertRegex(original['candidate_id'], r'^IC1N13Ckb0fb182PDP3pairRCsampleLAgaussTDpdpISO0h[0-9a-f]{12}$')
        reordered = copy.deepcopy(REPORT)
        reordered['factor_base'].reverse()
        # Set is identical, but relation indices now address another ordered base.
        self.assertNotEqual(original['candidate_id'], candidate_manifest(FIXTURE, reordered, self.method)['candidate_id'])
        for section, key, value in [('point_decomposition', 'solver', 'f4'),
                                    ('implementation', 'source_manifest_sha256', '3'*64),
                                    ('factor_base', 'nominal_bound', 65)]:
            changed = copy.deepcopy(self.method)
            changed[section][key] = value
            if key == 'source_manifest_sha256':
                for stage in STAGE_FIELDS:
                    changed[stage]['source_sha256'] = value
            self.assertNotEqual(original, candidate_manifest(FIXTURE, REPORT, changed))

    def test_unknown_stage_and_isogeny_are_not_results(self):
        for section, key, value in [('point_decomposition', 'encoding', None),
                                    ('point_decomposition', 'limits', {'max_degree': None}),
                                    ('relation_linear_algebra', 'modulus', 2)]:
            changed = copy.deepcopy(self.method)
            changed[section][key] = value
            with self.assertRaises(InvalidEvidence):
                candidate_manifest(FIXTURE, REPORT, changed)
        self.method['isogeny'] = {'unverified_neighbor': True}
        with self.assertRaises(InvalidEvidence):
            candidate_manifest(FIXTURE, REPORT, self.method)

    def test_no_floats_paths_or_run_seeds_in_method_hash(self):
        for key, value in [('run_seed', 2), ('source_path', '/tmp/source'),
                            ('threshold', .5), ('candidate_id', 'arbitrary')]:
            changed = copy.deepcopy(self.method)
            changed['implementation']['flags'][key] = value
            with self.subTest(key=key), self.assertRaises(InvalidEvidence):
                candidate_manifest(FIXTURE, REPORT, changed)
        changed = copy.deepcopy(self.method)
        changed['implementation']['flags']['threshold'] = '0.5'
        candidate_manifest(FIXTURE, REPORT, changed)

    def test_unknown_conductor_is_not_zero(self):
        self.method['endomorphism']['order_conductor'] = {'value': 0, 'proof_sha256': '1'*64}
        with self.assertRaises(InvalidEvidence):
            candidate_manifest(FIXTURE, REPORT, self.method)

    def test_base_count_is_not_nominal_bound_or_column_count(self):
        inventory = factor_base_inventory(REPORT, FIXTURE)
        self.assertEqual(inventory['usable_point_count'], 182)
        self.assertEqual(inventory['effective_columns'], 7)
        self.assertEqual(inventory['geometric_point_count'], 182)
        changed = copy.deepcopy(REPORT)
        changed['columns'] = 8
        with self.assertRaises(InvalidEvidence):
            factor_base_inventory(changed, FIXTURE)

    def test_cofactor_images_deduplicate_and_remove_torsion(self):
        c = Curve(FIXTURE)
        p = c.g
        torsion = (0, 1)
        q = c.add(p, torsion)
        self.assertIsNone(c.mul(torsion, c.h))
        self.assertEqual(c.mul(p, c.h), c.mul(q, c.h))
        inventory = factor_base_inventory({'factor_base': [list(p), list(q), list(torsion)],
                                           'columns': 1}, FIXTURE)
        self.assertEqual(inventory['geometric_point_count'], 3)
        self.assertEqual(inventory['usable_point_count'], 1)
        self.assertEqual(inventory['identity_images'], 1)
        self.assertEqual(inventory['duplicate_nonidentity_images'], 1)
        with self.assertRaises(InvalidEvidence):
            factor_base_inventory({'factor_base': [list(q)], 'columns': 1,
                                   'column_convention': 'representative'}, FIXTURE)

    def test_duplicate_or_identity_raw_points_rejected(self):
        for point in [None, REPORT['factor_base'][0]]:
            report = copy.deepcopy(REPORT)
            report['factor_base'].append(point)
            with self.assertRaises(InvalidEvidence):
                factor_base_inventory(report, FIXTURE)

    def test_workload_changes_independently_of_method(self):
        kwargs = dict(input_law='public-hash', algorithm_seed=9, resource_envelope={'memory_bytes': 1024})
        workload = workload_manifest(FIXTURE, **kwargs)
        changed = workload_manifest(FIXTURE, **dict(kwargs, algorithm_seed=10))
        self.assertNotEqual(workload['workload_id'], changed['workload_id'])
        candidate = candidate_manifest(FIXTURE, REPORT, self.method)
        self.assertEqual(run_id(candidate['candidate_id'], workload['workload_id'], 0),
                         candidate['candidate_id'] + 'W' + workload['workload_id'] + 'R0')
        with self.assertRaises(InvalidEvidence):
            run_id(candidate['candidate_id'], workload['workload_id'], True)

    def test_canonical_unicode_and_immutable_collision(self):
        self.assertEqual(canonical({'z': 1, 'a': 'λ'}), '{"a":"λ","z":1}'.encode())
        with tempfile.TemporaryDirectory() as d:
            path = Path(d)/'identity.json'
            write_immutable(path, {'value': 1})
            write_immutable(path, {'value': 1})
            with self.assertRaisesRegex(InvalidEvidence, 'immutable'):
                write_immutable(path, {'value': 2})
            self.assertEqual(path.read_bytes(), b'{"value":1}\n')


class AccountingTests(unittest.TestCase):
    def ledger(self):
        costs = {p: 10 for p in PHASES}
        costs['isogeny'] = 0
        return exclusive_ledger(costs, unit='synthetic-test-instructions', process_operations=100,
                                zero_reasons={'isogeny': 'no transport in this method'})

    def arguments(self):
        method = method_fixture()
        candidate = candidate_manifest(FIXTURE, REPORT, method)
        resources = {'memory_bytes': 1024}
        workload = workload_manifest(FIXTURE, input_law='public-hash', algorithm_seed=9,
                                     resource_envelope=resources)
        return dict(candidate=candidate, workload=workload, report=REPORT, fixture=FIXTURE,
                    method=method, number=0, ledger=self.ledger(), native_wall_ns=1000,
                    status='complete', provenance={'source_manifest_sha256': '1'*64,
                    'worker_sha256': '2'*64, 'host_id': 'synthetic', 'calibration_id': 'synthetic',
                    'resource_envelope_id': sha256(resources), 'report_sha256': report_sha256(REPORT)})

    def test_complete_cost_requires_closure_but_does_not_promote(self):
        result = measured_run(**self.arguments())
        self.assertEqual(result['total_operations'], 100)
        self.assertEqual(result['certificate']['verified_targets'], 2)
        self.assertFalse(result['promotion_eligible'])

    def test_unknown_cost_does_not_turn_into_zero_or_total(self):
        kwargs = self.arguments()
        costs = dict(kwargs['ledger']['operations'], matrix_build=None)
        kwargs['ledger'] = exclusive_ledger(costs, unit='test', process_operations=100,
                                            zero_reasons={'isogeny': 'absent'})
        result = measured_run(**kwargs)
        self.assertIsNone(result['total_operations'])
        self.assertEqual(result['certificate']['verified_targets'], 2)
        self.assertEqual(result['phase_ledger']['missing_phases'], ['matrix_build'])

    def test_legacy_intervals_are_diagnostics_only(self):
        ledger = legacy_ledger({'factor_base_and_tables': 60, 'everything_else': 40},
                               unit='test', process_operations=100)
        self.assertIsNone(ledger['cold_operations'])
        self.assertTrue(all(v is None for v in ledger['operations'].values()))
        self.assertEqual(ledger['whole_process_operations'], 100)

    def test_cost_gap_overcount_bool_negative_and_unexplained_zero_fail(self):
        for value in (11, -1, True, 0):
            costs = dict(self.ledger()['operations'], matrix_build=value)
            with self.subTest(value=value), self.assertRaises(InvalidEvidence):
                exclusive_ledger(costs, unit='test', process_operations=100,
                                 zero_reasons={'isogeny': 'absent'})
        with self.assertRaises(InvalidEvidence):
            exclusive_ledger(self.ledger()['operations'], unit='test', process_operations=101,
                             zero_reasons={'isogeny': 'absent'})

    def test_failures_keep_keys_and_null_totals(self):
        for status in ('timeout', 'oom', 'budget', 'error', 'insufficient_relations'):
            kwargs = self.arguments()
            kwargs.update(status=status, report=None, admission_report=REPORT)
            kwargs['provenance']['report_sha256'] = None
            result = measured_run(**kwargs)
            self.assertIsNone(result['total_operations'])
            self.assertIsNone(result['certificate'])
            self.assertEqual(result['status'], status)
            self.assertTrue(result['run_id'].endswith('R0'))

    def test_forged_labels_and_source_are_rejected_even_on_failure(self):
        for changed in ('candidate_id', 'record_sha256'):
            kwargs = self.arguments()
            kwargs['candidate'][changed] = 'forged'
            kwargs.update(status='timeout', report=None, admission_report=REPORT)
            kwargs['provenance']['report_sha256'] = None
            with self.assertRaises(InvalidEvidence):
                measured_run(**kwargs)
        kwargs = self.arguments()
        kwargs['provenance']['source_manifest_sha256'] = '4'*64
        with self.assertRaises(InvalidEvidence):
            measured_run(**kwargs)

    def test_wrong_scalar_and_changed_report_cannot_claim_cost(self):
        kwargs = self.arguments()
        kwargs['report'] = copy.deepcopy(REPORT)
        kwargs['report']['solutions'][0]['recovered'] = '0'
        kwargs['provenance']['report_sha256'] = report_sha256(kwargs['report'])
        with self.assertRaises(InvalidEvidence):
            measured_run(**kwargs)
        kwargs = self.arguments()
        kwargs['provenance']['report_sha256'] = '0'*64
        with self.assertRaises(InvalidEvidence):
            measured_run(**kwargs)

    def test_warm_workload_never_exports_a_cold_total(self):
        kwargs = self.arguments()
        kwargs['workload'] = workload_manifest(FIXTURE, input_law='public-hash', algorithm_seed=9,
            resource_envelope={'memory_bytes': 1024}, cache_policy='warm')
        result = measured_run(**kwargs)
        self.assertIsNone(result['total_operations'])
        self.assertFalse(result['complete_cold_cost'])

    def test_mismatched_query_diagnostics_rejected(self):
        kwargs = self.arguments()
        outcomes = dict.fromkeys(PDP_OUTCOMES, 0)
        outcomes['verified'] = 7
        kwargs['diagnostics'] = query_diagnostics(attempts=7, outcomes=outcomes, ordinary_queries=7,
            verified_relations=7, novel_rows=7, final_rank=7, effective_columns=7)
        with self.assertRaises(InvalidEvidence):
            measured_run(**kwargs)

    def test_zero_yield_retained_and_unresolved_is_not_unsat(self):
        outcomes = dict.fromkeys(PDP_OUTCOMES, 0)
        outcomes['unresolved'] = 20
        result = query_diagnostics(attempts=20, outcomes=outcomes, ordinary_queries=20,
                                   verified_relations=0, novel_rows=0, final_rank=0, effective_columns=7)
        self.assertEqual(result['outcomes']['proved_unsat'], 0)
        self.assertEqual(result['outcomes']['unresolved'], 20)
        outcomes['timeout'] = 1
        with self.assertRaises(InvalidEvidence):
            query_diagnostics(attempts=20, outcomes=outcomes, ordinary_queries=20,
                              verified_relations=0, novel_rows=0, final_rank=0, effective_columns=7)


if __name__ == '__main__':
    unittest.main()
