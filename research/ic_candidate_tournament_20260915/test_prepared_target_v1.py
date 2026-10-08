"""Target adapter math/accounting controls; callbacks are not new SAT results."""
import copy
import json
from pathlib import Path
import unittest

from oracle import InvalidEvidence
from prepared_ic_state_v1 import accepted_files
from prepared_target_v1 import (CERTIFICATE_SEALS, audit_native_target, load_preparation, native_job,
                               preparation_exposures, solve_sat_target)
from generic_stages import DEFAULTS

HERE = Path(__file__).resolve().parent
POINT = [52411, 72106]


class Clock:
    def __init__(self):
        self.value = 0

    def __call__(self):
        self.value += 100
        return self.value


def failure(item, status='CONFLICT_BUDGET_INCONCLUSIVE'):
    return dict(trial=item['trial'], probe_scalar=item['probe_scalar'],
        public_point=item['point'], verification_wall_ns=1,
        source_model_valid=None, point_witness=None, status=status)


class PreparedTargetTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.docs = {family:json.loads((HERE/'goal_20260924/prepared-ic-state-v1'/
                        (family+'-preparation.json')).read_text()) for family in ('f5', 'sat')}

    def sat(self, query, *, seed=2026093012, cap=8, progress=None, point=POINT):
        return solve_sat_target(self.docs['sat'], CERTIFICATE_SEALS['sat'], point=point,
            panel=dict(descent_query_seed=seed, max_descent_queries=cap),
            query=query, clock=Clock(), progress=progress)

    def test_import_reconstructs_rank_logs_and_killed_torsion_without_query_execution(self):
        for family in self.docs:
            curve, base, matrix, logs, receipt = load_preparation(
                self.docs[family], CERTIFICATE_SEALS[family], family=family)
            self.assertEqual((len(base), matrix.rank, len(logs)), (63, 29, 29))
            self.assertIn((0,1), base)
            self.assertIsNone(curve.mul((0,1), curve.h))
            self.assertEqual(receipt['state_id'], 'ICP1hedbff76da644')

    def test_native_job_contains_only_the_mathematical_preparation_and_new_supplied_point(self):
        job = native_job(self.docs['f5'], CERTIFICATE_SEALS['f5'], point=POINT,
                         algorithm_seed=2026093032)
        self.assertEqual(job['public_targets'], [['52411','72106']])
        self.assertEqual(job['target_seeds'], [])
        self.assertEqual(len(job['prepared']['factor_base']), 63)
        self.assertEqual(len(job['prepared']['columns']), 29)
        self.assertNotIn('certificate', job['prepared'])
        self.assertNotIn('provenance', job['prepared'])
        self.assertNotIn('recovered_scalar', job)
        self.assertEqual(job['config']['max_trials'], 8)

    def test_family_swaps_and_wrong_external_certificate_seals_fail_before_queries(self):
        for document, seal in [(self.docs['f5'], CERTIFICATE_SEALS['f5']),
                               (self.docs['sat'], '0'*64)]:
            with self.assertRaises(InvalidEvidence):
                solve_sat_target(document, seal, point=POINT,
                    panel=dict(descent_query_seed=2026093012, max_descent_queries=8),
                    query=lambda *args:self.fail('query executed before rejected import'))

    def native_math_control(self):
        # Copy only for an offline auditor control. The accepted archive and
        # its consumed registration remain byte-for-byte unchanged.
        report = copy.deepcopy(accepted_files('f5')['execution/entry-output/pipeline.stdout'])
        job = native_job(self.docs['f5'], CERTIFICATE_SEALS['f5'], point=POINT,
                         algorithm_seed=2026093032)
        report.update(preparation_mode='imported-certified-log-table-v1',
            preparation_mathematical_state_sha256=job['prepared']['mathematical_state_sha256'],
            reusable_symbolic_template_prepared=True,
            trials=0, relations=[], collection_reports=[], solve_attempts=0,
            effective_config=dict(copy.deepcopy(DEFAULTS), **job['config']), scalar_verified=True)
        return report, job

    def test_native_audit_replays_real_failed_queries_and_relation_recovery_without_source_claim(self):
        report, job = self.native_math_control()
        result = audit_native_target(report, job, self.docs['f5'], CERTIFICATE_SEALS['f5'])
        self.assertEqual(result['recovered_scalar'], 24886)
        self.assertEqual(result['query_audit']['collection_queries'], 0)
        self.assertEqual(result['query_audit']['descent_queries'], 3)
        self.assertEqual(sum(result['online_phases_ns'].values()), result['online_wall_ns'])
        self.assertFalse(result['source_bound_execution_admitted'] or result['fresh_paired_qualification'])

    def test_native_audit_rejects_engine_substitution_timing_tampering_and_ordinary_queries(self):
        for change in ('engine', 'timing', 'ordinary', 'target'):
            report, job = self.native_math_control()
            if change == 'engine':
                report['solutions'][0]['attempts'][0]['pdp']['stats']['engine'] = {'MatrixF4':{'max_degree':3}}
            elif change == 'timing':
                report['online_wall_ns'] += 1
            elif change == 'ordinary':
                report['trials'] = 1
            else:
                report['fixture']['targets'] = [['471','57570']]
            with self.assertRaises(InvalidEvidence):
                audit_native_target(report, job, self.docs['f5'], CERTIFICATE_SEALS['f5'])

    def test_replay_of_the_retained_sat_witness_recovers_through_the_real_relation(self):
        report = accepted_files('sat')['execution/entry-output/summary.json']
        source_row = report['target_attempts'][0]
        calls = []

        def replay(panel, item, curve, base):
            calls.append(item)
            self.assertEqual(item['point'], source_row['public_point'])
            self.assertEqual(item['probe_scalar'], source_row['a'])
            row = copy.deepcopy(source_row)
            row.update(trial=item['trial'], verification_wall_ns=1)
            return row

        result = self.sat(replay)
        self.assertEqual(result['status'], 'COMPLETE')
        self.assertEqual(result['recovered_scalar'], 24886)
        self.assertEqual(len(calls), 1)
        self.assertEqual(sum(result['online_phases_ns'].values()), result['online_wall_ns'])
        self.assertFalse(result['source_bound_execution_admitted'])
        self.assertFalse(result['fresh_paired_qualification'] or result['promotion_eligible'])
        self.assertIsNone(result['online_speedup'])

    def test_failed_and_timed_out_attempts_are_charged_before_relation_recovery(self):
        calls, progress = [], []

        def controlled_query(panel, item, curve, base):
            calls.append(item)
            if item['trial'] < 2:
                return failure(item, 'TIMEOUT' if item['trial'] else 'CONFLICT_BUDGET_INCONCLUSIVE')
            row = failure(item)
            row.update(status='VALID_POINT_WITNESS', source_model_valid=True,
                       point_witness=dict(group_replay=True, point_indices=[0,4,50]))
            return row

        result = self.sat(controlled_query, seed=2026093032, cap=3, progress=progress.append)
        self.assertEqual(result['recovered_scalar'], 24886)
        self.assertEqual(len(result['target_attempts']), 3)
        self.assertEqual(len(progress), 2)
        self.assertEqual(result['target_status_mix'], {'CONFLICT_BUDGET_INCONCLUSIVE':1,
                          'TIMEOUT':1, 'VALID_POINT_WITNESS':1})
        self.assertEqual(sum(result['online_phases_ns'].values()), result['online_wall_ns'])
        self.assertGreater(result['online_phases_ns']['target_pdp'], 2*99)

    def test_exhausted_attempts_remain_incomplete_with_no_verified_online_time(self):
        rows = []
        result = self.sat(lambda panel,item,*args:failure(item), progress=rows.append)
        self.assertEqual(len(rows), 8)
        self.assertEqual(len(result['target_attempts']), 8)
        self.assertEqual(result['status'], 'INCOMPLETE_TARGET')
        self.assertFalse(result['scalar_verified'])
        self.assertIsNone(result['online_wall_ns'])
        self.assertIsNone(result['recovered_scalar'])
        self.assertGreater(result['online_attempt_wall_ns'], 0)
        self.assertEqual(sum(result['online_phases_ns'].values()), result['online_attempt_wall_ns'])

    def test_false_group_relation_and_missing_source_model_are_rejected(self):
        for model, indices in [(False,[0,4,50]), (True,[0,0,0]), (True,[True,4,50])]:
            def bad_query(panel,item,*args):
                row = failure(item)
                row.update(status='VALID_POINT_WITNESS', source_model_valid=model,
                           point_witness=dict(group_replay=True, point_indices=indices))
                return row
            with self.assertRaises(InvalidEvidence):
                self.sat(bad_query)

    def test_wrong_query_receipt_and_overlapping_clock_are_rejected(self):
        for changed in [dict(public_point=POINT), dict(trial=9), dict(verification_wall_ns=101)]:
            def bad_query(panel,item,*args):
                return dict(failure(item), **changed)
            with self.assertRaises(InvalidEvidence):
                self.sat(bad_query)

    def test_malformed_targets_and_caps_reject_before_query(self):
        for point in [[True,72106], ['52411','72106'], [0,1], [2**17,0], [1,2,3]]:
            with self.assertRaises(InvalidEvidence):
                self.sat(lambda *args:self.fail('malformed target dispatched'), point=point)
        for cap in [True,0,9,1.0]:
            with self.assertRaises(InvalidEvidence):
                self.sat(lambda *args:self.fail('invalid cap dispatched'), cap=cap)

    def test_preparation_exclusions_cover_all_queries_and_known_orbits(self):
        exposure = preparation_exposures(self.docs)
        points = {tuple(p) for p in exposure['record']['points']}
        for family, document in self.docs.items():
            curve, _, matrix, _, _ = load_preparation(document, CERTIFICATE_SEALS[family], family=family)
            for attempt in document['certificate']['inputs']['attempts']:
                self.assertIn(curve.mul(curve.g, attempt['scalar']), points)
            for column in matrix.columns:
                current = column
                for _ in range(curve.n):
                    self.assertIn(current, points)
                    self.assertIn(curve.neg(current), points)
                    current = curve.frob(current)
        self.assertFalse(exposure['fresh_sampling_authorized'])
        retained = json.loads((HERE/'goal_20260924/prepared-target-runtime-v1/preparation-exclusions.json').read_text())
        self.assertEqual(exposure,retained)
        self.assertEqual(exposure['point_count'],1340)
        with self.assertRaises(InvalidEvidence):
            preparation_exposures({'f5':self.docs['f5']})


if __name__ == '__main__':
    unittest.main()
