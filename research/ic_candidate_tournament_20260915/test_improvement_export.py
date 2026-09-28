"""Reporting guards: target weighting, unknown rates and retained failures."""
import copy
import hashlib
import json
import math
from pathlib import Path
import unittest

from goal_20260924.improvement.export_report import aggregate, rate_interval, stage_diagnostics, reference_ratios


class ImprovementExportTests(unittest.TestCase):
    def test_correction_receipt_preserves_the_exact_previous_export(self):
        root = Path(__file__).parent/'goal_20260924/improvement/round2'
        raw = (root/'RESULTS.json').read_bytes()
        receipt = json.loads((root/'REPORTING-CORRECTION.json').read_text())
        self.assertEqual(hashlib.sha256(raw).hexdigest(), receipt['corrected_RESULTS_sha256'])
        previous = json.loads(raw)
        for change in receipt['rows']:
            row = next(r for r in previous['stages'][change['stage']]['table'] if r['alias'] == change['alias'])
            for key, values in change['fields'].items():
                self.assertEqual(row[key], values['after'])
                if key == 'online_ic_reference':
                    del row[key]
                else:
                    row[key] = values['before']
        reconstructed = (json.dumps(previous, indent=2, sort_keys=True, allow_nan=False)+'\n').encode()
        self.assertEqual(hashlib.sha256(reconstructed).hexdigest(), receipt['original_RESULTS_sha256'])

    def test_rho_speedup_does_not_divide_ratios_to_different_ic_references(self):
        # incumbent=100, ic_online=50, challenger=40, rho_online=200.
        # (rho/incumbent)/(challenger/ic_online) is 2.5, but rho/challenger is 5.
        rows = [dict(alias='incumbent', online_ms=100, online_over_ic=1,
                     cold_Ir_over_ic=1, complete=True, comparison=None),
                dict(alias='candidate', online_ms=40, online_over_ic=.8,
                     cold_Ir_over_ic=.7, complete=True,
                     comparison=dict(metric_references=dict(online_ns='ic_online'))),
                dict(alias='rho_online', online_ms=200, online_over_ic=2,
                     cold_Ir_over_ic=3, complete=True, comparison=None)]
        reference_ratios(rows, versioned=True)
        self.assertEqual(rows[1]['rho_online_over_IC_online'], 5)
        self.assertEqual(rows[1]['online_ic_reference'], 'ic_online')
        self.assertEqual(rows[0]['online_ic_reference'], 'incumbent')
        self.assertEqual(rows[2]['rho_online_over_IC_online'], 1)
        rows[2]['complete'] = False
        reference_ratios(rows, versioned=True)
        self.assertIsNone(rows[1]['rho_online_over_IC_online'])

    def test_legacy_common_denominator_export_is_unchanged(self):
        rows = [dict(alias='candidate', online_over_ic=.5, cold_Ir_over_ic=.75),
                dict(alias='rho_online', online_over_ic=2, cold_Ir_over_ic=3),
                dict(alias='rho', online_over_ic=3, cold_Ir_over_ic=2)]
        reference_ratios(rows, versioned=False)
        self.assertEqual(rows[0]['rho_online_over_IC_online'], 4)
        self.assertEqual(rows[0]['cold_Ir_over_cold_rho'], .375)
        self.assertNotIn('online_ic_reference', rows[0])

    def test_process_repetitions_are_not_targets_and_cells_have_equal_weight(self):
        rows = []
        for cell, case, values in (('a', 'a1', [1, 1, 100]), ('a', 'a2', [4, 4, 100]),
                                   ('b', 'b1', [16, 16, 1600])):
            rows.extend(dict(cell=cell, case=case, value=v) for v in values)
        self.assertAlmostEqual(aggregate(rows, lambda r:r['value']), math.sqrt(2*16))
        rows[-1]['value'] = None
        self.assertIsNone(aggregate(rows, lambda r:r['value']))

    def test_zero_yield_is_zero_but_an_absent_denominator_is_unknown(self):
        points = [dict(success=0, attempts=3), dict(success=0, attempts=7)]
        result = rate_interval(points, 'success', 'attempts', 7)
        self.assertEqual(result['value'], 0)
        self.assertEqual(result['ci95'], [0, 0])
        self.assertEqual(result['targets'], 2)
        self.assertIsNone(rate_interval([dict(success=0, attempts=0)], 'success', 'attempts', 7)['value'])
        self.assertIsNone(rate_interval(points[:1], 'success', 'attempts', 7)['ci95'])

    @staticmethod
    def rows():
        phases = dict(setup=11, isogeny=0, factor_base=2, precompute=3, queries=4, pdp=5,
            relation_check=6, matrix_build=7, relation_la=8, target_descent=9, recovery_check=10)
        rows = []
        for case in ('c1','c2'):
            for rep in range(3):
                rows.append(dict(cell='c', case=case, repetition=rep, status='VERIFIED',
                    total_operations=sum(phases.values()), phase_costs=phases,
                    native_process={'peak_rss_bytes':1024},
                    measurement=dict(phase_ledger={'operations':phases},
                        diagnostics=dict(attempts=7, verified_relations=4, novel_rows=3,
                                         final_rank=3, outcomes={'verified':4,'unresolved':3}))))
        return rows

    def test_exclusive_means_close_and_one_failed_process_blocks_cell_rates(self):
        rows = self.rows(); cell = stage_diagnostics(rows)['c']
        self.assertEqual(sum(cell['arithmetic_mean_phase_Ir'].values()), cell['arithmetic_mean_complete_cold_Ir'])
        self.assertEqual(cell['novel_rows_per_attempt']['targets'], 2)
        self.assertAlmostEqual(cell['collection_Ir_per_novel_row']['value'], (4+5+6+7+8)/3)
        self.assertIsNone(cell['base_only_peak_bytes'])
        rows[0]['status']='TIMEOUT'
        cell = stage_diagnostics(rows)['c']
        self.assertFalse(cell['complete'])
        self.assertIsNone(cell['verified_relations_per_attempt'])
        self.assertIsNone(cell['arithmetic_mean_complete_cold_Ir'])

    def test_repetitions_cannot_silently_become_new_walk_samples(self):
        rows=copy.deepcopy(self.rows())
        rows[0]['measurement']['diagnostics']['attempts'] += 1
        with self.assertRaisesRegex(ValueError, 'deterministic walk'):
            stage_diagnostics(rows)


if __name__ == '__main__':
    unittest.main()
