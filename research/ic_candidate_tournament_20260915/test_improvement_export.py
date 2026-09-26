"""Reporting guards: target weighting, unknown rates and retained failures."""
import copy
import math
import unittest

from goal_20260924.improvement.export_report import aggregate, rate_interval, stage_diagnostics


class ImprovementExportTests(unittest.TestCase):
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
