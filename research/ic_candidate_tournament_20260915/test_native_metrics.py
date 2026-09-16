import tempfile
import unittest
import os
import sys
from pathlib import Path

from tournament import comparison, execute, gate, rho_gate
from test_tournament import rows


class NativeMetricsTests(unittest.TestCase):
    def test_instruction_gain_without_time_gain_does_not_promote(self):
        c = {'confirmation_ratio': .8, 'max_cell_ratio': 1.1, 'require_native_progress': True}
        result = comparison(rows(.4), 'candidate', draws=200)
        self.assertFalse(gate(result, c))
        self.assertEqual(result['native_wall_ci95'], [1, 1])

    def test_two_metric_gain_passes_and_time_regression_blocks(self):
        c = {'confirmation_ratio': .8, 'max_cell_ratio': 1.1, 'require_native_progress': True}
        data = rows(.4)
        for row in data:
            if row['arm'] == 'candidate':
                row['native_process']['process_wall_seconds'] = .04
        self.assertTrue(gate(comparison(data, 'candidate', draws=200), c))
        for row in data:
            if row['arm'] == 'candidate' and row['cell'] == '0':
                row['native_process']['process_wall_seconds'] = .12
        self.assertFalse(gate(comparison(data, 'candidate', draws=200), c))

    def test_rho_objective_blocks_aa_and_requires_strict_rho_win(self):
        c = {'confirmation_ratio': .8, 'max_cell_ratio': 1.1, 'require_native_progress': True,
             'objective': 'rho', 'no_regression_ratio': .98}
        # An A/A control has ratio one: never a promotion under either objective.
        self.assertFalse(gate(comparison(rows(), 'candidate', draws=200), c))
        # A small but measurable gain over the incumbent passes the no-regression gate
        # under the rho objective and fails the default 20% gate.
        data = rows(.9)
        for row in data:
            if row['arm'] == 'candidate':
                row['native_process']['process_wall_seconds'] = .095
        result = comparison(data, 'candidate', draws=200)
        self.assertTrue(gate(result, c))
        self.assertFalse(gate(result, dict(c, objective='incumbent')))
        # The rho gate needs every upper limit and every cell strictly below one.
        self.assertTrue(rho_gate(result))
        for row in data:
            if row['arm'] == 'candidate' and row['cell'] == '0':
                row['total_operations'] = 1000
        self.assertFalse(rho_gate(comparison(data, 'candidate', draws=200)))
        self.assertFalse(rho_gate({'eligible': False}))

    def test_watchdog_still_retains_timeout(self):
        with tempfile.TemporaryDirectory() as tmp:
            result = execute([sys.executable, '-c', 'import time; time.sleep(10)'], {},
                             Path(tmp), .05, 1024**3, min(os.sched_getaffinity(0)))
            self.assertEqual(result['process_status'], 'TIMEOUT')
            self.assertNotEqual(result['exit_code'], 0)
            self.assertLess(result['process_wall_seconds'], 3)


if __name__ == '__main__':
    unittest.main()
