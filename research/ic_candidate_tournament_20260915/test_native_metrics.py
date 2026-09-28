import json
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

    @unittest.skipUnless(hasattr(os, 'sched_getaffinity'), 'Linux instruction-runner capability')
    def test_watchdog_still_retains_timeout(self):
        with tempfile.TemporaryDirectory() as tmp:
            result = execute([sys.executable, '-c', 'import time; time.sleep(10)'], {},
                             Path(tmp), .05, 1024**3, min(os.sched_getaffinity(0)))
            self.assertEqual(result['process_status'], 'TIMEOUT')
            self.assertNotEqual(result['exit_code'], 0)
            self.assertLess(result['process_wall_seconds'], 3)

    @unittest.skipUnless(hasattr(os, 'sched_getaffinity'), 'Linux instruction-runner capability')
    def test_execute_caps_and_pins_the_child_without_a_preexec_fork(self):
        from tournament import SPAWN
        probe = ('import json,os,resource,sys; sys.stdin.read(); '
                 'print(json.dumps({"as":resource.getrlimit(resource.RLIMIT_AS)[0],'
                 '"core":resource.getrlimit(resource.RLIMIT_CORE)[0],'
                 '"cpus":sorted(os.sched_getaffinity(0)),"sid":os.getsid(0)==os.getpid()}))')
        cpu = min(os.sched_getaffinity(0))
        before = os.sched_getaffinity(0)
        with tempfile.TemporaryDirectory() as tmp:
            result = execute([sys.executable, '-c', probe], {'k': 1}, Path(tmp), 30, 4 * 1024**3, cpu)
            seen = json.loads(Path(tmp, 'stdout.json').read_text())
        self.assertEqual(result['process_status'], 'EXITED')
        self.assertEqual(result['exit_code'], 0)
        self.assertEqual(result['spawn'], SPAWN)
        self.assertGreater(result['peak_rss_bytes'], 0)
        self.assertEqual(result['process_wall_seconds'], result['process_wall_ns']/1_000_000_000)
        self.assertEqual(seen, {'as': 4 * 1024**3, 'core': 0, 'cpus': [cpu], 'sid': True})
        self.assertEqual(os.sched_getaffinity(0), before)

    def test_built_worker_accepts_host_or_triple_layout(self):
        from oracle import InvalidEvidence
        from tournament import built_worker
        with tempfile.TemporaryDirectory() as tmp:
            build = Path(tmp)
            with self.assertRaises(InvalidEvidence):
                built_worker(build)
            triple = build / 'x86_64-unknown-linux-gnu/release/examples'
            triple.mkdir(parents=True)
            (triple / 'ic_tournament_worker').write_bytes(b'1')
            self.assertEqual(built_worker(build), triple / 'ic_tournament_worker')
            host = build / 'release/examples'
            host.mkdir(parents=True)
            (host / 'ic_tournament_worker').write_bytes(b'2')
            with self.assertRaises(InvalidEvidence):
                built_worker(build)
if __name__ == '__main__':
    unittest.main()
