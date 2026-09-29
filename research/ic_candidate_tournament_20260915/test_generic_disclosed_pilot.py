"""Controls for the disclosed-input diagnostic pilot's frozen schedule."""
import json
from pathlib import Path
import tempfile
import unittest

from run_generic_backend_disclosed_pilot import (PANEL, PANEL_SHA256, digest,
                                                  fixture_map, run_bounded_worker,
                                                  static_preflight)
from tournament import read


class DisclosedPilotTests(unittest.TestCase):
    def test_panel_uses_only_five_disclosed_points_and_feasible_layouts(self):
        self.assertEqual(digest(PANEL), PANEL_SHA256)
        panel = read(PANEL)
        fixtures = fixture_map(panel)
        self.assertEqual(set(fixtures), {'n17a1', 'n19a0', 'n23a0', 'n23a1', 'n31a0'})
        self.assertTrue(all(len(fixture['targets']) == 1 for fixture in fixtures.values()))
        layout = static_preflight(panel)
        self.assertEqual(layout['status'], 'PASS_STATIC_LAYOUT_ONLY')
        self.assertEqual(len(layout['rows']), 10)

    def test_rss_watch_preserves_raw_output_and_timeout(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            worker = root/'worker.py'
            worker.write_text('#!/usr/bin/env python3\n'
                              'import json,sys,time\n'
                              'job=json.load(sys.stdin)\n'
                              'time.sleep(job["sleep"])\n'
                              'print(json.dumps(job))\n')
            worker.chmod(0o755)
            panel = dict(timeout_seconds=2, memory_bytes=8*1024**3,
                         memory_poll_ms=10)
            success = root/'success'
            success.mkdir()
            stdout, stderr, code, disposition, wall_ns, peak = run_bounded_worker(
                worker, {'sleep': 0.02}, success, panel)
            self.assertEqual((code, disposition, stderr), (0, 'EXITED', ''))
            self.assertEqual(json.loads(stdout), {'sleep': 0.02})
            self.assertEqual(stdout, (success/'stdout.json').read_text())
            self.assertGreater(wall_ns, 0)
            self.assertGreaterEqual(peak, 0)
            timeout = root/'timeout'
            timeout.mkdir()
            panel['timeout_seconds'] = 0.05
            _, _, code, disposition, _, _ = run_bounded_worker(
                worker, {'sleep': 1}, timeout, panel)
            self.assertEqual(disposition, 'TIMEOUT')
            self.assertNotEqual(code, 0)
            self.assertTrue((timeout/'stdout.json').exists())
            self.assertTrue((timeout/'stderr.txt').exists())


if __name__ == '__main__':
    unittest.main()
