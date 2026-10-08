"""Raw diagnostic serialization must neither drop data nor invent IC costs."""
import copy
import json
from pathlib import Path
import tempfile
import unittest

from f5_runtime_receipts_v2 import publication_process, SECONDS
from identity import canonical, write_immutable
from oracle import InvalidEvidence

HERE = Path(__file__).resolve().parent


class F5RuntimeReceiptsV2Tests(unittest.TestCase):
    def process(self):
        # A retained actual receipt, read only. No worker runs in this test.
        record = json.loads((HERE/'goal_20260924/f5-source-bound-runtime-v1/'
            'results-20260930/NATIVE-PARSE-REPLAY.json').read_text())
        return record['native']

    def test_actual_raw_float_metrics_publish_and_preserve_exact_clock(self):
        process = self.process()
        original = copy.deepcopy(process)
        result = publication_process(process)
        self.assertEqual(process, original)
        self.assertEqual(result['native_wall_ns'], 6275834)
        self.assertEqual(result['returncode'], 2)
        self.assertNotIn('metrics', result)
        diagnostics = result['resource_diagnostics']
        self.assertEqual(diagnostics['peak_rss_bytes'], process['metrics']['peak_rss_bytes'])
        for key in SECONDS:
            self.assertEqual(float(diagnostics['seconds_decimal'][key]), process['metrics'][key])
        self.assertEqual(json.loads(canonical(result)), result)
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary)/'diagnostics.json'
            write_immutable(path, result)
            self.assertEqual(json.loads(path.read_text()), result)
        # Publication representation cannot add a mathematical success or cost.
        for key in ('online_wall_ns', 'online_speedup', 'complete_ic_admitted'):
            self.assertNotIn(key, result)

    def test_nonfinite_negative_boolean_and_unknown_diagnostics_fail_closed(self):
        for key, value in [('wall_seconds', float('nan')), ('wall_seconds', float('inf')),
                           ('user_seconds', -1.0), ('system_seconds', True),
                           ('peak_rss_bytes', 1.5), ('peak_rss_bytes', False),
                           ('unknown_seconds', 0.0)]:
            with self.subTest(key=key, value=value):
                process = self.process()
                process['metrics'][key] = value
                with self.assertRaises(InvalidEvidence):
                    publication_process(process)


if __name__ == '__main__':
    unittest.main()
