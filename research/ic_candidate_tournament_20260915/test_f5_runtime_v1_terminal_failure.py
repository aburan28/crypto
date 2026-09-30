"""Replay the closed F5 v1 schema failure without starting a native worker."""
import json
from pathlib import Path
import tempfile
import unittest

from audit_f5_runtime_v1 import audit
from oracle import InvalidEvidence
from publish_sat_runtime_failure_v3 import replay
from static_sat_native_v3 import audit_meter

BUNDLE = (Path(__file__).resolve().parent /
          'goal_20260924/f5-source-bound-runtime-v1/results-20260930')
EXPECTED = '5969837ddc69d25a3e44954f203077c0e6d682bbce2d1ac8eb0f9c3a10bb826b'


class F5RuntimeV1TerminalFailureTests(unittest.TestCase):
    def test_failure_transport_native_parse_and_admission_rejection(self):
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / 'transport'
            result = replay(BUNDLE, output, EXPECTED)
            self.assertEqual(result['status'], 'EXECUTION_FAILURE')
            self.assertTrue(result['python_source_audit']['complete_source_gates'])
            self.assertFalse(result['python_source_audit']['entrypoint_succeeded'])
            for key in ('mathematical_audit', 'verified_target_count', 'final_rank',
                        'online_wall_ns', 'online_phases_ns', 'online_speedup'):
                self.assertIsNone(result[key])
            self.assertFalse(result['complete_ic_source_bound'])
            self.assertFalse(result['headline_online_admissible'])
            self.assertFalse(result['promotion_eligible'])

            execution = output / 'execution'
            root = execution / 'entry-output'
            spec = json.loads((execution / 'execution.json').read_text())
            preflight = audit_meter(execution, root, 'build_identity',
                asset_role='bin/worker', arguments=['--build-identity'], seconds=10)
            native = audit_meter(execution, root, 'pipeline', asset_role='bin/worker',
                arguments=[], seconds=7170, stdin_argument='job')
            recorded = json.loads((BUNDLE / 'NATIVE-PARSE-REPLAY.json').read_text())
            self.assertEqual(preflight, recorded['preflight'])
            self.assertEqual(native, recorded['native'])
            self.assertEqual(native['returncode'], 2)
            self.assertFalse(native['timed_out'])
            self.assertEqual(spec['arguments']['job']['target_seeds'], [None])
            self.assertEqual(json.loads((root / 'pipeline.stdout').read_text()),
                             recorded['native_error'])
            self.assertEqual(recorded['native_error']['status'], 'error')
            self.assertFalse((root / 'summary.json').exists())
            with self.assertRaisesRegex(InvalidEvidence, 'failed F5 entrypoint'):
                audit(execution, spec)

            # These are alterations of an extracted evidence copy. No native
            # command, query stream, or historical registration is rerun.
            stdin = root / 'pipeline.stdin.json'
            original = stdin.read_bytes()
            self.assertIn(b'"target_seeds":[null]', original)
            stdin.chmod(0o644)
            stdin.write_bytes(original.replace(b'"target_seeds":[null]', b'"target_seeds":[]'))
            with self.assertRaisesRegex(InvalidEvidence, 'retained native stdin'):
                audit_meter(execution, root, 'pipeline', asset_role='bin/worker',
                    arguments=[], seconds=7170, stdin_argument='job')
            stdin.write_bytes(original)
            stderr = root / 'pipeline.stderr'
            stderr.write_bytes(b'changed')
            with self.assertRaisesRegex(InvalidEvidence, 'output bytes changed'):
                audit_meter(execution, root, 'pipeline', asset_role='bin/worker',
                    arguments=[], seconds=7170, stdin_argument='job')

    def test_external_invocation_seal_is_required_before_extraction(self):
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / 'wrong-seal'
            with self.assertRaisesRegex(InvalidEvidence, 'externally frozen'):
                replay(BUNDLE, output, '0' * 64)
            self.assertFalse(output.exists())


if __name__ == '__main__':
    unittest.main()
