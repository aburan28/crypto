"""Real synthetic failures transport losslessly and cannot become IC wins."""
import json
from pathlib import Path
import shutil
import tempfile
import unittest

from identity import sha256
from oracle import InvalidEvidence
from publish_sat_runtime_failure_v3 import assess, publish, replay
from sat_runtime_bundle import CHILD_SCRIPTS, DIRECTORY
from sat_runtime_execution_v3 import execute, register

HERE = Path(__file__).resolve().parent


class PublishSatRuntimeFailureV3Tests(unittest.TestCase):
    def execution(self, temporary, body, seconds=30):
        root = Path(temporary)
        source = root/'repository'
        directory = source/DIRECTORY
        (directory/'producer').mkdir(parents=True)
        for name in ('identity.py', 'oracle.py', 'sat_runtime_bundle.py', 'sat_runtime_execution_v3.py'):
            shutil.copyfile(HERE/name, directory/name)
        for name in ('evidence.py', 'timing.py'):
            (directory/'producer'/name).write_text('CONTROL = True\n')
        (directory/'control.py').write_text('def run(arguments, output):\n'+body)
        (source/'scripts').mkdir()
        for role in CHILD_SCRIPTS:
            (source/role).write_text('CONTROL = True\n')
        registration, execution = root/'registration', root/'execution'
        spec = register(source, registration, module='control', action='run', arguments={},
                        timeout_seconds=seconds)
        process = execute(registration, execution, expected_spec=spec, timeout_seconds=seconds)
        return root, registration, execution, spec, process

    def test_error_retains_truncated_progress_and_complete_python_failure_gates(self):
        with tempfile.TemporaryDirectory() as temporary:
            root, registration, execution, spec, process = self.execution(temporary,
                '    (output/"collection.progress.jsonl").write_bytes(b\'{"trial":0}\\n{"trial":\')\n'
                '    raise RuntimeError("synthetic registered failure")\n')
            self.assertEqual(process['exit_code'], 1)
            result = publish(registration, execution, root/'bundle', sha256(spec))
            self.assertEqual(result['status'], 'EXECUTION_FAILURE')
            self.assertTrue(result['python_source_audit']['complete_source_gates'])
            self.assertFalse(result['python_source_audit']['entrypoint_succeeded'])
            self.assertEqual(replay(root/'bundle', root/'transport', sha256(spec)), result)
            self.assertEqual((root/'transport/execution/entry-output/collection.progress.jsonl').read_bytes(),
                             (execution/'entry-output/collection.progress.jsonl').read_bytes())
            self.assertIsNone(result['final_rank'])
            self.assertIsNone(result['verified_target_count'])
            self.assertIsNone(result['online_wall_ns'])
            self.assertFalse(result['complete_ic_source_bound'])
            self.assertFalse(result['headline_online_admissible'])
            self.assertFalse(result['promotion_eligible'])
            self.assertIsNone(result['online_speedup'])
            with self.assertRaisesRegex(InvalidEvidence, 'never overwrite'):
                publish(registration, execution, root/'bundle', sha256(spec))

    def test_actual_watchdog_without_terminal_gate_is_transportable_but_incomplete(self):
        with tempfile.TemporaryDirectory() as temporary:
            root, registration, execution, spec, process = self.execution(
                temporary, '    import time\n    time.sleep(20)\n', seconds=1)
            self.assertTrue(process['timed_out'])
            self.assertFalse((execution/'after.json').exists())
            result = publish(registration, execution, root/'bundle', sha256(spec))
            self.assertEqual(result['status'], 'CONTROLLER_TIMEOUT')
            self.assertIsNone(result['python_source_audit'])
            self.assertIsNone(result['mathematical_audit'])
            self.assertEqual(replay(root/'bundle', root/'transport', sha256(spec)), result)
            self.assertEqual((root/'transport/execution/runtime/runtime.tar.gz').read_bytes(),
                             (execution/'runtime/runtime.tar.gz').read_bytes())

    def test_missing_process_receipt_cannot_infer_timeout_or_zero_cost(self):
        with tempfile.TemporaryDirectory() as temporary:
            root, registration, execution, spec, _ = self.execution(
                temporary, '    raise RuntimeError("synthetic registered failure")\n')
            # Artifact fault injection, not another measured invocation.
            (execution/'process.json').unlink()
            result = publish(registration, execution, root/'bundle', sha256(spec))
            self.assertEqual(result['status'], 'PROCESS_RECORD_MISSING')
            self.assertIsNone(result['process'])
            self.assertIsNone(result['online_wall_ns'])
            self.assertEqual(replay(root/'bundle', root/'transport', sha256(spec)), result)

    def test_wrong_external_seal_or_changed_archive_fails_before_extraction(self):
        with tempfile.TemporaryDirectory() as temporary:
            root, registration, execution, spec, _ = self.execution(
                temporary, '    raise RuntimeError("synthetic registered failure")\n')
            with self.assertRaisesRegex(InvalidEvidence, 'externally frozen'):
                publish(registration, execution, root/'wrong-publication', '0'*64)
            self.assertFalse((root/'wrong-publication').exists())
            publish(registration, execution, root/'bundle', sha256(spec))
            with self.assertRaisesRegex(InvalidEvidence, 'externally frozen'):
                replay(root/'bundle', root/'wrong-transport', '0'*64)
            self.assertFalse((root/'wrong-transport').exists())
            with (root/'bundle/evidence.tar.gz').open('ab') as stream:
                stream.write(b'changed')
            with self.assertRaisesRegex(InvalidEvidence, 'archive changed'):
                replay(root/'bundle', root/'changed-transport', sha256(spec))
            self.assertFalse((root/'changed-transport').exists())

    def test_success_or_relabelled_failed_process_is_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            root, registration, execution, spec, process = self.execution(temporary, '    return {}\n')
            self.assertEqual(process['exit_code'], 0)
            with self.assertRaisesRegex(InvalidEvidence, 'successful execution'):
                publish(registration, execution, root/'bundle', sha256(spec))
        with tempfile.TemporaryDirectory() as temporary:
            root, registration, execution, spec, process = self.execution(
                temporary, '    raise RuntimeError("synthetic registered failure")\n')
            original = (execution/'process.json').read_text()
            process['watchdog_seconds'] += 1
            (execution/'process.json').write_text(json.dumps(process))
            with self.assertRaisesRegex(InvalidEvidence, 'process receipt changed'):
                assess(registration, execution, sha256(spec))
            (execution/'process.json').write_text(original)
            with (execution/'stderr.txt').open('ab') as stream:
                stream.write(b'changed')
            with self.assertRaisesRegex(InvalidEvidence, 'process receipt changed'):
                assess(registration, execution, sha256(spec))


if __name__ == '__main__':
    unittest.main()
