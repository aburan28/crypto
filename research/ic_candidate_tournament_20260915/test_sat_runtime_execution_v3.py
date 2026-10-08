"""Real isolated imports and incomplete executions retain their source boundary."""
import copy
from concurrent.futures import ThreadPoolExecutor
import json
import os
from pathlib import Path
import shutil
import tempfile
import threading
import unittest
from unittest.mock import patch

from oracle import InvalidEvidence
from sat_runtime_bundle import CHILD_SCRIPTS, DIRECTORY
from sat_runtime_execution_v3 import (
    EXECUTION_POLICY, audit_execution, claim_execution, execute, register,
)

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


class SatRuntimeExecutionV3Tests(unittest.TestCase):
    def fixture(self, temporary, contents):
        root = Path(temporary)/'repository'
        directory = root/DIRECTORY
        (directory/'producer').mkdir(parents=True)
        for name in ('identity.py', 'oracle.py', 'sat_runtime_bundle.py',
                     'sat_runtime_execution_v3.py'):
            shutil.copyfile(HERE/name, directory/name)
        (directory/'producer/evidence.py').write_text('VALUE = 41\n')
        (directory/'producer/timing.py').write_text('VALUE = 1\n')
        (directory/'control.py').write_text(contents)
        (root/'scripts').mkdir()
        for role in CHILD_SCRIPTS:
            (root/role).write_text('CONTROL = True\n')
        return root, Path(temporary)/'registration', Path(temporary)/'execution'

    def test_actual_sat_imports_pass_in_fresh_isolated_snapshot(self):
        with tempfile.TemporaryDirectory() as temporary:
            registration = Path(temporary)/'registration'
            output = Path(temporary)/'execution'
            spec = register(ROOT, registration, module='sat_runtime_execution_v3',
                            action='import_probe', arguments=['run_static_sat_full_v2'],
                            timeout_seconds=60)
            process = execute(registration, output, expected_spec=spec, timeout_seconds=60)
            self.assertEqual(process['exit_code'], 0,
                             (output/'stderr.txt').read_text())
            audited = audit_execution(output, spec)
            self.assertTrue(audited['complete_source_gates'])
            self.assertTrue(audited['entrypoint_succeeded'])
            self.assertFalse(audited['promotion_eligible'])
            terminal = json.loads((output/'after.json').read_text())
            self.assertEqual(terminal['loaded_modules']['producer.evidence'],
                             (DIRECTORY/'producer/evidence.py').as_posix())
            self.assertEqual(terminal['loaded_modules']['producer.timing'],
                             (DIRECTORY/'producer/timing.py').as_posix())
            self.assertEqual(terminal['result']['measured_solver_executed'], False)

    def test_execution_uses_frozen_package_after_live_source_changes(self):
        with tempfile.TemporaryDirectory() as temporary:
            root, registration, output = self.fixture(temporary, '''def run(arguments, output):
    from producer.evidence import VALUE
    from producer.timing import VALUE as SECOND
    return {"value": VALUE+SECOND}
''')
            spec = register(root, registration, module='control', action='run', arguments={},
                            timeout_seconds=30)
            (root/DIRECTORY/'producer/evidence.py').write_text('VALUE = 10000\n')
            process = execute(registration, output, expected_spec=spec, timeout_seconds=30)
            self.assertEqual(process['exit_code'], 0,
                             (output/'stderr.txt').read_text())
            self.assertEqual(json.loads((output/'after.json').read_text())['result'],
                             {'value': 42})
            self.assertTrue(audit_execution(output, spec)['entrypoint_succeeded'])
            with self.assertRaisesRegex(InvalidEvidence, 'no retries'):
                execute(registration, output, expected_spec=spec, timeout_seconds=30)
            with self.assertRaisesRegex(InvalidEvidence, 'registration consumed'):
                execute(registration, Path(temporary)/'second-execution',
                        expected_spec=spec, timeout_seconds=30)
            with self.assertRaisesRegex(InvalidEvidence, 'already exists'):
                register(root, registration, module='control', action='run', arguments={},
                         timeout_seconds=30)

    def test_thrown_entrypoint_remains_failed_with_complete_source_gates(self):
        with tempfile.TemporaryDirectory() as temporary:
            root, registration, output = self.fixture(temporary,
                'def run(arguments, output):\n    raise RuntimeError("registered failure")\n')
            spec = register(root, registration, module='control', action='run', arguments={},
                            timeout_seconds=30)
            process = execute(registration, output, expected_spec=spec, timeout_seconds=30)
            self.assertEqual(process['exit_code'], 1)
            audited = audit_execution(output, spec)
            self.assertTrue(audited['complete_source_gates'])
            self.assertFalse(audited['entrypoint_succeeded'])
            self.assertIsNone(audited['online_speedup'])
            self.assertIn('registered failure', (output/'stderr.txt').read_text())
            with self.assertRaisesRegex(InvalidEvidence, 'registration consumed'):
                execute(registration, Path(temporary)/'retry-failure',
                        expected_spec=spec, timeout_seconds=30)

    def test_late_live_checkout_import_rejects_terminal_source_admission(self):
        with tempfile.TemporaryDirectory() as temporary:
            external = Path(temporary)/'external'
            external.mkdir()
            (external/'late.py').write_text('VALUE = True\n')
            root, registration, output = self.fixture(temporary, '''def run(arguments, output):
    import sys
    sys.path.insert(0, arguments["external"])
    import late
    return {"value": late.VALUE}
''')
            spec = register(root, registration, module='control', action='run',
                            arguments={'external': str(external)}, timeout_seconds=30)
            process = execute(registration, output, expected_spec=spec, timeout_seconds=30)
            self.assertEqual(process['exit_code'], 2)
            self.assertTrue((output/'before.json').is_file())
            self.assertFalse((output/'after.json').exists())
            self.assertIn('unregistered external SAT Python import',
                          (output/'stderr.txt').read_text())
            with self.assertRaises(FileNotFoundError):
                audit_execution(output, spec)

    def test_changed_interpreter_registration_fails_before_launch(self):
        with tempfile.TemporaryDirectory() as temporary:
            root, registration, output = self.fixture(temporary,
                'def run(arguments, output):\n    return {}\n')
            spec = register(root, registration, module='control', action='run', arguments={},
                            timeout_seconds=30)
            changed = copy.deepcopy(spec)
            changed['interpreter']['executable_sha256'] = '0'*64
            with patch('sat_runtime_execution_v3.interpreter_record',
                       return_value=changed['interpreter']):
                with self.assertRaisesRegex(InvalidEvidence, 'interpreter changed'):
                    execute(registration, output, expected_spec=spec, timeout_seconds=30)
            self.assertFalse(output.exists())
            self.assertFalse((registration/EXECUTION_POLICY['claim_file']).exists())

    def test_retained_receipt_cannot_change_invocation_or_promote_cost(self):
        with tempfile.TemporaryDirectory() as temporary:
            root, registration, output = self.fixture(temporary,
                'def run(arguments, output):\n    return {"count": arguments["count"]}\n')
            spec = register(root, registration, module='control', action='run',
                            arguments={'count': 1}, timeout_seconds=30)
            with self.assertRaisesRegex(InvalidEvidence, 'watchdog differs'):
                execute(registration, output, expected_spec=spec, timeout_seconds=60)
            self.assertFalse(output.exists())
            changed = copy.deepcopy(spec)
            changed['arguments']['count'] = 2
            (registration/'execution.json').write_text(json.dumps(changed))
            with self.assertRaisesRegex(InvalidEvidence, 'independently sealed'):
                execute(registration, output, expected_spec=spec, timeout_seconds=30)
            self.assertFalse(output.exists())
            (registration/'execution.json').write_text(json.dumps(spec))
            process = execute(registration, output, expected_spec=spec, timeout_seconds=30)
            self.assertEqual(process['exit_code'], 0,
                             (output/'stderr.txt').read_text())
            other = copy.deepcopy(spec)
            other['arguments']['count'] = 2
            with self.assertRaisesRegex(InvalidEvidence, 'sealed registration'):
                audit_execution(output, other)
            process['online_speedup'] = '2'
            (output/'process.json').write_text(json.dumps(process))
            with self.assertRaisesRegex(InvalidEvidence, 'evidence is incomplete or changed'):
                audit_execution(output, spec)

    def test_watchdog_retains_partial_sources_and_never_admits_timeout(self):
        with tempfile.TemporaryDirectory() as temporary:
            root, registration, output = self.fixture(temporary,
                'def run(arguments, output):\n    import time\n    time.sleep(10)\n')
            spec = register(root, registration, module='control', action='run', arguments={},
                            timeout_seconds=1)
            process = execute(registration, output, expected_spec=spec, timeout_seconds=1)
            self.assertTrue(process['timed_out'])
            self.assertTrue((output/'runtime/runtime.tar.gz').is_file())
            self.assertFalse((output/'after.json').exists())
            with self.assertRaises((FileNotFoundError, InvalidEvidence)):
                audit_execution(output, spec)
            with self.assertRaisesRegex(InvalidEvidence, 'registration consumed'):
                execute(registration, Path(temporary)/'retry-timeout',
                        expected_spec=spec, timeout_seconds=1)

    def test_claim_arbitrates_concurrent_distinct_outputs_and_cannot_replay(self):
        with tempfile.TemporaryDirectory() as temporary:
            root, registration, output = self.fixture(temporary,
                'def run(arguments, output):\n    return {}\n')
            spec = register(root, registration, module='control', action='run', arguments={},
                            timeout_seconds=30)
            barrier = threading.Barrier(2)

            def attempt(index):
                barrier.wait()
                try:
                    return claim_execution(registration, output.with_name(str(index)), spec)
                except InvalidEvidence:
                    return None

            with ThreadPoolExecutor(max_workers=2) as pool:
                claims = list(pool.map(attempt, (1, 2)))
            winner, = [claim for claim in claims if claim is not None]
            retained = json.loads((registration/EXECUTION_POLICY['claim_file']).read_text())
            self.assertEqual(retained, winner)
            # Identical bytes still cannot acquire an already consumed claim.
            with self.assertRaisesRegex(InvalidEvidence, 'registration consumed'):
                claim_execution(registration, Path(winner['output_directory']), spec)
            self.assertFalse(output.exists())

    def test_setup_failure_after_claim_is_consumed_and_claim_tamper_rejects(self):
        with tempfile.TemporaryDirectory() as temporary:
            root, registration, output = self.fixture(temporary,
                'def run(arguments, output):\n    return {}\n')
            spec = register(root, registration, module='control', action='run', arguments={},
                            timeout_seconds=30)
            with patch('sat_runtime_execution_v3.extract', side_effect=OSError('control setup failed')):
                with self.assertRaisesRegex(OSError, 'control setup failed'):
                    execute(registration, output, expected_spec=spec, timeout_seconds=30)
            self.assertTrue((registration/EXECUTION_POLICY['claim_file']).exists())
            with self.assertRaisesRegex(InvalidEvidence, 'registration consumed'):
                execute(registration, Path(temporary)/'new-output', expected_spec=spec, timeout_seconds=30)
        with tempfile.TemporaryDirectory() as temporary:
            root, registration, output = self.fixture(temporary,
                'def run(arguments, output):\n    return {}\n')
            spec = register(root, registration, module='control', action='run', arguments={},
                            timeout_seconds=30)
            execute(registration, output, expected_spec=spec, timeout_seconds=30)
            claim = json.loads((output/EXECUTION_POLICY['claim_file']).read_text())
            claim['output_directory'] = '/changed-retained-path'
            (output/EXECUTION_POLICY['claim_file']).write_text(json.dumps(claim))
            with self.assertRaisesRegex(InvalidEvidence, 'claim differs from process/source gates'):
                audit_execution(output, spec)

    def test_legacy_registration_cannot_launch_through_new_code(self):
        with tempfile.TemporaryDirectory() as temporary:
            root, registration, output = self.fixture(temporary,
                'def run(arguments, output):\n    return {}\n')
            spec = register(root, registration, module='control', action='run', arguments={},
                            timeout_seconds=30)
            spec.pop('execution_policy')
            (registration/'execution.json').write_text(json.dumps(spec))
            with self.assertRaisesRegex(InvalidEvidence, 'legacy registration is audit-only'):
                execute(registration, output, expected_spec=spec, timeout_seconds=30)
            self.assertFalse(output.exists())

    def test_moved_archive_audits_original_claim_without_reopening_registration(self):
        with tempfile.TemporaryDirectory() as temporary:
            root, registration, output = self.fixture(temporary,
                'def run(arguments, output):\n    return {}\n')
            spec = register(root, registration, module='control', action='run', arguments={},
                            timeout_seconds=30)
            execute(registration, output, expected_spec=spec, timeout_seconds=30)
            moved = Path(temporary)/'transported-archive'
            output.rename(moved)
            shutil.rmtree(registration)
            self.assertTrue(audit_execution(moved, spec)['entrypoint_succeeded'])
            claim = json.loads((moved/EXECUTION_POLICY['claim_file']).read_text())
            self.assertEqual(claim['output_directory'], str(output))

    def test_watchdog_does_not_resignal_a_group_after_successful_kill_and_reap(self):
        with tempfile.TemporaryDirectory() as temporary:
            root, registration, output = self.fixture(temporary,
                'def run(arguments, output):\n    import time\n    time.sleep(10)\n')
            spec = register(root, registration, module='control', action='run',
                            arguments={}, timeout_seconds=1)
            actual_killpg = os.killpg
            calls = []

            def killed_group_can_become_inaccessible(pgid, signal):
                calls.append(pgid)
                if len(calls) > 1:
                    raise PermissionError('simulated macOS post-reap orphan-zombie group')
                return actual_killpg(pgid, signal)

            with patch('sat_runtime_execution_v3.os.killpg',
                       side_effect=killed_group_can_become_inaccessible):
                process = execute(registration, output, expected_spec=spec, timeout_seconds=1)
            self.assertTrue(process['timed_out'])
            self.assertEqual(len(calls), 1)
            self.assertLess(process['exit_code'], 0)
            self.assertTrue((output/'process.json').exists())
            self.assertFalse((output/'after.json').exists())


if __name__ == '__main__':
    unittest.main()
