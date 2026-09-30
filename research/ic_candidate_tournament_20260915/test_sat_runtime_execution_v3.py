"""Real isolated imports and incomplete executions retain their source boundary."""
import copy
import json
from pathlib import Path
import shutil
import tempfile
import unittest
from unittest.mock import patch

from oracle import InvalidEvidence
from sat_runtime_bundle import CHILD_SCRIPTS, DIRECTORY
from sat_runtime_execution_v3 import audit_execution, execute, register

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


if __name__ == '__main__':
    unittest.main()
