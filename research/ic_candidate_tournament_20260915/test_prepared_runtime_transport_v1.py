"""Real isolated transport controls; no native solver or IC measurement.

The tiny known-name auditor is synthetic. These tests establish source and
artifact transport, not mathematical admission of an actual SAT/F5 run.
"""
import json
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest
from unittest.mock import patch

from identity import sha256
from oracle import InvalidEvidence
from prepared_runtime_transport_v1 import transport
from sat_runtime_bundle import CHILD_SCRIPTS, DIRECTORY
from sat_runtime_execution_v3 import audit_execution, execute, register

HERE = Path(__file__).resolve().parent
SEAL = dict(candidate_id='synthetic-transport-control', workload_id='synthetic', run_id='synthetic')
CONTROL = '''
from pathlib import Path
from identity import write_immutable
from oracle import require
from sat_runtime_execution_v3 import audit_execution, read
def run(arguments, output):
    write_immutable(output/'control.json', {'scope': 'synthetic; no native or IC solve'})
    return {'native_solvers_executed': 0}
def audit(execution, spec):
    require(audit_execution(execution, spec)['entrypoint_succeeded'], 'control execution failed')
    require(read(Path(execution)/'entry-output/control.json') ==
            {'scope': 'synthetic; no native or IC solve'}, 'synthetic retained artifact changed')
    return dict(spec['arguments']['seal'], status=STATUS,
                source_bound_execution_admitted=True, scalar_verified=False,
                online_wall_ns=None, headline_online_admissible=False,
                fresh_paired_qualification=False, promotion_eligible=False, online_speedup=None,
                scope='synthetic source-transport control; no mathematical or performance evidence')
'''


class PreparedRuntimeTransportV1Tests(unittest.TestCase):
    def fixture(self, temporary, name='prepared_sat_runtime_v1', extra=''):
        root = Path(temporary)/'repository'
        directory = root/DIRECTORY
        (directory/'producer').mkdir(parents=True)
        (directory/'producer/control.py').write_text('SYNTHETIC = True\n')
        for file in ('identity.py', 'oracle.py', 'sat_runtime_bundle.py',
                     'sat_runtime_execution_v3.py', 'prepared_runtime_transport_v1.py'):
            shutil.copyfile(HERE/file, directory/file)
        status = ('ADMITTED_INCOMPLETE_PREPARED_SAT_CONTROL' if name == 'prepared_sat_runtime_v1'
                  else 'ADMITTED_INCOMPLETE_TARGET_PREPARED_F5_CONTROL')
        (directory/(name+'.py')).write_text('STATUS = '+repr(status)+'\n'+CONTROL+extra)
        (root/'scripts').mkdir()
        for role in CHILD_SCRIPTS:
            (root/role).write_text('CONTROL = True\n')
        registration, execution = Path(temporary)/'registration', Path(temporary)/'execution'
        spec = register(root, registration, module=name, action='run', arguments={'seal': SEAL},
                        timeout_seconds=30)
        process = execute(registration, execution, expected_spec=spec, timeout_seconds=30)
        self.assertEqual(process['exit_code'], 0, (execution/'stderr.txt').read_text())
        return root, execution, spec, Path(temporary)/'transport'

    def test_real_cli_replays_frozen_sat_auditor_after_live_checkout_changes(self):
        with tempfile.TemporaryDirectory() as temporary:
            root, execution, spec, out = self.fixture(temporary)
            (root/DIRECTORY/'prepared_sat_runtime_v1.py').write_text(
                'raise AssertionError("live checkout must not supply auditor")\n')
            # Invoke the real public CLI, which in turn launches -I -S -B from
            # the preexecution-retained tree. This is not a mocked child.
            import sys
            cp = subprocess.run([sys.executable, '-I', '-S', '-B',
                                 str(HERE/'prepared_runtime_transport_v1.py'), 'audit',
                                 '--execution', str(execution), '--out', str(out),
                                 '--expected-execution-sha256', sha256(spec)],
                                capture_output=True, text=True, timeout=30)
            self.assertEqual(cp.returncode, 0, cp.stderr)
            receipt = json.loads(cp.stdout)
            self.assertEqual(receipt['status'], 'PASS_FROZEN_PREPARED_TRANSPORT')
            self.assertEqual(receipt['native_solvers_executed'], 0)
            self.assertFalse(receipt['promotion_eligible'])
            self.assertIsNone(receipt['online_speedup'])
            before, after = [json.loads((out/file).read_text()) for file in ('before.json', 'after.json')]
            self.assertEqual(before['flags'], dict(isolated=True, no_site=True, bytecode_writes=False))
            self.assertIn('prepared_sat_runtime_v1', after['loaded_modules'])
            self.assertLessEqual(before['loaded_modules'].items(), after['loaded_modules'].items())
            admission = json.loads((out/'admission.json').read_text())
            self.assertFalse(admission['scalar_verified'])
            self.assertIsNone(admission['online_wall_ns'])
            with self.assertRaisesRegex(InvalidEvidence, 'one-use output'):
                transport(execution, sha256(spec), out)

    def test_f5_auditor_uses_the_same_transport_and_incomplete_remains_incomplete(self):
        with tempfile.TemporaryDirectory() as temporary:
            _, execution, spec, out = self.fixture(temporary, name='prepared_f5_runtime_v2')
            result = transport(execution, sha256(spec), out)
            self.assertEqual(result['status'], 'PASS_FROZEN_PREPARED_TRANSPORT')
            admission = json.loads((out/'admission.json').read_text())
            self.assertEqual(admission['status'], 'ADMITTED_INCOMPLETE_TARGET_PREPARED_F5_CONTROL')
            self.assertFalse(admission['scalar_verified'])

    def test_family_identity_or_claim_boundary_tamper_rejects_in_frozen_child(self):
        for change in ("{'run_id': 'changed'}", "{'promotion_eligible': True}",
                       "{'online_speedup': '2'}", "{'online_wall_ns': 123}"):
            with self.subTest(change=change), tempfile.TemporaryDirectory() as temporary:
                extra = '''
original_audit = audit
def audit(execution, spec):
    result = original_audit(execution, spec)
    result.update(CHANGE)
    return result
'''.replace('CHANGE', change)
                _, execution, spec, out = self.fixture(temporary, extra=extra)
                with self.assertRaisesRegex(InvalidEvidence, 'transport rejected'):
                    transport(execution, sha256(spec), out)
                self.assertFalse((out/'admission.json').exists())
                self.assertEqual(json.loads((out/'transport.json').read_text())['exit_code'], 1)

    def test_tampered_artifact_rejects_in_frozen_child_and_retains_failure(self):
        with tempfile.TemporaryDirectory() as temporary:
            _, execution, spec, out = self.fixture(temporary)
            (execution/'entry-output/control.json').write_text('{}\n')
            with self.assertRaisesRegex(InvalidEvidence, 'transport rejected'):
                transport(execution, sha256(spec), out)
            receipt = json.loads((out/'transport.json').read_text())
            self.assertEqual(receipt['status'], 'REJECTED_PREPARED_TRANSPORT')
            self.assertEqual(receipt['exit_code'], 1)
            self.assertTrue((out/'before.json').is_file())
            self.assertFalse((out/'admission.json').exists())
            self.assertIn('synthetic retained artifact changed', (out/'stderr.txt').read_text())

    def test_live_import_during_audit_cannot_pass_terminal_gate(self):
        with tempfile.TemporaryDirectory() as temporary:
            external = Path(temporary)/'external'
            external.mkdir()
            (external/'unbound.py').write_text('VALUE = True\n')
            extra = '''
original_audit = audit
def audit(execution, spec):
    import sys
    sys.path.insert(0, EXTERNAL)
    import unbound
    return original_audit(execution, spec)
'''.replace('EXTERNAL', repr(str(external)))
            _, execution, spec, out = self.fixture(temporary, extra=extra)
            with self.assertRaisesRegex(InvalidEvidence, 'transport rejected'):
                transport(execution, sha256(spec), out)
            self.assertTrue((out/'before.json').exists())
            self.assertFalse((out/'after.json').exists())
            self.assertIn('unregistered external SAT Python import', (out/'stderr.txt').read_text())

    def test_seal_or_retained_source_change_rejects_before_audit_launch(self):
        with tempfile.TemporaryDirectory() as temporary:
            _, execution, spec, out = self.fixture(temporary)
            with patch('prepared_runtime_transport_v1.subprocess.run',
                       side_effect=AssertionError('audit child must not launch')):
                with self.assertRaisesRegex(InvalidEvidence, 'transport seal'):
                    transport(execution, '0'*64, out)
                file = execution/'extracted'/DIRECTORY/'prepared_sat_runtime_v1.py'
                file.chmod(0o644)
                file.write_text('raise RuntimeError("changed source")\n')
                with self.assertRaisesRegex(InvalidEvidence, 'extraction differs'):
                    transport(execution, sha256(spec), out)
            self.assertFalse(out.exists())

    def test_timeout_is_retained_and_does_not_consume_another_solver_attempt(self):
        with tempfile.TemporaryDirectory() as temporary:
            extra = '''
original_audit = audit
def audit(execution, spec):
    import time
    time.sleep(10)
    return original_audit(execution, spec)
'''
            _, execution, spec, out = self.fixture(temporary, extra=extra)
            with patch('prepared_runtime_transport_v1.TIMEOUT_SECONDS', 1):
                with self.assertRaisesRegex(InvalidEvidence, 'transport rejected'):
                    transport(execution, sha256(spec), out)
            receipt = json.loads((out/'transport.json').read_text())
            self.assertTrue(receipt['timed_out'])
            self.assertIsNone(receipt['exit_code'])
            self.assertEqual(receipt['native_solvers_executed'], 0)
            self.assertTrue(audit_execution(execution, spec)['entrypoint_succeeded'])


if __name__ == '__main__':
    unittest.main()
