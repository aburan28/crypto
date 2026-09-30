"""Transport real native evidence without running a producer; fail on tampering."""
import json
from pathlib import Path
import shutil
import tempfile
import unittest

from audit_static_sat_full_v3 import audit
from oracle import InvalidEvidence
from publish_static_sat_v3_control import replay

BUNDLE=Path(__file__).resolve().parent/'goal_20260924/static-sat-runtime-v3/native-smoke-20260929/results-20260929'
EXPECTED='3ce0af74afe448840fffc28f80eece0e8e3c04681f074a347b030bb87d2853f5'


class PublishStaticSatV3ControlTests(unittest.TestCase):
    def test_real_native_control_replays_after_transport(self):
        with tempfile.TemporaryDirectory() as temporary:
            output=Path(temporary)/'transport'
            result=replay(BUNDLE,output,EXPECTED)
            self.assertTrue(result['complete_source_bound'])
            self.assertEqual(result['statuses'],{'SOURCE_UNSAT':1})
            self.assertEqual(result['final_rank'],0)
            self.assertFalse(result['headline_online_admissible'])
            self.assertFalse(result['promotion_eligible'])
            self.assertIsNone(result['online_speedup'])
            self.assertIsNone(result['online_wall_ns'])

    def test_external_seal_or_changed_archive_fails_before_extraction(self):
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary)
            with self.assertRaisesRegex(InvalidEvidence,'externally frozen'):
                replay(BUNDLE,root/'wrong-seal','0'*64)
            self.assertFalse((root/'wrong-seal').exists())
            shutil.copytree(BUNDLE,root/'changed')
            archive=root/'changed/evidence.tar.gz'
            with archive.open('ab') as stream:
                stream.write(b'changed')
            with self.assertRaisesRegex(InvalidEvidence,'archive changed'):
                replay(root/'changed',root/'changed-output',EXPECTED)
            self.assertFalse((root/'changed-output').exists())

    def test_native_command_or_prime_matrix_label_cannot_be_rewritten(self):
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary)/'transport'
            replay(BUNDLE,root,EXPECTED)
            spec=json.loads((root/'registration/execution.json').read_text())
            intent=root/'execution/entry-output/collection/trial-00/cms.intent.json'
            old=intent.read_text()
            call=json.loads(old)
            call['arguments'][call['arguments'].index('--maxconfl')+1]='1000001'
            intent.write_text(json.dumps(call))
            with self.assertRaisesRegex(InvalidEvidence,'command or watchdog'):
                audit(root/'execution',spec)
            intent.write_text(old)
            matrix=root/'execution/entry-output/relation-matrix.json'
            summary=root/'execution/entry-output/summary.json'
            snapshot=json.loads(matrix.read_text())
            result=json.loads(summary.read_text())
            snapshot['modulus']='2'
            result['matrix']=snapshot
            matrix.write_text(json.dumps(snapshot))
            summary.write_text(json.dumps(result))
            with self.assertRaisesRegex(InvalidEvidence,'final matrix'):
                audit(root/'execution',spec)

    def test_missing_native_terminal_gate_never_receives_source_admission(self):
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary)/'transport'
            replay(BUNDLE,root,EXPECTED)
            spec=json.loads((root/'registration/execution.json').read_text())
            (root/'execution/entry-output/collection/trial-00/cms.after.json').unlink()
            with self.assertRaises(FileNotFoundError):
                audit(root/'execution',spec)


if __name__=='__main__':
    unittest.main()
