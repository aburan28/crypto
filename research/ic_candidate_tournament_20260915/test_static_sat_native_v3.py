"""Real non-SAT binaries test child source binding, relocation and watchdogs."""
import copy
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import time
import unittest

from oracle import InvalidEvidence
from sat_runtime_execution_v3 import audit_execution, execute, register
from static_sat_assets_v3 import freeze_assets
from static_sat_native_v3 import audit_meter

ROOT=Path(__file__).resolve().parents[2]


class StaticSatNativeV3Tests(unittest.TestCase):
    def registered(self,temporary,tool,argv,seconds,outer_seconds=30):
        root=Path(temporary)
        # macOS rejects relocated Apple system binaries through code signing.
        # Build a tiny portable control; it has no SAT or IC implementation.
        source=root/'control.c'
        source.write_text('#include <stdio.h>\n#include <stdlib.h>\n#include <unistd.h>\n'
                          'int main(int n,char **v){if(n<2)return 1;'
                          'if(n>2){sleep((unsigned)atoi(v[2]));return 0;}'
                          'puts(v[1]);return 0;}\n')
        tool=root/'control'
        compiler=shutil.which('cc')
        if compiler is None:
            self.skipTest('portable native process control requires a C compiler')
        subprocess.run([compiler,'-O0',str(source),'-o',str(tool)],
                       capture_output=True,check=True)
        assets=root/'assets'
        freeze_assets({'bin/control':Path(tool).read_bytes()},{'bin/control'},assets)
        spec=register(ROOT,root/'registration',module='static_sat_native_v3',
                      action='process_control',arguments=dict(role='bin/control',
                      argv=argv,seconds=seconds),timeout_seconds=outer_seconds,
                      asset_snapshot=assets)
        return root,spec

    def test_child_receipts_survive_transport_and_reject_command_changes(self):
        with tempfile.TemporaryDirectory() as temporary:
            root,spec=self.registered(temporary,'/bin/echo',['source control'],5)
            process=execute(root/'registration',root/'execution',expected_spec=spec,
                            timeout_seconds=30)
            self.assertEqual(process['exit_code'],0,(root/'execution/stderr.txt').read_text())
            self.assertTrue(audit_execution(root/'execution',spec)['entrypoint_succeeded'])
            shutil.copytree(root/'execution',root/'transported')
            transported=root/'transported'
            receipt=audit_meter(transported,transported/'entry-output','control',
                                asset_role='bin/control',arguments=['source control'],seconds=5)
            self.assertEqual(receipt['returncode'],0)
            self.assertEqual((transported/'entry-output/control.stdout').read_text(),'source control\n')
            with self.assertRaisesRegex(InvalidEvidence,'command or watchdog'):
                audit_meter(transported,transported/'entry-output','control',
                            asset_role='bin/control',arguments=['changed'],seconds=5)
            path=transported/'entry-output/control.after.json'
            path.chmod(0o644)
            gate=json.loads(path.read_text())
            gate['loaded_modules']={}
            path.write_text(json.dumps(gate))
            with self.assertRaisesRegex(InvalidEvidence,'source gate'):
                audit_meter(transported,transported/'entry-output','control',
                            asset_role='bin/control',arguments=['source control'],seconds=5)

    def test_native_timeout_is_retained_with_complete_wrapper_gates(self):
        with tempfile.TemporaryDirectory() as temporary:
            root,spec=self.registered(temporary,'/bin/sleep',['sleep','10'],1)
            process=execute(root/'registration',root/'execution',expected_spec=spec,
                            timeout_seconds=30)
            self.assertEqual(process['exit_code'],0,(root/'execution/stderr.txt').read_text())
            receipt=audit_meter(root/'execution',root/'execution/entry-output','control',
                                asset_role='bin/control',arguments=['sleep','10'],seconds=1)
            self.assertTrue(receipt['timed_out'])
            self.assertLess(receipt['returncode'],0)
            self.assertTrue(audit_execution(root/'execution',spec)['complete_source_gates'])

    def test_outer_watchdog_kills_inherited_native_group(self):
        with tempfile.TemporaryDirectory() as temporary:
            # Source/interpreter preflight varies with the installed stdlib and
            # hosted-runner I/O. This is a cleanup control, not a startup-speed
            # assertion. Leave ample time for launch; the native sleep remains
            # strictly longer than both registered watchdogs.
            root,spec=self.registered(temporary,'/bin/sleep',['sleep','120'],60,outer_seconds=30)
            process=execute(root/'registration',root/'execution',expected_spec=spec,
                            timeout_seconds=30)
            self.assertTrue(process['timed_out'])
            spawned=json.loads((root/'execution/entry-output/control.spawned.json').read_text())
            pid=spawned['native_pid']
            deadline=time.monotonic()+3
            while time.monotonic()<deadline:
                observed=subprocess.run(['/bin/ps','-p',str(pid),'-o','stat='],
                                        capture_output=True,text=True,check=False)
                if observed.returncode!=0 or not observed.stdout.strip() or observed.stdout.strip().startswith('Z'):
                    break
                time.sleep(0.05)
            else:
                self.fail('native process survived its enclosing watchdog group')
            self.assertFalse((root/'execution/after.json').exists())
            with self.assertRaises((FileNotFoundError,InvalidEvidence)):
                audit_execution(root/'execution',spec)


if __name__=='__main__':
    unittest.main()
