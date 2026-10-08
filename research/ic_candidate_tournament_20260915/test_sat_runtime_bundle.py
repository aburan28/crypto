"""Package imports, late imports and replay remain bound to frozen bytes."""
import copy
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
from types import SimpleNamespace
import unittest

from oracle import InvalidEvidence
from sat_runtime_bundle import (CHILD_SCRIPTS, DIRECTORY, check_loaded_modules,
                                extract, freeze, source_manifest, verified_files)


class SatRuntimeBundleTests(unittest.TestCase):
    def test_current_sat_runtime_package_imports_are_covered(self):
        root = Path(__file__).resolve().parents[2]
        script = '''import json,sys
from pathlib import Path
import run_static_sat_full_v2
from sat_runtime_bundle import check_loaded_modules,source_manifest
root=Path(sys.argv[1])
print(json.dumps(check_loaded_modules(root,source_manifest(root))))
'''
        environment = {key: value for key, value in os.environ.items()
                       if key not in ('PYTHONPATH', 'PYTHONSTARTUP')}
        loaded = json.loads(subprocess.check_output(
            [sys.executable, '-c', script, str(root)],
            cwd=root/DIRECTORY, env=environment, text=True, timeout=30))
        self.assertEqual(loaded['producer.evidence'],
                         (DIRECTORY/'producer/evidence.py').as_posix())
        self.assertEqual(loaded['producer.timing'],
                         (DIRECTORY/'producer/timing.py').as_posix())
        self.assertEqual(loaded['run_static_sat_full_v2'],
                         (DIRECTORY/'run_static_sat_full_v2.py').as_posix())

    def fixture(self, root):
        directory = root/DIRECTORY
        (directory/'producer').mkdir(parents=True)
        for name, contents in (
                ('runner.py', 'from producer.evidence import verify\n'),
                ('producer/evidence.py', 'def verify(): return True\n'),
                ('producer/timing.py', 'CLOCK = "exclusive"\n')):
            (directory/name).write_text(contents)
        (root/'scripts').mkdir()
        for role in CHILD_SCRIPTS:
            (root/role).write_text('CONTROL = True\n')
        return directory

    def test_namespace_package_is_bound_even_without_init_file(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            directory = self.fixture(root)
            manifest = source_manifest(root)
            roles = {item['role'] for item in manifest['components']}
            self.assertIn((DIRECTORY/'producer/evidence.py').as_posix(), roles)
            self.assertIn((DIRECTORY/'producer/timing.py').as_posix(), roles)
            modules = {'producer.evidence': SimpleNamespace(
                __file__=str(directory/'producer/evidence.py'))}
            self.assertEqual(check_loaded_modules(root, manifest, modules),
                             {'producer.evidence':
                              (DIRECTORY/'producer/evidence.py').as_posix()})

    def test_package_mutation_changes_identity_and_fails_loaded_gate(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            directory = self.fixture(root)
            original = source_manifest(root)
            module = directory/'producer/evidence.py'
            module.write_text('def verify(): return False\n')
            self.assertNotEqual(original, source_manifest(root))
            with self.assertRaisesRegex(InvalidEvidence, 'source changed'):
                check_loaded_modules(root, original, {'producer.evidence':
                    SimpleNamespace(__file__=str(module))})

    def test_late_repository_import_outside_surface_is_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self.fixture(root)
            manifest = source_manifest(root)
            module = root/'unregistered.py'
            module.write_text('VALUE = 1\n')
            with self.assertRaisesRegex(InvalidEvidence, 'unregistered local'):
                check_loaded_modules(root, manifest, {'late':
                    SimpleNamespace(__file__=str(module))})

    def test_live_checkout_cannot_supply_frozen_runtime_import(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)/'repository'
            directory = self.fixture(root)
            snapshot = Path(temporary)/'snapshot'
            manifest, seal = freeze(root, snapshot)
            restored = extract(snapshot, Path(temporary)/'restored', manifest, seal)
            with self.assertRaisesRegex(InvalidEvidence,
                                        'unregistered external'):
                check_loaded_modules(restored, manifest, {'producer.evidence':
                    SimpleNamespace(__file__=str(directory/'producer/evidence.py'))})

    def test_snapshot_replays_original_after_live_source_changes(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)/'repository'
            directory = self.fixture(root)
            snapshot = Path(temporary)/'snapshot'
            manifest, seal = freeze(root, snapshot)
            (directory/'producer/evidence.py').write_text('changed = True\n')
            restored = extract(snapshot, Path(temporary)/'restored', manifest, seal)
            self.assertEqual(source_manifest(restored), manifest)
            with self.assertRaisesRegex(InvalidEvidence, 'already exists'):
                freeze(root, snapshot)
            with (snapshot/'runtime.tar.gz').open('ab') as archive:
                archive.write(b'changed')
            with self.assertRaisesRegex(InvalidEvidence, 'archive changed'):
                verified_files(snapshot, manifest, seal)

    def test_registered_manifest_cannot_be_replaced_with_current_source(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)/'repository'
            self.fixture(root)
            snapshot = Path(temporary)/'snapshot'
            manifest, seal = freeze(root, snapshot)
            other = copy.deepcopy(manifest)
            other['components'][0]['sha256'] = '0'*64
            with self.assertRaisesRegex(InvalidEvidence,
                                        'candidate registration'):
                verified_files(snapshot, other, seal)

    def test_package_symlink_cannot_supply_unbound_source(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)/'repository'
            directory = self.fixture(root)
            external = Path(temporary)/'external.py'
            external.write_text('VALUE = 1\n')
            (directory/'producer/extra.py').symlink_to(external)
            with self.assertRaisesRegex(InvalidEvidence, 'symlinked'):
                source_manifest(root)


if __name__ == '__main__':
    unittest.main()
