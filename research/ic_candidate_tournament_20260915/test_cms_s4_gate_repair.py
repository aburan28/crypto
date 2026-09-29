"""Pre-continuation proof that only the manifest reader's key needs repair."""
from pathlib import Path
import tempfile
import unittest

from continue_cms_s4_controls import (
    PANEL, STAGE_A_SHA256, stage_a_files, verify_stage_a,
)
from run_cms_s4_controls import preflight, sha_bytes, validate_export
from tournament import read


class CmsS4GateRepairTests(unittest.TestCase):
    def test_sealed_exports_and_no_original_sat_dispatch(self):
        panel = read(PANEL)
        _, curve, base = preflight(panel)
        self.assertEqual(sha_bytes(
            (PANEL.parent/'stage-a-evidence.tar.gz').read_bytes()), STAGE_A_SHA256)
        files = stage_a_files()
        original, registered = verify_stage_a(files, panel, curve, base)
        self.assertEqual(len(registered), 4)
        self.assertTrue(all(row['status'] == 'INVALID_EXPORT'
                            for row in original['rows']))
        self.assertFalse(any('/cms.stdout' in name or '/cms.metrics.json' in name
                             for name in files))
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            for item, _, manifest in registered:
                instance = root/f"trial-{item['trial']:02d}"/'instance'
                instance.mkdir(parents=True)
                prefix = f"trial-{item['trial']:02d}/instance/"
                for name, payload in files.items():
                    if name.startswith(prefix):
                        (instance/name.removeprefix(prefix)).write_bytes(payload)
                exports = validate_export(manifest, item, panel, base, instance)
                self.assertEqual(set(exports),
                                 {'wdsat_anf', 'cryptominisat_xor_dimacs',
                                  'magma_boolean_f4'})


if __name__ == '__main__':
    unittest.main()
