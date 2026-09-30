"""Native inputs remain exact, read-only, and outside the Python source surface."""
import json
import os
from pathlib import Path
import tempfile
import unittest

from oracle import InvalidEvidence
from static_sat_assets_v3 import (check_extracted_assets, extract_assets, freeze_assets,
                                   verified_assets)


class StaticSatAssetsV3Tests(unittest.TestCase):
    def test_exact_roundtrip_modes_mutation_and_symlink_rejection(self):
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary)
            files={'bin/control':b'executable-control','fixture.json':b'{}\n'}
            manifest,seal=freeze_assets(files,{'bin/control'},root/'bundle')
            self.assertEqual(verified_assets(root/'bundle',manifest,seal),files)
            extracted=extract_assets(root/'bundle',root/'files',manifest,seal)
            self.assertEqual(check_extracted_assets(extracted,manifest),files)
            self.assertTrue(os.access(extracted/'bin/control',os.X_OK))
            self.assertFalse(os.access(extracted/'fixture.json',os.X_OK))
            (extracted/'fixture.json').chmod(0o644)
            (extracted/'fixture.json').write_bytes(b'{"changed":true}\n')
            with self.assertRaisesRegex(InvalidEvidence,'differ'):
                check_extracted_assets(extracted,manifest)
            (extracted/'fixture.json').unlink()
            (extracted/'fixture.json').symlink_to(extracted/'bin/control')
            with self.assertRaisesRegex(InvalidEvidence,'symlinked'):
                check_extracted_assets(extracted,manifest)

    def test_python_and_unsafe_members_never_enter_asset_bundle(self):
        for role in ('code.py','nested/cache.pyc','../escape','/absolute','a//b','.'):
            with self.subTest(role=role),tempfile.TemporaryDirectory() as temporary:
                with self.assertRaises(InvalidEvidence):
                    freeze_assets({role:b'x'},set(),Path(temporary)/'bundle')

    def test_archive_or_seal_change_cannot_rebind_assets(self):
        with tempfile.TemporaryDirectory() as temporary:
            bundle=Path(temporary)/'bundle'
            manifest,seal=freeze_assets({'fixture.json':b'{}'},set(),bundle)
            archive=bundle/'assets.tar.gz'
            archive.write_bytes(archive.read_bytes()+b'changed')
            with self.assertRaisesRegex(InvalidEvidence,'archive changed'):
                verified_assets(bundle,manifest,seal)
            with self.assertRaisesRegex(InvalidEvidence,'already exists'):
                freeze_assets({'fixture.json':b'{}'},set(),bundle)


if __name__=='__main__':
    unittest.main()
