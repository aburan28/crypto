"""Static contract for the campaign bootstrap / publish gate.

The 2026-09-17 fleet died because binaryKey existed at
bin/72dcb624cffd7de6/ecc2k130 while test-packed-storage-cuda 404'd, and
bootstrap.sh treated that as fatal instead of rebuilding from sourceKey.
These checks pin the scripts so that path cannot return.
"""
from pathlib import Path
import subprocess
import unittest

HERE = Path(__file__).resolve().parent
BOOTSTRAP = (HERE / "bootstrap.sh").read_text()
BUILD = (HERE / "build.sh").read_text()
INFRA = (HERE / "infra.sh").read_text()

FIXTURES = (
    "test-packed-cuda",
    "test-packed-storage-cuda",
    "test-shared-sigma-cuda",
)


class BootstrapContract(unittest.TestCase):
    def test_fixtures_named(self):
        for name in FIXTURES:
            self.assertIn(name, BOOTSTRAP)
            self.assertIn(name, BUILD)

    def test_incomplete_prefix_rebuilds(self):
        self.assertIn("prefix_complete", BOOTSTRAP)
        self.assertIn("published prefix incomplete or missing sm_", BOOTSTRAP)
        self.assertIn("prefix_usable", BOOTSTRAP)
        self.assertNotIn("fixture $f missing from $PREFIX; not starting workers", BOOTSTRAP)

    def test_rebuild_requires_source_key(self):
        self.assertIn("no sourceKey in campaign.json", BOOTSTRAP)
        self.assertIn("aws/build.sh", BOOTSTRAP)

    def test_build_refuses_incomplete_publish(self):
        self.assertIn("refusing to point campaign.json at", BUILD)
        self.assertIn("build did not produce", BUILD)
        gate = BUILD.split("Point the campaign at this build")[0]
        self.assertIn("head-object", gate)

    def test_infra_syncs_build_script(self):
        self.assertIn('aws s3 cp build.sh "s3://$BUCKET/aws/build.sh"', INFRA)

    def test_unversioned_campaign_allows_legacy_storage(self):
        self.assertIn("ECC_ALLOW_LEGACY_STORAGE=1", BOOTSTRAP)
        self.assertIn("field storageProtocol", BOOTSTRAP)
        # Must not stamp the strict protocol onto the live corpus.
        self.assertNotIn("storageProtocol=", BOOTSTRAP)

    def test_fat_client_covers_ada_and_blackwell(self):
        self.assertIn("local_cc", BOOTSTRAP)
        self.assertIn('ARCHES="${ARCHES:-89 120}"', BOOTSTRAP)
        self.assertIn("manifest.json", BOOTSTRAP)
        # A thin sm_89 publish must not become the live binaryKey.
        self.assertIn("g6/g6e (sm_89) and g7e (sm_120) share one binaryKey", BOOTSTRAP)

    def test_launch_g6_is_spot_ada_only(self):
        script = (HERE / "launch_g6.sh").read_text()
        self.assertIn("g6e.2xlarge,g6.2xlarge", script)
        self.assertIn('MarketType":"spot"', script)
        self.assertNotIn("g7e.2xlarge", script)
        self.assertNotIn("on-demand", script)
        subprocess.run(["bash", "-n", str(HERE / "launch_g6.sh")], check=True)


if __name__ == "__main__":
    unittest.main()
