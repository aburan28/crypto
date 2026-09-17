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
        gate = BUILD.split("Geometry in campaign.json must match")[0]
        self.assertIn("head-object", gate)

    def test_stage_does_not_point_campaign(self):
        self.assertIn("POINT_CAMPAIGN", BUILD)
        self.assertIn("--stage", BUILD)
        self.assertIn('if [ "$POINT_CAMPAIGN" = 1 ]', BUILD)
        self.assertIn("campaign.json unchanged", BUILD)
        # Bootstrap rebuilds still point: first Ada fat publish.
        self.assertIn("POINT_CAMPAIGN=1", BOOTSTRAP)
        self.assertNotIn("--stage", BOOTSTRAP)

    def test_kernel_versions_are_immutable(self):
        self.assertIn("ecc2k-kernel-v1", BUILD)
        self.assertIn("kernels/$kver.json", BUILD)
        self.assertIn("--if-none-match", BUILD)
        # Must not stamp storageProtocol onto the live corpus.
        self.assertNotIn("storageProtocol=", BUILD)

    def test_infra_syncs_build_script(self):
        self.assertIn('aws s3 cp build.sh "s3://$BUCKET/aws/build.sh"', INFRA)
        self.assertIn('aws s3 cp rollout.sh "s3://$BUCKET/aws/rollout.sh"', INFRA)
        self.assertIn('aws s3 cp rollout.py "s3://$BUCKET/aws/rollout.py"', INFRA)

    def test_unversioned_campaign_allows_legacy_storage(self):
        self.assertIn("ECC_ALLOW_LEGACY_STORAGE=1", BOOTSTRAP)
        self.assertIn("field storageProtocol", BOOTSTRAP)
        # Must not stamp the strict protocol onto the live corpus.
        self.assertNotIn("storageProtocol=", BOOTSTRAP)

    def test_fat_client_covers_ada_and_blackwell(self):
        self.assertIn("local_cc", BOOTSTRAP)
        self.assertIn('ARCHES="${ARCHES:-75 89 120}"', BOOTSTRAP)
        self.assertIn('CLMAD="${CLMAD:-1}"', BOOTSTRAP)
        self.assertIn("manifest.json", BOOTSTRAP)
        # A thin sm_89/sm_75 publish must not become the live binaryKey.
        self.assertIn("g4dn (sm_75), g6/g6e (sm_89) and g7e (sm_120) share one", BOOTSTRAP)
        self.assertIn("binaryKey. A thin Ada/T4 rebuild", BOOTSTRAP)
        self.assertIn("ECC_INSTANCE_TYPE", BOOTSTRAP)
        self.assertIn("ECC_DEVICE_NAME", BOOTSTRAP)

    def test_launch_g6_is_spot_ada_only(self):
        script = (HERE / "launch_g6.sh").read_text()
        self.assertIn("g6e.2xlarge,g6.2xlarge", script)
        self.assertIn('MarketType":"spot"', script)
        self.assertNotIn("g7e.2xlarge", script)
        self.assertNotIn("on-demand", script)
        self.assertIn("bash ./infra.sh", script)
        self.assertIn("bash ./push_source.sh", script)
        subprocess.run(["bash", "-n", str(HERE / "launch_g6.sh")], check=True)

    def test_launch_spot_all_is_discovered_spot_only(self):
        script = (HERE / "launch_spot_all.sh").read_text()
        self.assertIn("g7e.2xlarge,g7.2xlarge,g6e.2xlarge,g6.2xlarge,g4dn.2xlarge", script)
        self.assertIn('MarketType":"spot"', script)
        self.assertIn("describe-regions --all-regions", script)
        self.assertIn("opted-in", script)
        self.assertNotIn("on-demand", script)
        self.assertNotIn("TerminateInstances", script)
        self.assertIn("bash ./infra.sh", script)
        self.assertIn("bash ./push_source.sh", script)
        self.assertIn("meow34", script)
        self.assertIn('REGIONS:-', script)
        self.assertIn("describe-instance-type-offerings", script)
        self.assertIn("g4dn.xlarge", script)
        self.assertIn("g4dn.4xlarge", script)
        self.assertIn("create-fleet", script)
        self.assertIn("capacity-optimized", script)
        self.assertIn("--type", script)
        self.assertIn("instant", script)
        self.assertNotIn("create-fleet --type maintain", script)
        subprocess.run(["bash", "-n", str(HERE / "launch_spot_all.sh")], check=True)


if __name__ == "__main__":
    unittest.main()
