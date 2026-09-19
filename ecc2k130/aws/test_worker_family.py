"""Ada / Blackwell slot pinning and auto thread count.

g6/g6e must not resume a 385,024-worker Blackwell checkpoint (that retires
the run id) and must not pass --threads from campaign.json.
"""
import unittest
from unittest import mock

import worker
from worker import (gpuFamily, gpuName, idleSlotClaimable, isCpuSlotRecord,
                    slotFamilyCompatible, usesCampaignWorkers,
                    BLACKWELL_FAMILIES, ADA_FAMILIES)


class GpuFamily(unittest.TestCase):
    def test_instance_type_wins(self):
        self.assertEqual(gpuFamily("NVIDIA L4", "g6.2xlarge"), "g6")
        self.assertEqual(gpuFamily("NVIDIA L40S", "g6e.2xlarge"), "g6e")
        self.assertEqual(gpuFamily("", "g7e.48xlarge"), "g7e")
        self.assertEqual(gpuFamily("", "g7.2xlarge"), "g7")

    def test_name_fallback(self):
        self.assertEqual(gpuFamily("NVIDIA L40S"), "g6e")
        self.assertEqual(gpuFamily("NVIDIA L4"), "g6")
        self.assertEqual(gpuFamily("NVIDIA RTX PRO 6000 Blackwell Server Edition"), "g7e")
        self.assertEqual(gpuFamily("NVIDIA RTX PRO 4500 Blackwell"), "g7")
        self.assertEqual(gpuFamily("Tesla T4"), "g4dn")

    def test_l40s_not_classified_as_l4(self):
        self.assertEqual(gpuFamily("NVIDIA L40S"), "g6e")


class UnreadableGpu(unittest.TestCase):
    """A GPU host that cannot read its device must not claim to be a CPU one.

    gpuName goes into the slot record; isCpuSlotRecord reads it back, and
    idleSlotClaimable refuses a CPU-shaped slot to every GPU claimant. So
    returning "cpu" when nvidia-smi fails once strands that slot's checkpoint
    behind a claimant that can only be a CPU worker -- which would refuse the
    packed shape anyway.
    """

    def setUp(self):
        self.env = mock.patch.dict("os.environ", {}, clear=False)
        self.env.start()
        self.addCleanup(self.env.stop)
        for key in ("ECC_DEVICE_NAME", "ECC_DEVICE"):
            worker.os.environ.pop(key, None)

    def test_a_failing_nvidia_smi_is_an_unknown_gpu_not_a_cpu(self):
        with mock.patch.object(worker.subprocess, "run", side_effect=OSError("no nvidia-smi")):
            name = gpuName(0)
        self.assertFalse(isCpuSlotRecord({"gpuName": name}), name)
        self.assertTrue(idleSlotClaimable({"gpuName": name}, {"gpuName": name}, False))

    def test_it_retries_before_giving_up(self):
        replies = [mock.Mock(returncode=1, stdout=""),
                   mock.Mock(returncode=0, stdout="NVIDIA RTX PRO 6000 Blackwell\n")]
        with mock.patch.object(worker.subprocess, "run", side_effect=replies) as run, \
                mock.patch.object(worker.time, "sleep"):
            self.assertEqual(gpuName(0), "NVIDIA RTX PRO 6000 Blackwell")
        self.assertEqual(run.call_count, 2)

    def test_a_declared_cpu_worker_is_still_a_cpu_worker(self):
        worker.os.environ["ECC_DEVICE"] = "cpu"
        worker.os.environ["ECC_THREADS"] = "4"
        self.assertEqual(gpuName(0), "cpu/4")
        self.assertTrue(isCpuSlotRecord({"gpuName": "cpu/4"}))

    def test_an_unknown_gpu_still_pins_its_family_by_instance_type(self):
        self.assertEqual(gpuFamily("gpu0-unknown", "g7e.2xlarge"), "g7e")


class SlotFamily(unittest.TestCase):
    def test_ada_skips_legacy_blackwell_slots(self):
        for family in ADA_FAMILIES:
            self.assertFalse(slotFamilyCompatible(None, family))
            self.assertFalse(slotFamilyCompatible("", family))
            self.assertFalse(slotFamilyCompatible("g7e", family))

    def test_blackwell_may_resume_legacy(self):
        for family in BLACKWELL_FAMILIES:
            self.assertTrue(slotFamilyCompatible(None, family))
            self.assertTrue(slotFamilyCompatible("", family))
            self.assertTrue(slotFamilyCompatible(family, family))

    def test_same_family_only(self):
        self.assertTrue(slotFamilyCompatible("g6", "g6"))
        self.assertTrue(slotFamilyCompatible("g6e", "g6e"))
        self.assertFalse(slotFamilyCompatible("g6", "g6e"))
        self.assertFalse(slotFamilyCompatible("g6e", "g6"))
        self.assertFalse(slotFamilyCompatible("g7e", "g7"))

    def test_rehearsal_families_are_open(self):
        self.assertTrue(slotFamilyCompatible("g7e", "cpu"))
        self.assertTrue(slotFamilyCompatible(None, "local"))
        self.assertTrue(slotFamilyCompatible("g6", ""))


class CampaignWorkers(unittest.TestCase):
    def test_ada_uses_autothreads(self):
        self.assertFalse(usesCampaignWorkers("g6"))
        self.assertFalse(usesCampaignWorkers("g6e"))
        self.assertFalse(usesCampaignWorkers("g4dn"))
        self.assertFalse(usesCampaignWorkers("g5"))

    def test_blackwell_keeps_preset(self):
        self.assertTrue(usesCampaignWorkers("g7"))
        self.assertTrue(usesCampaignWorkers("g7e"))
        self.assertTrue(usesCampaignWorkers("cpu"))
        self.assertTrue(usesCampaignWorkers(""))


if __name__ == "__main__":
    unittest.main()
