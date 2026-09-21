"""Ada / Blackwell slot pinning and auto thread count.

g6/g6e must not resume a 385,024-worker Blackwell checkpoint (that retires
the run id) and must not pass --threads from campaign.json.
"""
import unittest
from unittest import mock

from worker import (gpuFamily, slotFamilyCompatible, usesCampaignWorkers,
                    visibleGpuCount, runAllGpus, BLACKWELL_FAMILIES, ADA_FAMILIES,
                    AUTO_FAMILIES, LOCAL_FAMILIES)


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
        self.assertEqual(gpuFamily("NVIDIA B200"), "b200")
        self.assertEqual(gpuFamily("NVIDIA B200 MIG 24GB"), "mig")
        self.assertEqual(gpuFamily("NVIDIA H100 80GB HBM3"), "h100")

    def test_l40s_not_classified_as_l4(self):
        self.assertEqual(gpuFamily("NVIDIA L40S"), "g6e")

    def test_mig_wins_over_parent_sku(self):
        self.assertEqual(gpuFamily("NVIDIA RTX PRO 6000 MIG 24GB"), "mig")


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

    def test_b200_skips_legacy_blackwell_slots(self):
        # Same trap as Ada: 385,024-worker checkpoints are unloadable on an
        # auto-sized B200/MIG grid, so those workers mint new slots.
        for family in ("b200", "mig", "h100"):
            self.assertFalse(slotFamilyCompatible(None, family))
            self.assertFalse(slotFamilyCompatible("g7e", family))
            self.assertTrue(slotFamilyCompatible(family, family))

    def test_unclassified_names_pin_to_their_own_slots(self):
        # An unmatched name auto-sizes, so it is not the family-less rehearsal
        # claimant: it must not resume the untagged Blackwell corpus or another
        # model's grid (exit 6), and Blackwell must not resume its slots.
        family = gpuFamily("NVIDIA GeForce RTX 4090")
        self.assertEqual(family, "nvidia-geforce-rtx-4090")
        self.assertNotIn(family, LOCAL_FAMILIES)
        self.assertFalse(usesCampaignWorkers(family))
        self.assertFalse(slotFamilyCompatible(None, family))
        self.assertFalse(slotFamilyCompatible("", family))
        self.assertFalse(slotFamilyCompatible("g7e", family))
        self.assertFalse(slotFamilyCompatible(gpuFamily("NVIDIA A40"), family))
        self.assertFalse(slotFamilyCompatible(family, "g7e"))
        self.assertTrue(slotFamilyCompatible(family, family))
        self.assertEqual(gpuFamily(""), "")

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

    def test_unknown_and_datacenter_auto_size(self):
        # 385,024 workers is the 188-SM RTX PRO 6000 grid. B200, MIG slices and
        # unclassified names omit --threads so autoThreads fills the device.
        self.assertFalse(usesCampaignWorkers("cpu"))
        self.assertFalse(usesCampaignWorkers(""))
        self.assertFalse(usesCampaignWorkers("b200"))
        self.assertFalse(usesCampaignWorkers("mig"))
        self.assertFalse(usesCampaignWorkers("h100"))
        self.assertFalse(usesCampaignWorkers(gpuFamily("NVIDIA B200")))
        self.assertFalse(usesCampaignWorkers(gpuFamily("NVIDIA B200 MIG 24GB")))
        self.assertFalse(BLACKWELL_FAMILIES & AUTO_FAMILIES)


class VisibleGpus(unittest.TestCase):
    def test_parses_nvidia_smi_list(self):
        listing = ("GPU 0: NVIDIA B200 (UUID: GPU-aaa)\n"
                   "GPU 1: NVIDIA B200 (UUID: GPU-bbb)\n"
                   "GPU 7: NVIDIA B200 (UUID: GPU-hhh)\n")
        result = mock.Mock(returncode=0, stdout=listing)
        with mock.patch("worker.subprocess.run", return_value=result):
            self.assertEqual(visibleGpuCount(), 3)

    def test_missing_binary_is_zero(self):
        with mock.patch("worker.subprocess.run", side_effect=OSError("no smi")):
            self.assertEqual(visibleGpuCount(), 0)


class AllGpus(unittest.TestCase):
    def test_single_device_runs_one_worker(self):
        with mock.patch("worker.visibleGpuCount", return_value=1), \
             mock.patch("worker.Worker") as workerCls:
            workerCls.return_value.run.return_value = 0
            self.assertEqual(runAllGpus(), 0)
            workerCls.return_value.run.assert_called_once()

    def test_eight_devices_spawn_eight_children(self):
        children = [mock.Mock() for _ in range(8)]
        for child in children:
            child.wait.return_value = 0
            child.poll.return_value = 0
        with mock.patch("worker.visibleGpuCount", return_value=8), \
             mock.patch("worker.subprocess.Popen", side_effect=children) as popen, \
             mock.patch("worker.signal.signal"):
            self.assertEqual(runAllGpus(), 0)
            self.assertEqual(popen.call_count, 8)
            gpus = [call.kwargs["env"]["ECC_GPU"] for call in popen.call_args_list]
            self.assertEqual(gpus, [str(i) for i in range(8)])
            self.assertTrue(all("ECC_ALL_GPUS" not in call.kwargs["env"]
                                for call in popen.call_args_list))


if __name__ == "__main__":
    unittest.main()
