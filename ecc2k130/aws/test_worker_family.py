"""Ada / Blackwell slot pinning and auto thread count.

g6/g6e must not resume a 385,024-worker Blackwell checkpoint (that retires
the run id) and must not pass --threads from campaign.json.
"""
import unittest

from worker import (gpuFamily, slotFamilyCompatible, usesCampaignWorkers,
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

    def test_blackwell_keeps_preset(self):
        self.assertTrue(usesCampaignWorkers("g7"))
        self.assertTrue(usesCampaignWorkers("g7e"))
        self.assertTrue(usesCampaignWorkers("cpu"))
        self.assertTrue(usesCampaignWorkers(""))


if __name__ == "__main__":
    unittest.main()
