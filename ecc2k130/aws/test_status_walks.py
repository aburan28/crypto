"""The dashboard's iteration count must use each slot's own walk count.

A checkpoint stores the per-walk iteration base.  Ada slots let the client
size the grid instead of taking campaign.json's 385,024 x 16 (worker
usesCampaignWorkers), so they run far fewer walks and their base climbs
proportionally faster.  Counting every slot at the Blackwell preset turned
that into tens of times the group operations actually walked, while the
points those slots produced stayed exactly the same.

No type hints, camelCase identifiers (project convention).
"""
import contextlib
import io
import json
import time
import unittest

import status
from worker import BANNER_RE

PRESET = 385024 * 16
ADA_WALKS = 14848 * 16          # an L4-sized grid: 25.9x fewer walks than the preset


def slot(n, ckptIter, walks=None, rate=0.0, dp=0):
    it = {"slot": n, "ckptIter": ckptIter, "rate": rate, "dpUploaded": dp, "state": "active",
          "leaseUntil": int(time.time()) + 180}
    if walks is not None:
        it["walks"] = walks
    return it


def summary(items, walks=PRESET):
    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        status.report(items, walks, True)
    return json.loads(out.getvalue())


class Walks(unittest.TestCase):
    def test_reported_walks_beat_the_preset(self):
        items = [slot(0, 1000, PRESET), slot(1, 1000, ADA_WALKS)]
        self.assertEqual(summary(items)["checkpointedIterations"],
                         1000 * PRESET + 1000 * ADA_WALKS)

    def test_ada_slot_is_not_counted_at_the_blackwell_preset(self):
        honest = summary([slot(1, 1000, ADA_WALKS)])["checkpointedIterations"]
        inflated = 1000 * PRESET
        self.assertEqual(honest, 1000 * ADA_WALKS)
        self.assertGreater(inflated / honest, 25)

    def test_same_group_operations_whatever_the_grid(self):
        # Two slots that walked the same number of group operations must count
        # the same, however differently their iteration bases climbed.
        small = ADA_WALKS
        big = summary([slot(0, small, PRESET)])["checkpointedIterations"]
        self.assertEqual(big, summary([slot(1, PRESET, small)])["checkpointedIterations"])

    def test_a_slot_without_walks_is_assumed_and_said_so(self):
        got = summary([slot(0, 1000), slot(1, 1000, ADA_WALKS)])
        self.assertEqual(got["walksAssumedSlots"], 1)
        self.assertEqual(got["assumedWalksPerSlot"], PRESET)
        self.assertEqual(got["checkpointedIterationsAssumed"], 1000 * PRESET)
        self.assertEqual(got["checkpointedIterations"], 1000 * PRESET + 1000 * ADA_WALKS)

    def test_idle_slot_is_neither_counted_nor_flagged(self):
        got = summary([slot(0, -1)])
        self.assertEqual(got["checkpointedIterations"], 0)
        self.assertEqual(got["walksAssumedSlots"], 0)

    def test_eta_and_fraction_follow_the_honest_count(self):
        got = summary([slot(0, 1000, ADA_WALKS, rate=1e9)])
        iters = 1000 * ADA_WALKS
        self.assertAlmostEqual(got["fractionOfExpected"], iters / status.EXPECTED_ITERS)
        self.assertAlmostEqual(got["etaSecondsAtCurrentRate"],
                               (status.EXPECTED_ITERS - iters) / 1e9)


class Banner(unittest.TestCase):
    """The worker learns the walk count from the client's opening line."""

    def test_packed_banner(self):
        line = ("backend packed-cuda: 385024 threads x 16 slots x 1 lanes = 6160384 walks, "
                "dp weight 34, 1024 steps per launch")
        self.assertEqual(int(BANNER_RE.search(line).group(1)), 6160384)

    def test_bitsliced_banner(self):
        line = ("backend cuda: 14848 threads x 16 slots x 64 lanes = 15204352 walks, "
                "dp weight 32, 1024 steps per launch")
        self.assertEqual(int(BANNER_RE.search(line).group(1)), 15204352)

    def test_progress_line_is_not_a_banner(self):
        self.assertIsNone(BANNER_RE.search(
            "  120.0 s  14110.000 M it/s  1693200000000 iterations  9 dp  9 stored  0 dropped"))


if __name__ == "__main__":
    unittest.main()
