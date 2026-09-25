"""The campaign's guard and the tools that replay its trails agree.

A walk that goes campaign.json's maxIters steps without a report is restarted
and its trail discarded.  The guard is checked every ECC_GUARD_PERIOD steps,
so a reported trail can run ECC_GUARD_PERIOD - 1 steps past it, and
build/witness and build/trailforest refuse anything longer than
ECC_REPLAY_MAX_ITERS by default: include/kernel.h has to follow the campaign.
For the sigma walk, which has no fruitless cycles, the guard only ever cuts
honest trails, so its value is pinned by what it discards
(../benchmarks/max-iters/).
"""
import json
import math
import re
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CAMPAIGN = json.loads((ROOT / "aws" / "campaign.json").read_text())
KERNEL = (ROOT / "include" / "kernel.h").read_text()


def define(name):
    m = re.search(r"^#define %s (.+?)\s*$" % name, KERNEL, re.M)
    if not m:
        raise AssertionError("include/kernel.h does not define %s" % name)
    return m.group(1)


def reportRate(curve, weight):
    """P(HW(x) <= weight) for a point of the order-l subgroup, whose x has
    trace zero and so even normal-basis weight."""
    return sum(math.comb(curve, k) for k in range(0, weight + 1, 2)) / 2 ** (curve - 1)


def discarded(cap, theta):
    """Share of trail steps a guard at `cap` discards when trail lengths are
    geometric on {0, 1, ...}: a trail longer than `cap` is cut after `cap`."""
    log = math.log1p(-theta)
    ac = math.exp(cap * log)
    cut = cap * ac * (1.0 - theta)
    kept = (1.0 - theta) / theta * (1.0 - ac * (1.0 + cap * theta))
    return cut / (kept + cut)


class CampaignGuard(unittest.TestCase):
    def test_header_follows_the_campaign(self):
        self.assertEqual(define("ECC_CAMPAIGN_MAX_ITERS"), "%dull" % CAMPAIGN["maxIters"])
        self.assertEqual(define("ECC_REPLAY_MAX_ITERS"),
                         "(ECC_CAMPAIGN_MAX_ITERS + ECC_GUARD_PERIOD - 1)")

    def test_replay_tools_default_to_the_replay_bound(self):
        for tool in ("witness.cpp", "trailforest.cpp"):
            src = (ROOT / "src" / tool).read_text()
            self.assertIn("unsigned long long maxIters = ECC_REPLAY_MAX_ITERS;", src, tool)

    def test_report_rate_is_the_fleets(self):
        # benchmarks/dp-interval: one report per 2^28.41 iterations at weight 32
        self.assertAlmostEqual(math.log2(reportRate(131, 32)), -28.41, delta=0.005)

    def test_sigma_guard_keeps_the_long_trails(self):
        if CAMPAIGN.get("walk", "sigma") != "sigma":
            self.skipTest("a walk with fruitless cycles needs its own guard")
        theta = reportRate(CAMPAIGN["curve"], CAMPAIGN["dpWeight"])
        self.assertLess(discarded(CAMPAIGN["maxIters"], theta), 1e-4)
        self.assertAlmostEqual(discarded(1 << 30, theta), 0.156, delta=0.001)


if __name__ == "__main__":
    unittest.main()
