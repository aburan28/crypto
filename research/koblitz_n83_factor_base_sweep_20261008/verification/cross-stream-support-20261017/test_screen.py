"""Independent partition and tamper controls for the N83 support census."""

import copy
import itertools
import json
from pathlib import Path
import sys
import unittest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import screen  # noqa: E402


class SupportScreenTests(unittest.TestCase):
    def test_small_multiset_counts_match_explicit_enumeration(self):
        for points in range(4):
            for count in range(5):
                explicit = sum(1 for _ in itertools.combinations_with_replacement(
                    range(points), count))
                self.assertEqual(screen.multiset(points, count), explicit)

    def test_disjoint_pool_partition_is_complete(self):
        for base, residual in itertools.product(range(3), repeat=2):
            for total in range(2, 7):
                partition = sum(screen.support(base, residual, total, large)
                                for large in range(total + 1))
                whole = screen.multiset(
                    (base + residual) * screen.POINTS_PER_SIGNED_ORBIT, total)
                self.assertEqual(partition, whole)

    def test_saved_census_and_hosted_pair_match(self):
        saved = json.loads(screen.SAVED.read_text())
        self.assertEqual(saved, screen.generate())
        self.assertEqual(saved["pair_count"], 648)
        self.assertEqual(saved["current_worker_admissible_pairs"], 558)
        base = ["public_x_sequential", 2026100802, 600]
        pair = {tuple(row["pool"]): row for row in saved["cases"]
                if row["base"] == base and row["pool"][2] == 600}
        gray = pair[("public_x_gray_prefix", 2026100801, 600)]
        hashed = pair[("public_x_hash", 2026100801, 600)]
        self.assertEqual(gray["residual_orbits"], 97)
        self.assertEqual(hashed["residual_orbits"], 600)
        self.assertEqual(gray["shapes"]["5/2"]["unordered_multisets"],
                         "21349930494938690639600")
        self.assertEqual(hashed["shapes"]["5/2"]["unordered_multisets"],
                         "816832271335643493360000")
        plan_path = (screen.STUDY / "verification"
                     / "k600-gray-hash-paired-protocol-20261014" / "plan.json")
        plan = json.loads(plan_path.read_text())
        by_role = {arm["role"]: arm for arm in plan["arms"]}
        self.assertEqual(gray["shapes"]["5/2"]["unordered_multisets"],
                         by_role["gray_pool"]["point_multisets_3_plus_2"])
        self.assertEqual(hashed["shapes"]["5/2"]["unordered_multisets"],
                         by_role["hash_pool"]["point_multisets_3_plus_2"])
        self.assertEqual(saved["admissible_bucket_counts"]["5/2"], {
            "under_1_percent": 432, "1_to_10_percent": 84,
            "10_to_50_percent": 42, "50_to_100_percent": 0,
            "vacuous_at_100_percent": 0,
        })

    def test_corrupt_residual_is_rejected(self):
        overlap = json.loads(screen.OVERLAP.read_text())
        audit = json.loads(screen.AUDIT.read_text())
        comparison = json.loads(screen.COMPARISON.read_text())
        corrupt = copy.deepcopy(overlap)
        corrupt["pairs"][0]["residual_orbits"] += 1
        with self.assertRaisesRegex(ValueError, "residual partition"):
            screen.build(corrupt, audit, comparison)


if __name__ == "__main__":
    unittest.main()
