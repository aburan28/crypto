"""Small-group and receipt checks for the rank-query necessary condition."""

import json
import unittest

from rank_query_screen import (
    PANEL,
    half_probability_query_floor,
    rank_probability_ceiling,
    summarize,
)


class RankQueryScreenTests(unittest.TestCase):
    def test_small_group_exhaustive_queries(self) -> None:
        # In Z/7Z the one-summand base {1,2} hits exactly two targets.
        # Every ordered pair of independent uniform targets is enumerated;
        # receiving two rows needs two hits, whether or not the rows span.
        group = range(7)
        base = {1, 2}
        two_hits = sum(a in base and b in base for a in group for b in group)
        self.assertEqual(two_hits, 4)
        upper_num, upper_den = rank_probability_ceiling(2, 2, 2, 7)
        self.assertLessEqual(two_hits * upper_den, upper_num * 49)
        self.assertEqual(rank_probability_ceiling(1, 2, 2, 7), (0, 1))
        # Correlated targets are also covered: the two targets may be equal.
        correlated_hits = sum(a in base for a in group)
        self.assertLessEqual(correlated_hits * upper_den, upper_num * 7)

    def test_half_probability_floor_is_exact(self) -> None:
        for rank in (1, 2, 64, 600):
            for order, multisets in ((7, 2), (7, 70), (2417851639230796216685689, 41013806413)):
                floor = half_probability_query_floor(rank, multisets, order)
                before = rank_probability_ceiling(floor - 1, rank, multisets, order)
                at = rank_probability_ceiling(floor, rank, multisets, order)
                self.assertLess(2 * before[0], before[1])
                self.assertGreaterEqual(2 * at[0], at[1])
                self.assertGreaterEqual(floor, rank)

    def test_receipt_is_deterministic_and_bound_to_support_screen(self) -> None:
        first = summarize()
        self.assertEqual(first, summarize())
        self.assertEqual(len(first["cases"]), 30)
        self.assertEqual({row["curve_a"] for row in first["cases"]}, {0, 1})
        saved = json.loads((PANEL / "rank-query-screen.json").read_text())
        self.assertEqual(saved, first)


if __name__ == "__main__":
    unittest.main()
