"""The worker learns each slot's walk count from the client's opening line.

The dashboard counts a slot's iterations as its checkpointed per-walk base
times the walks it reports, so a banner the worker misreads turns into a
wrong headline.  The dashboard's own walk-count tests moved to Rust with it
(`ecc2k130_status` in the crate).

No type hints, camelCase identifiers (project convention).
"""
import unittest

from worker import BANNER_RE


class Banner(unittest.TestCase):
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
