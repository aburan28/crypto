"""Exact-count and provenance controls for the paired N83 cap screen."""

import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

import screen


class PairedCapacityScreenTests(unittest.TestCase):
    def test_saved_screen_and_full_signed_pair_population(self):
        row = screen.build()
        self.assertEqual(screen.canonical(row), screen.OUTPUT.read_bytes())
        self.assertEqual(row["signed_pair_spec_count"], 990)
        self.assertEqual(row["unique_signed_size_pairs"], 55)
        self.assertEqual(len(row["size_pair_shape_rows"]), 220)
        self.assertEqual(
            [(s["summands"], s["residual_summands"],
              s["size_pairs_within_current_worker_caps"])
             for s in row["shape_summaries"]],
            [(5, 1, 54), (5, 2, 45), (6, 1, 48), (6, 2, 42)],
        )
        self.assertTrue(all(s["size_pairs_over_domain_clause_cap"] == 0
                            for s in row["shape_summaries"]))
        self.assertFalse(row["sat_search_executed"])

    def test_retained_counters_and_cap_boundary(self):
        row = json.loads(screen.OUTPUT.read_text())
        observed = row["retained_counter_replay"]
        self.assertEqual((observed["K600"]["observed_total_variables"],
                          observed["K600"]["observed_domain_clauses"],
                          observed["K600"]["observed_xor_rows"]),
                         (100472, 415686, 4980))
        self.assertEqual((observed["K1182"]["observed_total_variables"],
                          observed["K1182"]["observed_domain_clauses"],
                          observed["K1182"]["observed_xor_rows"]),
                         (101636, 465738, 4980))
        self.assertTrue(screen.requirements(1182, 8192, 5, 2)
                        ["current_worker_cap_pass"])
        self.assertFalse(screen.requirements(1182, 8192, 6, 2)
                         ["current_worker_cap_pass"])
        self.assertEqual(screen.requirements(1182, 8192, 6, 2)
                         ["total_variables"], 150992)

    def test_changed_retained_event_is_rejected(self):
        with TemporaryDirectory() as directory:
            changed = Path(directory) / "event.jsonl"
            original = screen.STUDY / screen.EVIDENCE["K600"]["event"][0]
            changed.write_bytes(original.read_bytes() + b" ")
            tampered = dict(screen.EVIDENCE["K600"])
            tampered["event"] = (str(changed),
                                 screen.EVIDENCE["K600"]["event"][1])
            # The receipt reader takes study-relative paths; the absolute
            # temporary path also resolves correctly under pathlib's /.
            with self.assertRaisesRegex(ValueError, "K600 event changed"):
                screen.retained_check("K600", tampered)


if __name__ == "__main__":
    unittest.main()
