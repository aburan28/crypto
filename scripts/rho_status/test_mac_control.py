#!/usr/bin/env python3
"""Protect the public Mac snapshot boundary from private worker fields."""

import unittest

from mac_control import CAMPAIGN_ID, PUBLIC_SCHEMA, SOURCE_SCHEMA, snapshot


class MacControlTests(unittest.TestCase):
    def worker(self, run="1" * 64, walk="a" * 64, dps=396):
        return {
            "schema": SOURCE_SCHEMA,
            "sequence": 16767,
            "updatedAt": "2026-10-09T04:24:13Z",
            "runIdentity": run,
            "walkIdentity": walk,
            "manifestKey": "private/manifest-key",
            "state": {"key": "private/state-key", "bytes": 500},
            "knownScalar": "private scalar",
            "cumulative": {"dpRecords": dps, "walkUpdates": 140620332532},
        }

    def test_only_allowlisted_counts_leave_private_manifest(self):
        public = snapshot([("m4pro-01", self.worker())], generated_at="2026-10-09T04:25:00Z")
        self.assertEqual(public["schema"], PUBLIC_SCHEMA)
        self.assertEqual(public["campaign_id"], CAMPAIGN_ID)
        self.assertEqual(public["source_state"], "available")
        self.assertEqual(public["dp_records"], 396)
        self.assertEqual(public["worker_count"], 1)
        self.assertEqual(set(public["workers"][0]),
                         {"worker_id", "updated_at", "sequence", "dp_records", "walk_updates"})
        self.assertNotIn("private", str(public))
        self.assertNotIn("knownScalar", str(public))
        self.assertNotIn("runIdentity", str(public))

    def test_never_sum_distinct_walks_or_reused_run_identity(self):
        for other in (self.worker(run="2" * 64, walk="b" * 64),
                      self.worker(run="1" * 64)):
            public = snapshot([("m4pro-01", self.worker()), ("m4max-02", other)],
                              generated_at="2026-10-09T04:25:00Z")
            self.assertIn(public["source_state"], ("mixed_walks", "ambiguous_runs"))
            self.assertIsNone(public["dp_records"])
            self.assertEqual(public["worker_count"], 2)

    def test_unavailable_and_bad_worker_are_visible_as_missing(self):
        down = snapshot([], generated_at="2026-10-09T04:25:00Z", source_available=False)
        self.assertEqual(down["source_state"], "unavailable")
        self.assertIsNone(down["dp_records"])
        bad = self.worker()
        bad["cumulative"]["dpRecords"] = -1
        partial = snapshot([("m4pro-01", self.worker()), ("m4max-02", bad)],
                           generated_at="2026-10-09T04:25:00Z")
        self.assertEqual(partial["source_state"], "partial")
        self.assertEqual(partial["read_errors"], 1)
        self.assertEqual(partial["worker_count"], 1)
        all_bad = snapshot([("m4max-02", bad)], generated_at="2026-10-09T04:25:00Z")
        self.assertIsNone(all_bad["dp_records"])


if __name__ == "__main__":
    unittest.main()
