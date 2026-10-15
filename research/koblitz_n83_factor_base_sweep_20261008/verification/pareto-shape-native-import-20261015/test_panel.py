"""No-launch gates for the frozen N83 source-Pareto import panel."""

import json
import unittest
from unittest import mock

import panel


class ParetoImportPlanTest(unittest.TestCase):
    def test_saved_plan_replays_all_distinct_shape_representatives(self):
        saved = json.loads(panel.PLAN.read_text())
        self.assertEqual(panel.make_plan(), saved)
        self.assertEqual(len(saved["arms"]), 24)
        self.assertEqual(len(saved["objects"]), 18)
        self.assertEqual(len({(arm["conditional_half_K_row_query_floor"],
                               arm["source_domain_clauses"],
                               arm["source_total_variables"])
                              for arm in saved["arms"]}), 24)
        self.assertFalse(saved["solver_search_executed"])
        self.assertIsNone(saved["selected_best_total_runtime"])

    def test_changed_pareto_receipt_is_rejected(self):
        with mock.patch.dict(panel.PINS, {panel.PARETO: "0" * 64}):
            with self.assertRaisesRegex(ValueError, "pinned source changed"):
                panel.make_plan()

    def test_s3_key_rejects_different_bucket_or_object(self):
        row = json.loads(panel.PLAN.read_text())["objects"][0]
        digest = row["compressed_blake3"]
        self.assertEqual(panel.s3_key(row["s3_uri"], digest),
                         f"{panel.PREFIX}{digest}.jsonl.gz")
        for changed in (row["s3_uri"].replace(panel.BUCKET, "other-bucket"),
                        row["s3_uri"].replace(digest, "0" * 64)):
            with self.subTest(uri=changed), self.assertRaisesRegex(ValueError, "S3 object URI"):
                panel.s3_key(changed, digest)


if __name__ == "__main__":
    unittest.main()
