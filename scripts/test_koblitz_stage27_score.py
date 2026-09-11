#!/usr/bin/env python3

from __future__ import annotations

import json
from pathlib import Path
import sys
import tempfile
import unittest


sys.path.insert(0, str(Path(__file__).resolve().parent))
import score_koblitz_stage27_direct_mitm as score


CELLS = {
    "n31-l5-m3-standard-a1-f0",
    "n31-l5-m3-ggmp-a0-f0",
    "n41-l5-m3-standard-a1-f0",
    "n59-l9-m3-standard-a1-f0",
}


class Stage27ScoreTests(unittest.TestCase):
    def test_workflow_requires_exact_run_and_cells(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "workflow.json"
            jobs = [
                {
                    "name": f"direct-mitm-cell ({cell})",
                    "databaseId": index + 1,
                    "status": "completed",
                    "conclusion": "success",
                    "startedAt": "2026-09-11T18:33:34Z",
                    "completedAt": "2026-09-11T18:38:41Z",
                    "url": f"https://example.test/{index}",
                }
                for index, cell in enumerate(sorted(CELLS))
            ]
            value = {
                "databaseId": score.RUN_ID,
                "headSha": score.RUN_COMMIT,
                "status": "completed",
                "conclusion": "success",
                "createdAt": "2026-09-11T18:33:20Z",
                "updatedAt": "2026-09-11T18:38:42Z",
                "url": "https://example.test/run",
                "jobs": jobs,
            }
            path.write_text(json.dumps(value))
            receipt = score.workflow_receipt(path, CELLS)
            self.assertEqual(receipt["workflow_wall_seconds"], 322)
            self.assertEqual(receipt["parallel_cell_job_span_seconds"], 307)
            value["databaseId"] += 1
            path.write_text(json.dumps(value))
            with self.assertRaises(score.Stage27ScoreError):
                score.workflow_receipt(path, CELLS)

    def test_artifacts_require_exact_four_cell_bijection(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "artifacts.json"
            artifacts = [
                {
                    "id": index + 1,
                    "name": f"koblitz-stage27-{cell}-{score.RUN_ID}",
                    "size_in_bytes": 100,
                    "digest": "sha256:" + f"{index + 1:064x}",
                    "expires_at": "2026-12-10T00:00:00Z",
                    "expired": False,
                }
                for index, cell in enumerate(sorted(CELLS))
            ]
            path.write_text(json.dumps({"total_count": 4, "artifacts": artifacts}))
            self.assertEqual(set(score.artifact_receipt(path, CELLS)), CELLS)
            artifacts[0]["expired"] = True
            path.write_text(json.dumps({"total_count": 4, "artifacts": artifacts}))
            with self.assertRaises(score.Stage27ScoreError):
                score.artifact_receipt(path, CELLS)


if __name__ == "__main__":
    unittest.main()
