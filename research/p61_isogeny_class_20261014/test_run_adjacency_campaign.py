#!/usr/bin/env python3
"""Unit tests for the p61 adjacency campaign supervisor."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest


SCRIPT = Path(__file__).with_name("run_adjacency_campaign.py")
SPEC = importlib.util.spec_from_file_location("run_adjacency_campaign", SCRIPT)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class AdjacencyCampaignTests(unittest.TestCase):
    def admission(self, root: Path, start: int = 10, end: int = 14) -> Path:
        segment = root / "segment.gz"
        verification = root / "verify.json"
        segment.write_bytes(b"segment")
        verification.write_bytes(b"verification")
        path = root / f"campaign-v1-admission-{start:012d}-{end:012d}.json"
        path.write_text(
            json.dumps(
                {
                    "schema": "ec-hardness/p61-campaign-segment-admission/v1",
                    "status": "pass",
                    "range": {"x_start": start, "x_end": end, "records": (end - start) * 18},
                    "segment": {"path": str(segment), "sha256": MODULE.sha256_file(segment)},
                    "verification_receipt": {
                        "path": str(verification),
                        "sha256": MODULE.sha256_file(verification),
                    },
                }
            ),
            encoding="utf-8",
        )
        return path

    def test_admission_digest_replay(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = self.admission(Path(directory))
            value = MODULE.load_vertex_admission(path, True)
            self.assertEqual(value.records, 72)
            value.segment.write_bytes(b"changed")
            with self.assertRaisesRegex(RuntimeError, "digest replay failed"):
                MODULE.load_vertex_admission(path, True)

    def test_windows_are_complete_and_nonoverlapping(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            admission = MODULE.load_vertex_admission(self.admission(root), True)
            windows = MODULE.windows_for_admissions([admission], root / "out", "p", 32)
            self.assertEqual([(item.offset, item.count) for item in windows], [(0, 32), (32, 32), (64, 8)])
            self.assertEqual(sum(item.count for item in windows), admission.records)
            self.assertEqual(len({item.paths.shard for item in windows}), 3)

    def test_summary_requires_all_seven_edges(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            admission = MODULE.load_vertex_admission(self.admission(root), True)
            window = MODULE.windows_for_admissions([admission], root / "out", "p", 32)[0]
            summary = {
                "source_record_offset": 0,
                "source_record_count": 32,
                "directed_edges": 224,
                "kernel_coefficients": 3680,
                "generation_failures": 0,
                "record_chain_sha256": "a" * 64,
                "degree_counts": [
                    {"degree": 5, "directed_edges": 32},
                    {"degree": 13, "directed_edges": 64},
                    {"degree": 37, "directed_edges": 64},
                    {"degree": 59, "directed_edges": 64},
                ],
            }
            MODULE.validate_summary(summary, window)
            summary["directed_edges"] = 223
            with self.assertRaisesRegex(RuntimeError, "directed-edge count"):
                MODULE.validate_summary(summary, window)

    def test_partial_generation_is_refused(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            paths = MODULE.WindowPaths(
                root / "shard.gz", root / "generate.json", root / "verify.json", root / "admit.json"
            )
            paths.shard.write_bytes(b"partial")
            with self.assertRaisesRegex(RuntimeError, "partial generation evidence"):
                MODULE.exact_generation_state(paths)

    def test_follow_refuses_bounded_plan(self) -> None:
        args = MODULE.parser().parse_args(
            [
                "--producer",
                "/producer",
                "--verifier",
                "/verifier",
                "--vertex-manifest-dir",
                "/manifests",
                "--root",
                "/root",
                "--mount-root",
                "/mount",
                "--producer-source-commit",
                "a" * 40,
                "--verifier-source-commit",
                "b" * 40,
                "--class-number-receipt-sha256",
                "c" * 64,
                "--follow",
                "--max-windows",
                "1",
            ]
        )
        with self.assertRaisesRegex(RuntimeError, "cannot be combined"):
            MODULE.validate_args(args)


if __name__ == "__main__":
    unittest.main()
