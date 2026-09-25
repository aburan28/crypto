#!/usr/bin/env python3
"""Fail-closed controls for the historical/current replay separation."""
from __future__ import annotations

import copy
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import replay


class ReplayControls(unittest.TestCase):
    def test_historical_reference_mismatch_is_rejected(self) -> None:
        data = copy.deepcopy(replay.frozen_input())
        data["reference_files"]["koblitz_builder"]["sha256"] = "0" * 64
        with self.assertRaisesRegex(AssertionError, "koblitz_builder"):
            replay.historical_replay(data)

    def test_current_width_guard_and_unsupported_status_are_required(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            builder_path = root / replay.HISTORICAL_BUILDER
            builder_path.parent.mkdir(parents=True)
            current_builder = (replay.ROOT / replay.HISTORICAL_BUILDER).read_text()
            frontend_path = root / "src/cryptanalysis/koblitz_index_calculus.rs"
            current_frontend = (replay.ROOT / frontend_path.relative_to(root)).read_text()
            frontend_path.write_text(current_frontend)
            builder_path.write_text(current_builder.replace(
                "if n_vars > MAX_VARS {", "if n_vars >= MAX_VARS {", 1
            ))
            with patch.object(replay, "ROOT", root):
                with self.assertRaises(AssertionError):
                    replay.current_contract()

            builder_path.write_text(current_builder)
            frontend_path.write_text(current_frontend.replace(
                "let unsupported = || SolveStats {\n        exhausted: true,\n        unsupported: true,",
                "let unsupported = || SolveStats {\n        exhausted: true,\n        unsupported: false,", 1
            ))
            with patch.object(replay, "ROOT", root):
                with self.assertRaises(AssertionError):
                    replay.current_contract()

    def test_reference_path_traversal_is_rejected(self) -> None:
        for path in ("../outside", "/absolute", "a/../../outside"):
            with self.subTest(path=path):
                with self.assertRaises(AssertionError):
                    replay.safe_reference(path)


if __name__ == "__main__":
    unittest.main()
