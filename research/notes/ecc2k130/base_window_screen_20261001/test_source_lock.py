#!/usr/bin/env python3
"""Source-only regression tests; no held-out Q or measured arm is run here."""
from __future__ import annotations

import math
from pathlib import Path
import statistics
import unittest

from analyze_screen import analyze
from prepare_inputs import candidate, canon_bytes, sha_bytes
from run_screen import config, schedule
from verify_screen import HERE, isolation_command_ok, paired


def child(total: float) -> dict:
    return {
        "exit_code": 0,
        "stopped_for": None,
        "child_user_cpu_seconds": total,
        "child_system_cpu_seconds": 0.0,
    }


def synthetic(values: list[list[float]], uncontended: bool = True) -> tuple[dict, dict]:
    cfg = config()
    assert len(values) == cfg["blocks"] and all(len(row) == 4 for row in values)
    report = {
        "schema": "ecc2k130-base-window-screen-run-v1",
        "status": "PASS",
        "config": cfg,
        "plan": [{"block": b, "arm": arm} for b, arm in schedule(cfg)],
        "runs": [],
    }
    for block, arm in schedule(cfg):
        if arm == "rho":
            stages = {"rho": child(1.0)}
        else:
            value = values[block][int(arm[1])]
            assert value > 0.1
            stages = {"generator": child(0.1), "compact": child(value - 0.1)}
        report["runs"].append({"block": block, "arm": arm, "children": stages})
    aa = {str(w): paired([1.0] * 5, cfg["log_t_critical_df4"]) for w in range(4)}
    fixed = {
        str(w): paired([values[b][w] for b in range(5)], cfg["log_t_critical_df4"])
        for w in range(4)
    }
    replay = {
        "schema": "ecc2k130-base-window-screen-replay-v1",
        "status": "PASS",
        "uncontended": uncontended,
        "aa_valid": {str(w): True for w in range(4)},
        "timing_eligible": uncontended,
        "aa": aa,
        "fixed_window_over_rho": fixed,
        "independent_rank_replays": 40,
        "verified_compact_target_logs": 40960,
        "verified_rho_target_logs": 5120,
    }
    return report, replay


class SourceLockTests(unittest.TestCase):
    def test_schedule_rotates_every_same_q_nine_arm_block(self) -> None:
        cfg = config()
        plan = schedule(cfg)
        self.assertEqual(len(plan), 45)
        original = tuple(cfg["arm_order_before_rotation"])
        for block in range(5):
            observed = tuple(arm for b, arm in plan if b == block)
            self.assertEqual(observed, original[block:] + original[:block])
            self.assertEqual(sum(arm == "rho" for arm in observed), 1)

    def test_candidate_domain_and_inventory_newline_are_fixed(self) -> None:
        domain = config()["target_domain"]
        self.assertEqual(candidate(domain, 0, 0), (
            19001090012,
            "dfddc9c8a52c402f87bea8ab0dd1c973b9f91da8263910834bfb69f3915b8407",
        ))
        self.assertEqual(candidate(domain, 4, 13), (
            497953983002,
            "13dfb086d045a8803bd15ca01ccead76c2eb7c7de2910f45a65335a6c636d251",
        ))
        inventory = [{"path": "a.points.jsonl", "rows": 1, "sha256": "ab"}]
        self.assertTrue(canon_bytes(inventory).endswith(b"\n"))
        self.assertEqual(
            sha_bytes(canon_bytes(inventory)),
            "5c38a6c3f7b1a67e39ef2e985d829cdeb9e45d806c632786b9b991de4dde9a95",
        )

    def test_no_window_opportunity_requires_every_h_above_one(self) -> None:
        report, replay = synthetic([[1.2, 1.4, 1.6, 1.8] for _ in range(5)])
        result = analyze(report, replay)
        self.assertEqual(result["decision"], "NO_WINDOW_OPPORTUNITY")
        self.assertEqual(result["hindsight_lower_envelope"]["values"], [1.2] * 5)
        self.assertGreater(result["hindsight_lower_envelope"]["interval_95pct"][0], 1)
        self.assertIsNone(result["method_speedup_claim"])

    def test_one_predeclared_fixed_window_is_only_a_candidate(self) -> None:
        report, replay = synthetic([[0.8, 1.4, 1.6, 1.8] for _ in range(5)])
        result = analyze(report, replay)
        self.assertEqual(result["decision"], "FIXED_WINDOW_CROSSOVER_CANDIDATE")
        self.assertEqual(result["crossing_fixed_windows_requiring_confirmation"], [0])
        self.assertIsNone(result["second_host_and_new_Q_confirmation"])
        self.assertIsNone(result["common_operation_equivalent_S"])

    def test_hindsight_wins_do_not_create_a_selector(self) -> None:
        values = [[1.5] * 4 for _ in range(5)]
        for block, winner in enumerate((0, 0, 1, 1, 2)):
            values[block][winner] = 0.8
        report, replay = synthetic(values)
        result = analyze(report, replay)
        self.assertEqual(result["decision"], "HINDSIGHT_VARIATION_ONLY")
        self.assertEqual(result["hindsight_lower_envelope"]["median"], 0.8)
        self.assertEqual(result["crossing_fixed_windows_requiring_confirmation"], [])
        self.assertIsNone(result["selection_cpu_charged"])

    def test_interval_overlap_is_inconclusive(self) -> None:
        values = [[1.05, 1.4, 1.6, 1.8] for _ in range(5)]
        values[0][0] = 0.95
        report, replay = synthetic(values)
        result = analyze(report, replay)
        self.assertEqual(result["decision"], "INCONCLUSIVE")
        self.assertGreater(statistics.median(result["hindsight_lower_envelope"]["values"]), 1)
        self.assertLess(result["hindsight_lower_envelope"]["interval_95pct"][0], 1)

    def test_isolation_failure_censors_even_if_ratio_looks_good(self) -> None:
        report, replay = synthetic([[0.8, 1.4, 1.6, 1.8] for _ in range(5)],
                                   uncontended=False)
        result = analyze(report, replay)
        self.assertEqual(result["decision"], "CENSORED")
        self.assertFalse(result["timing_eligible"])
        self.assertNotIn("hindsight_lower_envelope", result)

    def test_failure_and_wrong_cpu_receipt_are_not_wins(self) -> None:
        report, replay = synthetic([[0.8, 1.4, 1.6, 1.8] for _ in range(5)])
        report["status"] = "CENSORED"
        replay["status"] = "PRODUCER_FAILURE"
        self.assertEqual(analyze(report, replay)["decision"], "CENSORED")
        report, replay = synthetic([[0.8, 1.4, 1.6, 1.8] for _ in range(5)])
        report["runs"][0]["children"]["generator"]["child_user_cpu_seconds"] = math.nan
        with self.assertRaises(AssertionError):
            analyze(report, replay)

    def test_isolation_receipt_must_name_this_run_and_binaries(self) -> None:
        build_dir, run_dir = Path("/tmp/window-build"), Path("/tmp/window-run")
        binary_dir = build_dir / "target/release/examples"
        build = {"binaries": {
            name: {"path": str(binary_dir / name)} for name in (
                "koblitz_base_window", "koblitz_orbit_dlp_s3_batch",
                "koblitz_rho_batch_ks_v3",
            )
        }}
        command = [
            "python3", str(HERE / "run_screen.py"),
            "--source-root", str(build_dir / "source"),
            "--generator", build["binaries"]["koblitz_base_window"]["path"],
            "--compact", build["binaries"]["koblitz_orbit_dlp_s3_batch"]["path"],
            "--rho", build["binaries"]["koblitz_rho_batch_ks_v3"]["path"],
            "--materialization", str(build_dir / "materialization.json"),
            "--build-receipt", str(build_dir / "BUILD_RECEIPT.json"),
            "--out", str(run_dir), "--cpu", "3",
        ]
        self.assertTrue(isolation_command_ok({"command": command}, run_dir, 3, build, False))
        command[-1] = "4"
        self.assertFalse(isolation_command_ok({"command": command}, run_dir, 3, build, False))
        command[-1] = "3"
        command[command.index("--compact") + 1] = str(binary_dir / "unreviewed")
        self.assertFalse(isolation_command_ok({"command": command}, run_dir, 3, build, False))


if __name__ == "__main__":
    unittest.main()
