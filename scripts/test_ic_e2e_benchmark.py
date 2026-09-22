"""Tests for scripts/ic_e2e_benchmark.py: gates fail closed, freeze never overwrites."""
from __future__ import annotations

import copy
import json
import math
import os
import stat
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import ic_e2e_benchmark as bench  # noqa: E402

PARAMS = "docs/ic/params/k0n-test.json"
R = 1_009_001  # any odd "subgroup order" works for the arithmetic; not a real curve
N = 31


def synthetic_report(*, targets: int = 4, rho_seconds: float = 2.0, ic_precompute: float = 0.5, descent: float = 0.1) -> dict:
    """A minimal ic workflow report with the fields the checker reads."""
    expected_steps = math.sqrt(math.pi * R / 2) / math.sqrt(2 * N)
    per_target = int(expected_steps)
    items = [
        {"index": i, "expected": str(1000 + i) if i % 2 == 0 else "not_constructed",
         "recovered": str(1000 + i) if i % 2 == 0 else str(77 + i), "verified": True,
         "descent_trials": 3, "elapsed_seconds": descent / targets}
        for i in range(targets)
    ]
    rows = [{"index": i, "verified": True, "iterations": per_target + 1, "walk_group_additions": per_target,
             "seconds": rho_seconds / targets} for i in range(targets)]
    ic_whole = ic_precompute + descent
    return {
        "schema_version": 1, "operation": "workflow", "status": "complete", "failure": None,
        "name": "test-rung", "degree": N, "params_digest": "ab" * 32,
        "factor_base": {"points": 200, "columns": 35, "pair_table_stored_pairs": 20100,
                        "pair_table_tier": "compact"},
        "solutions": {"count": targets, "verified": targets, "items": items},
        "stages": [
            {"stage": "select", "status": "complete", "ran": True},
            {"stage": "collect", "status": "complete", "ran": True, "trials_total": 128,
             "summands_scanned_total": 25600, "relations_total": 78},
            {"stage": "logs", "status": "complete", "ran": True},
            {"stage": "solve", "status": "complete", "ran": True},
            {"stage": "baseline", "status": "complete", "ran": True, "vs_rho": {
                "n": N, "subgroup_order": str(R), "targets": targets, "claim_boundary": "synthetic_known_answer",
                "ic": {"descent_seconds_total": descent, "descent_seconds_per_target": descent / targets,
                       "descent_trials_total": 3 * targets, "precompute_seconds": ic_precompute,
                       "pair_table_seconds": 0.0, "amortised_seconds_per_target": ic_whole / targets,
                       "verified": targets},
                "rho": {"seconds_total": rho_seconds, "seconds_per_target": rho_seconds / targets,
                        "iterations_total": (per_target + 1) * targets, "verified": targets},
                "ratio": {"charged": (rho_seconds / targets) / (descent / targets), "amortised": (rho_seconds / targets) / (ic_whole / targets)},
                "verdict": {"charged_crossover": True, "amortised_crossover": rho_seconds > ic_whole,
                            "whole_process_crossover": rho_seconds > ic_whole},
                "targets_detail": rows,
            }},
        ],
    }


def write_run(root: Path, report: dict, *, exit_code: int = 0, params: str = PARAMS) -> Path:
    out = root / "run"
    rung = out / Path(params).stem
    rung.mkdir(parents=True)
    (rung / bench.REPORT_FILE).write_text(json.dumps(report))
    (rung / bench.STDERR_FILE).write_text("")
    manifest = {"schema_version": 1, "git": {"commit": "deadbeef", "dirty": False},
                "ic_binary": {"path": "ic", "sha256": "00"}, "host": {"cpu_count": 4},
                "rungs": [{"name": Path(params).stem, "params": params, "exit_code": exit_code,
                           "report": f"{Path(params).stem}/{bench.REPORT_FILE}", "stderr": f"{Path(params).stem}/{bench.STDERR_FILE}"}]}
    (out / bench.MANIFEST_FILE).write_text(json.dumps(manifest))
    return out


class MeasureTests(unittest.TestCase):
    def test_counters_and_wall(self) -> None:
        m = bench.measure(synthetic_report())
        self.assertEqual(m["counters"]["collection_trials"], 128)
        self.assertEqual(m["counters"]["rho_group_additions"], 4 * int(math.sqrt(math.pi * R / 2) / math.sqrt(2 * N)))
        self.assertAlmostEqual(m["wall"]["whole_process_ratio"], 2.0 / 0.6, places=9)
        self.assertTrue(m["wall"]["whole_process_crossover"])
        self.assertIsNone(m["ic_S"])
        self.assertAlmostEqual(m["rho_S"]["measured_over_expected"], 1.0, delta=0.01)

    def test_unverified_descent_is_refused(self) -> None:
        rep = synthetic_report()
        rep["solutions"]["items"][1]["verified"] = False
        with self.assertRaisesRegex(bench.CheckFailure, "target 1 not verified by the descent"):
            bench.measure(rep)

    def test_recovered_must_equal_planted_scalar(self) -> None:
        rep = synthetic_report()
        rep["solutions"]["items"][0]["recovered"] = "999"
        with self.assertRaisesRegex(bench.CheckFailure, "planted"):
            bench.measure(rep)

    def test_unverified_rho_is_refused(self) -> None:
        rep = synthetic_report()
        vs = rep["stages"][4]["vs_rho"]
        vs["targets_detail"][2]["verified"] = False
        with self.assertRaisesRegex(bench.CheckFailure, "target 2 not verified by rho"):
            bench.measure(rep)

    def test_incomplete_status_is_refused(self) -> None:
        rep = synthetic_report()
        rep["status"] = "failed"
        rep["failure"] = "3 of 4 targets unsolved"
        with self.assertRaisesRegex(bench.CheckFailure, "status 'failed'"):
            bench.measure(rep)

    def test_missing_baseline_is_refused(self) -> None:
        rep = synthetic_report()
        rep["stages"] = [s for s in rep["stages"] if s["stage"] != "baseline"]
        with self.assertRaisesRegex(bench.CheckFailure, "baseline.rho = true"):
            bench.measure(rep)


class CheckRungTests(unittest.TestCase):
    def setUp(self) -> None:
        self.ref = bench.measure(synthetic_report())

    def test_identical_run_passes(self) -> None:
        r = bench.check_rung(PARAMS, self.ref, bench.measure(synthetic_report()), 0.5, 3.0)
        self.assertTrue(r["ok"], r["problems"])
        self.assertTrue(r["counters_identical"])

    def test_counter_drift_fails_even_when_faster(self) -> None:
        rep = synthetic_report(ic_precompute=0.1)  # faster end to end …
        rep["stages"][1]["trials_total"] = 64  # … but the algorithm changed
        r = bench.check_rung(PARAMS, self.ref, bench.measure(rep), 0.5, 3.0)
        self.assertFalse(r["ok"])
        self.assertIn("collection_trials", r["counter_drift"])
        self.assertEqual(r["counter_drift"]["collection_trials"]["ratio"], 0.5)
        self.assertTrue(any("freeze a new reference" in p for p in r["problems"]))

    def test_wall_ratio_within_tolerance_passes(self) -> None:
        # ratio 2/0.6 = 3.33 frozen; 2/0.9 = 2.22 is a 33% fall, inside 0.5
        r = bench.check_rung(PARAMS, self.ref, bench.measure(synthetic_report(ic_precompute=0.8)), 0.5, 3.0)
        self.assertTrue(r["ok"], r["problems"])

    def test_wall_ratio_beyond_tolerance_fails(self) -> None:
        # 2/1.6 = 1.25 is a 62% fall
        r = bench.check_rung(PARAMS, self.ref, bench.measure(synthetic_report(ic_precompute=1.5)), 0.5, 3.0)
        self.assertFalse(r["ok"])
        self.assertTrue(any("wall ratio" in p for p in r["problems"]))

    def test_losing_the_e2e_crossover_fails(self) -> None:
        # frozen crosses (2.0 > 0.6); now 2.0 < 2.4 — even at a permissive tolerance
        r = bench.check_rung(PARAMS, self.ref, bench.measure(synthetic_report(ic_precompute=2.3)), 0.99, 3.0)
        self.assertFalse(r["ok"])
        self.assertTrue(any("no longer does" in p for p in r["problems"]))

    def test_faster_run_passes_and_is_reported(self) -> None:
        r = bench.check_rung(PARAMS, self.ref, bench.measure(synthetic_report(ic_precompute=0.2)), 0.5, 3.0)
        self.assertTrue(r["ok"], r["problems"])
        self.assertGreater(r["wall"]["now"]["whole_process_ratio"], r["wall"]["frozen"]["whole_process_ratio"])

    def test_verdict_inconsistent_with_timings_fails(self) -> None:
        rep = synthetic_report()
        rep["stages"][4]["vs_rho"]["verdict"]["whole_process_crossover"] = False  # timings say True
        r = bench.check_rung(PARAMS, self.ref, bench.measure(rep), 0.5, 3.0)
        self.assertFalse(r["ok"])
        self.assertTrue(any("inconsistent" in p for p in r["problems"]))

    def test_changed_parameter_file_fails(self) -> None:
        rep = synthetic_report()
        rep["params_digest"] = "cd" * 32
        r = bench.check_rung(PARAMS, self.ref, bench.measure(rep), 0.5, 3.0)
        self.assertFalse(r["ok"])
        self.assertTrue(any("parameter file changed" in p for p in r["problems"]))

    def test_broken_rho_baseline_fails(self) -> None:
        rep = synthetic_report()
        for row in rep["stages"][4]["vs_rho"]["targets_detail"]:
            row["walk_group_additions"] *= 10
        r = bench.check_rung(PARAMS, self.ref, bench.measure(rep), 0.5, 3.0)
        self.assertFalse(r["ok"])
        self.assertTrue(any("broken baseline" in p for p in r["problems"]))

    def test_markdown_has_one_row_per_rung(self) -> None:
        r = bench.check_rung(PARAMS, self.ref, bench.measure(synthetic_report()), 0.5, 3.0)
        summary = {"ok": True, "rungs_ok": 1, "rungs_total": 1, "rungs": [r], "errors": []}
        text = bench.render_markdown(summary)
        self.assertIn("**PASS**", text)
        self.assertEqual(text.count(f"| `{PARAMS}` |"), 1)
        self.assertIn("IC `S` is null", text)


class TierPinTests(unittest.TestCase):
    """The tier is pinned because no counter can see it.

    ``pair_table_stored_pairs`` is ``|F|(|F|+1)/2`` for the full tier and
    the compact one alike, so a run that switched between them used to
    read as "counters identical" while the algorithm had changed.  This
    pins the fix, not the bug.
    """

    def test_a_tier_change_is_caught_even_when_every_counter_matches(self):
        report = synthetic_report()
        m_before = bench.measure(report)
        switched = synthetic_report()
        switched["factor_base"]["pair_table_tier"] = "folded"
        m_after = bench.measure(switched)
        # The thing that made this slip through: nothing else moved.
        self.assertEqual(m_before["counters"], m_after["counters"])
        ref = dict(m_before)
        row = bench.check_rung("p.json", ref, m_after, 0.5, 3.0)
        self.assertFalse(row["ok"])
        self.assertTrue(
            any("tier changed" in p for p in row["problems"]),
            row["problems"],
        )

    def test_an_unchanged_tier_passes(self):
        m = bench.measure(synthetic_report())
        row = bench.check_rung("p.json", dict(m), m, 0.5, 3.0)
        self.assertTrue(row["ok"], row["problems"])

    def test_a_report_without_a_tier_is_refused(self):
        report = synthetic_report()
        del report["factor_base"]["pair_table_tier"]
        with self.assertRaises(bench.CheckFailure):
            bench.measure(report)


class CliTests(unittest.TestCase):
    def test_freeze_then_check_roundtrip_and_no_overwrite(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out = write_run(root, synthetic_report())
            ref = root / "ref.json"
            self.assertEqual(bench.main(["freeze", "--output", str(out), "--reference-out", str(ref), "--note", "t"]), 0)
            frozen = json.loads(ref.read_text())
            self.assertEqual(frozen["schema_version"], bench.SCHEMA_VERSION)
            self.assertIn(PARAMS, frozen["rungs"])
            self.assertEqual(frozen["frozen_from"]["git"]["commit"], "deadbeef")
            # never overwrite
            self.assertEqual(bench.main(["freeze", "--output", str(out), "--reference-out", str(ref)]), 2)
            summary = root / "summary.json"
            md = root / "summary.md"
            code = bench.main(["check", "--output", str(out), "--reference", str(ref),
                               "--summary-json", str(summary), "--summary-markdown", str(md)])
            self.assertEqual(code, 0)
            s = json.loads(summary.read_text())
            self.assertTrue(s["ok"])
            self.assertEqual(s["rungs_ok"], 1)
            self.assertIn("**PASS**", md.read_text())

    def test_check_fails_when_ic_exited_nonzero(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            good = write_run(root, synthetic_report())
            ref = root / "ref.json"
            self.assertEqual(bench.main(["freeze", "--output", str(good), "--reference-out", str(ref)]), 0)
            bad_root = root / "bad"
            bad_root.mkdir()
            bad = write_run(bad_root, synthetic_report(), exit_code=1)
            summary = root / "s.json"
            self.assertEqual(bench.main(["check", "--output", str(bad), "--reference", str(ref), "--summary-json", str(summary)]), 1)
            s = json.loads(summary.read_text())
            self.assertFalse(s["ok"])
            self.assertTrue(any("exited 1" in e for e in s["errors"]))

    def test_check_refuses_missing_or_extra_rungs(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out = write_run(root, synthetic_report())
            ref = root / "ref.json"
            self.assertEqual(bench.main(["freeze", "--output", str(out), "--reference-out", str(ref)]), 0)
            other_root = root / "other"
            other_root.mkdir()
            other = write_run(other_root, synthetic_report(), params="docs/ic/params/k0n-other.json")
            summary = root / "s.json"
            self.assertEqual(bench.main(["check", "--output", str(other), "--reference", str(ref), "--summary-json", str(summary)]), 1)
            errors = json.loads(summary.read_text())["errors"]
            self.assertTrue(any("not run" in e for e in errors))

    def test_freeze_refuses_unverified_run(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            rep = synthetic_report()
            rep["solutions"]["items"][0]["verified"] = False
            out = write_run(root, rep)
            self.assertEqual(bench.main(["freeze", "--output", str(out), "--reference-out", str(root / "ref.json")]), 1)
            self.assertFalse((root / "ref.json").exists())

    def test_run_invokes_ic_per_params_and_writes_manifest(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            stub = root / "ic"
            fixture = root / "fixture.json"
            fixture.write_text(json.dumps(synthetic_report()))
            # A stand-in ic: copies the fixture report to --out and echoes its argv to stderr.
            stub.write_text(
                "#!/usr/bin/env python3\nimport shutil,sys\nargs=sys.argv[1:]\n"
                f"shutil.copyfile({str(fixture)!r}, args[args.index('--out')+1])\n"
                "sys.stderr.write(' '.join(args))\n"
            )
            stub.chmod(stub.stat().st_mode | stat.S_IXUSR)
            p1 = root / "a.json"
            p2 = root / "b.json"
            p1.write_text("{}")
            p2.write_text("{}")
            out = root / "out"
            self.assertEqual(bench.main(["run", "--ic", str(stub), "--output", str(out), "--params", str(p1), str(p2)]), 0)
            manifest = json.loads((out / bench.MANIFEST_FILE).read_text())
            self.assertEqual([r["name"] for r in manifest["rungs"]], ["a", "b"])
            self.assertEqual(manifest["ic_binary"]["sha256"], bench._blake3_or_sha256(stub)["sha256"])
            for r in manifest["rungs"]:
                self.assertEqual(r["exit_code"], 0)
                self.assertEqual(r["command"][1], "workflow")
                self.assertTrue((out / r["report"]).is_file())
                self.assertIn("--dir", (out / r["stderr"]).read_text())
            # refuses a non-empty output directory
            self.assertEqual(bench.main(["run", "--ic", str(stub), "--output", str(out), "--params", str(p1)]), 2)

    def test_run_records_nonzero_exit(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            stub = root / "ic"
            stub.write_text("#!/bin/sh\necho boom >&2\nexit 3\n")
            stub.chmod(stub.stat().st_mode | stat.S_IXUSR)
            p1 = root / "a.json"
            p1.write_text("{}")
            out = root / "out"
            self.assertEqual(bench.main(["run", "--ic", str(stub), "--output", str(out), "--params", str(p1)]), 1)
            manifest = json.loads((out / bench.MANIFEST_FILE).read_text())
            self.assertEqual(manifest["rungs"][0]["exit_code"], 3)
            self.assertIn("boom", (out / manifest["rungs"][0]["stderr"]).read_text())


if __name__ == "__main__":
    unittest.main()
