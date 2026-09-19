#!/usr/bin/env python3
"""Unit tests for the Crossbred α AutoLab control plane."""
from __future__ import annotations

import os
import sys
import tempfile
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import crossbred_autolab as lab  # noqa: E402


X1_TABLE = """
=== Oracle cost per target, three engines, same verdicts ===

| n | m | ℓ | v | |F| | Q_enum | Q_word | Q/Q_enum | deg | targets | reference | agree | brute (bit ops) | F4 (bit ops) | crossbred (bit ops) | xb/brute | xb/F4 | D | k | kernel | filters | xb wall |
|--:|--:|--:|--:|----:|-------:|-------:|---------:|----:|--------:|:----------|:-----:|----------------:|-------------:|--------------------:|---------:|------:|--:|--:|-------:|--------:|--------:|
T4 pick n=5 m=3 v=17 D=3 k=8
| 5 | 3 | 4 | 17 | 21 | 210 | 17482 | 83.248 | 3 | 11 | exhaustive, == | yes | 43599313 | 17119202 | 1118853 | 0.026 | 0.065 | 3 | 8 | 52 | 0 | 0.8 ms |
| 7 | 3 | 3 | 16 | 1 | 0 | 5006 | — | 3 | 11 | exhaustive, == | yes | 20000395 | 25174603 | 320436 | 0.016 | 0.013 | 3 | 6 | 15 | 0 | 0.2 ms |
| 13 | 3 | 12 | 49 | — | — | — | — | — | — | — | — | — | — | — | — | — | — | — | — | — | no determining space |

X1/X3 unit: Q_word is 64-bit word XORs per target (extraction + search).
"""

FFD_TABLE = """
| n | ℓ | m | vars | eqs | deg | FFD min | FFD max | no fall | mean syz D=2 |
|--:|--:|--:|-----:|----:|----:|--------:|--------:|--------:|-------------:|
| 9 | 6 | 3 | 27 | 18 | 3 | 3 | 3 | 0/4 | 0.00 |
| 9 | 6 | 4 | 42 | 27 | 3 | 3 | 3 | 0/4 | 0.00 |
| 15 | 4 | 3 | 27 | 30 | 3 | 3 | 3 | 0/4 | 0.00 |
| 15 | 4 | 4 | 46 | 45 | 3 | 3 | 3 | 0/4 | 0.00 |
"""

GROWING_FFD = """
| n | ℓ | m | vars | eqs | deg | FFD min | FFD max | no fall | mean syz D=2 |
|--:|--:|--:|-----:|----:|----:|--------:|--------:|--------:|-------------:|
| 9 | 6 | 4 | 42 | 27 | 3 | 3 | 3 | 0/4 | 0.00 |
| 15 | 4 | 4 | 46 | 45 | 3 | 3 | 4 | 0/4 | 0.00 |
"""

ELL1_THEN_CUBIC = """
| n | ℓ | m | vars | eqs | deg | FFD min | FFD max | no fall | mean syz D=2 |
|--:|--:|--:|-----:|----:|----:|--------:|--------:|--------:|-------------:|
| 7 | 1 | 4 | 17 | 21 | 1 | 2 | 2 | 0/4 | 223.50 |
| 9 | 3 | 4 | 29 | 27 | 3 | 4 | 4 | 0/4 | 0.00 |
| 15 | 3 | 4 | 41 | 45 | 3 | 4 | 4 | 0/4 | 0.00 |
"""


class ProtocolTests(unittest.TestCase):
    def test_protocol_loads(self) -> None:
        protocol = lab.load_protocol()
        self.assertEqual(protocol["task_id"], lab.TASK_ID)
        self.assertIn("smoke.x1_n5", protocol["beats"])
        self.assertIn("replay.x3_k1_n7", protocol["beats"])
        self.assertIn("fit.alpha", protocol["beats"])
        self.assertIn("x5.ffd_chained_m4", protocol["beats"])
        self.assertIn("x5.ffd_chained_m4_16", protocol["beats"])
        self.assertIn("x5.ffd_chained_sym_m4_smoke", protocol["beats"])
        self.assertIn("x5.ffd_chained_sym_m4", protocol["beats"])
        self.assertIsNone(protocol["beats"]["fit.alpha"]["expect"]["fit"])
        self.assertEqual(protocol["beats"]["smoke.x1_n5"]["expect"]["Q_over_Q_enum"], 83.248)
        self.assertEqual(protocol["beats"]["replay.x3_k1_n7"]["expect"]["Q_over_Q_enum"], 4.653)

    def test_plan_lists_beats_and_refuses_a_fit(self) -> None:
        report = lab.plan(lab.load_protocol())
        beat_ids = [row["beat_id"] for row in report["beats"]]
        self.assertIn("x5.ffd_chained_sym_m4_smoke", beat_ids)
        self.assertIn("x5.ffd_chained_m4_16", beat_ids)
        self.assertIsNone(report["next"])
        self.assertIn("X5 FFD", report["next_note"])
        self.assertIsNone(report["incumbent"]["fit"]["fit"])
        self.assertEqual(report["incumbent"]["fit"]["frames"]["x1"]["n_rungs"], 2)
        self.assertEqual(report["incumbent"]["fit"]["frames"]["x3"]["n_rungs"], 2)
        rungs = report["incumbent"]["usable_rungs"]
        self.assertEqual(len(rungs), 4)
        frames = {r["frame"] for r in rungs}
        self.assertEqual(frames, {"x1", "x3"})


class ParserTests(unittest.TestCase):
    def test_pipe_in_header_does_not_shift_columns(self) -> None:
        rows = lab.parse_crossbred_output(X1_TABLE)
        self.assertEqual(len(rows), 3)
        n5 = rows[0]
        self.assertEqual(n5["n"], 5)
        self.assertEqual(n5["m"], 3)
        self.assertEqual(n5["ell"], 4)
        self.assertEqual(n5["F"], 21)
        self.assertEqual(n5["Q_enum"], 210)
        self.assertEqual(n5["Q_word"], 17482)
        self.assertAlmostEqual(n5["Q_over_Q_enum"], 83.248)
        self.assertTrue(n5["agree"])
        self.assertEqual(n5["filters"], 0)
        self.assertIsNone(rows[1]["Q_over_Q_enum"])
        self.assertEqual(rows[2]["reason"], "no determining space")

    def test_ffd_table(self) -> None:
        rows = lab.parse_ffd_table(FFD_TABLE)
        self.assertEqual(len(rows), 4)
        m4 = [r for r in rows if r["m"] == 4]
        self.assertEqual([r["n"] for r in m4], [9, 15])
        self.assertEqual([r["ffd_max"] for r in m4], [3, 3])
        self.assertEqual(m4[0]["trials"], 4)

    def test_table_excerpt_drops_cargo_preamble(self) -> None:
        text = "warning: unused\n" + X1_TABLE
        excerpt = lab.table_excerpt(text)
        self.assertTrue(excerpt.startswith("| n |"))
        self.assertNotIn("warning:", excerpt)


class FitHygieneTests(unittest.TestCase):
    def test_fewer_than_four_is_not_a_fit(self) -> None:
        rungs = [
            {"ell": 4, "Q_over_Q_enum": 83.248},
            {"ell": 6, "Q_over_Q_enum": 319.222},
        ]
        fit = lab.fit_alpha(rungs)
        self.assertIsNone(fit["fit"])
        self.assertEqual(fit["n_rungs"], 2)

    def test_four_rungs_with_falling_ratio_is_success(self) -> None:
        rungs = [
            {"ell": 4, "Q_over_Q_enum": 16.0},
            {"ell": 5, "Q_over_Q_enum": 8.0},
            {"ell": 6, "Q_over_Q_enum": 4.0},
            {"ell": 7, "Q_over_Q_enum": 2.0},
        ]
        fit = lab.fit_alpha(rungs)
        self.assertTrue(fit["fit"])
        self.assertAlmostEqual(fit["slope"], -1.0)
        self.assertAlmostEqual(fit["alpha"], 1.0)
        self.assertEqual(fit["verdict"], "success")
        self.assertEqual(fit["class"], "advance")

    def test_flat_ratio_is_falsified(self) -> None:
        rungs = [{"ell": ell, "Q_over_Q_enum": 4.0} for ell in (4, 5, 6, 7)]
        fit = lab.fit_alpha(rungs)
        self.assertEqual(fit["verdict"], "falsified")
        self.assertAlmostEqual(fit["alpha"], 2.0)

    def test_claim_fit_passes_on_frozen_no_fit_per_frame(self) -> None:
        protocol = lab.load_protocol()
        # Concatenating X1+X3 yields four rungs and would OLS-fit if we
        # mixed frames. The claim must not do that: α is per frame.
        x1, x3 = lab.load_frozen(protocol)
        self.assertEqual(len(lab.usable_rungs(x1)), 2)
        self.assertEqual(len(lab.usable_rungs(x3)), 2)
        claim = lab.claim_fit(protocol)
        self.assertEqual(claim["status"], "PASS")
        self.assertIsNone(claim["fit"]["fit"])

    def test_two_point_sketch_numbers_are_not_a_result(self) -> None:
        x1, _ = lab.load_frozen(lab.load_protocol())
        sketch = x1["two_point_sketch_not_a_fit"]
        self.assertEqual(sketch["class"], "not a result")
        delta_ell = sketch["delta_ell"]
        self.assertEqual(delta_ell, 2)
        rebuilt = sketch["delta_log2_Q_over_C"] / delta_ell
        self.assertAlmostEqual(rebuilt, sketch["sketch_slope"], places=3)


class ClaimTests(unittest.TestCase):
    def test_smoke_matches_frozen_q_over_c(self) -> None:
        protocol = lab.load_protocol()
        beat = protocol["beats"]["smoke.x1_n5"]
        claim = lab.claim_smoke(beat, X1_TABLE)
        self.assertEqual(claim["status"], "PASS")
        self.assertEqual(claim["class"], "accounting")

    def test_smoke_fails_on_missing_row(self) -> None:
        protocol = lab.load_protocol()
        beat = protocol["beats"]["replay.x3_k1_n7"]
        claim = lab.claim_smoke(beat, X1_TABLE)
        self.assertEqual(claim["status"], "FAIL")

    def test_ffd_incumbent_does_not_grow(self) -> None:
        protocol = lab.load_protocol()
        beat = protocol["beats"]["x5.ffd_chained_m4"]
        claim = lab.claim_ffd(FFD_TABLE, beat)
        self.assertEqual(claim["status"], "PASS")
        self.assertFalse(claim["m4_ffd_grows_with_n"])
        self.assertEqual(claim["class"], "measurement")
        self.assertIn("chained x-system", claim["note"])
        self.assertNotIn("unchained S4", claim["note"])

    def test_ffd_sym_smoke_note_is_not_the_x_arm(self) -> None:
        protocol = lab.load_protocol()
        beat = protocol["beats"]["x5.ffd_chained_sym_m4_smoke"]
        claim = lab.claim_ffd(FFD_TABLE, beat)
        self.assertEqual(claim["status"], "PASS")
        self.assertEqual(claim["class"], "measurement")
        self.assertIn("Not unchained S4", claim["note"])
        self.assertNotIn("chained x-system", claim["note"])

    def test_ffd_growth_is_measurement_not_advance(self) -> None:
        protocol = lab.load_protocol()
        beat = protocol["beats"]["x5.ffd_chained_m4"]
        claim = lab.claim_ffd(GROWING_FFD, beat)
        self.assertTrue(claim["m4_ffd_grows_with_n"])
        self.assertEqual(claim["class"], "measurement")

    def test_ell_1_row_is_not_h1(self) -> None:
        protocol = lab.load_protocol()
        beat = protocol["beats"]["x5.ffd_chained_sym_m4_smoke"]
        claim = lab.claim_ffd(ELL1_THEN_CUBIC, beat)
        self.assertEqual(claim["m4_ffd_maxima_including_ell_1"], [2, 4, 4])
        self.assertTrue(claim["m4_ffd_grows_with_n_including_ell_1"])
        self.assertEqual(claim["m4_ffd_maxima"], [4, 4])
        self.assertFalse(claim["m4_ffd_grows_with_n"])
        self.assertEqual(len(claim["m4_ell_eq_1"]), 1)
        self.assertEqual(claim["class"], "measurement")


class LockTests(unittest.TestCase):
    def test_stale_lock_clears(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            os.environ["CROSSBRED_AUTOLAB_RUNS"] = tmp
            try:
                path = lab.lock_path()
                path.write_text("1\n")  # PID 1 is not this process; may be alive.
                # Write a definitely-dead pid.
                path.write_text("99999999\n")
                lab.clear_stale_lock()
                self.assertFalse(path.exists())
            finally:
                os.environ.pop("CROSSBRED_AUTOLAB_RUNS", None)

    def test_lock_blocks_a_second_holder(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            os.environ["CROSSBRED_AUTOLAB_RUNS"] = tmp
            try:
                with lab.RunLock():
                    with self.assertRaises(SystemExit):
                        with lab.RunLock():
                            pass
            finally:
                os.environ.pop("CROSSBRED_AUTOLAB_RUNS", None)


class FitLaunchHygieneTests(unittest.TestCase):
    def test_launch_fit_alpha_under_tmp_runs(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            os.environ["CROSSBRED_AUTOLAB_RUNS"] = tmp
            try:
                claim = lab.launch("fit.alpha")
                self.assertEqual(claim["status"], "PASS")
                self.assertIsNone(claim["fit"]["fit"])
                st = lab.status()
                self.assertEqual(st["status"], "PASS")
                verified = lab.verify()
                self.assertEqual(verified["status"], "PASS")
            finally:
                os.environ.pop("CROSSBRED_AUTOLAB_RUNS", None)


if __name__ == "__main__":
    unittest.main()
