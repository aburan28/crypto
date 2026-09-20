#!/usr/bin/env python3
import json
import subprocess
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
FREEZE = HERE / "freeze.py"


def payload(gpu="NVIDIA RTX PRO 6000 Blackwell Server Edition", sms=188, auto=96256):
    rates = {
        1: (14470.0, 14480.0, 14460.0),
        4: (15110.0, 15120.0, 15100.0),
        6: (15100.0, 15090.0, 15105.0),
        8: (14900.0, 14880.0, 14920.0),
    }
    rows = []
    for wave, rs in rates.items():
        med = sorted(rs)[1]
        rows.append({
            "wave": wave,
            "workers": auto * wave,
            "valid": True,
            "rate": med,
            "minRate": min(rs),
            "maxRate": max(rs),
            "ratesM": list(rs),
            "perSmM": med / sms,
            "samples": [{"valid": True, "rate": r} for r in rs],
        })
    return {
        "gpu": gpu,
        "cc": "120",
        "valid": True,
        "sms": sms,
        "automaticThreads": auto,
        "residentBlocks": 2,
        "blockThreads": 256,
        "rows": rows,
    }


class FreezeTests(unittest.TestCase):
    def freeze(self, obj, modal="RTX-PRO-6000"):
        with tempfile.TemporaryDirectory() as tmp:
            src = Path(tmp) / "raw.json"
            dst = Path(tmp) / "out.json"
            src.write_text(json.dumps(obj))
            r = subprocess.run(
                ["python3", str(FREEZE), modal, str(src), str(dst)],
                cwd=ROOT, capture_output=True, text=True,
            )
            out = json.loads(dst.read_text()) if r.returncode == 0 and dst.exists() else None
            return r, out

    def test_records_per_sm_and_does_not_beat_four(self):
        r, out = self.freeze(payload())
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual(out["sms"], 188)
        self.assertEqual(out["bestWave"], 4)
        self.assertFalse(out["fourBeaten"])
        four = next(row for row in out["rows"] if row["wave"] == 4)
        self.assertAlmostEqual(four["perSmM"], 15110.0 / 188)

    def test_four_beaten_when_eight_is_2pct_faster(self):
        obj = payload()
        obj["rows"][-1]["rate"] = 15500.0
        obj["rows"][-1]["perSmM"] = 15500.0 / 188
        r, out = self.freeze(obj)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertTrue(out["fourBeaten"])
        self.assertEqual(out["bestWave"], 8)

    def test_refuses_wrong_gpu(self):
        r, out = self.freeze(payload(gpu="NVIDIA L40S"))
        self.assertNotEqual(r.returncode, 0)
        self.assertIsNone(out)

    def test_make_n_binds_waves_entrypoint(self):
        dry = subprocess.run(
            ["make", "-n", "bench-waves-modal"],
            cwd=ROOT, capture_output=True, text=True,
        )
        self.assertEqual(dry.returncode, 0, dry.stderr)
        self.assertIn("ECC_GPU=RTX-PRO-6000", dry.stdout)
        self.assertIn("modal_app.py::waves", dry.stdout)
        self.assertIn("--wave-list", dry.stdout)
        self.assertNotIn("--workers", dry.stdout)

    def test_parses_spinner_between_wave_and_repeats(self):
        mangled = """
"gpu": "NVIDIA RTX PRO 6000 Blackwell Server Edition",
"cc": "120",
automatic occupancy: 96256 threads on 188 SMs (2 x 256 resident)
=== wave 1: 96256 workers ===
Running (1/1 containers active)... View app at https://modal.com/apps/x
  repeat 1/3: 14470.000 M it/s (complete)
  repeat 2/3: 14480.000 M it/s (complete)
  repeat 3/3: 14460.000 M it/s (complete)
=== wave 4: 385024 workers ===
  repeat 1/3: 15110.000 M it/s (complete)
  repeat 2/3: 15120.000 M it/s (complete)
  repeat 3/3: 15100.000 M it/s (complete)
=== wave 6: 577536 workers ===
  repeat 1/3: 15090.000 M it/s (complete)
  repeat 2/3: 15100.000 M it/s (complete)
  repeat 3/3: 15105.000 M it/s (complete)
=== wave 8: 770048 workers ===
  repeat 1/3: 14900.000 M it/s (complete)
  repeat 2/3: 14880.000 M it/s (complete)
  repeat 3/3: 14920.000 M it/s (complete)
"""
        with tempfile.TemporaryDirectory() as tmp:
            src = Path(tmp) / "m.log"
            dst = Path(tmp) / "out.json"
            src.write_text(mangled)
            r = subprocess.run(
                ["python3", str(FREEZE), "RTX-PRO-6000", str(src), str(dst)],
                cwd=ROOT, capture_output=True, text=True,
            )
            self.assertEqual(r.returncode, 0, r.stderr)
            out = json.loads(dst.read_text())
            self.assertEqual(out["sms"], 188)
            self.assertEqual(out["residentBlocks"], 2)
            self.assertEqual(out["blockThreads"], 256)
            self.assertEqual([row["wave"] for row in out["rows"]], [1, 4, 6, 8])
            self.assertEqual(out["bestWave"], 4)
            self.assertFalse(out["fourBeaten"])

    def test_run_sh_tees_outside_the_image_tree(self):
        script = (HERE / "run.sh").read_text()
        self.assertIn("/tmp/ecc2k130-per-sm", script)
        self.assertNotIn('tee "$OUTDIR/', script)

    def test_committed_logs_freeze_to_receipts(self):
        for modal, slug, best, beaten in (
                ("RTX-PRO-6000", "rtx-pro-6000", 4, False),
                ("L40S", "l40s", 1, True)):
            log = HERE / f"{slug}-waves.log"
            receipt = json.loads((HERE / f"{slug}-waves.json").read_text())
            with tempfile.TemporaryDirectory() as tmp:
                dst = Path(tmp) / "out.json"
                r = subprocess.run(
                    ["python3", str(FREEZE), modal, str(log), str(dst)],
                    cwd=ROOT, capture_output=True, text=True,
                )
                self.assertEqual(r.returncode, 0, r.stderr)
                out = json.loads(dst.read_text())
            self.assertEqual(out["bestWave"], best)
            self.assertEqual(out["fourBeaten"], beaten)
            self.assertEqual(out["sms"], receipt["sms"])
            self.assertEqual(
                [row["medianB"] for row in out["rows"]],
                [row["medianB"] for row in receipt["rows"]])


if __name__ == "__main__":
    unittest.main()
