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

    def test_run_sh_tees_outside_the_image_tree(self):
        script = (HERE / "run.sh").read_text()
        self.assertIn("/tmp/ecc2k130-per-sm", script)
        self.assertNotIn('tee "$OUTDIR/', script)


if __name__ == "__main__":
    unittest.main()
