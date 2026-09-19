#!/usr/bin/env python3
import json
import subprocess
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
FREEZE = Path(__file__).resolve().parent / "freeze.py"


def sample(rate, clmad=1, sms=148, threads=75776, regs=94):
    raw = (
        f"NVIDIA B200, {sms} SMs, 2 block(s) of 256 packed threads resident per SM\n"
        f"packed kernel: {regs} registers/thread, 0 local bytes/thread\n"
        f"packed native carryless multiply: {clmad}\n"
        f"backend cuda-packed131: {threads} threads x 16 slots x 1 lanes\n"
        f"finished: {rate:.3f} M it/s, 0 distinguished points "
        "(0 verified against the reference, 0 dropped)\n"
    )
    return {
        "valid": True,
        "rate": rate,
        "packedClmad": bool(clmad),
        "expectedPackedClmad": bool(clmad),
        "raw": raw,
    }


def payload(gpu="NVIDIA B200", cc="100", rates=(9000.0, 8900.0, 8800.0)):
    return {
        "gpu": gpu,
        "cc": cc,
        "valid": True,
        "rate": sorted(rates)[1],
        "minRate": min(rates),
        "maxRate": max(rates),
        "packedClmad": True,
        "batch": 16,
        "threads": 256,
        "minBlocks": 2,
        "workers": 0,
        "steps": 1024,
        "launches": 32,
        "repeats": 3,
        "identity": {
            "sourceSha256": "abc",
            "binarySha256": "def",
            "activeBackend": "packed-poly131",
            "cudaImageVersion": "13.3.1",
            "compiler": "nvcc 13.3.73",
            "gpuState": "NVIDIA B200, ...",
        },
        "samples": [sample(r) for r in rates],
    }


class FreezeTests(unittest.TestCase):
    def freeze(self, obj, arm="clmad"):
        with tempfile.TemporaryDirectory() as tmp:
            src = Path(tmp) / "raw.json"
            dst = Path(tmp) / "out.json"
            src.write_text(json.dumps(obj))
            r = subprocess.run(
                ["python3", str(FREEZE), arm, str(src), str(dst)],
                cwd=ROOT, capture_output=True, text=True,
            )
            return r, json.loads(dst.read_text()) if r.returncode == 0 and dst.exists() else None

    def test_accepts_b200_and_records_ratios(self):
        r, out = self.freeze(payload())
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertTrue(out["valid"])
        self.assertEqual(out["gpu"], "NVIDIA B200")
        self.assertEqual(out["medianB"], 8.9)
        self.assertEqual(out["sms"], 148)
        self.assertEqual(out["automaticThreads"], 75776)
        self.assertAlmostEqual(out["ratioTo6000"], 8.9 / 15.115792)
        self.assertEqual(out["workers"], "automatic")
        self.assertEqual(out["requestedWorkers"], 0)

    def test_refuses_h100(self):
        r, out = self.freeze(payload(gpu="NVIDIA H100 80GB HBM3"))
        self.assertNotEqual(r.returncode, 0)
        self.assertIsNone(out)
        self.assertIn("non-B200", r.stderr)

    def test_refuses_wrong_cc(self):
        r, out = self.freeze(payload(cc="90"))
        self.assertNotEqual(r.returncode, 0)
        self.assertIsNone(out)

    def test_refuses_short_repeat_set(self):
        obj = payload()
        obj["samples"] = obj["samples"][:2]
        r, out = self.freeze(obj)
        self.assertNotEqual(r.returncode, 0)
        self.assertIsNone(out)

    def test_make_n_binds_b200_and_omits_fixed_worker_count(self):
        dry = subprocess.run(
            ["make", "-n", "bench-b200-modal"],
            cwd=ROOT, capture_output=True, text=True,
        )
        self.assertEqual(dry.returncode, 0, dry.stderr)
        self.assertIn("ECC_GPU=B200", dry.stdout)
        self.assertIn("ECC_PACKED_CLMAD=1", dry.stdout)
        self.assertIn("--batch 16", dry.stdout)
        self.assertNotIn("--workers", dry.stdout)
        val = subprocess.run(
            ["make", "-n", "validate-b200-modal"],
            cwd=ROOT, capture_output=True, text=True,
        )
        self.assertEqual(val.returncode, 0, val.stderr)
        self.assertIn("modal_app.py::validate", val.stdout)
        self.assertIn("ECC_GPU=B200", val.stdout)

    def test_run_modal_tees_outside_the_image_tree(self):
        script = (Path(__file__).resolve().parent / "run-modal.sh").read_text()
        self.assertIn("/tmp/ecc2k130-b200", script)
        self.assertNotIn('tee "$OUTDIR/', script)


if __name__ == "__main__":
    unittest.main()
