#!/usr/bin/env python3
import json
import subprocess
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
FREEZE = HERE / "freeze.py"
CATALOG = json.loads((HERE / "catalog.json").read_text())


def sample(rate, clmad=1, sms=142, threads=72704, regs=94, device="NVIDIA L40S"):
    raw = (
        f"{device}, {sms} SMs, 2 block(s) of 256 packed threads resident per SM\n"
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


def payload(gpu="NVIDIA L40S", cc="89", rates=(8900.0, 8800.0, 8700.0), clmad=True,
            sms=142, threads=72704):
    return {
        "gpu": gpu,
        "cc": cc,
        "valid": True,
        "rate": sorted(rates)[1],
        "minRate": min(rates),
        "maxRate": max(rates),
        "packedClmad": clmad,
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
            "gpuState": gpu + ", ...",
        },
        "samples": [sample(r, clmad=int(clmad), sms=sms, threads=threads, device=gpu)
                    for r in rates],
    }


class FreezeTests(unittest.TestCase):
    def freeze(self, obj, modal="L40S"):
        with tempfile.TemporaryDirectory() as tmp:
            src = Path(tmp) / "raw.json"
            dst = Path(tmp) / "out.json"
            src.write_text(json.dumps(obj))
            r = subprocess.run(
                ["python3", str(FREEZE), modal, str(src), str(dst)],
                cwd=ROOT, capture_output=True, text=True,
            )
            return r, json.loads(dst.read_text()) if r.returncode == 0 and dst.exists() else None

    def test_accepts_l40s_and_records_ratios(self):
        r, out = self.freeze(payload())
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertTrue(out["valid"])
        self.assertEqual(out["modal"], "L40S")
        self.assertEqual(out["medianB"], 8.8)
        self.assertEqual(out["sms"], 142)
        self.assertEqual(out["automaticThreads"], 72704)
        self.assertAlmostEqual(out["ratioTo6000"], 8.8 / 15.115792)
        self.assertEqual(out["workers"], "automatic")
        self.assertEqual(out["requestedWorkers"], 0)
        self.assertTrue(out["packedClmad"])

    def test_refuses_h200_when_h100_pinned(self):
        r, out = self.freeze(payload(gpu="NVIDIA H200", cc="90"), modal="H100!")
        self.assertNotEqual(r.returncode, 0)
        self.assertIsNone(out)
        self.assertIn("H100", r.stderr + r.stdout)

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

    def test_t4_wants_software_product(self):
        obj = payload(gpu="Tesla T4", cc="75", clmad=False, sms=40, threads=20480,
                      rates=(540.0, 533.0, 529.0))
        r, out = self.freeze(obj, modal="T4")
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertFalse(out["packedClmad"])
        self.assertAlmostEqual(out["medianB"], 0.533)

    def test_t4_refuses_clmad(self):
        obj = payload(gpu="Tesla T4", cc="75", clmad=True, sms=40, threads=20480)
        r, out = self.freeze(obj, modal="T4")
        self.assertNotEqual(r.returncode, 0)
        self.assertIsNone(out)

    def test_a100_40_refuses_80gb(self):
        r, out = self.freeze(
            payload(gpu="NVIDIA A100-SXM4-80GB", cc="80"), modal="A100-40GB")
        self.assertNotEqual(r.returncode, 0)
        self.assertIsNone(out)

    def test_gpu_arch_covers_every_catalog_run(self):
        text = (ROOT / "modal_app.py").read_text()
        start = text.index("GPU_ARCH = {")
        end = text.index("}", start)
        blob = text[start:end]
        for gpu in CATALOG["gpus"]:
            if gpu.get("run"):
                self.assertIn('"%s"' % gpu["modal"], blob, gpu["modal"])

    def test_make_n_binds_survey_gpu_and_omits_fixed_worker_count(self):
        dry = subprocess.run(
            ["make", "-n", "bench-modal-gpu", "SURVEY_GPU=B200"],
            cwd=ROOT, capture_output=True, text=True,
        )
        self.assertEqual(dry.returncode, 0, dry.stderr)
        self.assertIn("ECC_GPU=B200", dry.stdout)
        self.assertIn("ECC_PACKED_CLMAD=1", dry.stdout)
        self.assertIn("--batch 16", dry.stdout)
        self.assertNotIn("--workers", dry.stdout)

    def test_run_one_tees_outside_the_image_tree(self):
        script = (HERE / "run-one.sh").read_text()
        self.assertIn("/tmp/ecc2k130-survey", script)
        self.assertNotIn('tee "$OUTDIR/', script)

    def test_catalog_skips_upgrade_aliases(self):
        self.assertIn("A100", CATALOG["skipAliases"])
        self.assertIn("H100", CATALOG["skipAliases"])
        self.assertIn("B200+", CATALOG["skipAliases"])
        run = {g["modal"] for g in CATALOG["gpus"] if g.get("run")}
        self.assertNotIn("A100", run)
        self.assertNotIn("H100", run)
        self.assertNotIn("B200+", run)
        self.assertIn("H100!", run)
        self.assertIn("A100-40GB", run)
        self.assertIn("B300", run)

    def test_parses_spinner_split_identity(self):
        mangled = """
repeat 1/3: 8838.416 M it/s (complete)
repeat 2/3: 8872.941 M it/s (complete)
repeat 3/3: 8727.023 M it/s (complete)
"gpu": "NVIDIA L40S",
"cc": "89",
"packedClmad": true,
device: NVIDIA L40S, 142 SMs, 2 block(s) of 256 packed threads
packed native
carryless multiply: 1
packed kernel: 94 registers/thread, 0 local bytes/thread
backend cuda-packed131: 72704 threads x 16 slots
"sourceSha256":
"aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa",
"binarySha256":
"bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb",
"cudaImageVersion": "13.3.1"
"""
        with tempfile.TemporaryDirectory() as tmp:
            src = Path(tmp) / "m.log"
            dst = Path(tmp) / "out.json"
            src.write_text(mangled)
            r = subprocess.run(
                ["python3", str(FREEZE), "L40S", str(src), str(dst)],
                cwd=ROOT, capture_output=True, text=True,
            )
            self.assertEqual(r.returncode, 0, r.stderr)
            out = json.loads(dst.read_text())
            self.assertEqual(out["sms"], 142)
            self.assertTrue(out["packedClmad"])
            self.assertAlmostEqual(out["medianB"], 8.838416)
            self.assertEqual(out["automaticThreads"], 72704)
            self.assertEqual(out["parsedFrom"], "modal-log")


if __name__ == "__main__":
    unittest.main()
