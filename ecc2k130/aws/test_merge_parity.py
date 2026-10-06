"""The Rust merge against merge.py, byte for byte.

Both run on the same corpora into separate work directories, and every
pass must end with the same return code, the same stdout and the same
files, bit for bit -- buckets, state.json, campaign.lock.json, pair files,
solution.json -- except the solve's own wall-clock stamp ("when").

The corpora are what a campaign holds: real client output (this build's
corpus layout, v2 when the witness is compiled in), the 32-byte and
table-v3 layouts, nested directories, re-reports of the same walk, a
cross-corpus collision that solves, and a strict campaign's committed
chunks.

Skipped unless ECC_MERGE_BIN names the Rust binary
(`cargo build --release --bin ecc2k-merge`); CI's merge-parity job sets it.
This file goes when merge.py does.
"""
import json
import os
import re
import shutil
import struct
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import protocol

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
CLIENT = ROOT / "ecc2k130-cpu"
FIXTURE = ROOT / "build/test-production"
RUST = os.environ.get("ECC_MERGE_BIN")
PYTHON = [sys.executable, str(HERE / "merge.py")]
RECORD = struct.Struct("<4Q")
WHEN = re.compile(rb'"when": "[^"]*"')


@unittest.skipUnless(RUST, "set ECC_MERGE_BIN to the Rust merge binary")
class MergeParity(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.fixture = json.loads(subprocess.check_output([str(FIXTURE), "--fixture"], text=True))

    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="ecc-merge-parity-")
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.corpus = self.root / "corpus"
        self.corpus.mkdir()
        self.works = {"python": self.root / "work-python", "rust": self.root / "work-rust"}

    # ---- corpora -----------------------------------------------------------
    def workers(self, count=12, steps=3, launches=3, base=200, name="slot-%02d.bin"):
        """Real client corpora: one file per worker, disjoint run ids."""
        for i in range(count):
            subprocess.run([str(CLIENT), "--curve", "41", "--threads", "1", "--steps", str(steps),
                            "--launches", str(launches), "--dp-weight", "18", "--verify", "0",
                            "--run-id", str(base + i), "--dp-file", str(self.corpus / (name % i))],
                           check=True, capture_output=True, timeout=60)

    def v1(self, rel, records):
        path = self.corpus / rel
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(b"".join(RECORD.pack(*r) for r in records))
        return path

    def table3(self, rel, records):
        path = self.corpus / rel
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(b"ECC2KDT3" + struct.pack("<II", 3, 32) + b"".join(RECORD.pack(*r) for r in records))
        return path

    def fixturePair(self):
        f = self.fixture
        self.v1("worker-a.bin", [(f["seedA"], *f["key"]), (f["seedA"], *f["key"])])
        self.v1("nested/worker-b.bin", [(f["seedB"], *f["key"])])

    # ---- running and comparing --------------------------------------------
    def normal(self, name, data):
        """The two runs differ in their work directory and the solve's clock;
        the client prints the pair file's path, so its tail names the first."""
        return WHEN.sub(b'"when": "-"', data.replace(str(self.works[name]).encode(), b"WORK"))

    def both(self, *args):
        out = {}
        for name, cmd in (("python", PYTHON), ("rust", [RUST])):
            out[name] = subprocess.run([*cmd, "--work", str(self.works[name]), "--client", str(CLIENT), *args],
                                       capture_output=True, timeout=300)
        py, rs = out["python"], out["rust"]
        self.assertEqual(py.returncode, rs.returncode,
                         "python %d, rust %d\n%s\n%s" % (py.returncode, rs.returncode, py.stderr[-2000:], rs.stderr[-2000:]))
        self.assertEqual(self.normal("python", py.stdout), self.normal("rust", rs.stdout))
        self.assertSameTree()
        return py

    def tree(self, name):
        root = self.works[name]
        if not root.exists():
            return {}
        out = {}
        for path in sorted(root.rglob("*")):
            if path.is_file():
                out[str(path.relative_to(root))] = self.normal(name, path.read_bytes())
        return out

    def assertSameTree(self):
        py, rs = self.tree("python"), self.tree("rust")
        self.assertEqual(sorted(py), sorted(rs))
        for name in py:
            self.assertEqual(py[name], rs[name], name)

    def legacy(self, *extra):
        return self.both("--legacy", "--local", str(self.corpus), "--curve", "41",
                         "--dp-weight", str(self.fixture["weight"]), "--buckets", "64", *extra)

    # ---- the cases -----------------------------------------------------------
    def test_client_corpora_and_every_layout_detect_identically(self):
        self.workers()
        self.v1("legacy/one.bin", [(1, 2, 3, 4), (1, 2, 3, 4), (5, 2, 3, 4)])
        self.table3("table/three.bin", [(9, 2, 3, 4), (6, 7, 8, 9)])
        p = self.legacy("--detect-only")
        summary = json.loads(p.stdout)
        self.assertGreater(summary["corpus"], 1000)
        # (1|5|9, key 2,3,4) is a three-way meeting: two adjacent pairs.
        self.assertGreaterEqual(summary["collisions"], 2)

    def test_incremental_passes_stay_identical(self):
        self.workers(count=6)
        self.legacy("--detect-only")
        self.workers(count=4, base=400, name="late-%02d.bin")
        # A worker's corpus grows past its committed offset; the new tail
        # re-reports two of its own points, which the sort must drop.
        grown = self.corpus / "slot-00.bin"
        data = grown.read_bytes()
        head, stride = (16, 72) if data[:8] == b"ECC2KDP2" else (0, 32)
        self.assertGreaterEqual(len(data), head + 2 * stride)
        grown.write_bytes(data + data[head:head + 2 * stride])
        self.v1("late.bin", [(77, 1, 1, 1)])
        p = self.legacy("--detect-only")
        self.assertGreater(json.loads(p.stdout)["added"], 0)
        p = self.legacy("--detect-only")                     # nothing new
        self.assertEqual(json.loads(p.stdout)["added"], 0)

    def test_a_cross_corpus_collision_solves_identically(self):
        self.fixturePair()
        self.legacy("--detect-only")
        p = self.legacy()
        self.assertEqual(json.loads(p.stdout)["solution"]["k"], self.fixture["k"])
        p = self.legacy()                                    # already solved
        self.assertEqual(json.loads(p.stdout)["k"], self.fixture["k"])

    def test_an_unsolvable_pair_records_the_same_failure(self):
        self.v1("a.bin", [(11, 21, 31, 41)])
        self.v1("b.bin", [(12, 21, 31, 41)])
        p = self.legacy("--solve-timeout", "0.000001")
        self.assertEqual(p.returncode, 2)

    def test_a_strict_campaign_merges_identically(self):
        binary = protocol.sha256File(CLIENT)
        config = dict(storageProtocol=protocol.PROTOCOL, curve=41, dpWeight=self.fixture["weight"],
                      maxIters=100000, packed=False, workers=1, batch=32, extraArgs=[],
                      binarySha256=binary, hostBinarySha256=binary, sourceSha256="a" * 64)
        campaign = protocol.campaignContract(config)
        configPath = self.root / "campaign.json"
        protocol.atomicJson(configPath, config)
        self.fixturePair()
        uncommitted = self.v1("pending.bin", [(3, 3, 3, 3)])
        for path in self.corpus.rglob("*.bin"):
            if path != uncommitted:
                protocol.atomicJson(str(path) + ".json", protocol.envelope(path, campaign, "dp"))
        args = ("--campaign", str(configPath), "--local", str(self.corpus), "--buckets", "16")
        self.both(*args, "--detect-only")
        protocol.atomicJson(str(uncommitted) + ".json", protocol.envelope(uncommitted, campaign, "dp"))
        p = self.both(*args)
        self.assertEqual(json.loads(p.stdout)["solution"]["k"], self.fixture["k"])

    def test_failures_are_the_same_failures(self):
        self.v1("a.bin", [(1, 2, 3, 4)])
        self.legacy("--detect-only")
        for work in self.works.values():
            bucket = next((work / "buckets").glob("*.bin"))
            data = bytearray(bucket.read_bytes()); data[3] ^= 1; bucket.write_bytes(bytes(data))
        p = self.legacy("--detect-only")                     # bucket integrity
        self.assertEqual(p.returncode, 1)
        for name, work in self.works.items():
            shutil.rmtree(work)
        self.legacy("--detect-only")
        (self.corpus / "a.bin").write_bytes(b"")             # source shrank
        self.assertEqual(self.legacy("--detect-only").returncode, 1)
        for flag in (("--buckets", "3"), ("--solve-timeout", "0")):
            self.assertEqual(self.both("--legacy", "--local", str(self.corpus), *flag).returncode, 2)


if __name__ == "__main__":
    unittest.main()
