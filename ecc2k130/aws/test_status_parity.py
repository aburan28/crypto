"""ecc2k-status against status.py, on random slot tables from every source.

Both read the same table -- a rehearsal store's slots.json, a DynamoDB scan
or the bucket's slots/ objects, the last two through a fake `aws` CLI on
PATH -- and must print the same thing.  The Rust adds what status.py did not
show, and that is taken off before comparing: the `stepsPerWalk*` and
`offWeightSlots` keys after status.py's own, the steps-per-walk and
off-weight lines, and the steps/walk column at the end of each row.  Leases are drawn at least 30 s
from now, so only the countdown can differ between the two runs, and it is
masked.

Python exits 0 when a snapshot fails; ecc2k-status exits 1.  The message on
stderr and whatever was printed before it must still agree.

Skipped unless ECC_STATUS_BIN names the Rust binary.  This file goes with
status.py.
"""
import json
import os
import random
import re
import subprocess
import sys
import tempfile
import time
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
RUST = os.environ.get("ECC_STATUS_BIN")
TABLES = int(os.environ.get("ECC_STATUS_PARITY_TABLES", "300"))
COUNTDOWN = re.compile(r"(?<= )\d+s(?= |$)", re.MULTILINE)

FAKE_AWS = r'''#!/usr/bin/env python3
import json, os, shutil, sys
d = os.environ["FAKE_AWS_DIR"]
a = sys.argv[1:]
if a[:2] == ["dynamodb", "scan"]:
    sys.stdout.write(open(os.path.join(d, "dynamo.json")).read())
elif a[:2] == ["s3api", "list-objects-v2"]:
    sys.stdout.write(open(os.path.join(d, "list.json")).read())
elif a[:2] == ["s3api", "get-object"]:
    i = a.index("--key")
    src = os.path.join(d, "objects", a[i + 1].replace("/", "_"))
    if not os.path.exists(src):
        sys.stderr.write("An error occurred (NoSuchKey) when calling the GetObject operation: "
                         "The specified key does not exist.\n")
        sys.exit(254)
    shutil.copy(src, a[i + 2])
    print(json.dumps({"AcceptRanges": "bytes", "ETag": '"%s"' % a[i + 1]}))
else:
    sys.stderr.write("fake aws: unexpected %r\n" % a)
    sys.exit(2)
'''


def randomItem(rng, now):
    it = {}
    keys = ["ckptIter", "walks", "rate", "dpUploaded", "state", "leaseUntil", "owner",
            "gpuName", "solution", "spoolBytes", "iters", "kernelVersion"]
    rng.shuffle(keys)
    for k in keys:
        if rng.random() < 0.15:
            continue
        if k == "ckptIter":
            v = rng.choice([-1, 0, rng.randrange(1, 1 << 34), rng.randrange(1 << 30, 1 << 42), None, 1500.75])
        elif k == "walks":
            v = rng.choice([0, 6160384, 237568, 15204352, None, rng.randrange(1, 1 << 24)])
        elif k == "rate":
            v = rng.choice([0, 0.0, rng.uniform(0, 3e10), 14110000000.0, None, rng.randrange(0, 10 ** 11), 1e-300])
        elif k == "dpUploaded":
            v = rng.choice([0, rng.randrange(0, 100000), rng.randrange(50000, 10 ** 8), None])
        elif k == "state":
            v = rng.choice(["active", "active", "idle", "retired", "error", "solved", None, "weird"])
        elif k == "leaseUntil":
            v = rng.choice([0, int(now) + rng.randrange(30, 100000), int(now) - rng.randrange(30, 100000), None])
        elif k == "owner":
            v = rng.choice(["i-0abc123", "modal-run-3", "x" * 40, "\u00fcn\u00efc\u00f6d\u00e9-owner-\u20ac", None, 17, ""])
        elif k == "gpuName":
            v = rng.choice(["NVIDIA RTX PRO 6000 Blackwell Server Edition", "NVIDIA L4", None, "", "A" * 30])
        elif k == "solution":
            v = rng.choice(["k = 0x1234 (verified)", None, {"k": "0x1234", "verified": True}])
        else:
            v = rng.choice([1, "text", 2.5, True, [1, "a"], {"n": 1}])
        it[k] = v
    return it


def dynamoValue(v):
    if isinstance(v, bool):
        return {"BOOL": v}
    if isinstance(v, int):
        return {"N": str(v)}
    if isinstance(v, float):
        return {"N": repr(v)}
    if isinstance(v, str):
        return {"S": v}
    if v is None:
        return {"NULL": True}
    if isinstance(v, list):
        return {"L": [dynamoValue(x) for x in v]}
    return {"M": {k: dynamoValue(x) for k, x in v.items()}}


@unittest.skipUnless(RUST, "set ECC_STATUS_BIN to the Rust status binary")
class StatusParity(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="ecc-status-parity-")
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        bin_ = self.root / "bin"
        bin_.mkdir()
        (bin_ / "aws").write_text(FAKE_AWS)
        (bin_ / "aws").chmod(0o755)
        (self.root / "objects").mkdir()
        self.env = dict(os.environ, PATH="%s:%s" % (bin_, os.environ["PATH"]), FAKE_AWS_DIR=str(self.root),
                        TMPDIR=str(self.root))
        for var in ("ECC_BUCKET", "ECC_TABLE", "ECC_LOCAL_STORE"):
            self.env.pop(var, None)
        self.compared = 0

    def store(self, slots):
        """Write one table to every source."""
        (self.root / "slots.json").write_text(json.dumps({str(k): v for k, v in slots}))
        items = [dict({"slot": {"N": str(k)}}, **{f: dynamoValue(x) for f, x in v.items()}) for k, v in slots]
        (self.root / "dynamo.json").write_text(json.dumps({"Items": items, "Count": len(items)}))
        keys = ["slots/slot-%05d.json" % k for k, _ in slots] + ["slots/README.txt", "slots/slot-99999.json"]
        (self.root / "list.json").write_text(json.dumps(keys or None))
        for old in (self.root / "objects").iterdir():
            old.unlink()
        for k, v in slots:
            (self.root / "objects" / ("slots_slot-%05d.json" % k)).write_text(json.dumps(v))

    def execute(self, cmd):
        return subprocess.run(cmd, capture_output=True, text=True, env=self.env, timeout=60)

    def both(self, *args):
        py = self.execute([sys.executable, str(HERE / "status.py"), *args])
        rs = self.execute([RUST, *args])
        return py, rs

    def assertSameSnapshot(self, py, rs, label):
        failed = "status failed" in py.stderr
        self.assertEqual(py.returncode, 0, label)
        self.assertEqual(rs.returncode, 1 if failed else 0, "%s\n%s" % (label, rs.stderr))
        self.assertEqual(py.stderr, rs.stderr, label)
        if "--json" in label:
            if failed:
                self.assertEqual(py.stdout, rs.stdout, label)
                return
            self.assertTrue(py.stdout.endswith("\n}\n"), label)
            head = py.stdout[:-3]
            self.assertEqual(rs.stdout[:len(head)], head, label)
            self.assertTrue(rs.stdout[len(head):].startswith(',\n "stepsPerWalk": '), label)
        else:
            self.assertEqual(self.pythonText(rs.stdout), COUNTDOWN.sub("Ns", py.stdout), label)
        self.compared += 1

    @staticmethod
    def pythonText(text):
        """The Rust snapshot without what status.py did not print."""
        out, inTable = [], False
        for line in text.splitlines(keepends=True):
            if line.startswith(("steps per walk", "  off weight", "  without them")):
                continue
            if line.startswith(" slot owner"):
                inTable = True
            if inTable:
                line = line[:-12] + "\n"
            out.append(COUNTDOWN.sub("Ns", line))
        return "".join(out)

    def test_random_tables_from_every_source(self):
        for seed in range(TABLES):
            rng = random.Random(seed)
            now = time.time()
            numbers = rng.sample(range(0, 10000), rng.randrange(0, 13))
            self.store([(n, randomItem(rng, now)) for n in numbers])
            walks = rng.choice([[], ["--walks", "237568"], ["--walks", "1_000"]])
            for source in (["--local", str(self.root)], ["--table", "slots"], ["--bucket", "b"]):
                for mode in ([], ["--json"]):
                    args = [*source, *walks, *mode]
                    py, rs = self.both(*args)
                    self.assertSameSnapshot(py, rs, "seed %d %s" % (seed, " ".join(args)))
        print("\n%d snapshots identical over %d tables" % (self.compared, TABLES), file=sys.stderr)

    def test_failures_say_the_same_thing(self):
        now = time.time()
        cases = [
            [(3, {"leaseUntil": "soon", "state": "active"})],
            [(4, {"ckptIter": [1], "leaseUntil": 0})],
            [(5, {"rate": "fast", "leaseUntil": int(now) + 600})],
            [(6, {"rate": "fast", "leaseUntil": 0}), (1, {"walks": 16, "ckptIter": 3})],
            [(7, {"walks": "many", "leaseUntil": 0})],
        ]
        failures = 0
        for slots in cases:
            self.store(slots)
            for source in (["--local", str(self.root)], ["--table", "slots"], ["--bucket", "b"]):
                for mode in ([], ["--json"]):
                    py, rs = self.both(*source, *mode)
                    failures += "status failed" in py.stderr
                    self.assertSameSnapshot(py, rs, "%r %s" % (slots, " ".join(source + mode)))
        # 30 snapshots, 5 of which succeed: DynamoDB carries the list as L,
        # which fromDv reads as None, and a dead slot's rate is read only for
        # its row, so the JSON never touches it.
        self.assertEqual(failures, 25)
        (self.root / "slots.json").unlink()
        py, rs = self.both("--local", str(self.root))
        self.assertSameSnapshot(py, rs, "missing slots.json")
        py, rs = self.both()
        self.assertEqual((py.returncode, py.stderr), (rs.returncode, rs.stderr))

    def test_dynamo_numbers_and_s3_quirks(self):
        now = time.time()
        self.store([(2, {"leaseUntil": int(now) + 900, "rate": 1.5e10, "walks": 16, "ckptIter": 10, "state": "active"})])
        # A float slot number from DynamoDB, a key that names a missing object,
        # a nested key the scan reads by its slot number, and no listing at all.
        dynamo = json.loads((self.root / "dynamo.json").read_text())
        dynamo["Items"].append({"slot": {"N": "9.0"}, "state": {"S": "idle"}, "ckptIter": {"N": "1e3"}})
        (self.root / "dynamo.json").write_text(json.dumps(dynamo))
        (self.root / "list.json").write_text(json.dumps(["slots/x/slot-00002.json", "slots/slot-00002.json",
                                                         "slots/slot-00008.json"]))
        for args in (["--table", "slots"], ["--table", "slots", "--json"], ["--bucket", "b"], ["--bucket", "b", "--json"]):
            py, rs = self.both(*args)
            self.assertSameSnapshot(py, rs, " ".join(args))
        (self.root / "list.json").write_text("null")
        py, rs = self.both("--bucket", "b", "--json")
        self.assertSameSnapshot(py, rs, "empty bucket --json")


if __name__ == "__main__":
    unittest.main()
