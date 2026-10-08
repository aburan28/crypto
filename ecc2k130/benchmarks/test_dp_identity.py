import importlib.util
import pathlib
import struct
import subprocess
import sys
import tempfile
import unittest


HERE = pathlib.Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location("dp_identity", HERE / "dp_identity.py")
DP = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(DP)


def record(seed: int, stride: int) -> bytes:
    return seed.to_bytes(8, "little") + bytes([seed & 255]) * (stride - 8)


class CorpusIdentityTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = pathlib.Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def write(self, name, body):
        path = self.root / f"dp-{name}.bin"
        path.write_bytes(body)
        return path

    def test_headerless_v1_and_record_order(self):
        a, b = record(1, 32), record(2, 32)
        self.write("ref", a + b)
        self.write("candidate", b + a)
        rows, valid = DP.compare(self.root, ["ref", "candidate"])
        self.assertTrue(valid)
        self.assertIn("2 records", rows[1])

    def test_table_v3_header_is_not_part_of_first_record(self):
        header = struct.pack("<8sII", b"ECC2KDT3", 3, 32)
        a, b = record(3, 32), record(4, 32)
        self.write("ref", header + a + b)
        self.write("candidate", header + b + a)
        identity = DP.corpus_identity(self.root / "dp-ref.bin")
        self.assertEqual(identity["format"], "table-v3")
        self.assertEqual(identity["headerBytes"], 16)
        self.assertTrue(DP.compare(self.root, ["ref", "candidate"])[1])

    def test_witness_v2_uses_declared_72_byte_records(self):
        header = struct.pack("<8sII", b"ECC2KDP2", 2, 72)
        self.write("ref", header + record(5, 72))
        identity = DP.corpus_identity(self.root / "dp-ref.bin")
        self.assertEqual(identity["records"], 1)
        self.assertEqual(identity["format"], "v2")

    def test_bad_header_and_trailing_bytes_fail(self):
        bad = struct.pack("<8sII", b"ECC2KDT3", 2, 32) + record(1, 32)
        self.write("bad", bad)
        with self.assertRaisesRegex(ValueError, "invalid table-v3 header"):
            DP.corpus_identity(self.root / "dp-bad.bin")
        self.write("tail", record(1, 32) + b"x")
        with self.assertRaisesRegex(ValueError, "not a multiple"):
            DP.corpus_identity(self.root / "dp-tail.bin")

    def test_cli_exit_status_tracks_identity(self):
        self.write("ref", record(1, 32))
        self.write("same", record(1, 32))
        self.write("different", record(2, 32))
        command = [sys.executable, str(HERE / "dp_identity.py"), str(self.root)]
        self.assertEqual(subprocess.run(command + ["ref", "same"], check=False).returncode, 0)
        self.assertEqual(subprocess.run(command + ["ref", "different"], check=False).returncode, 1)


if __name__ == "__main__":
    unittest.main()
