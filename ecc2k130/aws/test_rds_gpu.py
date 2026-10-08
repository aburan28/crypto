#!/usr/bin/env python3
"""Packing tests for GPU → RDS records (no database)."""

import struct
import unittest

from rds_gpu import RECORD_BYTES, coeffBytes, iterGpuRecords, parseGpuRecord


class GpuRecordTests(unittest.TestCase):
    def test_parse_keeps_starting_seed_and_orbit_key(self):
        seed = 0x0123456789ABCDEF
        key = (0x1111111111111111, 0x2222222222222222, 0x3333333333333333)
        blob = struct.pack("<4Q", seed, *key)
        got_seed, point_key = parseGpuRecord(blob)
        self.assertEqual(got_seed, seed)
        self.assertEqual(point_key, struct.pack("<3Q", *key))
        self.assertEqual(len(point_key), 24)

    def test_iter_skips_trailing_partial_record(self):
        rec = struct.pack("<4Q", 7, 1, 2, 3)
        rows = list(iterGpuRecords(rec + b"\x00\x01"))
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0][0], 7)

    def test_walk_seed_encoding_is_big_endian_integer(self):
        self.assertEqual(coeffBytes(0), b"\x00")
        self.assertEqual(coeffBytes(1), b"\x01")
        self.assertEqual(coeffBytes(0x4F09), b"\x4f\x09")
        # a/b width used by the Type-II ingest (17 bytes mod 2^131)
        a = coeffBytes(0x0123456789ABCDEF, 1 << 131)
        self.assertEqual(len(a), 17)
        self.assertEqual(int.from_bytes(a, "big"), 0x0123456789ABCDEF)

    def test_record_size(self):
        self.assertEqual(RECORD_BYTES, 32)


if __name__ == "__main__":
    unittest.main()
