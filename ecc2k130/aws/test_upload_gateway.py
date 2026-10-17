"""Pure validation tests for the credential-free worker upload gateway."""

import hashlib
import json
import os
import struct
import tempfile
import unittest

from upload_gateway.common import (DP_HEADER_BYTES, DP_V2_MAGIC, ValidationError,
                                   campaignIdentity, reserveStreamRange,
                                   validateDpFile,
                                   validateTicket)
from protocol import PROTOCOL, campaignContract


def campaign(cutoff=32):
    return {"storageProtocol": PROTOCOL, "curve": 131, "dpWeight": cutoff,
            "maxIters": 1 << 32, "packed": True, "workers": 385024,
            "batch": 16, "binarySha256": "1" * 64,
            "hostBinarySha256": "2" * 64, "sourceSha256": "3" * 64,
            "walk": "sigma"}


def canonicalValue(bits):
    ring = 263
    orbit = []
    value = bits
    for _ in range(131):
        orbit.append(value)
        mapped = 0
        for source in range(1, 132):
            if not value & (1 << (source - 1)):
                continue
            target = (2 * source) % ring
            if target > 131:
                target = ring - target
            mapped |= 1 << (target - 1)
        value = mapped
    return min(orbit)


def record(slot, lane, point, restart=0):
    seed = ((slot + 1) << 48) | (lane << 16) | restart
    return struct.pack("<QQQQ", seed, point & ((1 << 64) - 1),
                       (point >> 64) & ((1 << 64) - 1), point >> 128)


class GatewayValidation(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.path = os.path.join(self.tmp.name, "body.bin")

    def write(self, body):
        with open(self.path, "wb") as output:
            output.write(body)

    def test_strict_campaign_identity_matches_the_existing_protocol(self):
        config = campaign()
        self.assertEqual(campaignIdentity(config), campaignContract(config)["id"])

    def test_valid_v1_records_bind_run_lane_weight_and_canonical_orbit(self):
        point = canonicalValue((1 << 130) | (1 << 73) | 1)
        self.assertEqual(point, 0xa040)
        self.write(record(246, 19, point) + record(246, 20, point, restart=1))
        got = validateDpFile(self.path, 246, campaign())
        self.assertEqual(got["records"], 2)
        self.assertEqual(got["recordBytes"], 64)

    def test_wrong_run_id_is_rejected(self):
        self.write(record(10, 1, 0))
        with self.assertRaisesRegex(ValidationError, "run id"):
            validateDpFile(self.path, 11, campaign())

    def test_noncanonical_orbit_representative_is_rejected(self):
        self.write(record(4, 1, 1 << 130))
        with self.assertRaisesRegex(ValidationError, "canonical"):
            validateDpFile(self.path, 4, campaign())

    def test_point_above_dp_cutoff_is_rejected(self):
        point = canonicalValue((1 << 33) - 1)
        self.write(record(4, 1, point))
        with self.assertRaisesRegex(ValidationError, "DP predicate"):
            validateDpFile(self.path, 4, campaign(cutoff=32))

    def test_v2_witness_counts_must_sum_to_iterations(self):
        seed = (8 << 48) | (5 << 16)
        header = DP_V2_MAGIC + struct.pack("<II", 2, 72)
        body = struct.pack("<QQQQQ8I", seed, 9, 0, 0, 0,
                           1, 1, 1, 1, 1, 1, 1, 1)
        self.assertEqual(len(body), 72)
        self.write(header + body)
        with self.assertRaisesRegex(ValidationError, "do not sum"):
            validateDpFile(self.path, 7, campaign())

    def test_valid_v2_reports_source_record_offset_after_header(self):
        seed = (8 << 48) | (5 << 16)
        header = DP_V2_MAGIC + struct.pack("<II", 2, 72)
        body = struct.pack("<QQQQQ8I", seed, 8, 0, 0, 0,
                           1, 1, 1, 1, 1, 1, 1, 1)
        self.write(header + body)
        got = validateDpFile(self.path, 7, campaign())
        self.assertEqual(got["headerBytes"], DP_HEADER_BYTES)
        self.assertEqual(got["recordBytes"], 72)
        self.assertEqual(reserveStreamRange(16, 72, 16), 88)
        self.assertEqual(reserveStreamRange(88, 144, 16, nextOffset=88), 232)
        with self.assertRaisesRegex(ValidationError, "first record"):
            reserveStreamRange(0, 72, 16)

    def test_ticket_hash_must_match_content_addressed_key(self):
        digest = hashlib.sha256(b"body").hexdigest()
        body = {"key": "dp/slot-00001/%s-0-%s.bin" % ("a" * 32, "b" * 64),
                "bytes": 32, "sha256": digest, "workerId": "worker-a",
                "owner": "owner-a"}
        with self.assertRaisesRegex(ValidationError, "immutable key"):
            validateTicket(body)

    def test_ticket_manifest_must_bind_the_body(self):
        digest = "a" * 64
        body = {"key": "dp/slot-00001/%s-0-%s.bin" % ("b" * 32, digest),
                "bytes": 32, "sha256": digest, "workerId": "worker-a",
                "owner": "owner-a", "manifest": {"bytes": 64, "sha256": digest}}
        with self.assertRaisesRegex(ValidationError, "manifest"):
            validateTicket(body)


if __name__ == "__main__":
    unittest.main()
