#!/usr/bin/env python3
"""Offline tests for the sanitized DP-weight status snapshot."""

from __future__ import annotations

import io
import json
import os
import struct
import sys
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import weight_status  # noqa: E402


class FakePaginator:
    def __init__(self, pages):
        self.pages = pages

    def paginate(self, **_kwargs):
        return self.pages


class FakeS3:
    def __init__(self, objects, bodies):
        self.objects = objects
        self.bodies = bodies

    def get_paginator(self, _name):
        return FakePaginator([{"Contents": self.objects}])

    def get_object(self, Bucket, Key, Range=None):  # noqa: N803 - boto shape
        value = self.bodies[Key]
        if Range:
            start, end = Range.removeprefix("bytes=").split("-")
            value = value[int(start):int(end) + 1]
        return {"Body": io.BytesIO(value)}


def checkpoint(iter_base=100):
    return struct.pack("<8s6IQ", b"ECC2KCK\0", 0x101, 131, 4, 16, 1, 77, iter_base) + b"payload"


class WeightStatusTests(unittest.TestCase):
    def test_build_status_is_counts_only_and_keeps_campaign_separate(self):
        objects = [
            {"Key": "campaigns/demo/campaign.json", "Size": 40, "LastModified": "2026-10-11T00:00:00Z"},
            {"Key": "campaigns/demo/dp/slot-00000/a.bin", "Size": 64, "LastModified": "2026-10-11T11:30:00Z"},
            {"Key": "campaigns/demo/dp/slot-00000/b.bin", "Size": 32, "LastModified": "2026-10-11T11:45:00Z"},
            {"Key": "campaigns/demo/slots/slot-00000.json", "Size": 20, "LastModified": "2026-10-11T11:50:00Z"},
            {"Key": "campaigns/demo/ckpt/slot-00000.ck", "Size": 48, "LastModified": "2026-10-11T11:50:00Z"},
        ]
        bodies = {
            "campaigns/demo/campaign.json": json.dumps({
                "campaignId": "demo-table-walk",
                "curveId": "certicom-ecc2k-130",
                "dpWeight": 34,
                "owner": "must-not-publish",
            }).encode(),
            "campaigns/demo/slots/slot-00000.json": json.dumps({
                "dpUploaded": 3,
                "createdAt": 1791718200,
                "updatedAt": 1791719400,
                "gpuName": "RTX fixture",
                "owner": "private-owner",
            }).encode(),
            "campaigns/demo/ckpt/slot-00000.ck": checkpoint(100),
        }
        # 2026-10-11T12:00:00Z, so both records are in the last hour.
        status = weight_status.build_status(
            FakeS3(objects, bodies),
            "bucket",
            "campaigns/demo",
            weight_status.as_datetime("2026-10-11T12:00:00Z"),
        )
        self.assertEqual(status["campaign_id"], "demo-table-walk")
        self.assertEqual(status["dp_mask_bits"], 34)
        self.assertEqual(status["dps"], 3)
        self.assertEqual(status["dps_last_hour"], 3)
        self.assertEqual(status["work"]["iterations"], 4 * 16 + 100)
        self.assertEqual(status["per_worker"][0]["dps"], 3)
        self.assertNotIn("owner", json.dumps(status))
        self.assertNotIn("campaigns/demo", json.dumps(status))
        self.assertIn("not merged", status["claim_boundary"])


if __name__ == "__main__":
    unittest.main()
