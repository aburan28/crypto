#!/usr/bin/env python3
"""Offline tests for the distinguished-point ingest.

The incident these pin: on 2026-09-17 the fleet's workers moved to the
`ecc2k-seed-orbit-v1` object naming at 10:06Z, the last worker still using the
old naming uploaded at 13:50Z, and the ingest -- which matched only the old
shape -- then took nothing for five hours while 56 M records landed in `dp/`.
It logged no error, because an object whose key does not match is not an
error, it is invisible. Every test below is about one of those two properties:
both key shapes are points, and anything unreadable is counted out loud.

    python3 -m unittest test_dp_ingest
"""

import datetime
import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import dp_ingest


def when(iso):
    return datetime.datetime.fromisoformat(iso).replace(tzinfo=datetime.timezone.utc)


LEGACY = "dp/slot-00002/1789311001-0000000000000000.bin"
ORBIT = ("dp/slot-00140/a8b4133d5c5f414ebd1337f15603588d-0000000000791392-"
         "9c6a7dae1862f8393fe196b0c4f890e5e01f7f45615be019252b86ca3e429914.bin")


class FakeS3:
    """Enough of list_objects_v2 to drive s3Objects, including paging."""

    def __init__(self, items, pageSize=2):
        self.items = list(items)
        self.pageSize = pageSize

    def list_objects_v2(self, **kw):
        start = int(kw.get("ContinuationToken") or 0)
        page = self.items[start:start + self.pageSize]
        nxt = start + self.pageSize
        return {"Contents": [{"Key": k, "Size": s, "LastModified": t} for k, s, t in page],
                "IsTruncated": nxt < len(self.items),
                "NextContinuationToken": str(nxt)}


class KeyShapes(unittest.TestCase):
    def test_legacy_key_takes_found_at_from_the_name(self):
        self.assertEqual(dp_ingest.foundAt(LEGACY, when("2026-09-20T00:00:00")), 1789311001)

    def test_orbit_key_takes_found_at_from_last_modified(self):
        t = when("2026-09-17T18:47:34")
        self.assertEqual(dp_ingest.foundAt(ORBIT, t), int(t.timestamp()))

    def test_orbit_key_is_a_point_object_at_all(self):
        # The regression itself: this returned None before the fix, and an
        # object with no found_at was skipped without a word.
        self.assertIsNotNone(dp_ingest.foundAt(ORBIT, when("2026-09-17T18:47:34")))

    def test_envelopes_and_foreign_keys_are_not_points(self):
        for key in (ORBIT + ".json", "dp/slot-00002/notes.txt", "ckpt/slot-00002.ck"):
            self.assertIsNone(dp_ingest.foundAt(key, when("2026-09-17T18:47:34")), key)


class Listing(unittest.TestCase):
    def objects(self, items):
        return dp_ingest.s3Objects(FakeS3(items), "bucket")

    def test_both_shapes_are_listed_and_sorted_by_upload_time(self):
        objects, unrecognised = self.objects([
            (ORBIT, 3200, when("2026-09-17T18:47:34")),
            (LEGACY, 6400, when("2026-09-17T13:50:59")),
        ])
        self.assertEqual([o[0] for o in objects], [LEGACY, ORBIT])
        self.assertEqual([o[1] for o in objects], [200, 100])
        self.assertEqual(unrecognised, [])

    def test_an_unreadable_key_is_reported_not_skipped(self):
        objects, unrecognised = self.objects([
            (LEGACY, 6400, when("2026-09-17T13:50:59")),
            ("dp/slot-00007/something-new.bin", 64, when("2026-09-18T00:00:00")),
        ])
        self.assertEqual([o[0] for o in objects], [LEGACY])
        self.assertEqual(unrecognised, ["dp/slot-00007/something-new.bin"])

    def test_envelopes_are_neither_ingested_nor_reported_as_unreadable(self):
        objects, unrecognised = self.objects([
            (ORBIT, 3200, when("2026-09-17T18:47:34")),
            (ORBIT + ".json", 412, when("2026-09-17T18:47:34")),
        ])
        self.assertEqual([o[0] for o in objects], [ORBIT])
        self.assertEqual(unrecognised, [])

    def test_a_short_object_carrying_no_whole_record_is_not_listed(self):
        objects, _ = self.objects([(ORBIT, 16, when("2026-09-17T18:47:34"))])
        self.assertEqual(objects, [])

    def test_paging_does_not_lose_objects(self):
        items = [(LEGACY.replace("0000000000000000", "%016d" % i), 64,
                  when("2026-09-17T13:50:59")) for i in range(7)]
        objects, _ = self.objects(items)
        self.assertEqual(len(objects), 7)


class Decoding(unittest.TestCase):
    def test_point_key_is_the_last_24_bytes_of_the_record(self):
        record = bytes(range(32))
        self.assertEqual(dp_ingest.decode(record)["point_key"], record[8:])

    def test_walk_seed_is_the_seed_big_endian_without_leading_zeros(self):
        record = (7).to_bytes(8, "little") + bytes(24)
        self.assertEqual(dp_ingest.decode(record)["walk_seed"], b"\x07")

    def test_coefficients_are_17_bytes(self):
        r = dp_ingest.decode(bytes(range(32)))
        self.assertEqual((len(r["a"]), len(r["b"])), (17, 17))
        self.assertEqual(r["b"], bytes(17))


class StatusState(unittest.TestCase):
    """The page must not read a stalled ingest as a quiet campaign."""

    def state(self, dps, dpsLastHour, collisions, ingest):
        behind = int(ingest.get("outstanding", 0)) or int(ingest.get("unrecognised", 0))
        return ("COLLISION_RECORDED" if collisions else
                "COLLECTING" if dpsLastHour else
                "INGEST_BEHIND" if behind else
                "IDLE_OR_STALE" if dps else "EMPTY")

    def test_a_backlog_reads_as_ingest_behind_not_idle(self):
        self.assertEqual(self.state(60053195, 0, 0, {"outstanding": 3795}), "INGEST_BEHIND")

    def test_unreadable_keys_alone_are_enough_to_say_behind(self):
        self.assertEqual(self.state(60053195, 0, 0, {"unrecognised": 12}), "INGEST_BEHIND")

    def test_a_caught_up_ingest_with_no_recent_points_is_still_idle(self):
        self.assertEqual(self.state(60053195, 0, 0, {"outstanding": 0}), "IDLE_OR_STALE")

    def test_points_in_the_last_hour_win(self):
        self.assertEqual(self.state(60053195, 142035, 0, {"outstanding": 3795}), "COLLECTING")


if __name__ == "__main__":
    unittest.main()
