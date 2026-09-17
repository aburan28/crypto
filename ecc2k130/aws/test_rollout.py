"""Stage / activate geometry gate and worker self-recycle.

A kernel compile must not flip live geometry. These checks pin that
activate refuses frozen knobs, allows CLMAD, and that workers reload the
client after binaryKey moves without releasing the slot.
"""
from pathlib import Path
import hashlib
import json
import os
import shutil
import subprocess
import tempfile
import time
import unittest

import rollout
from worker import (Worker, campaignPointerMoved, frozenCampaignMoved)


HERE = Path(__file__).resolve().parent

SHIPPING = (
    "BATCH=16 THREADS=256 MINBLOCKS=2 PACKED_SINGLE_PRODUCT=1 PACKED_CACHE_DENOM=1 "
    "PACKED_BY_VALUE=1 PACKED_PERM_SIGMA=3 PACKED_POLY_CHAIN=1 PACKED_UNROLL_INV=1 "
    "PACKED_PAIR_PRODUCTS=1 PACKED_POLY_STATE=1 PACKED_DIRECT_REDUCE=1 "
    "PACKED_GENERATED_PRODUCT=1 PACKED_CLMAD=1 PACKED_STATE_TILE=256 "
    "PACKED_WEIGHTED_PREFIX=2 PACKED_COMPACT_STATE=1 PACKED_SHARED_SIGMA=1 "
    "WALK_TABLE=0 TABLE_BRANCHES=8"
)


def campaign(**extra):
    c = {
        "curve": 131, "walk": "sigma", "dpWeight": 32, "packed": True,
        "batch": 16, "blockThreads": 256, "minBlocks": 2, "workers": 385024,
        "steps": 1024, "checkpointEvery": 600,
        "binaryKey": "bin/oldoldoldoldol/ecc2k130",
        "hostBinaryKey": "bin/oldoldoldoldol/ecc2k130-cpu",
        "binarySha256": "a" * 64,
        "hostBinarySha256": "b" * 64,
        "sourceSha256": "c" * 64,
    }
    c.update(extra)
    return c


def manifest(knobs=SHIPPING, arches=None, **extra):
    m = {
        "knobs": knobs,
        "arches": arches if arches is not None else ["89", "120"],
        "binarySha256": "d" * 64,
        "hostBinarySha256": "e" * 64,
        "sourceSha256": "f" * 64,
        "buildSha256": "1" * 64,
    }
    m.update(extra)
    return m


def knobsReplace(src, **changes):
    parts = dict(tok.split("=", 1) for tok in src.split())
    parts.update({k: str(v) for k, v in changes.items()})
    return " ".join("%s=%s" % (k, parts[k]) for k in parts)


class GeometryGate(unittest.TestCase):
    def test_clmad_may_move(self):
        live = campaign()
        staged = manifest(knobs=knobsReplace(SHIPPING, PACKED_CLMAD=0))
        self.assertEqual(rollout.geometryReasons(live, staged, manifest()), [])

    def test_product_knob_may_move(self):
        staged = manifest(knobs=knobsReplace(SHIPPING, PACKED_PAIR_PRODUCTS=0))
        self.assertEqual(rollout.geometryReasons(campaign(), staged, manifest()), [])

    def test_batch_blocked(self):
        staged = manifest(knobs=knobsReplace(SHIPPING, BATCH=32))
        reasons = rollout.geometryReasons(campaign(), staged, manifest())
        self.assertTrue(any("geometry" in r for r in reasons))

    def test_walk_blocked(self):
        staged = manifest(knobs=knobsReplace(SHIPPING, WALK_TABLE=1))
        reasons = rollout.geometryReasons(campaign(), staged, manifest())
        self.assertTrue(any("walk" in r for r in reasons))

    def test_state_tile_blocked(self):
        staged = manifest(knobs=knobsReplace(SHIPPING, PACKED_STATE_TILE=128))
        reasons = rollout.geometryReasons(campaign(), staged, manifest())
        self.assertTrue(any("PACKED_STATE_TILE" in r for r in reasons))

    def test_compact_state_blocked(self):
        staged = manifest(knobs=knobsReplace(SHIPPING, PACKED_COMPACT_STATE=0))
        reasons = rollout.geometryReasons(campaign(), staged, manifest())
        self.assertTrue(any("PACKED_COMPACT_STATE" in r for r in reasons))

    def test_dropping_arch_blocked(self):
        staged = manifest(arches=["120"])
        reasons = rollout.geometryReasons(campaign(), staged, manifest(arches=["89", "120"]))
        self.assertTrue(any("drop live" in r for r in reasons))

    def test_unmanifested_live_requires_fat_client(self):
        staged = manifest(arches=["120"])
        reasons = rollout.geometryReasons(campaign(), staged, None)
        self.assertTrue(any("must cover" in r for r in reasons))

    def test_empty_pointer_refuses_activate(self):
        reasons = rollout.geometryReasons(campaign(binaryKey=""), manifest(), None)
        self.assertTrue(any("no binaryKey" in r for r in reasons))

    def test_apply_pointer_only(self):
        live = campaign()
        staged = manifest()
        out = rollout.applyPointer(live, "bin/c68a7b6eb45dbbc3", staged)
        self.assertTrue(rollout.pointerOnly(live, out))
        self.assertEqual(out["binaryKey"], "bin/c68a7b6eb45dbbc3/ecc2k130")
        self.assertEqual(out["workers"], 385024)
        self.assertEqual(out["walk"], "sigma")
        self.assertEqual(out["batch"], 16)
        self.assertEqual(out["binarySha256"], "d" * 64)

    def test_cli_check_and_apply(self):
        tmp = tempfile.mkdtemp()
        try:
            camp = os.path.join(tmp, "campaign.json")
            staged = os.path.join(tmp, "manifest.json")
            json.dump(campaign(), open(camp, "w"))
            json.dump(manifest(), open(staged, "w"))
            r = subprocess.run(
                ["python3", str(HERE / "rollout.py"), "check",
                 "--campaign", camp, "--staged", staged],
                capture_output=True, text=True)
            self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
            json.dump(manifest(knobs=knobsReplace(SHIPPING, BATCH=32)), open(staged, "w"))
            r = subprocess.run(
                ["python3", str(HERE / "rollout.py"), "check",
                 "--campaign", camp, "--staged", staged],
                capture_output=True, text=True)
            self.assertEqual(r.returncode, 2)
        finally:
            shutil.rmtree(tmp)


class Adoption(unittest.TestCase):
    def test_walking_slots(self):
        now = int(time.time())
        slots = [
            {"state": "active", "leaseUntil": now + 120, "binary": "bin/new/ecc2k130"},
            {"state": "active", "leaseUntil": now + 120, "binary": "bin/old/ecc2k130"},
            {"state": "active", "leaseUntil": now + 120},
            {"state": "idle", "leaseUntil": 0, "binary": "bin/old/ecc2k130"},
            {"state": "active", "leaseUntil": now - 10, "binary": "bin/old/ecc2k130"},
        ]
        rep = rollout.walkingAdoption(slots, "bin/new/ecc2k130")
        self.assertEqual(rep["walking"], 3)
        self.assertEqual(rep["adopted"], 1)
        self.assertEqual(rep["stale"], 1)
        self.assertEqual(rep["unknown"], 1)
        self.assertFalse(rep["done"])

    def test_done_when_all_adopted(self):
        now = int(time.time())
        slots = [{"state": "active", "leaseUntil": now + 60, "binarySha256": "d" * 64}]
        rep = rollout.walkingAdoption(slots, "bin/new/ecc2k130", "d" * 64)
        self.assertTrue(rep["done"])


class PointerHelpers(unittest.TestCase):
    def test_worker_and_rollout_agree(self):
        cur = campaign()
        nxt = campaign(binaryKey="bin/new/ecc2k130", binarySha256="d" * 64)
        self.assertTrue(rollout.campaignPointerMoved(cur, nxt))
        self.assertTrue(campaignPointerMoved(cur, nxt))
        bad = campaign(workers=1)
        self.assertFalse(rollout.campaignPointerMoved(cur, bad))
        self.assertEqual(frozenCampaignMoved(cur, bad), "workers")

    def test_geometry_change_is_not_a_pointer_move(self):
        cur = campaign()
        nxt = campaign(walk="table", binaryKey="bin/new/ecc2k130")
        self.assertEqual(frozenCampaignMoved(cur, nxt), "walk")
        self.assertFalse(campaignPointerMoved(cur, nxt))


class WorkerReload(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.store = os.path.join(self.tmp, "store")
        self.root = os.path.join(self.tmp, "root")
        os.makedirs(self.store)
        os.makedirs(self.root)
        self.client = os.path.join(self.root, "ecc2k130")
        open(self.client, "wb").write(b"old-client")
        os.chmod(self.client, 0o755)
        self.env = {
            "ECC_LOCAL_STORE": self.store,
            "ECC_ROOT": self.root,
            "ECC_CLIENT": self.client,
            "ECC_ALLOW_LEGACY_STORAGE": "1",
            "ECC_GPU": "0",
            "ECC_DEVICE_NAME": "NVIDIA RTX PRO 6000 Blackwell Server Edition",
            "ECC_INSTANCE_TYPE": "g7e.2xlarge",
        }
        self.saved = {k: os.environ.get(k) for k in self.env}
        os.environ.update(self.env)
        self.w = Worker()

    def tearDown(self):
        for k, v in self.saved.items():
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v
        shutil.rmtree(self.tmp)

    def _put(self, key, data):
        path = os.path.join(self.store, key)
        os.makedirs(os.path.dirname(path), exist_ok=True)
        if isinstance(data, (dict, list)):
            json.dump(data, open(path, "w"))
        else:
            open(path, "wb").write(data)

    def test_pointer_change_reloads_client(self):
        new = b"new-client-bytes"
        digest = hashlib.sha256(new).hexdigest()
        live = campaign()
        nxt = campaign(binaryKey="bin/newnewnewnewne/ecc2k130",
                       hostBinaryKey="bin/newnewnewnewne/ecc2k130-cpu",
                       binarySha256=digest)
        self._put("campaign.json", live)
        self.w.loadConfig()
        self._put("campaign.json", nxt)
        self.assertTrue(self.w.campaignPointerChanged())
        self._put("bin/newnewnewnewne/ecc2k130", new)
        self.w.reloadClient()
        self.assertEqual(open(self.client, "rb").read(), new)
        self.assertEqual(self.w.cfg["binaryKey"], "bin/newnewnewnewne/ecc2k130")

    def test_frozen_geometry_keeps_old_client(self):
        self._put("campaign.json", campaign())
        self.w.loadConfig()
        self._put("campaign.json", campaign(workers=1, binaryKey="bin/new/ecc2k130"))
        self.assertFalse(self.w.campaignPointerChanged())
        self.assertEqual(open(self.client, "rb").read(), b"old-client")

    def test_local_key_does_not_fetch(self):
        self.w.cfg = campaign(binaryKey="local")
        self.assertEqual(self.w.clientStoreKey(), "")
        self.assertFalse(self.w.fetchClient())


class Scripts(unittest.TestCase):
    def test_rollout_sh_parses(self):
        subprocess.run(["bash", "-n", str(HERE / "rollout.sh")], check=True)

    def test_build_sh_parses(self):
        subprocess.run(["bash", "-n", str(HERE / "build.sh")], check=True)

    def test_bootstrap_still_points(self):
        text = (HERE / "bootstrap.sh").read_text()
        self.assertIn("POINT_CAMPAIGN=1", text)


if __name__ == "__main__":
    unittest.main()
