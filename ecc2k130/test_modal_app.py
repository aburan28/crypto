#!/usr/bin/env python3
"""The campaign rules modal_app.py enforces before a GPU does anything.

`modal` is not needed to test them, so it is stubbed: the decorators become
identity functions and the image/volume builders accept anything.
"""
import os
import sys
import tempfile
import types
import unittest
from pathlib import Path


class _Chain:
    def __getattr__(self, name):
        return lambda *a, **k: self


class _App:
    def __init__(self, *a, **k):
        pass

    def function(self, **k):
        return lambda fn: fn

    def local_entrypoint(self):
        return lambda fn: fn


if "modal" not in sys.modules:
    stub = types.ModuleType("modal")
    stub.is_local = lambda: True
    stub.Image = _Chain()
    stub.Volume = _Chain()
    stub.App = _App
    sys.modules["modal"] = stub

sys.path.insert(0, str(Path(__file__).parent))
import modal_app  # noqa: E402
import modal_sync  # noqa: E402


class CampaignRules(unittest.TestCase):
    def test_the_range_matches_modal_sync_and_stays_a_five_digit_slot(self):
        self.assertEqual(modal_app.MODAL_RUN_ID_MIN, modal_sync.MODAL_RUN_ID_MIN)
        self.assertEqual(modal_app.MODAL_RUN_ID_MAX, modal_sync.MODAL_RUN_ID_MAX)
        self.assertEqual(modal_sync.slot_for_run(modal_app.MODAL_RUN_ID_MAX), 99999)

    def test_run_ids_an_aws_slot_can_reach_are_refused(self):
        for bad in (0, 1, 3, 196, 4242, modal_app.MODAL_RUN_ID_MIN - 1, modal_app.MODAL_RUN_ID_MAX + 1):
            with self.assertRaises(ValueError, msg=str(bad)):
                modal_app.checkCampaignRunId(bad)
        self.assertEqual(modal_app.checkCampaignRunId(8000), 8000)
        self.assertEqual(modal_app.checkCampaignRunId("9999"), 9999)

    def test_the_weight_comes_from_campaign_json_and_nothing_else(self):
        self.assertEqual(modal_app.campaignDpWeight(), 32)
        # -1 used to mean "size it to the pass"; on the campaign it means 32.
        self.assertEqual(modal_app.campaignDpWeightFor(-1), 32)
        self.assertEqual(modal_app.campaignDpWeightFor(32), 32)
        with self.assertRaises(ValueError) as ctx:
            modal_app.campaignDpWeightFor(35)
        self.assertIn("stop at their own distinguished point", str(ctx.exception))

    def test_only_the_campaign_curve_is_a_campaign_run(self):
        self.assertTrue(modal_app.isCampaignRun(131))
        self.assertFalse(modal_app.isCampaignRun(131, offCampaign=True))
        self.assertFalse(modal_app.isCampaignRun(97))

    def test_off_campaign_runs_live_where_nothing_uploads(self):
        self.assertEqual(modal_app.dataRoot(131), "/data")
        self.assertEqual(modal_app.dataRoot(131, offCampaign=True), modal_app.OFF_CAMPAIGN_ROOT)
        self.assertNotEqual(modal_app.OFF_CAMPAIGN_ROOT, "/data")
        self.assertTrue(modal_app.OFF_CAMPAIGN_ROOT.startswith("/data/"))
        # modal_sync lists dp/ at the volume root; the off-campaign tree is
        # a subdirectory it never descends into.
        self.assertFalse(modal_sync.CORPUS_RE.match("dp/offcampaign/curve131-run1.bin"))
        self.assertEqual(modal_app.dataRoot(97, offCampaign=True), "/data")

    def test_next_free_run_id_reads_the_volume_and_starts_the_range(self):
        with tempfile.TemporaryDirectory() as root:
            self.assertEqual(modal_app.nextFreeRunId(131, root), 8000)
            self.assertEqual(modal_app.nextFreeRunId(97, root), 1)
            os.makedirs(os.path.join(root, "dp"))
            os.makedirs(os.path.join(root, "ckpt"))
            for name in ("dp/curve131-run3.bin", "dp/curve131-run8000.bin",
                         "ckpt/curve131-run8003.hdr", "dp/curve97-run5.bin",
                         "ckpt/curve131-run8001.ck"):
                Path(root, name).write_bytes(b"")
            # Legacy ids below the range are ignored; the highest in-range id
            # decides, whichever artefact carries it.
            self.assertEqual(modal_app.nextFreeRunId(131, root), 8004)
            self.assertEqual(sorted(modal_app.usedRunIds(131, root)), [3, 8000, 8001, 8003])
            self.assertEqual(modal_app.nextFreeRunId(97, root), 6)

    def test_campaign_entrypoints_demand_an_explicit_run_id(self):
        with self.assertRaises(SystemExit) as ctx:
            modal_app.requireExplicitRunId(131, 0, False)
        self.assertIn("next_run_id", str(ctx.exception))
        with self.assertRaises(ValueError):
            modal_app.requireExplicitRunId(131, 4242, False)
        self.assertEqual(modal_app.requireExplicitRunId(131, 8000, False), 8000)
        # Off-campaign and other curves keep their old freedom.
        self.assertEqual(modal_app.requireExplicitRunId(131, 1, True), 1)
        self.assertEqual(modal_app.requireExplicitRunId(97, 1, False), 1)

    def test_run_sh_routes_the_rules(self):
        script = Path(__file__).with_name("run.sh").read_text()
        self.assertIn("require_run_id", script)
        self.assertIn('--run-id-base "$RUNID"', script)
        self.assertIn("next-run-id", script)
        # No implicit run id at the top level any more; the default of 1 is
        # only reached on the other curves, inside require_run_id's else.
        self.assertNotRegex(script, r"(?m)^RUNID=\$\{RUNID:-1\}")
        self.assertRegex(script, r"(?m)^\s+RUNID=\$\{RUNID:-1\}")

    def test_run_search_puts_campaign_files_where_modal_sync_looks(self):
        src = Path(__file__).with_name("modal_app.py").read_text()
        self.assertIn('dpFile = f"{root}/dp/curve{curve}-run{runId}.bin"', src)
        self.assertIn('hdrFile = f"{root}/ckpt/curve{curve}-run{runId}.hdr"', src)
        self.assertIn("runId = checkCampaignRunId(runId)", src)
        self.assertIn("dpWeight = campaignDpWeightFor(dpWeight)", src)


if __name__ == "__main__":
    raise SystemExit(unittest.main())
