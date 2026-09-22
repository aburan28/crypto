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

    # --- the CPU walker beside the GPU ------------------------------------

    def test_cpu_sidecar_ids_stay_inside_the_range_and_off_the_gpu_ids(self):
        self.assertEqual(modal_app.cpuRunId(8000), 9000)
        self.assertEqual(modal_app.cpuRunId(modal_app.MODAL_GPU_RUN_ID_MAX), modal_app.MODAL_RUN_ID_MAX)
        self.assertEqual(modal_app.checkCampaignRunId(8999, withCpu=True), 8999)
        with self.assertRaises(ValueError) as ctx:
            modal_app.checkCampaignRunId(9000, withCpu=True)
        self.assertIn("sidecar", str(ctx.exception))
        # Without a sidecar the whole range is still a GPU run's to take.
        self.assertEqual(modal_app.checkCampaignRunId(9500), 9500)
        self.assertEqual(modal_sync.slot_for_run(modal_app.cpuRunId(8999)), 99999)

    def test_next_free_run_id_ignores_sidecar_ids(self):
        with tempfile.TemporaryDirectory() as root:
            os.makedirs(os.path.join(root, "dp"))
            for r in (8000, 8001, 9000, 9001):
                Path(root, "dp", "curve131-run%d.bin" % r).write_bytes(b"")
            # 9001 is run 8001's CPU walker, not the highest GPU run.
            self.assertEqual(modal_app.nextFreeRunId(131, root), 8002)

    def test_the_host_build_is_chosen_by_the_cpu_flags_it_needs(self):
        with tempfile.TemporaryDirectory() as root:
            for name in modal_app.CPU_BINARIES.values():
                Path(root, name).write_bytes(b"")
            Path(root, "ecc2k130-cpu").write_bytes(b"")
            full = set(modal_app.X86_64_V4_FLAGS) | {"avx2", "sse4_2"}
            self.assertEqual(modal_app.chooseCpuBinary(full, root)[0], "v4")
            # One missing flag is a SIGILL an hour in, not a slower walk.
            self.assertEqual(modal_app.chooseCpuBinary(full - {"avx512vl"}, root)[0], "v3")
            self.assertEqual(modal_app.chooseCpuBinary({"avx2"}, root)[0], "v3")
            os.unlink(os.path.join(root, modal_app.CPU_BINARIES["v3"]))
            level, path = modal_app.chooseCpuBinary({"avx2"}, root)
            self.assertEqual((level, os.path.basename(path)), ("native", "ecc2k130-cpu"))

    def test_cpu_flags_are_read_from_cpuinfo(self):
        with tempfile.NamedTemporaryFile("w", suffix=".cpuinfo", delete=False) as fh:
            fh.write("processor\t: 0\nflags\t\t: fpu vme avx2 avx512f\n")
            path = fh.name
        try:
            self.assertEqual(modal_app.cpuFlags(path), {"fpu", "vme", "avx2", "avx512f"})
        finally:
            os.unlink(path)
        self.assertEqual(modal_app.cpuFlags("/nonexistent/cpuinfo"), set())

    def test_the_host_client_is_run_as_the_aws_cpu_worker_runs_it(self):
        cmd = modal_app.cpuClientCommand("/root/ecc2k130/ecc2k130-cpu-v4", 131, 9000, 32,
                                         "/data/dp/curve131-run9000.bin",
                                         "/data/ckpt/curve131-run9000.ck", 32, 2000000,
                                         checkpointEvery=60, loads=["/data/dp/curve131-run8000.bin"])
        self.assertTrue(cmd.startswith("/root/ecc2k130/ecc2k130-cpu-v4 --curve 131 "))
        for part in ("--run-id 9000", "--threads 32", "--dp-weight 32", "--launches 0",
                     "--verify 0", "--checkpoint-every 60", "--load-max 2000000",
                     "--load /data/dp/curve131-run8000.bin"):
            self.assertIn(part, cmd)
        # The host binary is bitsliced; --packed is a GPU flag and exit 6 here.
        self.assertNotIn("--packed", cmd)

    def test_cpu_threads_default_to_the_image_and_zero_means_none(self):
        self.assertEqual(modal_app.resolveCpuThreads(-1), modal_app.CPU_THREADS)
        self.assertEqual(modal_app.resolveCpuThreads(0), 0)
        self.assertEqual(modal_app.resolveCpuThreads(16), 16)
        # Default off: the container's CPU request stays Modal's default. Read
        # from the source so an ECC_CPU_THREADS in the test's own environment
        # does not decide the verdict.
        src = Path(__file__).with_name("modal_app.py").read_text()
        self.assertIn('int(os.environ.get("ECC_CPU_THREADS", "0"))', src)
        if modal_app.CPU_THREADS == 0:
            self.assertIsNone(modal_app.CPU_CORES_REQUEST)
        else:
            self.assertEqual(modal_app.CPU_CORES_REQUEST, (modal_app.CPU_THREADS + 1) // 2)

    def test_clients_are_exec_d_so_sigterm_reaches_them(self):
        # /bin/sh in the image is dash, which keeps itself between Popen and
        # the client; SIGTERM then kills the shell and orphans the client.
        src = Path(__file__).with_name("modal_app.py").read_text()
        self.assertIn('subprocess.Popen("exec " + cmd, shell=True', src)
        self.assertIn('subprocess.Popen("exec " + self.cmd, shell=True', src)
        # shStream (builds) may keep the plain shell; nothing signals a build.
        self.assertEqual(src.count("subprocess.Popen(cmd, shell=True"), 1)

    def test_the_image_builds_both_host_binaries_and_keeps_v3_as_the_default(self):
        src = Path(__file__).with_name("modal_app.py").read_text()
        self.assertIn("make cpu MARCH=x86-64-v3 && cp ecc2k130-cpu {CPU_BINARIES['v3']}", src)
        self.assertIn("make -B cpu MARCH=x86-64-v4 && cp ecc2k130-cpu {CPU_BINARIES['v4']}", src)
        self.assertIn("cp {CPU_BINARIES['v3']} ecc2k130-cpu", src)
        self.assertIn("cpu=CPU_CORES_REQUEST", src)

    def test_run_sh_routes_the_rules(self):
        script = Path(__file__).with_name("run.sh").read_text()
        self.assertIn("require_run_id", script)
        self.assertIn('--run-id-base "$RUNID"', script)
        self.assertIn("next-run-id", script)
        self.assertIn("--cpu-threads $CPU_THREADS", script)
        self.assertIn('export ECC_CPU_THREADS="$CPU_THREADS"', script)
        # Passes run detached and bounded: a dropped local connection must not
        # stop the app, and a hung client must not stop the loop.
        self.assertIn("modal run --detach", script)
        self.assertIn('timeout --signal=INT --kill-after=60 "$PASS_TIMEOUT"', script)
        self.assertNotRegex(script, r"(?m)^\s+modal run modal_app\.py::(search|fanout)")
        # No implicit run id at the top level any more; the default of 1 is
        # only reached on the other curves, inside require_run_id's else.
        self.assertNotRegex(script, r"(?m)^RUNID=\$\{RUNID:-1\}")
        self.assertRegex(script, r"(?m)^\s+RUNID=\$\{RUNID:-1\}")

    def test_run_search_puts_campaign_files_where_modal_sync_looks(self):
        src = Path(__file__).with_name("modal_app.py").read_text()
        self.assertIn('dpFile = f"{root}/dp/curve{curve}-run{runId}.bin"', src)
        self.assertIn('hdrFile = f"{root}/ckpt/curve{curve}-run{runId}.hdr"', src)
        self.assertIn("runId = checkCampaignRunId(runId, withCpu=cpuThreads > 0)", src)
        self.assertIn("dpWeight = campaignDpWeightFor(dpWeight)", src)


if __name__ == "__main__":
    raise SystemExit(unittest.main())
