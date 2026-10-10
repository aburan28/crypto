"""Exercise the Modal analyzer with archived T205 reports and synthetic telemetry."""

from __future__ import annotations

import copy
import gzip
import hashlib
import json
from pathlib import Path
import tempfile
import unittest

import analyze
import remote_pilot


HERE = Path(__file__).resolve().parent
LOCAL = HERE.parent / "local/paired_confirm"


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def snapshot() -> dict:
    return {
        "unix_ns": 1_000_000_000,
        "cpu_ticks": {str(cpu): {"steal": 0} for cpu in range(16)},
        "interrupts": {str(cpu): 0 for cpu in range(16)},
        "cgroup": {"cpu.max": "max 100000", "cpu.stat": "nr_throttled 0\nthrottled_usec 0",
                   "cpu.pressure": "some avg10=0.00 total=0",
                   "memory.pressure": "some avg10=0.00 total=0",
                   "io.pressure": "some avg10=0.00 total=0",
                   "memory.events": "oom 0\noom_kill 0"},
    }


def synthetic_receipt() -> dict:
    plan_bytes = analyze.PLAN.read_bytes()
    plan = json.loads(plan_bytes)
    bundle_bytes = analyze.BUNDLE.read_bytes()
    bundle = json.loads(bundle_bytes)
    inventory = [{"remote": "/opt/s3/" + arm,
                  "sha256": plan["binary_sha256"][arm],
                  "bytes": next(item["bytes"] for item in json.loads(
                      (HERE / "TRANSFER_MANIFEST.json").read_bytes())["files"]
                      if item["remote_path"] == "/opt/s3/" + arm)}
                 for arm in ("baseline", "expanded")]
    inventory.append({"remote": "/opt/s3/public_bundle.json",
                      "sha256": sha(bundle_bytes), "bytes": len(bundle_bytes)})
    cycles = {"raw_count": 100, "time_enabled_ns": 1000,
              "time_running_ns": 1000, "running_fraction": 1.0,
              "scaled_count": 100.0}
    instructions = dict(cycles, raw_count=200, scaled_count=200.0)
    perf = remote_pilot.perf_summary(cycles, instructions, None, None)
    before = snapshot()
    after = copy.deepcopy(before)
    after["unix_ns"] += 2_000_000_000
    rows = []
    repetitions = {"baseline": 0, "expanded": 0}
    case = bundle["curves"][0]
    for _, arm in bundle["run_order"]:
        repetitions[arm] += 1
        repetition = repetitions[arm]
        local_arm = "baseline" if arm == "baseline" else "combined"
        path = LOCAL / f"T205_{local_arm}_R{repetition}/report.json.gz"
        report = json.loads(gzip.decompress(path.read_bytes()))
        rows.append({
            "run_id": f"{case['arms'][arm]['candidate_id']}W{case['workload_id']}R{repetition}",
            "candidate_id": case["arms"][arm]["candidate_id"],
            "workload_id": case["workload_id"],
            "native_verified": True, "phase_ok": True, "error": None,
            "outer_wall_ns": 2_000_000_000,
            "sampled_peak_rss_bytes": 25_000_000_000,
            "sampled_numa_pages": {"N0": 1000},
            "max_threads_seen": 14,
            "affinity": list(range(16)),
            "thread_affinity_masks_seen": ["0-15"],
            "before": before, "after": after,
            "perf": perf, "report": report,
        })
    return {
        "kind": "s3_modal_c1800_launch_receipt_v1",
        "freeze_sha256": sha(plan_bytes),
        "transfer_inventory": inventory,
        "error": None,
        "guest": {
            "bundle_sha256": sha(bundle_bytes),
            "binary_sha256": plan["binary_sha256"],
            "run_order": bundle["run_order"], "rows": rows,
            "stop_reason": None,
            "host": {
                "selected_cpus": list(range(16)),
                "topology": {"cpus": {str(cpu): {
                    "physical_package_id": "0", "core_id": str(cpu // 2),
                    "thread_siblings_list": f"{2 * (cpu // 2)}-{2 * (cpu // 2) + 1}"}
                    for cpu in range(16)},
                    "numa_nodes": {"node0": "0-15"}},
                "cpu_governor": {str(cpu): "performance" for cpu in range(16)},
                "cpu_frequency_limits_khz": {str(cpu): {
                    "scaling_min_freq": "2500000", "scaling_max_freq": "2500000"}
                    for cpu in range(16)},
                "idle_before": before, "idle_after": after,
            },
        },
    }


class ReceiptContractTest(unittest.TestCase):
    def call_analyzer(self, receipt: dict) -> dict:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "synthetic_receipt.json"
            path.write_text(json.dumps(receipt, separators=(",", ":")))
            return analyze.analyze(path)

    def test_four_archived_reports_reach_full_analysis(self) -> None:
        result = self.call_analyzer(synthetic_receipt())
        self.assertEqual(result["status"], "verified_complete")
        self.assertEqual(result["completed_rows"], 4)
        self.assertEqual(len(result["curves"]["53"]["baseline_online_ms"]), 2)
        self.assertTrue(all(row["guest_quality"]["checks"][
            "perf_counters_unmultiplexed"] for row in result["rows"]))
        self.assertIsNone(result["controlled_speedup"])

    def test_fixture_mismatch_invalidates_paired_result(self) -> None:
        receipt = synthetic_receipt()
        receipt["guest"]["rows"][0]["report"]["recovered_scalar"] += 1
        result = self.call_analyzer(receipt)
        self.assertEqual(result["status"], "incomplete_or_invalid")
        self.assertFalse(result["rows"][0]["verified_native_and_fixture"])
        self.assertEqual(result["curves"], {})

    def test_interleaved_smt_selection_passes_affinity_check(self) -> None:
        receipt = synthetic_receipt()
        host = receipt["guest"]["host"]
        for cpu in range(16):
            host["topology"]["cpus"][str(cpu)]["core_id"] = str(cpu % 8)
            host["topology"]["cpus"][str(cpu)]["thread_siblings_list"] = (
                f"{cpu % 8},{cpu % 8 + 8}")
        host["selected_cpus"] = remote_pilot.choose_cpus(
            list(range(16)), host["topology"], 16)
        self.assertNotEqual(host["selected_cpus"], sorted(host["selected_cpus"]))
        result = self.call_analyzer(receipt)
        self.assertTrue(all(row["guest_quality"]["checks"][
            "observed_thread_affinity"] for row in result["rows"]))

    def test_cpu_selection_keeps_whole_siblings_with_extra_cores(self) -> None:
        cpus = {str(cpu): {
            "physical_package_id": "0", "core_id": str(cpu % 16),
            "thread_siblings_list": f"{cpu % 16},{cpu % 16 + 16}"}
            for cpu in range(32)}
        selected = remote_pilot.choose_cpus(
            list(range(32)), {"cpus": cpus, "numa_nodes": {"node0": "0-31"}}, 16)
        self.assertEqual(len(selected), 16)
        self.assertEqual(len({cpus[str(cpu)]["core_id"] for cpu in selected}), 8)
        self.assertEqual(set(selected), set(range(8)) | set(range(16, 24)))

    def test_cpu_selection_prefers_one_numa_node(self) -> None:
        cpus = {str(cpu): {
            "physical_package_id": "0", "core_id": str(cpu % 16),
            "thread_siblings_list": f"{cpu % 16},{cpu % 16 + 16}"}
            for cpu in range(32)}
        selected = remote_pilot.choose_cpus(list(range(32)), {
            "cpus": cpus,
            "numa_nodes": {"node0": "0-3,16-19", "node1": "4-11,20-27"},
        }, 16)
        self.assertEqual(set(selected), set(range(4, 12)) | set(range(20, 28)))

    def test_initial_affinity_mismatch_fails_quality_check(self) -> None:
        receipt = synthetic_receipt()
        receipt["guest"]["rows"][0]["affinity"] = list(range(15))
        result = self.call_analyzer(receipt)
        self.assertFalse(result["rows"][0]["guest_quality"]["checks"][
            "observed_thread_affinity"])

    def test_failed_producer_is_preserved(self) -> None:
        receipt = synthetic_receipt()
        receipt["guest"] = None
        receipt["error"] = "synthetic_failure"
        result = self.call_analyzer(receipt)
        self.assertEqual(result["status"], "producer_failure")
        self.assertEqual(result["error"], "synthetic_failure")


if __name__ == "__main__":
    unittest.main()
