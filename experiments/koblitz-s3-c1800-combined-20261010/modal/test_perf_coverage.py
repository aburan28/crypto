"""Synthetic perf scheduling and IPC fidelity checks for the Modal pilot."""

from __future__ import annotations

import os
import struct
import unittest

import analyze
import remote_pilot


def sample(raw_count: int, enabled: int, running: int) -> tuple[dict | None, str | None]:
    read_fd, write_fd = os.pipe()
    try:
        os.write(write_fd, struct.pack("QQQ", raw_count, enabled, running))
    finally:
        os.close(write_fd)
    return remote_pilot.perf_close(read_fd)


def quality(perf: dict) -> dict:
    selected = list(range(16))
    host = {
        "selected_cpus": selected,
        "topology": {"cpus": {str(cpu): {"physical_package_id": "0",
                                             "core_id": str(cpu // 2),
                                             "thread_siblings_list":
                                             f"{2 * (cpu // 2)}-{2 * (cpu // 2) + 1}"}
                               for cpu in selected},
                     "numa_nodes": {"node0": "0-15"}},
        "cpu_governor": {}, "cpu_frequency_limits_khz": {},
    }
    row = {
        "phase_ok": True, "sampled_numa_pages": {},
        "thread_affinity_masks_seen": [], "max_threads_seen": 14,
        "perf": perf,
        "before": {"cgroup": {}, "cpu_ticks": {}, "interrupts": {}},
        "after": {"cgroup": {}, "cpu_ticks": {}, "interrupts": {}},
        "outer_wall_ns": 1_000_000_000,
    }
    return analyze.guest_quality(row, host, True)["checks"]


class PerfCoverageTest(unittest.TestCase):
    def test_scaled_counts_retain_running_fraction(self) -> None:
        event, error = sample(500, 1000, 500)
        self.assertIsNone(error)
        self.assertEqual(event["raw_count"], 500)
        self.assertEqual(event["time_enabled_ns"], 1000)
        self.assertEqual(event["time_running_ns"], 500)
        self.assertEqual(event["running_fraction"], .5)
        self.assertEqual(event["scaled_count"], 1000)

    def test_short_read_is_unavailable(self) -> None:
        read_fd, write_fd = os.pipe()
        os.write(write_fd, struct.pack("Q", 100))
        os.close(write_fd)
        event, error = remote_pilot.perf_close(read_fd)
        self.assertIsNone(event)
        self.assertEqual(error, "short_perf_read:8")

    def test_guest_gate_distinguishes_multiplexing_and_unavailable(self) -> None:
        cycles, _ = sample(100, 1000, 1000)
        instructions, _ = sample(200, 1000, 1000)
        full = remote_pilot.perf_summary(cycles, instructions, None, None)
        self.assertEqual(full["ipc"], 2.0)
        self.assertIs(quality(full)["perf_counters_unmultiplexed"], True)
        multiplexed_cycles, _ = sample(98, 1000, 980)
        multiplexed = remote_pilot.perf_summary(
            multiplexed_cycles, instructions, None, None)
        self.assertIsNone(multiplexed["ipc"])
        self.assertIsNotNone(multiplexed["ipc_scaled_estimate"])
        self.assertIs(quality(multiplexed)["perf_counters_unmultiplexed"], False)
        unavailable = remote_pilot.perf_summary(
            None, instructions, "EACCES:Permission denied", None)
        self.assertIsNone(unavailable["ipc"])
        self.assertIsNone(quality(unavailable)["perf_counters_unmultiplexed"])


if __name__ == "__main__":
    unittest.main()
