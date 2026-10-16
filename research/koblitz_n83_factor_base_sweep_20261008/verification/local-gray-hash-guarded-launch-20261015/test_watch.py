"""No-launch controls for the N83 quiet-host watcher."""

from __future__ import annotations

import importlib.util
from pathlib import Path
from unittest import TestCase, main


spec = importlib.util.spec_from_file_location("n83_quiet_watch", Path(__file__).with_name("watch.py"))
assert spec is not None and spec.loader is not None
watch = importlib.util.module_from_spec(spec)
spec.loader.exec_module(watch)


def row(**updates: object) -> dict:
    value = {"idle_percent": 40.0, "load_1m": 4.0,
             "logical_cpus": 14, "other_containers": [],
             "system_free_bytes": 20 * (1 << 30)}
    value.update(updates)
    return value


class QuietHostControls(TestCase):
    def test_quiet_requires_cpu_load_container_and_disk_admission(self) -> None:
        self.assertTrue(watch.is_quiet(row()))
        self.assertFalse(watch.is_quiet(row(idle_percent=0.0)))
        self.assertFalse(watch.is_quiet(row(load_1m=20.0)))
        self.assertFalse(watch.is_quiet(row(other_containers=["another-job"])))
        self.assertFalse(watch.is_quiet(row(system_free_bytes=9 * (1 << 30))))

    def test_missing_probe_fields_fail_closed(self) -> None:
        self.assertFalse(watch.is_quiet({}))
        self.assertFalse(watch.is_quiet(row(logical_cpus=0)))
        self.assertFalse(watch.is_quiet(row(idle_percent=None)))


if __name__ == "__main__":
    main()
