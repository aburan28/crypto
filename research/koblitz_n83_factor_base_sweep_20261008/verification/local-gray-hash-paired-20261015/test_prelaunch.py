"""No-launch admission controls for the local paired gate."""

from __future__ import annotations

import importlib.util
from pathlib import Path
from unittest import TestCase, main


spec = importlib.util.spec_from_file_location("local_prelaunch", Path(__file__).with_name("prelaunch.py"))
assert spec is not None and spec.loader is not None
gate = importlib.util.module_from_spec(spec)
spec.loader.exec_module(gate)


def row(stdout: str, exit_code: int = 0) -> dict:
    return {"stdout": stdout, "exit_code": exit_code}


class PrelaunchControls(TestCase):
    def inputs(self, idle: str = "50.0", load: str = "2.0") -> tuple:
        return (row(f"CPU usage: 40.0% user, 10.0% sys, {idle}% idle"),
                row(f"{{ {load} 2.0 2.0 }}"), row("14\n"),
                row("conductor-db\n"), row(f"{gate.IMAGE} arm64 linux\n"),
                row("4294967296\n0\n"), 20 * (1 << 30))

    def test_busy_host_cannot_launch(self) -> None:
        status, values = gate.classify(*self.inputs(idle="0.0", load="61.0"))
        self.assertEqual(status, "HOST_BUSY_NO_SOLVER_LAUNCH")
        self.assertEqual(values["idle_percent"], 0.0)

    def test_identity_or_probe_failure_cannot_launch(self) -> None:
        inputs = list(self.inputs())
        inputs[5] = row("4294967296\n1073741824\n")
        self.assertEqual(gate.classify(*inputs)[0], "PRODUCER_FAILURE_PRELAUNCH_IDENTITY")
        inputs[5] = row("", exit_code=1)
        self.assertEqual(gate.classify(*inputs)[0], "PRODUCER_FAILURE_PRELAUNCH_PROBE")

    def test_idle_is_admitted_only_without_other_containers(self) -> None:
        inputs = list(self.inputs())
        self.assertEqual(gate.classify(*inputs)[0], "PASS_LOCAL_LINUX_ADMISSION")
        inputs[3] = row("conductor-db\nunrelated-build\n")
        self.assertEqual(gate.classify(*inputs)[0], "HOST_BUSY_NO_SOLVER_LAUNCH")


if __name__ == "__main__":
    main()
