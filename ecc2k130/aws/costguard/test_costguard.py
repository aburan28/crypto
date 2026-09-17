#!/usr/bin/env python3
"""Unit tests for costguard policy helpers (no AWS calls)."""
from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
PRICES = ROOT / "prices.json"
HOST_SCRIPT = ROOT / "host-idle-stop.sh"


class PricesContract(unittest.TestCase):
    def test_required_keys(self):
        data = json.loads(PRICES.read_text())
        self.assertIn("on_demand_hourly", data)
        self.assertIn("gpu_count", data)
        west = data["on_demand_hourly"]["us-west-2"]
        # List on-demand, not spot (README: g7e.2xlarge $3.36/h, 48xlarge $33.14/h).
        self.assertGreaterEqual(west["g7e.2xlarge"], 3.36)
        self.assertGreaterEqual(west["g7e.48xlarge"], 33.14)
        self.assertGreaterEqual(west["g7.2xlarge"], 2.52)
        self.assertEqual(data["gpu_count"]["g7e.48xlarge"], 8)

    def test_shell_syntax(self):
        for name in ("costguard.sh", "install-host.sh", "host-idle-stop.sh"):
            subprocess.run(["bash", "-n", str(ROOT / name)], check=True)


class HostIdleLogic(unittest.TestCase):
    def test_exempt_exits_zero(self):
        with tempfile.TemporaryDirectory() as td:
            Path(td, "exempt").write_text("1")
            env = {
                "COSTGUARD_STATE_DIR": td,
                "COSTGUARD_LOG": str(Path(td) / "log"),
                "COSTGUARD_MAX_AGE_HOURS": "1",
                "COSTGUARD_IDLE_HOURS": "1",
                "PATH": "/usr/bin:/bin",
            }
            r = subprocess.run(
                ["bash", str(HOST_SCRIPT)], env=env, capture_output=True, text=True
            )
            self.assertEqual(r.returncode, 0)
            self.assertIn("exempt marker present", Path(td, "log").read_text())

    def test_idle_disabled_by_default_path(self):
        with tempfile.TemporaryDirectory() as td:
            env = {
                "COSTGUARD_STATE_DIR": td,
                "COSTGUARD_LOG": str(Path(td) / "log"),
                "COSTGUARD_MAX_AGE_HOURS": "0",
                "COSTGUARD_IDLE_HOURS": "0",
                "PATH": "/usr/bin:/bin",
            }
            r = subprocess.run(
                ["bash", str(HOST_SCRIPT)], env=env, capture_output=True, text=True
            )
            self.assertEqual(r.returncode, 0)
            self.assertIn("idle stop disabled", Path(td, "log").read_text())


class BudgetDefaults(unittest.TestCase):
    def test_header_documents_5k(self):
        text = (ROOT / "costguard.sh").read_text()
        self.assertIn("MONTHLY_BUDGET_USD=${MONTHLY_BUDGET_USD:-5000}", text)
        self.assertIn("MAX_HOURLY_USD=${MAX_HOURLY_USD:-0}", text)


if __name__ == "__main__":
    unittest.main()
