"""Optimization must stop before benchmark checks or remote setup can be skipped."""
from pathlib import Path
import os
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parent
SCRIPTS = (
    "capacity_gate.py", "check_mapping.py", "compile.py",
    "compile_capacity.py", "reference.py", "run_capacity.py", "validate.py",
)
MODES = ((("-O",), None), (("-OO",), None), ((), "1"), ((), "2"))


class AssertionModeTests(unittest.TestCase):
    def check_rejected(self, arguments):
        with tempfile.TemporaryDirectory() as directory:
            for options, optimize in MODES:
                with self.subTest(arguments=arguments, options=options, optimize=optimize):
                    environment = dict(os.environ)
                    environment.pop("PYTHONOPTIMIZE", None)
                    if optimize:
                        environment["PYTHONOPTIMIZE"] = optimize
                    # -S removes third-party packages. Reaching Modal or any
                    # other later import would produce the wrong failure.
                    result = subprocess.run(
                        [sys.executable, "-B", "-S", *options, *arguments],
                        cwd=directory, env=environment, capture_output=True,
                        text=True, timeout=10,
                    )
                    self.assertNotEqual(result.returncode, 0)
                    self.assertIn("RuntimeError: This benchmark requires Python assertions", result.stderr)
                    self.assertNotIn("ModuleNotFoundError", result.stderr)
                    self.assertEqual(result.stdout, "")

    def test_optimized_interpreters_stop_before_imports_or_execution(self):
        for script in SCRIPTS:
            self.check_rejected([str(ROOT / script)])

    def test_optimized_module_imports_are_rejected(self):
        for module in ("capacity_gate", "reference"):
            self.check_rejected([
                "-c", f"import sys; sys.path.insert(0, {str(ROOT)!r}); import {module}",
            ])


if __name__ == "__main__":
    unittest.main()
