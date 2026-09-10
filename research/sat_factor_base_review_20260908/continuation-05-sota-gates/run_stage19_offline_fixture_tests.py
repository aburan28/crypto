#!/usr/bin/env python3
"""Run pinned legacy tests with host-independent synthetic response fixtures.

Only synthetic response persistence gets the recorded CA identity. Real CA validation tests,
transport implementations and committed historical files remain unchanged.
"""
import argparse
import importlib.util
from pathlib import Path
import socket
import subprocess
import unittest
from unittest.mock import patch

SUITES = {
    "legacy": ("180f8e8d3b273009953974bac31b88e32e7d3a0f", "test_stage19_magma_calculator_panel.py"),
    "amendment03": ("5a09482aea0ce57239b9a06a7cf93c03121730a7", "test_stage19_amendment03_magma_calculator_panel.py"),
}
GATE = Path("research/sat_factor_base_review_20260908/continuation-05-sota-gates")


def deny_network(*_args, **_kwargs):
    raise AssertionError("historical fixture test attempted a network connection")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("suite", choices=SUITES)
    parser.add_argument("--case", action="append", default=[])
    args = parser.parse_args()
    commit, filename = SUITES[args.suite]
    observed = subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip()
    if observed != commit:
        parser.error("run this wrapper from the exact pinned preparation checkout")
    path = Path.cwd() / GATE / filename
    if path.read_bytes() != subprocess.check_output(["git", "show", f"{commit}:{GATE / filename}"]):
        parser.error("historical test source differs from its committed bytes")
    with patch.object(socket.socket, "connect", deny_network), patch.object(socket.socket, "connect_ex", deny_network):
        spec = importlib.util.spec_from_file_location("stage19_pinned_fixture_tests", path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        original = module.CHILD.persist_result

        def persist_synthetic_result(output, *arguments, **keywords):
            if Path(output).resolve().is_relative_to(module.CHILD.REPO.resolve()):
                raise AssertionError("fixture persistence must target an external temporary artifact")
            with patch.object(module.CHILD, "ca_bundle_record", module.CHILD.expected_ca_bundle_record):
                return original(output, *arguments, **keywords)

        with patch.object(module.CHILD, "persist_result", persist_synthetic_result):
            loader = unittest.TestLoader()
            suite = loader.loadTestsFromNames(args.case, module) if args.case else loader.loadTestsFromModule(module)
            result = unittest.TextTestRunner(verbosity=2).run(suite)
    raise SystemExit(0 if result.wasSuccessful() else 1)


if __name__ == "__main__":
    main()
