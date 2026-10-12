#!/usr/bin/env python3
"""Prospective measured-contamination gate on the archived v1 false censor."""
from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path
import tarfile
import unittest

from verify_cold import observed_isolation_ok

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
RAW = ROOT / ("research/notes/ecc2k130/disjoint_cold_outcome_20261001/"
              "evidence_run_36794339148/raw/n37_L1.tar.gz")
RAW_SHA = "dec433e570a28750d29ab344f7b5405cc97c16f9159f4202f149e46a20b3321f"


class IsolationRegression(unittest.TestCase):
    def test_idle_affinity_is_not_contention(self) -> None:
        self.assertEqual(hashlib.sha256(RAW.read_bytes()).hexdigest(), RAW_SHA)
        with tarfile.open(RAW, "r:gz") as bundle:
            stream = bundle.extractfile("disjoint-cold-n37_L1/n37_L1/isolation.jsonl")
            self.assertIsNotNone(stream)
            record = json.loads(stream.read())
        self.assertEqual(record["contended_samples"], 0)
        self.assertGreater(len(record["left_on_reserved"]["user_threads"]), 0)
        label, cpu = "disjoint-cold-n37_L1", 3
        self.assertTrue(observed_isolation_ok(record, label, cpu))
        self.assertFalse(observed_isolation_ok(record, label + "-wrong", cpu))
        self.assertFalse(observed_isolation_ok(record, label, 1))
        noisy = copy.deepcopy(record)
        noisy["samples"][0]["contended"] = True
        noisy["contended_samples"] = 1
        self.assertFalse(observed_isolation_ok(noisy, label, cpu))
        busy = copy.deepcopy(record)
        busy["preflight"]["settle"]["other_cpu_seconds"] = 0.201
        self.assertFalse(observed_isolation_ok(busy, label, cpu))
        missing = copy.deepcopy(record)
        missing["samples"] = []
        self.assertFalse(observed_isolation_ok(missing, label, cpu))
        no_psi = copy.deepcopy(record)
        no_psi["preflight"]["conditions"].pop("psi_cpu")
        self.assertFalse(observed_isolation_ok(no_psi, label, cpu))


if __name__ == "__main__":
    unittest.main()
