"""Adversarial controls for the independent v2 raw-receipt comparison."""

import importlib.util
import json
from pathlib import Path
import shutil
import tempfile
import unittest


HERE = Path(__file__).resolve().parent
REGISTRATION = HERE / "goal_20260924/generic-backend-qualification-v2"
SPEC = importlib.util.spec_from_file_location(
    "independent_pairs", REGISTRATION / "independent_pairs.py")
PAIRS = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(PAIRS)


def write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, sort_keys=True) + "\n")


def make_bundle(root):
    shutil.copy2(REGISTRATION / "panel.json", root / "registered-panel.json")
    shutil.copy2(REGISTRATION / "lost-campaign-exposures.json",
                 root / "prior-censored-exposures.json")
    tournament = root / "tournament"
    write(tournament / "contract.json", dict(seed=2026092902, repetitions=1,
        cells=list(PAIRS.CELLS), stages=list(PAIRS.STAGE_COUNTS)))
    fixtures = {}
    serial = 0
    for stage, per_cell in (("aa", 1), ("smoke", 1), ("development", 3)):
        fixtures[stage] = []
        for cell in PAIRS.CELLS:
            for index in range(per_cell):
                serial += 1
                fixtures[stage].append(dict(id=f"{cell}-{index:03d}", cell=cell,
                    fixture=dict(targets=[[str(serial), str(1000 + serial)]])))
    write(tournament / "fixtures.json", fixtures)
    for stage, cases in fixtures.items():
        online = []
        for case in cases:
            arms = (("incumbent", "aa_control") if stage == "aa"
                    else (*PAIRS.IC, *PAIRS.RHO))
            for arm in arms:
                mode = "rho" if arm in PAIRS.RHO else "ic"
                cert = dict(verified_targets=1, solutions=["7"])
                timing = dict(target_count=1, unit="native_monotonic_ns",
                    online=dict(wall_ns=2, phase_wall_ns=dict(target_query=1,
                        recovery_check=1), scalar_replay_included=True,
                        target_generation_included=False),
                    cold=dict(wall_ns=5, phase_wall_ns=dict(setup=3,
                        target_query=1, recovery_check=1)))
                measured = dict(status="complete", certificate=cert,
                    native_timing=timing, native_wall_ns=5, total_operations=20,
                    workload_id=f"workload-{stage}-{case['id']}",
                    run_id=f"run-{stage}-{case['id']}-{arm}",
                    provenance=dict(host_id="host", resource_envelope_id="resources"))
                if mode == "ic":
                    measured["candidate_id"] = f"candidate-{case['cell']}-{arm}"
                else:
                    measured["reference_id"] = f"reference-{case['cell']}-{arm}"
                row = dict(stage=stage, case=case["id"], cell=case["cell"], arm=arm,
                    repetition=0, status="VERIFIED", mode=mode, certificate=cert,
                    case_sha256=f"case-{stage}-{case['id']}", measurement=measured,
                    total_operations=20, phase_costs=dict(setup=10, solve=10),
                    native_process=dict(process_wall_ns=5, process_status="EXITED",
                                        exit_code=0),
                    profile_process=dict(process_status="EXITED", exit_code=0))
                directory = tournament / "runs" / stage / case["id"] / arm / "rep-0"
                write(directory / "receipt.json", row)
                write(directory / "native/stdout.json", dict(status="complete",
                    mode=mode, fixture=case["fixture"],
                    solutions=[dict(recovered="7")]))
            if stage != "aa":
                for arm in PAIRS.IC:
                    for rho in PAIRS.RHO:
                        online.append(dict(case=case["id"], arm=arm, rho_alias=rho,
                            public_target=case["fixture"]["targets"][0],
                            verified=True, workload_ids=[f"workload-{stage}-{case['id']}"],
                            run_ids=[f"run-{stage}-{case['id']}-{arm}",
                                     f"run-{stage}-{case['id']}-{rho}"],
                            candidate_ids=[f"candidate-{case['cell']}-{arm}"],
                            rho_reference_ids=[f"reference-{case['cell']}-{rho}"],
                            IC_online_ms=0.000002, rho_online_ms=0.000002,
                            online_speedup=1.0))
        if stage != "aa":
            write(tournament / "summaries" / (stage + ".json"),
                  dict(single_target_online=online))
    return tournament, fixtures


class IndependentPairsTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.bundle = Path(self.temp.name)
        self.round, self.fixtures = make_bundle(self.bundle)

    def receipt(self, stage="development", arm="generic_f4_dense"):
        case = self.fixtures[stage][0]
        return self.round / "runs" / stage / case["id"] / arm / "rep-0/receipt.json"

    def test_complete_schedule_recomputes_all_400_pairs(self):
        result = PAIRS.audit(self.bundle)
        self.assertEqual((result["status"], result["retained_receipts"],
                          result["paired_rows_checked"]),
                         ("AUDITED_COMPLETE_SCHEDULE", 250, 400))
        self.assertEqual(result["paired_online"]["generic_f4_dense"]["rho_online"]
                         ["equal_cell_geometric_speedup"], 1.0)

    def test_missing_receipt_censors_the_panel(self):
        self.receipt().unlink()
        result = PAIRS.audit(self.bundle)
        self.assertEqual(result["status"], "PARTIAL_CAMPAIGN")
        self.assertIsNone(result["paired_online"])

    def test_failed_arm_remains_a_row_without_aggregate_speedup(self):
        path = self.receipt()
        row = PAIRS.read(path)
        row["status"] = "ERROR"
        row["total_operations"] = None
        row["measurement"]["native_timing"] = None
        row["measurement"]["total_operations"] = None
        write(path, row)
        summary_path = self.round / "summaries/development.json"
        summary = PAIRS.read(summary_path)
        for item in summary["single_target_online"]:
            if item["case"] == row["case"] and item["arm"] == row["arm"]:
                item.update(verified=False, IC_online_ms=None, rho_online_ms=None,
                            online_speedup=None)
        write(summary_path, summary)
        result = PAIRS.audit(self.bundle)
        self.assertEqual(result["failure_statuses"], {"ERROR": 1})
        self.assertIsNone(result["paired_online"][row["arm"]]["rho_online"]
                          ["equal_cell_geometric_speedup"])

    def test_rejects_mismatched_resource_envelope(self):
        path = self.receipt()
        row = PAIRS.read(path)
        row["measurement"]["provenance"]["resource_envelope_id"] = "other"
        write(path, row)
        with self.assertRaisesRegex(ValueError, "hosts or resource"):
            PAIRS.audit(self.bundle)

    def test_rejects_unclosed_online_interval(self):
        path = self.receipt()
        row = PAIRS.read(path)
        row["measurement"]["native_timing"]["online"]["wall_ns"] += 1
        write(path, row)
        with self.assertRaisesRegex(ValueError, "interval does not close"):
            PAIRS.audit(self.bundle)

    def test_rejects_native_scalar_mismatch(self):
        path = self.receipt().parent / "native/stdout.json"
        report = PAIRS.read(path)
        report["solutions"][0]["recovered"] = "8"
        write(path, report)
        with self.assertRaisesRegex(ValueError, "returned another scalar"):
            PAIRS.audit(self.bundle)

    def test_rejects_optimistic_published_speedup(self):
        path = self.round / "summaries/development.json"
        summary = PAIRS.read(path)
        summary["single_target_online"][0]["online_speedup"] = 2.0
        write(path, summary)
        with self.assertRaisesRegex(ValueError, "differs from raw receipts"):
            PAIRS.audit(self.bundle)


if __name__ == "__main__":
    unittest.main()
