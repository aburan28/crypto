#!/usr/bin/env python3
"""Focused control-plane tests for the Stage-21 relation-yield bridge."""

from __future__ import annotations

from copy import deepcopy
import json
import os
from pathlib import Path
import sys
import subprocess
import tempfile
import unittest

import run_koblitz_relation_yield_bridge as bridge


def valid_result(protocol: dict, production: bool) -> dict:
    counts = bridge.PRODUCTION_COUNTS if production else bridge.SMOKE_COUNTS
    n, curve_a = (23, 0) if production else (7, 1)
    factor_points = 4281 if production else 15
    dimension = 12 if production else 3
    abscissae = 4096 if production else 8
    natural_hits = 0
    arm_hits = {"natural": natural_hits, "planted_sat": counts["planted_sat"], "proven_unsat": 0}
    rows = []
    ordinal = 0
    for name, count in counts.items():
        for arm_ordinal in range(count):
            hit = arm_ordinal < arm_hits[name]
            planted = name == "planted_sat"
            packed = ordinal + 1
            encoded = packed - 1
            x, y = encoded >> n, encoded & ((1 << n) - 1)
            rows.append({
                "global_ordinal": ordinal, "arm": name, "arm_ordinal": arm_ordinal,
                "packed_target": packed, "x_hex": f"{x:0{(n + 3) // 4}x}",
                "y_hex": f"{y:0{(n + 3) // 4}x}", "x_hamming_weight": x.bit_count(),
                "y_hamming_weight": y.bit_count(), "frobenius_orbit_length": 1,
                "candidate_attempts": arm_ordinal + 1,
                "candidate_kind": "complete_canonical_pair_priority_scan" if planted else "domain_separated_blake3_uniform_affine_then_cofactor_projection",
                "candidate_attempts_scope": "one-based canonical-pair scan position of the selected construction witness" if planted else "hash candidates since the preceding accepted target in this arm",
                "selection_counter": arm_ordinal,
                "selection_priority_blake3": "4" * 64 if planted else None,
                "construction_witness_indices": [0, 0] if planted else None,
                "construction_witness_verified": True if planted else None,
                "exact_pair_table_hit": hit,
                "verified_witness_indices": [0, 0] if hit else None,
                "point_witness_verified": True if hit else None,
            })
            ordinal += 1
    policy = {
        "schema": "koblitz_relation_yield_policy.v1", "n": n, "a": curve_a, "m": 2,
        "irreducible_low_terms": [1],
        "group_order": protocol["curve"]["group_order"] if production else "142",
        "subgroup_order": protocol["curve"]["subgroup_order"] if production else "71",
        "cofactor": protocol["curve"]["cofactor"] if production else "2",
        "divisor_indices": [0, 2], "target_mix": dict(counts),
        "seed_hex": {name: protocol["target_arms"][name]["seed_hex"] for name in counts},
        "point_encoding": "u64 little-endian; 0 for infinity; 1 + (x << n | y) for affine",
        "hash_to_curve": {
            "hash": "BLAKE3", "domain_hex": "00", "message": "domain",
            "x_draw": "uniform", "lift": "deterministic", "projection": "public cofactor",
            "rejections": [], "target_scalar_constructed_or_recorded": False,
        },
        "selection": {
            "natural": "first distinct nonidentity hash-to-curve/cofactor targets; pair oracle unavailable and unused",
            "planted_sat": "pair priority", "proven_unsat": "pair-table absence",
            "arms_disjoint": True, "hash_arm_attempt_cap_formula": "10000*requested + 10000",
        },
        "pair_priority": {"hash": "BLAKE3", "domain_hex": "00", "message": "domain"},
        "oracle": "all canonical i<=j factor-base pairs exactly once, with repetition allowed",
    }
    predicate = {
        "curve": {"n": n, "a": curve_a, "irreducible_low_terms": [1]},
        "group_order": policy["group_order"], "subgroup_order": policy["subgroup_order"],
        "cofactor": policy["cofactor"], "m": 2,
        "construction_method": "linearized kernel of a frozen divisor of T^n-1 over F2",
        "divisor_indices": [0, 2], "divisor_polynomial_bitmask": 5279 if production else 5,
        "linearised_exponents": protocol["factor_base"]["linearised_exponents"] if production else [0, 1, 3],
        "dimension": dimension, "abscissae": abscissae, "rational_points": factor_points,
        "signed_frobenius_orbits_before_projection": 95 if production else 3,
        "projected_signed_frobenius_columns": 93 if production else 3,
    }
    policy_hash, predicate_hash, rows_hash = map(
        bridge.json_blake3, (policy, predicate, rows)
    )
    factor_hash, pair_hash = "c" * 64, "d" * 64
    arms = {}
    arm_hashes = {}
    field_names = {
        "natural": ("natural_targets", "natural_witness_results"),
        "planted_sat": ("planted_targets", "planted_witness_results"),
        "proven_unsat": ("proven_unsat_targets", "proven_unsat_witness_results"),
    }
    canonical_pairs = factor_points * (factor_points + 1) // 2
    for name, count in counts.items():
        selected = [row for row in rows if row["arm"] == name]
        target_bytes = bytearray(bridge.ARM_TARGET_HASH_DOMAINS[name])
        witness_bytes = bytearray(bridge.WITNESS_HASH_DOMAIN)
        for row in selected:
            packed = row["packed_target"].to_bytes(8, "little")
            target_bytes.extend(packed); witness_bytes.extend(packed)
            witness = row["verified_witness_indices"]
            witness_bytes.extend((2**64 - 1).to_bytes(8, "little") if witness is None else witness[0].to_bytes(4, "little") + witness[1].to_bytes(4, "little"))
        target_hash, witness_hash = bridge.blake3_hex(bytes(target_bytes)), bridge.blake3_hex(bytes(witness_bytes))
        arm_hashes[field_names[name][0]], arm_hashes[field_names[name][1]] = target_hash, witness_hash
        hits = arm_hits[name]
        low, high = bridge.wilson_interval_95(hits, count)
        generation = {field: 0 for field in bridge.GENERATION_FIELDS}
        generation["candidates"] = canonical_pairs if name == "planted_sat" else count
        if name != "planted_sat":
            generation.update({"hash_candidates": count, "uniform_affine_decode_successes": count, "cofactor_projection_scalar_multiplications": count})
            if name == "proven_unsat": generation["selection_pair_table_lookups"] = count
        arms[name] = {
            "count": count, "hits": hits, "misses": count - hits, "hit_rate": hits / count,
            "wilson_95_ci": {"low": low, "high": high}, "targets_blake3": target_hash,
            "witness_results_blake3": witness_hash, "verified_point_witnesses": hits,
            "timing_ns": {"target_generation_and_selection": 1, "exact_pair_table_queries": 1, "point_witness_verification": 1},
            "high_level_operations": {"generation": generation, "final_pair_table_lookups": count, "point_witness_verification_additions": hits},
        }
    result_binding = {
        "policy_blake3": policy_hash, "predicate_blake3": predicate_hash,
        "factor_base_blake3": factor_hash, "canonical_pair_transcript_blake3": pair_hash,
        "ordered_covariate_rows_blake3": rows_hash, "arm_hashes": arm_hashes,
        "natural_hits": natural_hits, "planted_hits": counts["planted_sat"], "proven_unsat_hits": 0,
    }
    unique_targets, capacity = canonical_pairs - 1, canonical_pairs
    field_bytes, encoded_bytes = (n + 7) // 8, 1 + 2 * ((n + 7) // 8)
    factor_payload = abscissae * field_bytes + dimension * field_bytes + factor_points * encoded_bytes
    clone_payload, capacity_payload = factor_points * encoded_bytes, capacity * 16
    target_payload = sum(counts.values()) * encoded_bytes
    timers = {
        "clock": "std::time::Instant monotonic elapsed wall time",
        "curve_and_subgroup_construction": 1, "factor_base_predicate_and_materialization": 1,
        "factor_base_projected_column_census": 1, "natural_target_generation": 1,
        "cofactor_class_construction": 1, "canonical_pair_table_and_planted_priority_selection": 1,
        "planted_construction_witness_and_subgroup_verification": 1,
        "proven_unsat_target_screening": 1, "covariate_extraction": 1, "end_to_end": 10,
        "overlap_note": "planted target selection is performed inside the canonical pair-table scan and is not an additive stage",
    }
    top_arm = lambda name: {
        "generation": arms[name]["high_level_operations"]["generation"],
        "final_pair_table_lookups": counts[name],
        "witness_verification_group_additions": arm_hits[name],
    }
    return {
        "schema": bridge.RESULT_SCHEMA, "status": "complete",
        "evidence_class": "finite_public_synthetic_measurement_pending_independent_replay_and_external_resource_receipt" if production else "operational_smoke_only_ineligible_for_ledger_promotion",
        "production_defaults_used": production, "ledger_promotion_eligible_from_this_output_alone": False,
        "measurement_schema": {"stage": "relation_yield", "n_or_bits": n, "base_id_or_hash": factor_hash,
            "eta_or_coverage_policy": {"m": 2, "target_distribution": "synthetic", "oracle": "exact", "policy_blake3": policy_hash},
            "pr_decomposition_or_hit_rate_with_ci": {"natural_hit_rate": arms["natural"]["hit_rate"], "wilson_95_ci": arms["natural"]["wilson_95_ci"]},
            "trials_per_relation": "infinity", "target_mix": dict(counts)},
        "frozen_policy": policy,
        "factor_base": {"predicate": predicate, "predicate_blake3": predicate_hash, "factor_base_blake3": factor_hash,
            "point_order": "ascending exact packed point key", "selection_and_construction_boundaries": {
                "public_field_and_curve_parameters_only": True, "target_available_during_selection": False,
                "target_subgroup_enumerated_for_factor_base": False, "factor_base_discrete_log_labels_constructed": False,
                "target_discrete_log_labels_constructed": False, "relation_yield_used_for_selection": False,
                "solver_timing_used_for_selection": False}},
        "exact_pair_table": {"summands": 2, "repeated_factor_points_allowed": True, "canonical_pair_policy": "all indices i<=j exactly once",
            "canonical_pairs": canonical_pairs, "enumerated_pairs": canonical_pairs, "unique_target_entries": unique_targets,
            "hash_table_capacity": capacity, "duplicate_pair_sums": 1, "packed_key": "exact", "canonical_pair_transcript_blake3": pair_hash,
            "unsat_proof_boundary": "exact m=2 only"},
        "arms": arms, "ordered_covariate_rows": rows, "ordered_covariate_rows_blake3": rows_hash,
        "timing_ns": timers,
        "high_level_operation_counts": {"pair_table_group_additions": canonical_pairs,
            "factor_base_projection_scalar_multiplications": factor_points, "cofactor_class_scalar_multiplications": factor_points,
            "cofactor_class_negations": factor_points, "planted_pair_priority_hashes": 0,
            "planted_final_subgroup_scalar_multiplications": counts["planted_sat"],
            "planted_construction_witness_verification_additions": counts["planted_sat"],
            "target_frobenius_applications": sum(counts.values()), "natural": top_arm("natural"),
            "planted_sat": {"canonical_pair_candidates": canonical_pairs, "final_pair_table_lookups": counts["planted_sat"],
                "witness_verification_group_additions": counts["planted_sat"]}, "proven_unsat": top_arm("proven_unsat"), "scope_note": "synthetic"},
        "retained_size_lower_bounds_bytes": {"scope": "synthetic", "field_element_bytes": field_bytes,
            "encoded_point_bytes": encoded_bytes, "factor_base_payload": factor_payload,
            "canonical_factor_point_clone_payload": clone_payload, "pair_table_live_key_and_witness_payload": unique_targets * 16,
            "pair_table_capacity_key_and_witness_payload": capacity_payload, "selected_target_coordinate_payload": target_payload,
            "sum_using_pair_table_capacity_payload": factor_payload + clone_payload + capacity_payload + target_payload},
        "hashes": {**{key: value for key, value in result_binding.items() if key.endswith("blake3")}, "arms": arm_hashes,
            "result_binding_blake3": bridge.json_blake3(result_binding)},
        "resource_accounting_boundary": {"single_core_elapsed_seconds": None, "user_cpu_seconds": None,
            "system_cpu_seconds": None, "total_core_seconds": None, "peak_rss_bytes": None, "host_identity": None,
            "executable_hash": None, "source_revision": None, "external_process_receipt_required": True, "reason": "external"},
        "claim_boundary": "public synthetic only", "non_claims": ["no SOTA or cryptographic-size security conclusion"],
    }


def valid_smoke_result(protocol: dict) -> dict:
    return valid_result(protocol, False)


def valid_production_result(protocol: dict) -> dict:
    return valid_result(protocol, True)


def metrics(command: list[str], *, wall: float, core: float, rss: int, watchdog: float) -> dict:
    return {
        "command": command,
        "returncode": 0,
        "watchdog_seconds": watchdog,
        "timed_out": False,
        "orphan_group_terminated": False,
        "metrics": {
            "wall_seconds": wall,
            "user_seconds": core,
            "system_seconds": 0.0,
            "total_core_seconds": core,
            "single_core_seconds": core,
            "peak_rss_bytes": rss,
            "meter": "fresh-process getrusage(RUSAGE_CHILDREN)",
        },
    }


class Stage21ControlPlaneTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.protocol, _ = bridge.read_json(bridge.DEFAULT_PROTOCOL, "test protocol")

    def test_protocol_freezes_seeds_commands_and_build_policy(self) -> None:
        bridge.validate_protocol(deepcopy(self.protocol), True)
        for mutate in (
            lambda value: value["target_arms"]["natural"].update(seed_hex="00"),
            lambda value: value["producer"].update(production_arguments=["--secret", "1"]),
            lambda value: value["execution"].update(build_offline_and_locked=False),
        ):
            changed = deepcopy(self.protocol)
            mutate(changed)
            with self.assertRaises(bridge.Stage21Error):
                bridge.validate_protocol(changed, True)

    def test_duplicate_json_and_hardlinks_fail_closed(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            duplicate = root / "duplicate.json"
            duplicate.write_text('{"a":1,"a":2}\n')
            with self.assertRaisesRegex(bridge.Stage21Error, "duplicate JSON key"):
                bridge.read_json(duplicate, "duplicate fixture")
            original = root / "original"
            linked = root / "linked"
            original.write_text("x")
            os.link(original, linked)
            with self.assertRaisesRegex(bridge.Stage21Error, "hard-linked"):
                bridge.regular_bytes(original, "hard-link fixture")

    def test_tool_identity_preserves_multicall_shim_basename(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            shim = Path(directory) / "python-shim"
            shim.symlink_to(Path(sys.executable).resolve())
            identity = bridge.tool_identity(shim, "shim fixture")
            self.assertEqual(identity["path"], str(shim.absolute()))
            self.assertEqual(identity["sha256"], bridge.sha256_bytes(Path(sys.executable).resolve().read_bytes()))
            environment = bridge.safe_child_environment(
                cargo_home=Path(directory), rustc=shim
            )
            self.assertEqual(environment["RUSTC"], str(shim.absolute()))

    def test_blake3_vectors_and_derivable_payload_tampering(self) -> None:
        self.assertEqual(
            bridge.blake3_hex(b""),
            "af1349b9f5f9a1a6a0404dea36dcc9499bcb25c9adc112b7cc9a93cae41f3262",
        )
        self.assertEqual(
            bridge.blake3_hex(b"abc"),
            "6437b3ac38465133ffb63b75273a8db548c558465d79db03fd359c6cd5bd9d85",
        )
        value = valid_smoke_result(self.protocol)
        mutations = (
            lambda result: result["frozen_policy"].update(point_encoding="changed"),
            lambda result: result["arms"]["natural"].update(targets_blake3="0" * 64),
            lambda result: result["hashes"].update(result_binding_blake3="0" * 64),
            lambda result: result["ordered_covariate_rows"][0].pop("x_hex"),
            lambda result: result["high_level_operation_counts"].update(target_frobenius_applications=0),
            lambda result: result["retained_size_lower_bounds_bytes"].update(factor_base_payload=0),
        )
        for mutate in mutations:
            changed = deepcopy(value)
            mutate(changed)
            with self.assertRaises(bridge.Stage21Error):
                bridge.validate_yield_result(changed, False, self.protocol)

    def test_exclusive_meter_preserves_preexisting_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            stdout = root / "stdout"
            stderr = root / "stderr"
            metrics_path = root / "metrics.json"
            stdout.write_bytes(b"sentinel")
            completed = subprocess.run(
                [
                    sys.executable,
                    str(bridge.DEFAULT_METER),
                    "--cwd", str(root),
                    "--timeout", "10",
                    "--stdout", str(stdout),
                    "--stderr", str(stderr),
                    "--metrics", str(metrics_path),
                    "--exclusive-create",
                    "--", sys.executable, "-c", "print('must not run')",
                ],
                capture_output=True,
                check=False,
            )
            self.assertNotEqual(completed.returncode, 0)
            self.assertEqual(stdout.read_bytes(), b"sentinel")
            self.assertFalse(stderr.exists())
            self.assertFalse(metrics_path.exists())

    def test_run_and_verification_roots_must_be_outside_repository(self) -> None:
        with self.assertRaisesRegex(bridge.Stage21Error, "unsafe"):
            bridge.safe_new_directory(bridge.REPO / "forbidden-stage21-run", "run")

    def test_smoke_result_enforces_counts_policy_pairs_and_unique_targets(self) -> None:
        value = valid_smoke_result(self.protocol)
        bridge.validate_yield_result(value, False, self.protocol)
        mutations = (
            lambda result: result["arms"]["natural"].update(count=7),
            lambda result: result["frozen_policy"]["selection"].update(natural="oracle filtered"),
            lambda result: result["exact_pair_table"].update(enumerated_pairs=119),
            lambda result: result["ordered_covariate_rows"][1].update(packed_target=1),
            lambda result: result["ordered_covariate_rows"][0].update(target_scalar=7),
        )
        for mutate in mutations:
            changed = deepcopy(value)
            mutate(changed)
            with self.assertRaises(bridge.Stage21Error):
                bridge.validate_yield_result(changed, False, self.protocol)

    def test_production_result_binds_exact_factor_base_and_target_mix(self) -> None:
        value = valid_production_result(self.protocol)
        bridge.validate_yield_result(value, True, self.protocol)
        changed = deepcopy(value)
        changed["factor_base"]["predicate"]["rational_points"] = 4280
        with self.assertRaisesRegex(bridge.Stage21Error, "factor-base identity"):
            bridge.validate_yield_result(changed, True, self.protocol)

    def test_child_environment_excludes_ambient_scalar_and_closes_io(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            environment = bridge.safe_child_environment()
            self.assertNotIn("TARGET_SCALAR", environment)
            executable = bridge.tool_identity(Path(sys.executable), "test Python")
            command = [
                executable["path"],
                "-c",
                (
                    "import json,os,sys\n"
                    "fds=[]\n"
                    "for i in range(3,32):\n"
                    " try:\n"
                    "  os.fstat(i)\n"
                    " except OSError:\n"
                    "  continue\n"
                    " fds.append(i)\n"
                    "print(json.dumps({'scalar':os.environ.get('TARGET_SCALAR'),"
                    "'stdin':sys.stdin.buffer.read(1).hex(),'fds':fds}))\n"
                ),
            ]
            previous = os.environ.get("TARGET_SCALAR")
            os.environ["TARGET_SCALAR"] = "forbidden"
            try:
                receipt, stdout = bridge.run_metered(
                    root,
                    "probe",
                    command,
                    30,
                    [],
                    bridge.DEFAULT_METER,
                    environment,
                    executable,
                )
            finally:
                if previous is None:
                    os.environ.pop("TARGET_SCALAR", None)
                else:
                    os.environ["TARGET_SCALAR"] = previous
            observed = json.loads(stdout.read_text())
            self.assertIsNone(observed["scalar"])
            self.assertEqual(observed["stdin"], "")
            self.assertEqual(observed["fds"], [])
            self.assertEqual(receipt["environment"], environment)

    def test_outer_receipt_must_enclose_children_and_use_frozen_watchdog(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "outer.json"
            command = ["/python", "/runner", "run"]
            seal = {"outer_expected_command": command}
            summary = {"resources": {
                "total_core_seconds": 2.0,
                "summed_process_wall_seconds": 3.0,
                "peak_process_rss_bytes": 10,
            }}
            bridge.write_json_new(path, metrics(
                command, wall=2.0, core=1.0, rss=9,
                watchdog=self.protocol["execution"]["whole_driver_watchdog_seconds"],
            ))
            with self.assertRaises(bridge.Stage21Error):
                bridge.verify_outer(path, seal, summary, self.protocol)
            path.unlink()
            bridge.write_json_new(path, metrics(
                command, wall=4.0, core=3.0, rss=10,
                watchdog=self.protocol["execution"]["whole_driver_watchdog_seconds"],
            ))
            bridge.verify_outer(path, seal, summary, self.protocol)

    def test_driver_command_binds_inner_meter(self) -> None:
        command = bridge.expected_driver_command(
            bridge.DEFAULT_PROTOCOL,
            Path("/tmp/stage21-new-output"),
            bridge.DEFAULT_METER,
            False,
            False,
        )
        self.assertEqual(command.count("--meter"), 1)
        index = command.index("--meter")
        self.assertEqual(command[index + 1], str(bridge.DEFAULT_METER.resolve()))

    def test_cli_plan_exits_cleanly(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            completed = subprocess.run(
                [
                    sys.executable,
                    str(Path(bridge.__file__).resolve()),
                    "plan",
                    "--smoke",
                    "--output",
                    str(root / "run"),
                    "--outer-stdout",
                    str(root / "outer.stdout"),
                    "--outer-stderr",
                    str(root / "outer.stderr"),
                    "--outer-metrics",
                    str(root / "outer.metrics.json"),
                ],
                text=True,
                capture_output=True,
                check=False,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            rendered = json.loads(completed.stdout)
            self.assertEqual(rendered["schema"], "koblitz_relation_yield_execution_plan.v1")
            self.assertEqual(rendered["mode"], "smoke")


if __name__ == "__main__":
    unittest.main()
