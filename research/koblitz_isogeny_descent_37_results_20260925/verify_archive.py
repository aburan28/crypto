#!/usr/bin/env python3
"""Audit the committed n37 archive against its frozen inputs and summary.

This is a custody and internal-consistency check. It does not rerun Sage or
independently prove the algebraic claims reported by the original producer.
"""

import argparse
from collections import defaultdict
import gzip
import hashlib
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SOURCE = ROOT / "research/koblitz_isogeny_descent_37_20260925"
SOURCE_FILES = {
    "experiment.py": SOURCE / "experiment.py",
    "contract.json": SOURCE / "contract.json",
    "matrix.py": ROOT / "research/toy_f5_neighbors_20260924/matrix.py",
}


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def digest(value):
    return sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode())


def check(condition, message):
    if not condition:
        raise ValueError(message)


def summarize(result):
    """Reconstruct experiment.py's aggregate without importing Sage."""
    groups = defaultdict(list)
    for case in result["cases"]:
        groups[(case["split"], case["summands"], case["encoding"], case["model"])].append(case)
    rows = []
    for (split, summands, encoding, model), cases in sorted(groups.items()):
        rows.append({
            "split": split, "summands": summands, "encoding": encoding, "model": model,
            "systems": len(cases),
            "completion_degree_max": max(c["completion_degree"] for c in cases),
            "completion_degree_sum": sum(c["completion_degree"] for c in cases),
            "f4_xors": sum(c["f4_xors"] for c in cases),
            "f5_xors": sum(c["f5_xors"] for c in cases),
            "matrix_rows": sum(c["matrix_rows"] for c in cases),
            "criterion_rows": sum(c["criterion_rows"] for c in cases),
            "all_verified": all(c["status"] == "VERIFIED" for c in cases),
        })
    indexed = {(c["split"], c["summands"], c["target_index"], c["encoding"], c["model"]): c
               for c in result["cases"]}
    paired_higher = 0
    for case in result["cases"]:
        if case["model"] == "degree_73":
            key = (case["split"], case["summands"], case["target_index"], case["encoding"], "source")
            paired_higher += int(case["completion_degree"] > indexed[key]["completion_degree"])
    ratios = []
    for row in rows:
        if row["model"] == "degree_73":
            base = next(x for x in rows if (x["split"], x["summands"], x["encoding"], x["model"])
                        == (row["split"], row["summands"], row["encoding"], "source"))
            ratios.append({
                "split": row["split"], "summands": row["summands"], "encoding": row["encoding"],
                "f5_xor_ratio_to_source": row["f5_xors"] / base["f5_xors"],
            })
    passed = bool(ratios) and all(r["f5_xor_ratio_to_source"] <= 0.90 for r in ratios) and paired_higher == 0
    return {
        "rows": rows, "paired_completion_degree_increases": paired_higher,
        "f5_xor_ratios": ratios, "registered_success_rule_passed": passed,
        "classification": "engineering diagnostic" if passed else "registered criterion failed",
        "full_dlp_speedup": None, "systems": len(result["cases"]),
        "solver_replays": len(result["cases"]) * result["contract"]["repetitions"] * 2,
        "rejected_algebraic_roots": sum(c["rejected_algebraic_roots"] for c in result["cases"]),
    }


def audit():
    evidence = json.loads((HERE / "EVIDENCE.json").read_text())
    check(evidence["schema"] == "ecc2k130-n37-degree73-archive-repair-v1", "unknown evidence schema")
    check((evidence["hosted_run_id"], evidence["artifact_id"], evidence["artifact_name"])
          == (36087300313, 10843942663, "koblitz-isogeny-descent-37"), "unexpected hosted artifact")
    raw_bytes = (HERE / "raw.json.gz").read_bytes()
    summary_bytes = (HERE / "summary.json").read_bytes()
    check(len(raw_bytes) == evidence["restored_raw_bytes"], "raw byte count mismatch")
    check(sha256(raw_bytes) == evidence["restored_raw_sha256"], "raw hash mismatch")
    check(sha256(summary_bytes) == evidence["summary_sha256"], "summary hash mismatch")
    source_hashes = {name: sha256(path.read_bytes()) for name, path in SOURCE_FILES.items()}
    check(source_hashes == evidence["source_file_sha256"], "frozen source hashes changed")
    try:
        result = json.loads(gzip.decompress(raw_bytes))
    except (OSError, ValueError) as error:
        raise ValueError("raw file is not valid gzip JSON") from error
    summary = json.loads(summary_bytes)
    contract_bytes = SOURCE_FILES["contract.json"].read_bytes()
    contract = json.loads(contract_bytes)
    check(result["contract"] == contract, "raw contract differs from frozen contract")
    check(result["contract_sha256"] == sha256(contract_bytes), "contract hash mismatch")
    check(result["source_sha256"] == source_hashes, "raw source hashes differ from frozen source")

    cert = result["isogeny_certificate"]
    check(contract["field_exponent"] == 37 and contract["prime_degree"] == 73, "wrong registered family")
    check(cert["degree"] == 73 and cert["separable"] is True and cert["inventory_count"] == 74,
          "degree/inventory certificate mismatch")
    check(cert["kernel_polynomial_degree"] == 36 and
          cert["kernel_division_polynomial_quotient_degree"] == 2628, "kernel certificate mismatch")
    check(sha256(cert["kernel_polynomial"].encode()) == cert["kernel_polynomial_sha256"],
          "kernel string hash mismatch")
    check(cert["source_group_order"] == cert["target_group_order"] == contract["source_group_order"],
          "point-count certificate mismatch")
    check(cert["trace_over_gf_2_37"] == -534059 and
          cert["frobenius_order_discriminant"] == -7 * 194399**2 and
          cert["frobenius_order_conductor"] == 73 * 2663, "Frobenius-order certificate mismatch")
    check(cert["source_endomorphism_order_discriminant"] == -7 and
          cert["target_endomorphism_order_discriminant"] == -7 * 73**2 and
          cert["legendre_symbol_minus7_mod_73"] == -1, "endomorphism certificate mismatch")
    check(cert["homomorphism_checked_pairs"] == 128, "homomorphism replay count mismatch")

    cases = result["cases"]
    expected = {(split, summands, target_index, encoding, model)
                for split in contract["seeds"] for summands in contract["summands"]
                for target_index in range(contract["targets_per_summand_count"])
                for encoding in contract["encodings"] for model in ("source", "degree_73")}
    indexed = {}
    for case in cases:
        key = (case["split"], case["summands"], case["target_index"], case["encoding"], case["model"])
        check(key in expected and key not in indexed, f"duplicate or unexpected case: {key}")
        indexed[key] = case
        check(case["seed"] == contract["seeds"][case["split"]], f"wrong seed: {key}")
        check(case["target_pattern"] == contract["target_patterns"][str(case["summands"])][case["target_index"]],
              f"wrong target pattern: {key}")
        check(case["status"] == "VERIFIED" and case["true_selector_assignments"] > 0 and
              case["verified_signed_lifts"] >= case["true_selector_assignments"], f"unverified case: {key}")
        check(case["rejected_algebraic_roots"] == 0, f"unaccounted roots: {key}")
        check(len(case["support_scalars"]) == 4 and
              all(0 < k < contract["factor_base_subgroup_order"] for k in case["support_scalars"]),
              f"invalid support scalars: {key}")
        check(len(case["generators"]) == 37 + (case["encoding"] == "canonical") and
              digest(case["generators"]) == case["input_sha256"], f"input digest mismatch: {key}")
        check(len(case["repetition_sha256"]) == contract["repetitions"] and
              len(set(case["repetition_sha256"])) == 1, f"matrix replay mismatch: {key}")
    check(set(indexed) == expected and len(cases) == 128, "incomplete case grid")
    for (split, summands, target_index, encoding, model), case in indexed.items():
        if model == "degree_73":
            source = indexed[(split, summands, target_index, encoding, "source")]
            for field in ("support_scalars", "target_pattern", "target_is_infinity", "truth_sha256",
                          "true_selector_assignments", "verified_signed_lifts"):
                check(case[field] == source[field], f"transported workload differs in {field}")
    check(summary == summarize(result), "summary does not reconstruct from raw cases")
    check(result["accounting"]["matrix_stage_only"] is True and
          all(result["accounting"][key] is None for key in
              ("full_dlp_total_operations", "S", "rho_ratio", "floor_ratio", "speedup")) and
          result["accounting"]["matrix_solver_repetitions_are_independent_samples"] is False,
          "stage-only accounting boundary changed")
    check(summary["classification"] == "registered criterion failed" and
          summary["registered_success_rule_passed"] is False, "registered verdict changed")
    return {
        "schema": "ecc2k130-n37-degree73-archive-replay-v1",
        "status": "PASS",
        "raw_sha256": sha256(raw_bytes),
        "summary_sha256": sha256(summary_bytes),
        "source_file_sha256": source_hashes,
        "cases": len(cases),
        "paired_workloads": len(cases) // 2,
        "solver_replays": summary["solver_replays"],
        "rejected_algebraic_roots": summary["rejected_algebraic_roots"],
        "classification": summary["classification"],
        "full_dlp_speedup": None,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, help="write deterministic JSON receipt; refuses overwrite")
    args = parser.parse_args()
    receipt = audit()
    payload = json.dumps(receipt, indent=2) + "\n"
    if args.out:
        with args.out.open("x") as output:
            output.write(payload)
    else:
        print(payload, end="")


if __name__ == "__main__":
    main()
