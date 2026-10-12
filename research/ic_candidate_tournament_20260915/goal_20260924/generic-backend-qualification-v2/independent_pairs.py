"""Independently cross-check v2 online pairs from retained raw receipts.

This post-registration auditor uses only Python's standard library. It does not
replace the frozen group-arithmetic verifier or natural-query auditor. It
checks their retained inputs and recomputes the one-target timing table without
calling tournament.py, qualification.py, or driver_admission.py.
"""

import argparse
from collections import Counter
import hashlib
import json
import math
from pathlib import Path
import statistics


PANEL_SHA256 = "d283a869b0412228d1c66260fdfd8f387d7243bd15456c7febf3c46ee5da27a8"
EXPOSURES_SHA256 = "a728677b199eac02800d8338204d5306f391ec5da757c910bec1e51955fe7b41"
CELLS = ("n17a1", "n19a0", "n23a0", "n23a1", "n31a0")
GENERIC = (
    "generic_pair_dense", "generic_pair_sparse", "generic_f4_dense",
    "generic_f4_sparse", "generic_f5_dense", "generic_sat_xor_dense",
    "generic_sat_cnf_dense", "generic_inherited_f4_dense",
)
IC = ("incumbent", "ic_online", *GENERIC)
RHO = ("rho", "rho_online")
STAGE_COUNTS = {"aa": 5, "smoke": 5, "development": 15}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def read(path):
    return json.loads(path.read_text())


def sha256(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def verified_interval(row, fixture, raw_dir):
    """Return one native online interval only after independent closure checks."""
    if row["status"] != "VERIFIED":
        require(row["measurement"].get("native_timing") is None
                and row["total_operations"] is None,
                "failed receipt has competitive timing or cost")
        return None
    measured = row["measurement"]
    require(measured["status"] == "complete"
            and row["certificate"] == measured["certificate"]
            and measured["certificate"]["verified_targets"] == 1
            and len(measured["certificate"]["solutions"]) == 1,
            "verified receipt lacks one certified target")
    interval = measured["native_timing"]
    require(interval["target_count"] == 1 and interval["unit"] == "native_monotonic_ns",
            "changed single-target timing unit")
    for part in ("online", "cold"):
        phases = interval[part]["phase_wall_ns"]
        require(phases and all(type(v) is int and v >= 0 for v in phases.values())
                and type(interval[part]["wall_ns"]) is int
                and sum(phases.values()) == interval[part]["wall_ns"],
                "native phase interval does not close")
    online, cold = interval["online"], interval["cold"]
    require(online["wall_ns"] > 0 and cold["wall_ns"] >= online["wall_ns"]
            and online["scalar_replay_included"] is True
            and online["target_generation_included"] is False
            and measured["native_wall_ns"] == cold["wall_ns"]
            == row["native_process"]["process_wall_ns"],
            "online boundary or full process wall changed")
    require(row["native_process"]["process_status"] == "EXITED"
            and row["native_process"]["exit_code"] == 0
            and row["profile_process"]["process_status"] == "EXITED"
            and row["profile_process"]["exit_code"] == 0,
            "verified receipt has failed native or profiled process")
    require(all(v is None or (type(v) is int and v >= 0)
                for v in row["phase_costs"].values()),
            "negative or nonintegral instruction phase")
    require(row["total_operations"] == measured["total_operations"]
            == sum(v for v in row["phase_costs"].values() if v is not None)
            and row["total_operations"] > 0,
            "instruction costs do not close")
    native = read(raw_dir / "native/stdout.json")
    require(native["fixture"] == fixture and native["mode"] == row["mode"]
            and native["status"] == "complete" and len(native["solutions"]) == 1
            and native["solutions"][0]["recovered"]
                == measured["certificate"]["solutions"][0],
            "native worker solved a different point or returned another scalar")
    return online["wall_ns"]


def audit(bundle):
    bundle = Path(bundle).resolve()
    require(sha256(bundle / "registered-panel.json") == PANEL_SHA256
            and sha256(bundle / "prior-censored-exposures.json") == EXPOSURES_SHA256,
            "not the registered fresh panel and exclusion corpus")
    panel = read(bundle / "registered-panel.json")
    require(panel["seed"] == 2026092902
            and [c["id"] for c in panel["candidates"]] ==
                ["incumbent", "prepared_both", *GENERIC]
            and [r["id"] for r in panel["rho_arms"]] == list(RHO),
            "changed registered algorithm arms")
    root = bundle / "tournament"
    contract, fixtures = read(root / "contract.json"), read(root / "fixtures.json")
    require(contract["seed"] == panel["seed"] and contract["repetitions"] == 1
            and contract["cells"] == list(CELLS)
            and contract["stages"] == list(STAGE_COUNTS)
            and set(fixtures) == set(STAGE_COUNTS)
            and all(len(fixtures[s]) == count for s, count in STAGE_COUNTS.items()),
            "changed frozen one-process schedule")
    cases = {}
    public_points = set()
    for stage, case_rows in fixtures.items():
        expected_per_cell = 3 if stage == "development" else 1
        require(Counter(row["cell"] for row in case_rows)
                == Counter({cell: expected_per_cell for cell in CELLS}),
                "changed per-cell public-point allocation")
        for case in case_rows:
            key = stage, case["id"]
            require(key not in cases and len(case["fixture"]["targets"]) == 1,
                    "duplicate case or multi-target fixture")
            cases[key] = case
            point = case["cell"], json.dumps(case["fixture"]["targets"][0])
            require(point not in public_points, "reused public point")
            public_points.add(point)
    require(len(public_points) == 25, "not 25 distinct public points")

    expected = {(stage, case["id"], arm, 0)
                for stage, case_rows in fixtures.items() for case in case_rows
                for arm in (("incumbent", "aa_control") if stage == "aa"
                            else (*IC, *RHO))}
    require(len(expected) == 250, "changed registration size")
    receipts = {}
    for path in root.glob("runs/**/receipt.json"):
        parts = path.relative_to(root / "runs").parts
        require(len(parts) == 5 and parts[-1] == "receipt.json"
                and parts[3] == "rep-0", "unexpected receipt path")
        key = parts[0], parts[1], parts[2], 0
        require(key in expected and key not in receipts, "extra or duplicate trial receipt")
        row = read(path)
        require((row["stage"], row["case"], row["arm"], row["repetition"]) == key
                and row["cell"] == cases[key[:2]]["cell"],
                "receipt disagrees with its frozen case or path")
        receipts[key] = row
    missing = sorted(expected - receipts.keys())
    if missing:
        return dict(status="PARTIAL_CAMPAIGN", expected_slots=250,
                    retained_receipts=len(receipts), missing_slots=len(missing),
                    missing_by_stage=dict(sorted(Counter(key[0] for key in missing).items())),
                    verified_receipts=sum(r["status"] == "VERIFIED" for r in receipts.values()),
                    paired_online=None, reason="incomplete schedule; no competitive result")

    require(len({r["measurement"]["provenance"]["host_id"]
                 for r in receipts.values()}) == 1
            and len({r["measurement"]["provenance"]["resource_envelope_id"]
                     for r in receipts.values()}) == 1,
            "campaign mixed hosts or resource envelopes")

    intervals = {}
    for key, row in receipts.items():
        case = cases[key[:2]]
        intervals[key] = verified_interval(row, case["fixture"],
                                           root / "runs" / key[0] / key[1] / key[2] / "rep-0")
    run_ids = [r["measurement"]["run_id"] for r in receipts.values()]
    require(len(set(run_ids)) == len(run_ids), "colliding raw execution identities")
    pairs = []
    for stage in ("smoke", "development"):
        saved = read(root / "summaries" / (stage + ".json"))["single_target_online"]
        saved_by_key = {(r["case"], r["arm"], r["rho_alias"]): r for r in saved}
        require(len(saved_by_key) == len(saved)
                == STAGE_COUNTS[stage] * len(IC) * len(RHO),
                "published online table has missing or duplicate rows")
        for case in fixtures[stage]:
            source = [receipts[stage, case["id"], arm, 0] for arm in (*IC, *RHO)]
            require(len({r["measurement"]["workload_id"] for r in source}) == 1
                    and len({r["case_sha256"] for r in source}) == 1
                    and len({r["measurement"]["provenance"]["host_id"] for r in source}) == 1
                    and len({r["measurement"]["provenance"]["resource_envelope_id"]
                             for r in source}) == 1,
                    "paired arms have different point, host or resource envelope")
            for arm in IC:
                ic = intervals[stage, case["id"], arm, 0]
                for rho_arm in RHO:
                    rho = intervals[stage, case["id"], rho_arm, 0]
                    stored = saved_by_key[case["id"], arm, rho_arm]
                    verified = ic is not None and rho is not None
                    ic_record = receipts[stage, case["id"], arm, 0]["measurement"]
                    rho_record = receipts[stage, case["id"], rho_arm, 0]["measurement"]
                    require(stored["public_target"] == case["fixture"]["targets"][0]
                            and stored["verified"] is verified
                            and stored["workload_ids"] ==
                                [source[0]["measurement"]["workload_id"]]
                            and stored["run_ids"] ==
                                [ic_record["run_id"], rho_record["run_id"]]
                            and stored["candidate_ids"] == [ic_record["candidate_id"]]
                            and stored["rho_reference_ids"] ==
                                [rho_record["reference_id"]],
                            "published online pair changed point or verification")
                    expected_ic = ic / 1_000_000 if verified else None
                    expected_rho = rho / 1_000_000 if verified else None
                    expected_ratio = rho / ic if verified else None
                    require(stored["IC_online_ms"] == expected_ic
                            and stored["rho_online_ms"] == expected_rho
                            and stored["online_speedup"] == expected_ratio,
                            "published paired online cost differs from raw receipts")
                    pairs.append(dict(stage=stage, case=case["id"], cell=case["cell"],
                                      arm=arm, rho_alias=rho_arm, verified=verified,
                                      ic_online_ns=ic if verified else None,
                                      rho_online_ns=rho if verified else None,
                                      rho_over_ic=expected_ratio))
    require(len(pairs) == 400, "not all registered online pairs checked")
    complete = {}
    for arm in IC:
        complete[arm] = {}
        for rho_arm in RHO:
            selected = [p for p in pairs if p["arm"] == arm and p["rho_alias"] == rho_arm]
            require(len(selected) == 20, "lost a public point in paired comparison")
            ratios = [p["rho_over_ic"] for p in selected if p["verified"]]
            complete[arm][rho_arm] = dict(verified_points=len(ratios),
                all_points_verified=len(ratios) == 20,
                equal_cell_geometric_speedup=(math.exp(statistics.mean(
                    statistics.mean(math.log(p["rho_over_ic"]) for p in selected
                                    if p["cell"] == cell)
                    for cell in CELLS)) if len(ratios) == 20 else None))
    return dict(status="AUDITED_COMPLETE_SCHEDULE", expected_slots=250,
                retained_receipts=len(receipts), missing_slots=0,
                verified_receipts=sum(r["status"] == "VERIFIED" for r in receipts.values()),
                failure_statuses=dict(sorted(Counter(r["status"] for r in receipts.values()
                    if r["status"] != "VERIFIED").items())),
                distinct_public_points=25, paired_online=complete,
                paired_rows_checked=len(pairs),
                scope="Independent raw-receipt timing and pairing check; group replay, "
                      "natural-query audit and family qualification remain separate gates")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    result = audit(args.bundle)
    args.out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": result["status"],
                      "retained_receipts": result["retained_receipts"]}))


if __name__ == "__main__":
    main()
