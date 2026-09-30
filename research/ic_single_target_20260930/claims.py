#!/usr/bin/env python3
"""Ledger §23's claim rows: one per target, checked (PROTOCOL.md, "Rows").

For every row the run produced (the first clean, complete attempt of each
`T<i>-R<run>`), this

1. builds the canonical identities with the tournament's adapter
   (`research/ic_candidate_tournament_20260915/identity.py`): the IC1
   candidate from the recipe and the factor base's actual orbits, the
   workload from the curve, the one public target, the seeds and the
   resource envelope, and the run ID `<candidate>W<workload>R<run>`;
2. replays both arms' certificates outside the Rust process, in the
   independent Python checker (`oracle.py`): the digest is recomputed from
   the certificate's canonical JSON, and `[scalar]G = target` is checked
   in its own field arithmetic;
3. gives rho its own reference identity (`tools/curve_identity.py`,
   never an IC1 label), bound to the curve's EC1 record;
4. assembles the `vs_rho` claim and runs the repository's checker
   (`research/sat_factor_base_review_20260908/autolab`, `validate_claim`).

Identity sources are hashed at the commit the binary was built from
(`host.json`'s `ic_built_from`), with `git show`, so the working tree
does not enter them.

    IC_RUNS=runs python3 claims.py        # -> claims/ and manifests/ (under IC_OUT, default here)
"""
from __future__ import annotations

import hashlib
import json
import os
import re
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / "research" / "ic_candidate_tournament_20260915"))
sys.path.insert(0, str(ROOT / "research" / "sat_factor_base_review_20260908" / "autolab"))
sys.path.insert(0, str(ROOT / "tools"))
import boundary_autolab as lab  # noqa: E402
import identity  # noqa: E402
import make_params  # noqa: E402
import oracle  # noqa: E402
from curve_identity import candidate_identity  # noqa: E402

RUNS = (HERE / os.environ.get("IC_RUNS", "runs")).resolve()
BASE = (HERE / os.environ.get("IC_OUT", ".")).resolve()
OUT = BASE / "claims"
MANIFESTS = BASE / "manifests"

# The sources the two arms execute, hashed at the binary's commit.
IC_SOURCES = [
    "src/cryptanalysis/koblitz_index_calculus.rs",
    "src/cryptanalysis/koblitz_fast.rs",
    "src/cryptanalysis/koblitz_sparse_la.rs",
    "src/cryptanalysis/semaev_decomp.rs",
    "src/cryptanalysis/ic_measurement.rs",
    "src/bin/ic/price.rs",
    "src/bin/ic/workflow.rs",
    "src/bin/ic/experiment.rs",
]
RHO_SOURCES = [
    "src/cryptanalysis/ic_boundary.rs",
    "src/cryptanalysis/koblitz_fast.rs",
    "src/cryptanalysis/semaev_decomp.rs",
    "src/cryptanalysis/ic_measurement.rs",
    "src/bin/ic/price.rs",
]
INPUT_LAW = "public hash-to-curve, ic workflow's domain ic-workflow-public-target-v1, seed 23000 + i"
NON_CLAIMS = [
    "one public target per row; no multi-target or batch result",
    "the online interval excludes the reusable set-up, which is reported separately and dominates the cold cost",
    "rho here has no precomputed table of distinguished points; a generic walk with one (Bernstein and Lange) "
    "is a different reference and is not measured",
    "no statement about ECC2K-130 or any curve past n = 61; m = 83 is not run (AGENTS.md section 8a)",
    "wall time on one x86-64 container; no claim for other hardware classes",
]


def rel(path: Path) -> str:
    """A path relative to the repository where it lies in it."""
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def git_bytes(commit: str, path: str) -> bytes:
    return subprocess.run(["git", "-C", str(ROOT), "show", f"{commit}:{path}"], capture_output=True,
                          check=True).stdout


def source_hashes(commit: str, paths: list[str]) -> dict[str, str]:
    return {p: hashlib.sha256(git_bytes(commit, p)).hexdigest() for p in paths}


def rows(a: int, n: int):
    """The row of each (target, run): the first clean, complete attempt."""
    d = RUNS / f"k{a}n{n}"
    for base in sorted(d.glob("T*-R*.price.json")):
        if "-retry" in base.name:
            continue
        m = re.fullmatch(r"T(\d+)-R(\d+)\.price\.json", base.name)
        if not m:
            continue
        i, run = int(m[1]), int(m[2])
        attempts = [base] + sorted(d.glob(base.name.replace(".price.json", "-retry*.price.json")))
        chosen, tried = None, []
        for att in attempts:
            rec = att.with_name(att.name[: -len(".price.json")] + ".isolation.jsonl")
            iso = json.loads(rec.read_text().splitlines()[-1])["run"] if rec.exists() else None
            rep = json.loads(att.read_text()) if att.exists() else {"status": "no report"}
            ok = iso is not None and iso["exit_status"] == 0 and not iso["contended"] \
                and rep.get("status") == "complete"
            tried.append({"file": att.name, "status": rep.get("status"),
                          "contended": None if iso is None else iso["contended"],
                          "exit_status": None if iso is None else iso["exit_status"]})
            if ok and chosen is None:
                chosen = (att, rep, iso)
        yield i, run, chosen, tried


def method_record(params: dict, fixture: dict, hashes: dict[str, str]) -> dict:
    lib = hashes["src/cryptanalysis/koblitz_index_calculus.rs"]
    spec = params["factor_base"]["spec"]
    m = params.get("descent_summands", params["summands"])
    return {
        "isogeny": "none",
        "endomorphism": {"order_conductor": None, "frobenius_order_conductor": None, "volcano_levels": []},
        "factor_base": {"construction": spec, "nominal_bound": spec["points"]},
        "point_decomposition": {
            "summands": params["summands"], "solver": "pairtable",
            "summation_polynomial": "none: the pair-sum table is looked up",
            "encoding": "pair sums folded by sign and Frobenius, keyed by the normal-basis least rotation of the abscissa",
            "equation_order": "not applicable: no polynomial system",
            "monomial_order": "not applicable: no polynomial system",
            "internal_matrix_kernel": "not applicable: no polynomial system",
            "limits": {"max_trials": params["max_trials"], "table_tier": "auto, the workflow's probe budget"},
            "cache_policy": "rebuilt from nothing in every process", "source_sha256": lib,
        },
        "relation_collection": {
            "collector": "aimed", "query_distribution": "[a]G + [b]Q-free probes from seeded work units",
            "query_rule": "aimed at the least-mentioned projected columns",
            "filtering": f"collection window {params['collection_window']}",
            "verification": "every relation checked in the group before it is pushed",
            "duplicates": "dropped by the log solver", "dependencies": "rank of the sparse system",
            "stop_rule": "extend by units until every column is determined",
            "source_sha256": lib,
        },
        "relation_linear_algebra": {
            "solver": "sparse", "modulus": int(fixture["subgroup_order"]),
            "matrix_construction": "one column per projected signed Frobenius orbit",
            "orbit_quotient": "sign-and-Frobenius", "rank_criterion": "every column determined and verified",
            "block_parameters": "library defaults", "preconditioner": "library defaults", "source_sha256": lib,
        },
        "target_descent": {
            "method": "walk",
            "policy": f"{m}-summand pair-table walk of 64 lanes stepped by G, one start by two scalar multiplications",
            "recursive_solvers": "none", "success_rule": "[d]G = Q on the single-word ladder",
            "stop_rule": f"max_trials {params['max_trials']}", "source_sha256": lib,
        },
        "implementation": {
            "source_manifest_sha256": identity.sha256(hashes),
            "components": [{"role": path, "sha256": h} for path, h in sorted(hashes.items())],
            "flags": {"rayon_threads": 1, "single_target": True},
        },
    }


def rho_manifest(curve_ids: dict, a: int, n: int, rep: dict, hashes: dict[str, str]) -> dict:
    rec = next(c for c in curve_ids["curves"] if (c["a"], c["n"]) == (a, n))
    policy = rep["rho_policy"]
    return {
        "field": rec["field"],
        "curve": {**rec["curve"], "curve_id": rec["curve_id"]},
        "method": "pollard-rho",
        "factor_base": "none",
        "isogeny": "none",
        "configuration": {
            "walk": policy["walk_policy"], "collision": policy["collision_policy"],
            "lanes": policy["lanes"], "distinguished_point_bits": policy["distinguished_point_bits"],
            "automorphisms": 2 * n, "jumps": policy["counters"]["jumps"],
        },
        "implementation": {"source_manifest_sha256": identity.sha256(hashes),
                           "components": [{"role": p, "sha256": h} for p, h in sorted(hashes.items())]},
    }


def replay(fixture: dict, rep: dict, arm: str) -> dict:
    cert = rep["certificates"][arm]
    text = json.dumps(cert, sort_keys=True, separators=(",", ":"), ensure_ascii=False)
    digest = hashlib.sha256(text.encode()).hexdigest()
    curve = oracle.Curve(fixture)
    target = curve.decode(cert["target"])
    scalar = int(cert["scalar"])
    holds = 0 <= scalar < curve.r and curve.mul(curve.g, scalar) == target \
        and cert["target"] == fixture["targets"][0]
    return {"arm": arm, "digest_recomputed": digest, "digest_matches": digest == rep["certificates"][arm + "_sha256"],
            "statement_holds": holds, "scalar": cert["scalar"], "checker": "oracle.py (Python field arithmetic)"}


def main() -> None:
    host = json.loads((RUNS / "host.json").read_text())
    commit = host["ic_built_from"] or host["commit"]
    ic_hashes, rho_hashes = source_hashes(commit, IC_SOURCES), source_hashes(commit, RHO_SOURCES)
    curve_ids = json.loads((HERE / "curve_ids.json").read_text())
    ledger = lab.load_ledger(lab.load_protocol())
    summary = []
    for a, n in make_params.SIZES:
        for i, run, chosen, tried in rows(a, n):
            label = f"k{a}n{n}/T{i:02d}-R{run}"
            if chosen is None:
                summary.append({"row": label, "status": "no clean complete attempt", "attempts": tried})
                continue
            path, rep, iso = chosen
            params = json.loads((RUNS / f"k{a}n{n}" / f"T{i:02d}.params.json").read_text())
            ii = rep["identity_inputs"]
            fixture = ii["fixture"]
            cand = identity.candidate_manifest(
                fixture,
                {"factor_base_orbits": ii["factor_base_orbits"], "columns": ii["columns"],
                 "column_convention": ii["column_convention"]},
                method_record(params, fixture, ic_hashes))
            envelope = {
                "worker_count": 1, "rayon_threads": 1,
                "isolation": "tools/isolated_bench.py run --wait --cpus 2, uncontended",
                "cpus_allowed_list": rep["resource_envelope"]["cpus_allowed_list"],
                "process": "one process: the index calculus, then rho, both rebuilt each repetition",
            }
            work = identity.workload_manifest(fixture, input_law=INPUT_LAW, algorithm_seed=params["seed"],
                                              resource_envelope=envelope, cache_policy="cold")
            run_id = identity.run_id(cand["candidate_id"], work["workload_id"], run)
            rho_ref = rho_manifest(curve_ids, a, n, rep, rho_hashes)
            rho_uid = candidate_identity(rho_ref)
            identity.write_immutable(MANIFESTS / "candidates" / f"{cand['candidate_id']}.json", cand)
            identity.write_immutable(MANIFESTS / "workloads" / f"{work['workload_id']}.json", work)
            identity.write_immutable(MANIFESTS / "rho" / f"{rho_uid['candidate_sha256']}.json",
                                     {**rho_uid, "manifest": rho_ref})
            replays = {arm: replay(fixture, rep, arm) for arm in ("ic", "rho")}
            replay_path = OUT / f"k{a}n{n}" / f"T{i:02d}-R{run}.replay.json"
            replay_path.parent.mkdir(parents=True, exist_ok=True)
            replay_path.write_text(json.dumps(replays, indent=1) + "\n")
            m = rep["median"]
            reps = rep["repetitions"]
            ic_ok = all(r["ic_online"]["verified"] for r in reps)
            rho_ok = all(r["rho_online"]["verified"] for r in reps)
            speedup = m["rho_online_wall_ms"] / m["ic_online_wall_ms"]
            claim = {
                "n": n, "a": a, "log2_r": rep["log2_r"], "curve": rep["curve"],
                "curve_id": next(c["curve_id"] for c in curve_ids["curves"] if (c["a"], c["n"]) == (a, n)),
                "candidate_id": cand["candidate_id"], "candidate_manifest_sha256": cand["record_sha256"],
                "workload_id": work["workload_id"], "workload_manifest_sha256": work["record_sha256"],
                "run_id": run_id, "rho_reference_uid": rho_uid["candidate_uid"],
                "target_count": 1,
                "ic_target_hash": rep["target"]["sha256"], "rho_target_hash": rep["target"]["sha256"],
                "timing_class": "single_target_online_wall",
                "ic_online_wall_ms": m["ic_online_wall_ms"], "rho_online_wall_ms": m["rho_online_wall_ms"],
                "online_speedup": speedup,
                "ic_online_phase_ms": m["ic_online_phase_ms"],
                "ic_online_phases_not_entered": m["ic_online_phases_not_entered"],
                "online_interval": rep["online_interval"],
                "same_resource_envelope": True,
                "independent_validation": all(r["digest_matches"] and r["statement_holds"] for r in replays.values()),
                "ic_replay_certificate_sha256": rep["certificates"]["ic_sha256"],
                "rho_replay_certificate_sha256": rep["certificates"]["rho_sha256"],
                "ic_resource_envelope": envelope, "rho_resource_envelope": dict(envelope),
                "ic_scalar_verified": ic_ok, "rho_scalar_verified": rho_ok,
                "rho_policy": {k: rep["rho_policy"][k] for k in
                               ("worker_count", "walk_policy", "collision_policy",
                                "distinguished_point_memory_bytes", "lanes", "distinguished_point_bits")},
                "verdict": "index calculus online faster" if speedup > 1 else "rho online faster",
                "claim_boundary": (f"one public target on {rep['curve']}; the index calculus's reusable set-up "
                                   "is outside its online interval and is reported beside it"),
                "independent_replay_pointer": rel(replay_path),
                "fixture_hash": identity.sha256(fixture),
                "executable_or_source_hash": {"ic_binary_sha256": host["ic_binary_sha256"], "source_commit": commit,
                                              "ic_source_manifest_sha256": identity.sha256(ic_hashes),
                                              "rho_source_manifest_sha256": identity.sha256(rho_hashes)},
                "host_id": {"cpu_model": host["cpu_model"], "logical_cores": host["logical_cores"],
                            "os": host["os"], "host_manifest": "runs/host.json"},
                "resource_caps": {"rayon_threads": 1, "cpus": "2", "contended": iso["contended"]},
                "seeds": {"recipe_seed": params["seed"], "public_hash_seed": make_params.target_seed(i),
                          "rho_seed": rep["rho_policy"]["seed"]},
                "claim_boundary_non_claims": NON_CLAIMS,
                # Beside the claim: the same row in the thread's unit.
                "units": {k: m[k] for k in (
                    "ic_online_units", "rho_online_units", "setup_units", "rho_setup_units", "rho_model_units",
                    "s_ic_online", "s_rho_online", "s_setup", "s_rho_setup", "s_ic_cold", "s_rho_cold",
                    "s_rho_model", "cold_ratio_ic_over_rho", "online_speedup_rho_model", "rho_units_per_step",
                    "rho_step_over_model", "ic_replay_ns", "rho_replay_ns")},
                "report": rel(path),
                "attempts": tried,
            }
            result = lab.validate_claim(claim, stage="vs_rho", ledger=ledger)
            claim_path = OUT / f"k{a}n{n}" / f"T{i:02d}-R{run}.claim.json"
            claim_path.write_text(json.dumps({"claim": claim, "validation": result}, indent=1) + "\n")
            summary.append({"row": label, "run_id": run_id, "status": result["status"],
                            "errors": result["validation_errors"],
                            "missing": result["missing_stage_fields"] + result["missing_global_provenance"],
                            "independent_validation": claim["independent_validation"]})
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    passed = sum(1 for s in summary if s["status"] == "PASS")
    print(f"{passed} of {len(summary)} rows pass the vs_rho check")


if __name__ == "__main__":
    main()
