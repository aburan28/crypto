"""One-target compact-orbit vs strong-rho panel (`launch-single`).

W one-target workloads on public hash-to-curve points with no known scalar.
Every workload runs each arm in its own process: the compact-orbit IC
(`koblitz_orbit_dlp_fast_online`) rebuilds its reusable setup and then solves
the one point inside its online interval; strong rho rung 3
(`koblitz_rho_batch_ks_strong_online`) starts with an empty
distinguished-point table. K comes from a tune on disjoint one-target
workloads. Before the eval workloads, both online producers are run against
their frozen originals on one known-answer target and must emit identical
untimed records.

Each workload yields one `vs_rho` claim keyed by IC1 candidate / workload / run
identities (`research/ic_candidate_tournament_20260915/identity.py`), with both
arms' certificates replayed in pure Python (`oracle.py`). Progress is
checkpointed after every process, so `--resume <run-id>` continues a stalled
run without repeating finished work.
"""
from __future__ import annotations

import hashlib
import json
import math
import os
import random
import statistics
import subprocess
import sys
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
for extra in (REPO / "research/ic_candidate_tournament_20260915", REPO / "tools"):
    if str(extra) not in sys.path:
        sys.path.insert(0, str(extra))
import identity  # noqa: E402
import oracle  # noqa: E402
from curve_identity import candidate_identity  # noqa: E402

PHASE_FIELDS = {
    "T_target_query_ms": "target_query_ms",
    "T_target_PDP_ms": "target_pdp_ms",
    "T_target_relation_check_ms": "target_relation_check_ms",
    "T_target_descent_ms": "target_descent_ms",
    "T_target_recovery_check_ms": "target_recovery_check_ms",
}
CURVE_FIELDS = ("degree", "curve_a", "subgroup_order", "irreducible", "group_order", "cofactor",
                "generator", "lambda")

_curves: dict[str, oracle.Curve] = {}
_OracleCurve = oracle.Curve


def cached_curve(fixture: dict[str, Any]) -> oracle.Curve:
    """oracle.Curve proves r prime by trial division (seconds at n=61): build once per curve."""
    key = json.dumps({k: fixture[k] for k in CURVE_FIELDS}, sort_keys=True)
    if key not in _curves:
        _curves[key] = _OracleCurve(fixture)
    return _curves[key]


identity.Curve = cached_curve


def half_trace(curve: oracle.Curve, beta: int) -> int:
    """H(b) = sum b^(4^i), i <= (n-1)/2; solves z^2 + z = b when Tr(b) = 0 (n odd)."""
    total, term = 0, beta
    for _ in range((curve.n - 1) // 2 + 1):
        total ^= term
        term = curve.fm(curve.fm(term, term), curve.fm(term, term))
    return total


def public_point(curve: oracle.Curve, domain: str, role: str, seed: int) -> tuple[int, int]:
    """Hash to x, solve y by half-trace, project by the cofactor; no scalar is known."""
    counter = 0
    while True:
        digest = hashlib.sha256(f"{domain}|{role}|{seed}|{counter}".encode()).digest()
        counter += 1
        x = int.from_bytes(digest[:8], "little") & ((1 << curve.n) - 1)
        if x == 0:
            continue
        inverse = curve.inv(x)
        beta = x ^ curve.a ^ curve.fm(inverse, inverse)
        z = half_trace(curve, beta)
        if curve.fm(z, z) ^ z != beta:
            continue
        if digest[8] & 1:
            z ^= 1
        point = curve.mul((x, curve.fm(x, z)), curve.h)
        if point is None:
            continue
        assert curve.decode([point[0], point[1]]) == point and curve.mul(point, curve.r) is None
        return point


def input_law(beat: dict[str, Any]) -> str:
    law = beat["target_law"]
    return (f"public hash-to-curve: x = SHA-256('{law['domain']}|<role>|<seed>|<counter>')[0:8] "
            "little-endian mod 2^n for the first counter whose x lifts; y by half-trace, sign by bit 0 "
            "of digest byte 8; Q = [h]P. No scalar is known to either arm.")


def primary_claim_fields(ic_record: dict[str, Any], rho_record: dict[str, Any],
                         interval: dict[str, Any]) -> dict[str, Any]:
    """Build the paired one-target fields from the two verified producer rows."""
    ic_ms, rho_ms = float(ic_record["online_ms"]), float(rho_record["online_ms"])
    probes, steps = int(ic_record["probes"]), int(rho_record["walk_steps"])
    if min(ic_ms, rho_ms, probes, steps) <= 0:
        raise ValueError("a verified one-target claim needs positive online costs and counters")
    return {
        "timing_class": "single_target_online",
        "record_class": "verified_answer_exploratory_wall",
        "ic_cost": ic_ms, "rho_cost": rho_ms,
        "ic_online_ms": ic_ms, "rho_online_ms": rho_ms,
        "ic_online_wall_ms": ic_ms, "rho_online_wall_ms": rho_ms,
        "online_speedup": rho_ms / ic_ms,
        "controlled_online_speedup": None,
        "automorphism_discount": {
            "A": int(rho_record["automorphism_size"]),
            "formula": "sqrt(A) for signed Frobenius classes",
        },
        "all_stages_charged_same_series": True,
        "paired_target": {
            "ic_public_q": ic_record["target"],
            "rho_public_q": rho_record["published_q"],
            "same_public_point": ic_record["target"] == rho_record["published_q"],
        },
        "ic_verified": ic_record["group_verified"] is True,
        "rho_verified": rho_record["verified"] is True,
        "ic_online_interval": f"{interval['ic_start_event']} -> {interval['ic_stop_event']}",
        "rho_online_interval": f"{interval['rho_start_event']} -> {interval['rho_stop_event']}",
        "operation_accounting": {
            "unit_assumption": "IC root-index probes and rho walk steps are distinct native counters",
            "comparison_status": "native_counters_only",
            "operation_units": {"ic": "root-index probes", "rho": "walk steps"},
            "ic_online_operations": probes,
            "rho_online_operations": steps,
            "rho_per_ic_native_counter": steps / probes,
            "ops_speedup_online": None,
        },
    }


def untimed(record: dict[str, Any], extra_skip: tuple[str, ...] = ()) -> dict[str, Any]:
    """Drop timers (`*_ms*`, `*_ns`, `*_event`), as `boundary_autolab.untimed_digest` does."""
    return {k: v for k, v in record.items()
            if "_ms" not in k and not k.endswith(("_ns", "_event")) and k not in extra_skip}


def git_bytes(commit: str, path: str) -> bytes:
    return subprocess.run(["git", "-C", str(REPO), "show", f"{commit}:{path}"],
                          capture_output=True, check=True).stdout


def committed_source_hashes(commit: str, paths: list[str], require: Any) -> dict[str, str]:
    """Hash sources at the commit; refuse when the working tree differs from it."""
    dirty = subprocess.run(["git", "-C", str(REPO), "status", "--porcelain", "--", *paths],
                           capture_output=True, text=True, check=True).stdout.strip()
    require(not dirty, f"producer sources differ from {commit}; commit them first:\n{dirty}")
    return {p: hashlib.sha256(git_bytes(commit, p)).hexdigest() for p in paths}


def bootstrap_median_ci(values: list[float], seed: int = 0, resamples: int = 10_000) -> list[float] | None:
    if len(values) < 2:
        return None
    rng = random.Random(seed)
    medians = sorted(statistics.median(rng.choices(values, k=len(values))) for _ in range(resamples))
    return [medians[int(0.025 * resamples)], medians[int(0.975 * resamples) - 1]]


def method_record(beat: dict[str, Any], k: int, r: int, hashes: dict[str, str]) -> dict[str, Any]:
    source = hashes[beat["panel_producers"]["ic"]["source"]]
    return {
        "isogeny": "none",
        "endomorphism": {"order_conductor": None, "frobenius_order_conductor": None, "volcano_levels": []},
        "factor_base": {
            "construction": {"kind": "deterministic public x-scan, cofactor projection, whole signed "
                                     "Frobenius orbit per accepted x", "orbit_columns": k},
            "nominal_bound": k,
        },
        "point_decomposition": {
            "summands": 4, "solver": "s3rootindex",
            "summation_polynomial": "third summation polynomial: S3 roots give the x of each pair sum",
            "encoding": "Frobenius-quotiented regular-root index keyed by the least normal-basis rotation",
            "equation_order": "not applicable: no polynomial system",
            "monomial_order": "not applicable: no polynomial system",
            "internal_matrix_kernel": "not applicable: no polynomial system",
            "limits": {"probe_budget": "one full scan of the index states and rotations"},
            "cache_policy": "rebuilt from nothing in every process", "source_sha256": source,
        },
        "relation_collection": {
            "collector": "guided", "query_distribution": "[a]G - R_j, a from the seeded LCG",
            "query_rule": "the first pivotless column j",
            "filtering": "none", "verification": "every relation lifted and summed in the group",
            "duplicates": "rows without rank gain are dropped", "dependencies": "incremental echelon rank",
            "stop_rule": "rank equals the column count", "source_sha256": source,
        },
        "relation_linear_algebra": {
            "solver": "dense", "modulus": r,
            "matrix_construction": "one column per signed Frobenius orbit",
            "orbit_quotient": "sign-and-Frobenius", "rank_criterion": "full column rank",
            "block_parameters": "none: incremental row echelon", "preconditioner": "none",
            "source_sha256": source,
        },
        "target_descent": {
            "method": "single", "policy": "one 4-summand decomposition of Q, scan origin from a hash of Q",
            "recursive_solvers": "none", "success_rule": "[d]G = Q in the producer",
            "stop_rule": "first lifted relation or the end of the scan", "source_sha256": source,
        },
        "implementation": {
            "source_manifest_sha256": identity.sha256(hashes),
            "components": [{"role": path, "sha256": h} for path, h in sorted(hashes.items())],
            "flags": {"rayon_threads": 1, "single_target": True},
        },
    }


def rho_manifest(curve_record: dict[str, Any], policy: dict[str, Any], hashes: dict[str, str]) -> dict[str, Any]:
    return {
        "field": curve_record["field"],
        "curve": {**curve_record["curve"], "curve_id": curve_record["curve_id"]},
        "method": "pollard-rho",
        "factor_base": "none",
        "isogeny": "none",
        "configuration": {
            "walk": policy["walk_policy"], "collision": policy["collision_policy"],
            "lanes": policy["lanes"], "distinguished_point_bits": policy["distinguished_point_bits"],
            "rung": policy["rung"], "automorphisms": policy["automorphism_size"],
        },
        "implementation": {"source_manifest_sha256": identity.sha256(hashes),
                           "components": [{"role": p, "sha256": h} for p, h in sorted(hashes.items())]},
    }


def replay_certificate(path: Path, digest: str, fixture: dict[str, Any],
                       base_points: list[tuple[int, int]] | None) -> dict[str, Any]:
    """Recompute the certificate digest and check [d]G = Q (and the IC relation) in oracle.py."""
    cert = json.loads(path.read_text())
    recomputed = identity.sha256(cert)
    curve = cached_curve(fixture)
    target = curve.decode(cert["target"])
    scalar = int(cert["scalar"])
    checks = {
        "digest_matches": recomputed == digest,
        "target_is_workload_target": cert["target"] == [int(v) for v in fixture["targets"][0]],
        "target_in_subgroup": target is not None and curve.mul(target, curve.r) is None,
        "generator_matches_fixture": cert["generator"] == [int(v) for v in fixture["generator"]],
        "scalar_in_range": 0 <= scalar < curve.r,
        "scalar_times_g_is_target": curve.mul(curve.g, scalar) == target,
    }
    if cert["arm"] == "ic":
        indices = cert["relation"]["point_indices"]
        points = [base_points[i] if base_points and 0 <= i < len(base_points) else None for i in indices]
        total = None
        for point in points:
            total = curve.add(total, point)
        checks["relation_points_on_curve"] = all(
            p is not None and curve.decode([p[0], p[1]]) == p for p in points)
        checks["relation_x_codes_match_points"] = [p[0] if p else None for p in points] == cert["relation"]["x_codes"]
        checks["relation_sums_to_target"] = total == target
    return {"arm": cert["arm"], "digest_recomputed": recomputed, "checks": checks,
            "statement_holds": all(checks.values()), "checker": "oracle.py (pure-Python field arithmetic)"}


def launch_single(arguments: Any, lab: Any) -> dict[str, Any]:
    protocol = lab.load_protocol()
    ledger = lab.load_ledger(protocol)
    beat_id = arguments.beat
    lab.require(beat_id in protocol["beats"], f"unknown beat id: {beat_id}")
    beat = protocol["beats"][beat_id]
    lab.require(beat.get("launch_mode") == "single_target_panel", f"{beat_id} is not a single-target panel beat")
    n, a = int(beat["n"]), int(beat["a"])
    producers = beat["panel_producers"]
    frozen = beat["frozen_originals"]

    if arguments.resume:
        run = lab.RUNS_DIR / arguments.resume
        lab.require((run / "state.json").is_file(), f"no run to resume: {run}")
        state = lab.read_json(run / "state.json")
        lab.require(state.get("beat_id") == beat_id, "resume beat differs from the run's beat")
        workloads, tune_workloads = int(state["workloads"]), int(state["tune_workloads"])
        candidates = list(state["k_candidates"])
        state.setdefault("resumed_at", []).append(lab.now())
    else:
        workloads = arguments.workloads or int(beat["workloads_default"])
        lab.require(workloads >= int(beat["workloads_minimum"]),
                    f"workloads must be >= {beat['workloads_minimum']}")
        tune_workloads = arguments.tune_workloads or int(beat["tune_workloads_default"])
        lab.require(tune_workloads >= 1, "tune workloads must be positive")
        if arguments.k is not None:
            candidates = [arguments.k]
        elif arguments.k_candidates:
            candidates = [int(v) for v in arguments.k_candidates.split(",") if v.strip()]
        else:
            candidates = [int(v) for v in beat["k_candidates"]]
        lab.require(bool(candidates) and all(k > 0 for k in candidates), "no K candidates")
        run_id = arguments.run_id or lab.make_run_id(beat_id)
        run = lab.RUNS_DIR / run_id
        lab.require(not run.exists(), f"run already exists: {run}")
        state = {
            "schema_version": "1.0", "task_id": lab.TASK_ID, "run_id": run_id, "beat_id": beat_id,
            "launch_mode": "single_target_panel", "status": "ACTIVE", "phase": "preflight",
            "workloads": workloads, "tune_workloads": tune_workloads, "k_candidates": candidates,
            "pinned_cpu": arguments.cpu, "created_at": lab.now(), "updated_at": lab.now(),
        }
    with lab.RunnerLock():
        for name in ("artifacts", "inputs", "logs", "receipts", "artifacts/certificates", "artifacts/claims"):
            (run / name).mkdir(parents=True, exist_ok=True)

        def advance(phase: str, **extra: Any) -> None:
            state.update(phase=phase, updated_at=lab.now(), **extra)
            lab.write_json(run / "state.json", state)

        preflight_receipt = lab.preflight(protocol, ledger)
        for arm, producer in list(producers.items()) + [(f"frozen_{k}", v) for k, v in frozen.items()]:
            source = REPO / producer["source"]
            preflight_receipt["checks"].append(
                {"name": f"panel_producer_source_{arm}", "ok": source.is_file(), "detail": str(source)})
        preflight_receipt["ok"] = all(
            c["ok"] for c in preflight_receipt["checks"] if c["name"] != "cryptominisat5_optional")
        lab.require(preflight_receipt["ok"], "preflight failed")
        receipt_name = "preflight.json" if not arguments.resume else f"preflight_resume{len(state['resumed_at'])}.json"
        lab.write_json(run / "artifacts" / receipt_name, preflight_receipt)
        git_head = next(c["detail"] for c in preflight_receipt["checks"] if c["name"] == "git_head")
        if not arguments.resume:
            lab.write_json(lab.CURRENT_PATH, {"run_id": run.name, "beat_id": beat_id, "updated_at": lab.now()})
            lab.write_json(run / "inputs/protocol.json", protocol)
            lab.write_json(run / "inputs/ledger_pin.json",
                           {"path": protocol["ledger"]["path"], "sha256": lab.sha256(lab.ledger_path(protocol)),
                            "schema_version": ledger["schema_version"]})
            state["git_head"] = git_head
        else:
            lab.require(state.get("git_head") == git_head,
                        f"resume needs the run's commit {state.get('git_head')}, not {git_head}")
        ic_hashes = committed_source_hashes(git_head, beat["ic_sources"], lab.require)
        rho_hashes = committed_source_hashes(git_head, beat["rho_sources"], lab.require)
        advance("build")

        command = ["cargo", "build", "--release"]
        for producer in list(producers.values()) + list(frozen.values()):
            command += ["--example", producer["example"]]
        completed = subprocess.run(command, cwd=REPO, capture_output=True, text=True)
        lab.require(completed.returncode == 0, "cargo build failed:\n" + completed.stderr[-4000:])
        binaries = {arm: str((REPO / "target/release/examples" / p["example"]).resolve())
                    for arm, p in list(producers.items()) + [(f"frozen_{k}", v) for k, v in frozen.items()]}
        executables = {f"{arm}_binary_sha256": lab.sha256(path) for arm, path in binaries.items()}
        executables |= {"git_head": git_head,
                        "ic_source_manifest_sha256": identity.sha256(ic_hashes),
                        "rho_source_manifest_sha256": identity.sha256(rho_hashes),
                        "rustc": subprocess.run(["rustc", "--version"], capture_output=True,
                                                text=True).stdout.strip()}
        binaries_path = run / "artifacts/binaries.json"
        if binaries_path.is_file():
            previous = lab.read_json(binaries_path)["hashes"]
            lab.require(all(previous.get(k) == v for k, v in executables.items() if k.endswith("_sha256")),
                        "resume rebuilt different binaries")
        lab.write_json(binaries_path, {"paths": binaries, "hashes": executables,
                                       "ic_sources": ic_hashes, "rho_sources": rho_hashes})

        env = os.environ.copy()
        env["RAYON_NUM_THREADS"] = "1"
        rho_env = dict(env, **{k: str(v) for k, v in producers["rho"].get("env", {}).items()})
        fixture_base = dict(beat["curve_fixture"])
        curve = cached_curve(fixture_base)
        r = curve.r
        generator = [int(v) for v in fixture_base["generator"]]

        advance("targets")
        targets_path = run / "inputs/targets.json"
        if targets_path.is_file():
            targets = lab.read_json(targets_path)
        else:
            law = beat["target_law"]
            targets = {}
            for role, count in (("tune", tune_workloads), ("eval", workloads)):
                targets[role] = [
                    {"seed": seed, "point": list(public_point(curve, law["domain"], role, seed))}
                    for seed in range(count)
                ]
            tune_points = {tuple(t["point"]) for t in targets["tune"]}
            lab.require(not tune_points & {tuple(t["point"]) for t in targets["eval"]},
                        "tune and eval targets overlap")
            targets |= {"domain": law["domain"], "input_law": input_law(beat), "disjoint": True}
            lab.write_json(targets_path, targets)
        for role in ("tune", "eval"):
            for entry in targets[role]:
                path = run / f"inputs/{role}_w{entry['seed']:03d}.target"
                if not path.is_file():
                    path.write_text(json.dumps(entry["point"]) + "\n")

        def ic_command(k: int, target_path: Path, out_path: str, binary: str = "ic") -> list[str]:
            return [binaries[binary], f"construct:{n}:{a}:{k}", str(target_path),
                    str(beat["ic_rank_seed"]), out_path]

        def rss_estimate(k: int) -> tuple[bool, int, int | None]:
            return lab.k_fits(beat, k)

        advance("k_tune")
        tune_file = run / "artifacts/k_tune.json"
        tune = lab.read_json(tune_file) if tune_file.is_file() else {"rows": []}
        done_tune = {(row["K"], row.get("workload")) for row in tune["rows"]}
        for k in candidates:
            fits, need, available = rss_estimate(k)
            if not fits:
                if (k, None) not in done_tune:
                    tune["rows"].append({"K": k, "workload": None, "state": "skipped_memory",
                                         "estimated_ic_rss_bytes": need, "available_memory_bytes": available})
                    lab.write_json(tune_file, tune)
                continue
            if len(candidates) == 1:
                if (k, None) not in done_tune:
                    tune["rows"].append({"K": k, "workload": None, "state": "fixed",
                                         "estimated_ic_rss_bytes": need})
                    lab.write_json(tune_file, tune)
                break
            for entry in targets["tune"]:
                seed = entry["seed"]
                if (k, seed) in done_tune:
                    continue
                out = run / f"logs/tune_K{k}_w{seed:03d}.jsonl"
                measured = lab.run_panel_process(
                    ic_command(k, run / f"inputs/tune_w{seed:03d}.target", str(out)), env=env,
                    stdout_path=run / f"logs/tune_K{k}_w{seed:03d}.summary.json",
                    stderr_path=run / f"logs/tune_K{k}_w{seed:03d}.stderr.txt", cpu=arguments.cpu)
                records = [x for x in lab.read_jsonl(out) if x.get("kind") == "compact_orbit_dlp_target"]
                summary_rows = lab.read_jsonl(run / f"logs/tune_K{k}_w{seed:03d}.summary.json")
                summary = summary_rows[-1] if summary_rows else {}
                ok = (measured["exit_code"] == 0 and len(records) == 1
                      and records[0].get("group_verified") is True and records[0]["target"] == entry["point"])
                tune["rows"].append({
                    "K": k, "workload": seed, "state": "ok" if ok else "failed",
                    "estimated_ic_rss_bytes": need, "available_memory_bytes": available,
                    "run": measured, "timing_ms": summary.get("timing_ms"),
                    "peak_rss_bytes": summary.get("peak_rss_bytes"),
                    "online_ms": records[0].get("online_ms") if records else None,
                    "probes": records[0].get("probes") if records else None,
                })
                lab.write_json(tune_file, tune)
        by_k: dict[int, list[dict[str, Any]]] = {}
        for row in tune["rows"]:
            if row["state"] in ("ok", "failed") and row["K"] in candidates:
                by_k.setdefault(row["K"], []).append(row)
        fixed = [row for row in tune["rows"] if row["state"] == "fixed"]
        table = []
        for k, rows in sorted(by_k.items()):
            complete = len(rows) == tune_workloads and all(x["state"] == "ok" for x in rows)
            table.append({
                "K": k, "complete": complete,
                "median_process_wall_s": statistics.median(x["run"]["wall_s"] for x in rows),
                "median_user_s": statistics.median(x["run"]["user_s"] for x in rows),
                "median_online_ms": statistics.median(x["online_ms"] for x in rows if x["online_ms"] is not None)
                if any(x["online_ms"] is not None for x in rows) else None,
                "median_instructions": statistics.median(x["run"]["instructions_retired"] for x in rows)
                if all("instructions_retired" in x["run"] for x in rows) else None,
                "max_rss_bytes": max(x["run"]["max_rss_bytes"] for x in rows),
                "estimated_ic_rss_bytes": rows[0]["estimated_ic_rss_bytes"],
            })
        if fixed:
            k_choice, k_source = int(fixed[0]["K"]), "fixed"
        else:
            eligible = [row for row in table if row["complete"]]
            lab.require(bool(eligible), "no K candidate completed every tune workload")
            k_choice = int(min(eligible, key=lambda row: row["median_process_wall_s"])["K"])
            k_source = "tuned"
        tune |= {"table": table, "selection": beat["k_selection"], "chosen_K": k_choice, "k_source": k_source,
                 "tune_targets": [t["seed"] for t in targets["tune"]]}
        lab.write_json(tune_file, tune)

        advance("producer_identity", chosen_K=k_choice)
        identity_file = run / "artifacts/producer_identity.json"
        if not identity_file.is_file():
            scalar = int(beat["identity_known_scalar"])
            scalar_path = run / "inputs/identity_scalar.txt"
            scalar_path.write_text(f"{scalar}\n")
            outputs = {}
            for arm in ("ic", "frozen_ic"):
                out = run / f"logs/identity_{arm}.jsonl"
                measured = lab.run_panel_process(
                    ic_command(k_choice, scalar_path, str(out), binary=arm), env=env,
                    stdout_path=run / f"logs/identity_{arm}.summary.json",
                    stderr_path=run / f"logs/identity_{arm}.stderr.txt", cpu=arguments.cpu)
                outputs[arm] = (measured, lab.read_jsonl(out), lab.read_jsonl(run / f"logs/identity_{arm}.summary.json"))
            for arm in ("rho", "frozen_rho"):
                out = run / f"logs/identity_{arm}.jsonl"
                measured = lab.run_panel_process(
                    [binaries[arm], str(n), str(a), "signed_frobenius", "1", str(beat["rho_seed_base"])],
                    env=dict(rho_env, KIC_RHO_EXPLICIT_SCALAR=str(scalar)), stdout_path=out,
                    stderr_path=run / f"logs/identity_{arm}.stderr.txt", cpu=arguments.cpu)
                outputs[arm] = (measured, lab.read_jsonl(out), [])
            ic_new = [untimed(x, ("relation_checks",)) for x in outputs["ic"][1]]
            ic_old = [untimed(x) for x in outputs["frozen_ic"][1]]
            summary_keys = ("base_hash", "rank", "rank_attempts", "rank_relations", "rank_failures",
                            "regular_states", "root_table_entries", "targets_solved")
            ic_summaries = [o[2][-1] if o[2] else {} for o in (outputs["ic"], outputs["frozen_ic"])]
            rho_new = [untimed(x, ("producer_version",)) for x in outputs["rho"][1]]
            rho_old = [untimed(x, ("producer_version",)) for x in outputs["frozen_rho"][1]]
            report = {
                "known_scalar": scalar, "K": k_choice, "rho_seed": beat["rho_seed_base"],
                "exit_codes": {arm: o[0]["exit_code"] for arm, o in outputs.items()},
                "ic_untimed_records_identical": bool(ic_new) and ic_new == ic_old,
                "ic_setup_summary_identical": all(
                    ic_summaries[0].get(k) == ic_summaries[1].get(k) for k in summary_keys),
                "rho_untimed_records_identical": bool(rho_new) and rho_new == rho_old,
                "compared": "every record field except timers (*_ms, *_ns, *_event), the IC's new "
                            "relation_checks counter and rho's producer_version",
            }
            report["identical"] = all(code == 0 for code in report["exit_codes"].values()) and all(
                report[k] for k in ("ic_untimed_records_identical", "ic_setup_summary_identical",
                                    "rho_untimed_records_identical"))
            lab.write_json(identity_file, report)
        identity_report = lab.read_json(identity_file)
        lab.require(identity_report["identical"], "online producers differ from their frozen originals")

        advance("workloads")
        rows_path = run / "artifacts/workloads.jsonl"
        rows = lab.read_jsonl(rows_path)
        done = {row["workload"] for row in rows}
        base_path = run / f"logs/base_n{n}_K{k_choice}.jsonl"
        for entry in targets["eval"]:
            seed = entry["seed"]
            if seed in done:
                continue
            order = ("ic", "rho") if seed % 2 == 0 else ("rho", "ic")
            row: dict[str, Any] = {"workload": seed, "order": "_then_".join(order), "target": entry["point"]}
            rho_seed = int(beat["rho_seed_base"]) + seed
            for arm in order:
                if arm == "ic":
                    ic_env = dict(env, KIC_DUMP_BASE=str(base_path)) if not base_path.is_file() else env
                    measured = lab.run_panel_process(
                        ic_command(k_choice, run / f"inputs/eval_w{seed:03d}.target",
                                   str(run / f"logs/ic_w{seed:03d}.jsonl")), env=ic_env,
                        stdout_path=run / f"logs/ic_w{seed:03d}.summary.json",
                        stderr_path=run / f"logs/ic_w{seed:03d}.stderr.txt", cpu=arguments.cpu)
                else:
                    measured = lab.run_panel_process(
                        [binaries["rho"], str(n), str(a), "signed_frobenius", "1", str(rho_seed)],
                        env=dict(rho_env, KIC_RHO_TARGET_POINT=f"{entry['point'][0]},{entry['point'][1]}"),
                        stdout_path=run / f"logs/rho_w{seed:03d}.jsonl",
                        stderr_path=run / f"logs/rho_w{seed:03d}.stderr.txt", cpu=arguments.cpu)
                row[arm] = measured
                lab.write_json(run / f"receipts/{arm}_w{seed:03d}.resource.json", measured)
            ic_records = [x for x in lab.read_jsonl(run / f"logs/ic_w{seed:03d}.jsonl")
                          if x.get("kind") == "compact_orbit_dlp_target"]
            ic_summary = (lab.read_jsonl(run / f"logs/ic_w{seed:03d}.summary.json") or [{}])[-1]
            rho_rows = lab.read_jsonl(run / f"logs/rho_w{seed:03d}.jsonl")
            rho_records = [x for x in rho_rows if x.get("kind") == "rho_ks_batch_fixture"]
            rho_summary = ([x for x in rho_rows if x.get("kind") == "rho_ks_batch_summary"] or [{}])[-1]
            row["ic_record"] = ic_records[0] if len(ic_records) == 1 else None
            row["ic_summary"] = ic_summary
            row["rho_record"] = rho_records[0] if len(rho_records) == 1 else None
            row["rho_summary"] = rho_summary
            rows.append(row)
            with rows_path.open("a") as handle:
                handle.write(json.dumps(row, sort_keys=True) + "\n")
            advance("workloads", workloads_done=len(rows))

        advance("replay")
        base_rows = lab.read_jsonl(base_path)[:1]
        lab.require(bool(base_rows), f"base dump missing: {base_path}")
        base = base_rows[0]
        base_points = [tuple(p) if p else None for p in base["factor_base_point_coordinates"]]
        representatives = [list(p) for p in base["factor_base_representatives"]]
        curve_ids = lab.read_json(REPO / beat["curve_ids"])
        curve_rec = next(c for c in curve_ids["curves"] if (c["a"], c["n"]) == (a, n))
        inventory = {"factor_base_orbits": representatives, "columns": len(representatives),
                     "column_convention": "representative"}
        candidate = identity.candidate_manifest(fixture_base | {"targets": []}, inventory,
                                                method_record(beat, k_choice, r, ic_hashes))
        identity.write_immutable(run / f"artifacts/manifests/candidates/{candidate['candidate_id']}.json", candidate)
        isolation = (f"children pinned to CPU {arguments.cpu} with taskset; wrap in tools/isolated_bench.py "
                     "reserve for section-10 evidence" if arguments.cpu is not None
                     else "none: unpinned shared host; per-run load, memory and swap recorded")
        envelope = {
            "worker_count": 1, "rayon_threads": 1, "memory_cap": "none",
            "cpus_allowed_list": str(arguments.cpu) if arguments.cpu is not None else "unpinned",
            "isolation": isolation,
            "process": "one process per arm per workload; the IC rebuilds its reusable setup in every process "
                       "and rho starts with an empty distinguished-point table",
        }
        host = lab.host_record() | {"cpu": subprocess.run(
            ["sysctl", "-n", "machdep.cpu.brand_string"], capture_output=True, text=True).stdout.strip()
            if sys.platform == "darwin" else lab.platform.processor()}
        rho_policy_base = {
            "worker_count": 1,
            "walk_policy": "rung 3: 32 lanes in lockstep, one batched inversion per step; signed-Frobenius "
                           "class representative (least normal-basis rotation, sign by y); 32 r-adding jumps "
                           "in G; stride starts [c]G + Q",
            "collision_policy": "distinguished points in one table that starts empty in every process; a stored "
                                "point reached with a different coefficient of Q gives d, checked as [d]G = Q",
            "lanes": int(producers["rho"]["env"]["KIC_RHO_LANES"]),
            "distinguished_point_bits": int(producers["rho"]["env"]["KIC_RHO_DP_BITS"]),
            "rung": int(producers["rho"]["env"]["KIC_RHO_RUNG"]),
        }
        online_interval = {
            "ic_start_event": "target_query_begin: the first operation on the public point (its scan-origin "
                              "hash), after base, index, rank stage and linear algebra are complete",
            "ic_stop_event": "recovery_check_end: [d]G = Q evaluated in the producer",
            "rho_start_event": "after_target_built: the walk's first start [c]G + Q, after the curve, jump "
                               "table and target point are ready",
            "rho_stop_event": "recovery_check_true: a distinguished-point collision returned d with [d]G = Q",
            "ic_included_stages": ["target_query", "target_PDP", "target_relation_check", "target_descent",
                                   "target_recovery_check"],
            "rho_included_stages": ["walk", "collision", "recovery_check"],
            "ic_phase_map": {"target_query": "target_query_ms", "target_PDP": "target_pdp_ms",
                             "target_relation_check": "target_relation_check_ms",
                             "target_descent": "target_descent_ms",
                             "target_recovery_check": "target_recovery_check_ms"},
            "rho_phase_map": {"walk_and_collision_ms": ["walk", "collision"],
                              "recovery_check_ms": ["recovery_check"]},
            "replay": "both certificates replayed after the run in oracle.py, outside both producers",
        }
        non_claims = [
            "one public target per row; no multi-target or batch result",
            "the IC's reusable setup (base, index, rank stage, linear algebra) is outside its online interval, "
            "is paid again in every process, and is reported beside the claim with the cold ratio",
            "rho has no precomputed distinguished-point table; Bernstein-Lange is a different reference and "
            "is not measured",
            "no key recovery; public synthetic hash-to-curve points only; no asymptotic claim",
            "out-of-process replay on the same host; independent-host validation is outstanding",
            f"wall time on one host ({host['cpu'] or host['machine']}); no claim for other hardware classes",
        ]
        verify_rows = []
        claim_rows = []
        for row in sorted(rows, key=lambda x: x["workload"]):
            seed = row["workload"]
            ic_record, rho_record = row["ic_record"], row["rho_record"]
            fixture = fixture_base | {"targets": [[str(v) for v in row["target"]]], "target_seeds": [seed],
                                      "target_scalar_constructed": False}
            producers_ok = row["ic"]["exit_code"] == 0 and row["rho"]["exit_code"] == 0
            ic_ok = bool(ic_record) and ic_record.get("group_verified") is True and ic_record["target"] == row["target"]
            rho_ok = bool(rho_record) and rho_record.get("verified") is True and rho_record["published_q"] == row["target"]
            verify: dict[str, Any] = {"workload": seed, "producers_ok": producers_ok, "ic_ok": ic_ok, "rho_ok": rho_ok,
                                      "base_hash": (row["ic_summary"] or {}).get("base_hash")}
            if not (producers_ok and ic_ok and rho_ok):
                verify["status"] = "FAILED"
                verify_rows.append(verify)
                continue
            target_hash = identity.sha256({"curve_id": curve_rec["curve_id"], "target": row["target"]})
            certificates = {
                "ic": {"arm": "ic", "curve_id": curve_rec["curve_id"], "generator": ic_record["generator"],
                       "target": ic_record["target"], "scalar": ic_record["recovered_scalar"],
                       "relation": {"base_hash": row["ic_summary"]["base_hash"],
                                    "point_indices": ic_record["point_indices"], "x_codes": ic_record["x_codes"]}},
                "rho": {"arm": "rho", "curve_id": curve_rec["curve_id"], "generator": generator,
                        "target": rho_record["published_q"], "scalar": rho_record["recovered_fixture_scalar"],
                        "walk_steps": rho_record["walk_steps"]},
            }
            digests, replays = {}, {}
            for arm, cert in certificates.items():
                path = run / f"artifacts/certificates/w{seed:03d}.{arm}.json"
                path.write_bytes(identity.canonical(cert) + b"\n")
                digests[arm] = identity.sha256(cert)
                replays[arm] = replay_certificate(path, digests[arm], fixture, base_points if arm == "ic" else None)
            replay_path = run / f"artifacts/claims/w{seed:03d}.replay.json"
            lab.write_json(replay_path, replays)
            independent = all(rep["statement_holds"] for rep in replays.values())
            verify |= {"replay_ic": replays["ic"]["statement_holds"], "replay_rho": replays["rho"]["statement_holds"],
                       "scalars_agree": certificates["ic"]["scalar"] == certificates["rho"]["scalar"],
                       "base_hash_matches_dump": row["ic_summary"]["base_hash"] == base["base_hash"],
                       "setup_before_online": row["ic_summary"]["setup_complete_ns"] < ic_record["online_start_ns"]}
            verify["status"] = "PASS" if all(verify[k] for k in (
                "replay_ic", "replay_rho", "scalars_agree", "base_hash_matches_dump", "setup_before_online")) else "FAILED"
            verify_rows.append(verify)
            work = identity.workload_manifest(fixture, input_law=targets["input_law"],
                                              algorithm_seed=int(beat["ic_rank_seed"]),
                                              resource_envelope=envelope, cache_policy="cold")
            identity.write_immutable(run / f"artifacts/manifests/workloads/{work['workload_id']}.json", work)
            run_key = identity.run_id(candidate["candidate_id"], work["workload_id"], 1)
            policy = rho_policy_base | {
                "distinguished_point_memory_bytes": int(row["rho_summary"]["table_payload_lower_bound_bytes"]),
                "seed": int(beat["rho_seed_base"]) + seed,
                "automorphism_size": int(rho_record["automorphism_size"]),
            }
            rho_ref = rho_manifest(curve_rec, policy, rho_hashes)
            rho_uid = candidate_identity(rho_ref)
            identity.write_immutable(run / f"artifacts/manifests/rho/{rho_uid['candidate_sha256']}.json",
                                     {**rho_uid, "manifest": rho_ref})
            ic_ms, rho_ms = float(ic_record["online_ms"]), float(rho_record["online_ms"])
            speedup = rho_ms / ic_ms
            ic_cold_ms = ic_record["online_stop_ns"] / 1e6
            rho_cold_ms = rho_record["online_stop_ns"] / 1e6
            claim = {
                **primary_claim_fields(ic_record, rho_record, online_interval),
                "schema_version": 2, "task_id": lab.TASK_ID, "beat_id": beat_id, "autolab_run_id": run.name,
                "stage": "vs_rho", "status": "PENDING_INDEPENDENT_VALIDATION",
                "regime": beat["regime"], "result_class": beat["result_class"],
                "n_or_bits": {"n": n, "a": a, "subgroup_order_bits": round(math.log2(r), 2)},
                "curve_id": curve_rec["curve_id"], "K": k_choice, "workload_index": seed,
                "candidate_id": candidate["candidate_id"], "candidate_manifest_sha256": candidate["record_sha256"],
                "workload_id": work["workload_id"], "workload_manifest_sha256": work["record_sha256"],
                "run_id": run_key, "rho_reference_uid": rho_uid["candidate_uid"],
                "target_count": 1, "ic_target_hash": target_hash, "rho_target_hash": target_hash,
                "ic_online_phase_ms": {schema: float(ic_record[field]) for schema, field in PHASE_FIELDS.items()},
                "rho_online_phase_ms": {"walk_and_collision_ms": rho_record["walk_and_collision_ms"],
                                        "recovery_check_ms": rho_record["recovery_check_ms"]},
                "online_interval": online_interval,
                "same_resource_envelope": True,
                "independent_validation": independent,
                "independent_validation_scope": "certificates replayed out of process in oracle.py on the same "
                                                "host; independent-host replay outstanding",
                "ic_replay_certificate_sha256": digests["ic"], "rho_replay_certificate_sha256": digests["rho"],
                "ic_resource_envelope": envelope, "rho_resource_envelope": dict(envelope),
                "ic_scalar_verified": ic_ok and replays["ic"]["statement_holds"],
                "rho_scalar_verified": rho_ok and replays["rho"]["statement_holds"],
                "rho_policy": policy,
                "verdict": "index calculus online faster" if speedup > 1 else "rho online faster",
                "claim_boundary": (f"one public point on {curve_rec['curve_id']} (a={a}, n={n}); the IC's reusable "
                                   "setup is outside its online interval and is reported beside it"),
                "independent_replay_pointer": str(replay_path.relative_to(REPO)),
                "fixture_hash": identity.sha256(fixture),
                "executable_or_source_hash": executables,
                "host_id": host,
                "resource_caps": {"rayon_threads": 1, "cpus": envelope["cpus_allowed_list"], "memory": "none"},
                "seeds": {"target_domain": targets["domain"], "target_seed": seed,
                          "ic_rank_seed": int(beat["ic_rank_seed"]), "rho_seed": policy["seed"]},
                "claim_boundary_non_claims": non_claims,
                "cold": {
                    "ic_setup_ms": row["ic_summary"]["setup_complete_ns"] / 1e6,
                    "ic_in_process_to_solved_ms": ic_cold_ms,
                    "rho_in_process_to_solved_ms": rho_cold_ms,
                    "cold_ratio_ic_over_rho": ic_cold_ms / rho_cold_ms,
                    "ic_process_wall_s": row["ic"]["wall_s"], "rho_process_wall_s": row["rho"]["wall_s"],
                    "ic_timing_ms": row["ic_summary"].get("timing_ms"),
                },
                "counts": {"ic_probes": ic_record["probes"], "ic_relation_checks": ic_record["relation_checks"],
                           "rho_walk_steps": rho_record["walk_steps"]},
                "memory": {"ic_max_rss_bytes": row["ic"]["max_rss_bytes"],
                           "rho_max_rss_bytes": row["rho"]["max_rss_bytes"]},
            }
            validation = lab.validate_claim(claim, stage="vs_rho", ledger=ledger)
            lab.write_json(run / f"artifacts/claims/w{seed:03d}.claim.json", claim)
            lab.write_json(run / f"artifacts/claims/w{seed:03d}.claim_check.json", validation)
            claim_rows.append({"workload": seed, "run_id": run_key, "status": validation["status"],
                               "errors": validation["validation_errors"],
                               "missing": validation["missing_stage_fields"] + validation["missing_global_provenance"],
                               "online_speedup": speedup, "ic_online_ms": ic_ms, "rho_online_ms": rho_ms,
                               "cold_ratio_ic_over_rho": claim["cold"]["cold_ratio_ic_over_rho"]})

        advance("analysis")
        speedups = [c["online_speedup"] for c in claim_rows]
        colds = [c["cold_ratio_ic_over_rho"] for c in claim_rows]
        claimed = {c["workload"] for c in claim_rows}

        def column(path: tuple[str, ...]) -> list[float]:
            values = []
            for row in rows:
                if row["workload"] not in claimed:
                    continue
                value: Any = row
                for key in path:
                    value = value.get(key) if isinstance(value, dict) else None
                if isinstance(value, (int, float)):
                    values.append(float(value))
            return values

        def describe(values: list[float]) -> dict[str, Any] | None:
            if not values:
                return None
            return {"median": statistics.median(values), "min": min(values), "max": max(values),
                    "median_95ci_bootstrap": bootstrap_median_ci(values)}

        summary = {
            "schema_version": "1.0", "beat_id": beat_id, "run_id": run.name, "n": n, "a": a,
            "subgroup_order_bits": round(math.log2(r), 2), "K": k_choice, "k_source": k_source,
            "workloads": workloads, "workloads_measured": len(rows), "workloads_claimed": len(claim_rows),
            "candidate_id": candidate["candidate_id"], "comparator": beat["comparator_status"],
            "producer_identity": identity_report, "k_tune_table": table,
            "online_speedup_rho_over_ic": describe(speedups),
            "cold_ratio_ic_over_rho": describe(colds),
            "ic_online_ms": describe(column(("ic_record", "online_ms"))),
            "rho_online_ms": describe(column(("rho_record", "online_ms"))),
            "ic_setup_ms": describe([v / 1e6 for v in column(("ic_summary", "setup_complete_ns"))]),
            "ic_process_wall_s": describe(column(("ic", "wall_s"))),
            "rho_process_wall_s": describe(column(("rho", "wall_s"))),
            "ic_instructions": describe(column(("ic", "instructions_retired"))),
            "rho_instructions": describe(column(("rho", "instructions_retired"))),
            "ic_max_rss_bytes": describe(column(("ic", "max_rss_bytes"))),
            "rho_max_rss_bytes": describe(column(("rho", "max_rss_bytes"))),
            "ic_online_phase_median_ms": {
                schema: statistics.median(column(("ic_record", field))) for schema, field in PHASE_FIELDS.items()
            } if claim_rows else None,
            "rho_walk_steps": describe(column(("rho_record", "walk_steps"))),
            "ic_probes": describe(column(("ic_record", "probes"))),
            "verification": {"rows": verify_rows,
                             "all_pass": len(verify_rows) == workloads and all(v["status"] == "PASS" for v in verify_rows)},
            "claim_checks": {"rows": claim_rows,
                             "pass": sum(1 for c in claim_rows if c["status"] == "PASS"),
                             "fail": sum(1 for c in claim_rows if c["status"] != "PASS")},
            "envelope": envelope, "host": host,
        }
        lab.write_json(run / "artifacts/single_target_summary.json", summary)
        first = run / f"artifacts/claims/w{targets['eval'][0]['seed']:03d}.claim.json"
        if first.is_file():
            lab.write_json(run / "artifacts/claim_draft.json", lab.read_json(first))
            lab.write_json(run / "artifacts/claim_check.json",
                           lab.validate_claim(lab.read_json(first), stage="vs_rho", ledger=ledger))
        lab.write_json(run / "artifacts/claim_check_all.json", summary["claim_checks"])

        exit_codes = [row[arm]["exit_code"] for row in rows for arm in ("ic", "rho")]
        if any(code != 0 for code in exit_codes):
            status_value = "PRODUCER_FAILURE"
        elif not summary["verification"]["all_pass"]:
            status_value = "VERIFICATION_FAILURE"
        elif summary["claim_checks"]["fail"] or len(claim_rows) != workloads:
            status_value = "SCHEMA_INCOMPLETE"
        else:
            status_value = "PENDING_INDEPENDENT_VALIDATION"
        state.update(status=status_value, phase="done", updated_at=lab.now(), chosen_K=k_choice,
                     claim_check="PASS" if status_value == "PENDING_INDEPENDENT_VALIDATION" else "FAIL",
                     claim_check_counts={"pass": summary["claim_checks"]["pass"],
                                         "fail": summary["claim_checks"]["fail"]},
                     online_speedup=summary["online_speedup_rho_over_ic"],
                     cold_ratio=summary["cold_ratio_ic_over_rho"])
        lab.write_json(run / "state.json", state)
        lab.write_json(run / "artifacts/candidate.json", {
            "schema_version": "1.0", "task_id": lab.TASK_ID, "run_id": run.name, "beat_id": beat_id,
            "status": status_value, "claim_checks": summary["claim_checks"],
            "ledger_sha256": lab.sha256(lab.ledger_path(protocol)), "created_at": lab.now(),
            "note": "One-target online vs_rho rows; PENDING_INDEPENDENT_VALIDATION, never a ledger promotion.",
        })
        files = {str(p.relative_to(run)): lab.sha256(p) for p in sorted(run.rglob("*")) if p.is_file()
                 and p.name != "review_manifest.json"}
        lab.write_json(run / "artifacts/review_manifest.json",
                       {"schema_version": "1.0", "task_id": lab.TASK_ID, "files": files})
        return state
