"""Adapter for the cryptanalysis repository's IC end-to-end benchmark.

cryptanalysis ``experiments/ic-bench`` runs complete, verified index-calculus
DLPs on the toy Koblitz curve K_0: y^2 + xy = x^3 + 1 over F_2^n.  It meters
eleven exclusive phases with exact integer counters (opcount.Meter) and
prices them with a frozen per-degree weight table in "rps" (reference
picoseconds), so its totals are exact on any host for one calibration_id.
This adapter runs ONE cell, in the pinned process ICMS forks, through
``cryptanalysis_cell.py`` (bench.py itself only runs whole suites through an
unpinned process pool, and ``--record`` rewrites its baseline).

What it has to be honest about, each recorded in the output:

* The workflow collects relations to the achievable rank and then descends
  every target (``relations.stop: full_rank``, ``descent.method: pdp``); a
  spec asking for crypto's first-log stop is refused.  Relations needed at
  full rank are about twice those at a first-log stop for m = 2 (audit X05).
* Its rho reference is the analytic plain walk sqrt(pi r / 2) (``rho.plain``),
  even on Koblitz curves; the signed-Frobenius figure sqrt(pi r / (4n)) is
  reported beside it.  The spec must name which one it compares against, and
  ``rho.measured_matched`` (crypto's counted walk) is refused.
* The unit is ``cryptanalysis.rps``.  It is not crypto's GAE and the two are
  never compared; ``count.group_additions`` (the ec_add counter alone) is
  offered as the one count both repositories share.
* The receipt has no Boolean-system shape (no n_vars, degrees or matrix
  sizes): those PDP fields are recorded as null, not zero.
"""
from __future__ import annotations

import json
import math
import os
import subprocess
from typing import Any

from ..canonical import sha256_file
from ..registry import Registry
from . import Context

UNITS = ("cryptanalysis.rps", "count.group_additions", "time.wall_ns")
LAWS = ("prefix", "geometric", "geomtrace", "geomtraceu", "random", "normal", "kertrace", "invariant")
# The solver limits bench.py hard-codes (bench.py LIMITS at cryptanalysis 46ad7a8).  The runner refuses to
# run if the checkout's constants differ from what the spec pinned, so a changed bench is a changed spec.
SOLVER = "macaulay-xl"
RUNNER = os.path.join(os.path.dirname(os.path.abspath(__file__)), "cryptanalysis_cell.py")


def _root(spec: dict[str, Any], ctx: Context) -> str:
    root = (spec["execution"].get("args") or {}).get("root", "../cryptanalysis")
    return root if os.path.isabs(root) else os.path.normpath(os.path.join(ctx["repo_root"], root))


def _git(root: str, *args: str) -> str | None:
    try:
        out = subprocess.run(["git", "-C", root, *args], capture_output=True, text=True, timeout=30)
    except (OSError, subprocess.SubprocessError):
        return None
    return out.stdout.strip() if out.returncode == 0 else None


class CryptanalysisIcBench:
    name = "cryptanalysis.ic_bench"

    def __init__(self, registry: Registry | None = None):
        self.reg = registry or Registry.load()

    # ------------------------------------------------------------------ check
    def check(self, spec: dict[str, Any]) -> list[str]:
        p: list[str] = []
        curve, wl = spec["instance"]["curve"], spec["instance"]["workload"]
        if curve["regime"] != "koblitz" or curve.get("koblitz_a") != 0 or "degree" not in curve:
            p.append("cryptanalysis ic-bench runs K_0 only: instance.curve must be {regime: koblitz, koblitz_a: 0, degree: n}")
        if wl["law"] != "known_answer":
            p.append("ic-bench plants its targets from the workload seed: workload.law must be known_answer")
        if "ref" in curve or "explicit" in curve:
            p.append("ic-bench builds K_0 from the degree; it does not resolve instance.curve.ref or explicit")
        fb = spec["factor_base"]
        params = fb.get("params") or {}
        if "basis_seed" not in params:
            p.append("the basis sampler's seed defaults to 1 in bench.py: set factor_base.params.basis_seed")
        if fb["family"] != "binary_subspace":
            p.append("ic-bench bases are F_2-subspaces: factor_base.family must be binary_subspace")
        if params.get("basis_law") not in LAWS:
            p.append(f"factor_base.params.basis_law must be one of {list(LAWS)}")
        if sorted(fb.get("quotient") or []) != ["negation"]:
            p.append("ic-bench folds P and -P into one column (quotient_rule sign) and folds Frobenius orbits only for "
                     "the invariant law; declare quotient [negation]")
        if fb.get("large_primes", {}).get("mode", "none") != "none" or fb.get("partition"):
            p.append("ic-bench has no large-prime or partitioned base")
        dec = spec["decomposition"]
        if dec["method"] != "algebraic" or dec.get("encoding") != "expanded_semaev":
            p.append("ic-bench decomposes by Weil descent of the summation polynomial: "
                     "decomposition must be {method: algebraic, encoding: expanded_semaev}")
        solver = dec.get("solver") or {}
        if solver.get("name") != SOLVER:
            p.append(f"ic-bench's PDP solver is {SOLVER}: decomposition.solver.name must be {SOLVER}")
        else:
            opts = solver.get("options") or {}
            if opts.get("mode") not in ("xl", "mxl"):
                p.append("decomposition.solver.options.mode must be xl or mxl")
            if opts.get("formulation", "direct") != "direct":
                p.append("ic-bench cells use the direct formulation")
        if spec["relations"]["stop"] != "full_rank":
            p.append("ic-bench collects to the achievable rank, then descends: relations.stop must be full_rank")
        if spec["relations"]["collector"] != "sample":
            p.append("ic-bench samples queries R = [k]G: relations.collector must be sample")
        if "max_trials" not in spec["relations"]:
            p.append("set relations.max_trials (bench max_attempts; 200000 in the ci suite)")
        if spec["linear_algebra"]["method"] != "gauss":
            p.append("ic-bench solves the relation matrix by Gaussian elimination mod r: linear_algebra.method gauss")
        if (spec.get("descent") or {}).get("method") != "pdp":
            p.append("ic-bench descends each target with the PDP: descent.method must be pdp")
        ref = (spec.get("reference") or {}).get("rho")
        if ref not in ("rho.plain", "rho.signed_frobenius"):
            p.append("ic-bench prices rho analytically: reference.rho must be rho.plain (its ratio_to_rho) or "
                     "rho.signed_frobenius (its ratio_to_floor), never rho.measured_matched")
        if spec["accounting"]["unit"] not in UNITS:
            p.append(f"ic-bench reports {list(UNITS)}")
        win = spec["measurement"]["window"]
        if win == "online_one_target" and wl["targets"] != 1:
            p.append("online_one_target needs workload.targets = 1")
        if win == "stage_only":
            p.append("ic-bench runs whole DLPs, never a stage alone")
        if spec["execution"]["threads"] != 1:
            p.append("one cell is one single-threaded interpreter: execution.threads must be 1")
        return p

    # ---------------------------------------------------------------- command
    def cell(self, spec: dict[str, Any], workload: dict[str, Any]) -> dict[str, Any]:
        fbp = spec["factor_base"]["params"]
        return {"n": spec["instance"]["curve"]["degree"], "m": spec["decomposition"]["arity"],
                "l": int(fbp["dimension"]), "family": fbp["basis_law"], "seed": int(fbp["basis_seed"]),
                "mode": spec["decomposition"]["solver"]["options"]["mode"],
                "workload_seed": workload["record"]["seed"], "targets": workload["record"]["targets"],
                "max_attempts": spec["relations"]["max_trials"], "run": 1}

    def limits(self, spec: dict[str, Any]) -> dict[str, Any]:
        o = spec["decomposition"]["solver"]["options"]
        return {k: o[k] for k in ("d_max", "max_cols", "max_rows")}

    def command(self, spec: dict[str, Any], workload: dict[str, Any], ctx: Context) -> dict[str, Any]:
        root = _root(spec, ctx)
        argv = ["python3", RUNNER, "--root", root,
                "--cell", json.dumps(self.cell(spec, workload), sort_keys=True),
                "--limits", json.dumps(self.limits(spec), sort_keys=True)]
        return {"argv": argv, "cwd": ctx["repo_root"], "env": dict(spec["execution"].get("env") or {})}

    # ------------------------------------------------------------- provenance
    def provenance(self, spec: dict[str, Any], ctx: Context) -> dict[str, Any]:
        root = _root(spec, ctx)
        bench = os.path.join(root, "experiments", "ic-bench")
        status = _git(root, "status", "--porcelain", "--untracked-files=no")
        return {"cryptanalysis": {"commit": _git(root, "rev-parse", "HEAD"), "dirty": bool(status) if status is not None else None,
                                  "remote": _git(root, "config", "--get", "remote.origin.url")},
                "bench_py_sha256": sha256_file(os.path.join(bench, "bench.py")),
                "calibration_json_sha256": sha256_file(os.path.join(bench, "calibration.json")),
                "runner_sha256": sha256_file(RUNNER),
                "kernel_note": "pdp-degree-heuristics/kernel.py compiles pdpkernel.c with -march=native on first use; "
                               "the record's capsule carries the compiler and CPU flags"}

    # ------------------------------------------------------------------ parse
    def parse(self, spec: dict[str, Any], workload: dict[str, Any], ctx: Context, stdout_path: str) -> dict[str, Any]:
        try:
            with open(stdout_path, encoding="utf-8") as fh:
                doc = json.load(fh)
        except (OSError, ValueError) as exc:
            return {"outcome": {"status": "error", "verified": False, "reason": f"unparseable runner output: {exc}"}}
        if "error" in doc:
            return {"outcome": {"status": "error", "verified": False, "reason": doc["error"]}}
        rec, runner = doc["receipt"], doc.get("icms_runner") or {}
        status = rec["status"] if rec["status"] in ("complete", "insufficient_relations", "budget", "error") else "error"
        verified = bool(rec.get("verified_scalar")) and status == "complete"
        stage, counts = rec.get("stage") or {}, rec.get("counts") or {}
        phases_list = runner.get("phases") or list(rec["phase_operations"])
        phases = {p: {"ops": rec["phase_operations"].get(p), "ops_unit": "cryptanalysis.rps",
                      "wall_ns": rec["phase_wall_ns"].get(p),
                      "native": rec["phase_counters"].get(p) or {}} for p in phases_list}
        ec_adds = sum((rec["phase_counters"].get(p) or {}).get("ec_add", 0) for p in phases_list)
        cold_wall = sum(rec["phase_wall_ns"].get(p, 0) for p in phases_list)
        r = int(rec["subgroup_order"])
        n = spec["instance"]["curve"]["degree"]
        rho_plain, rho_frob = rec.get("rho_operations"), rec.get("rho_floor_operations")
        ref_id = (spec.get("reference") or {}).get("rho")
        total = rec.get("total_operations")
        units = {
            "cryptanalysis.rps": {"total": total, "S": rec.get("S_rps"), "sqrt_r": math.sqrt(r), "deterministic": True,
                                  "host_dependent_because": [], "calibration_id": runner.get("calibration_id"),
                                  "note": "exact counters times frozen per-degree weights; the weights were measured once "
                                          "on the calibration host, so totals pair only under one calibration_id"},
            "count.group_additions": {"total": ec_adds if status == "complete" else None, "deterministic": True,
                                      "host_dependent_because": [],
                                      "note": "the ec_add counter summed over the eleven charged phases; ec_lift, field and "
                                              "Macaulay work are not group additions and are not in this total"},
            "time.wall_ns": {"total": cold_wall, "deterministic": False,
                             "host_dependent_because": ["wall time of the eleven exclusive phase clocks on this host"]},
        }
        online_parts = ["target_descent", "recovery_check"]
        windows = {
            "cold_end_to_end": {"conformance": "exact", "composition": {p: f"ic-bench phase {p}" for p in phases_list},
                                "ops": total, "ops_unit": "cryptanalysis.rps", "wall_ns": cold_wall},
            "whole_process": {"conformance": "runner", "ops": None, "ops_unit": None, "wall_ns": None},
        }
        if workload["record"]["targets"] == 1:
            ops = [rec["phase_operations"].get(p) for p in online_parts]
            windows["online_one_target"] = {
                "conformance": "derived",
                "composition": {"target_descent": "target_query + target_pdp + target_relation_check + target_descent, "
                                                  "on one ic-bench clock", "recovery_check": "recovery_check"},
                "ops": None if status != "complete" or None in ops else float(sum(ops)), "ops_unit": "cryptanalysis.rps",
                "wall_ns": float(sum(rec["phase_wall_ns"].get(p, 0) for p in online_parts)),
                "why_derived": "ic-bench times the four target subphases on one clock (target_descent), so the five "
                               "exclusive online clocks AGENTS.md asks for are not separable"}
        fb_points, cols = stage.get("fb_points"), stage.get("effective_columns")
        l = stage.get("nominal_dimension")
        pdp = rec["phase_counters"].get("pdp") or {}
        return {
            "outcome": {"status": status, "verified": verified, "targets": counts.get("targets"),
                        "targets_verified": counts.get("targets_verified"), "run_id": rec.get("run_id"),
                        "candidate_id": rec.get("candidate_id"), "workload_id": rec.get("workload_id")},
            "units": units,
            "reference": {"id": ref_id,
                          "ops": rho_plain if ref_id == "rho.plain" else rho_frob,
                          "ratio": rec.get("ratio_to_rho") if ref_id == "rho.plain" else rec.get("ratio_to_floor"),
                          "unit": "cryptanalysis.rps",
                          "rho.plain": {"ops": rho_plain, "ratio": rec.get("ratio_to_rho"), "rule": "sqrt(pi r / 2) x w_ec_add"},
                          "rho.signed_frobenius": {"ops": rho_frob, "ratio": rec.get("ratio_to_floor"),
                                                   "rule": f"sqrt(pi r / (4n)) x w_ec_add, n = {n}"},
                          "method": "analytic expected count, not a run"},
            "metrics": {
                "instance": {"curve_id": rec.get("source_curve_ref"), "r": str(r), "log2_r": math.log2(r),
                             "degree": n},
                "factor_base": {"family": "binary_subspace", "producer_name": stage.get("family") or rec["cell"]["family"],
                                "nominal_dimension": l, "usable_points": fb_points, "signed_points": fb_points,
                                "abscissae_with_points": None, "abscissae_allowed": (2 ** l) if l is not None else None,
                                "columns": cols, "orbit_representatives": cols,
                                "points_per_column": (fb_points / cols) if fb_points and cols else None,
                                "native": {"geometric_points": stage.get("geometric_points"),
                                           "factor_base_sha256": stage.get("factor_base_sha256"),
                                           "exact_p_decomposable": stage.get("exact_p_decomposable"),
                                           "achievable_rank": stage.get("achievable_rank")}},
                "decomposition": {"oracle": f"{SOLVER} ({rec['cell']['mode']})", "arity": rec["cell"]["m"],
                                  "targets_tried": counts.get("ordinary_queries"),
                                  "relations_found": counts.get("verified_relations"),
                                  "hit_rate": stage.get("yield_per_query"),
                                  "native": {k: v for k, v in counts.items()}},
                "system": None,
                "solver": {"name": SOLVER, "calls": counts.get("pdp_attempts"), "ops": pdp.get("mac_op"),
                           "op_unit": "mac_op (Macaulay word operations)", "priced_by": "pinned",
                           "wall_ns": rec["phase_wall_ns"].get("pdp"),
                           "extra": {"pdp_counters": pdp, "limits": runner.get("limits"),
                                     "pdp_mac_ops_per_query": stage.get("pdp_mac_ops_per_query")}},
                "linear_algebra": {"method": "gauss", "rows": counts.get("novel_rows"), "rank": counts.get("final_rank"),
                                   "dependent": None, "work": None, "work_unit": None,
                                   "native": rec["phase_counters"].get("relation_la") or {}},
            },
            "phases": phases,
            "windows": windows,
            "consistency": [],
            "producer": {"operation": "ic-bench run_cell", "schema_version": rec.get("schema_version"),
                         "kind": rec.get("kind"), "operation_unit": rec.get("operation_unit"),
                         "provenance": rec.get("provenance"), "resource_envelope": rec.get("resource_envelope"),
                         "instrument": rec.get("instrument"), "warm": rec.get("warm"),
                         "producer_wall_ns": rec.get("wall_ns"), "peak_rss_bytes": rec.get("peak_rss_bytes")},
        }
