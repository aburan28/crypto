"""Adapter for crypto-autoresearcher's prime-field index calculus.

``python -m crypto_autoresearcher.index_calculus solve`` solves one planted
logarithm on a random prime-field curve of a given bit size and reports exact
counters: S3 root solves, group operations (target operations among them),
linear-algebra field operations, relations and attempts, with producer
clocks for the factor base, the tail table, the relation search and the
linear algebra.  ``--rho`` runs Pollard rho on the same instance.

What this adapter has to be honest about, each recorded in the output:

* The curve, the target and the algorithm's randomness all come from the one
  workload seed (curve_seed = target_seed = seed), so the curve is part of
  the workload identity (``instance.curve.selected_by_seed: true``).
* Its rho is a Teske r-adding walk with no negation map, one run, counted in
  walk plus setup operations: ``rho.plain`` (A = 1), never a matched or
  strong reference.  Against the negation walk it flatters IC by sqrt(2)
  (audit X03).
* ``count.s3_solves`` and ``count.group_additions`` are separate units and
  are never summed into one total; the field operations of the linear
  algebra are in neither.
* A factor base of ``size`` abscissae has 2 * size points on a curve of odd
  order (crypto and cryptanalysis count signed points, so their B is twice
  this producer's |F|).
* The numpy scan is turned off (``--no-accel``) so counts do not depend on
  the vectorised path; msolve (the algebraic engine) is refused until its
  version is pinned.
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

UNITS = ("count.s3_solves", "count.group_additions", "time.wall_ns")
ENGINES = {"enumerate": "enumerate", "mitm": "mitm"}
FB_KINDS = {"prime_small_abscissa": "small_x", "prime_random_abscissa": "random"}
PIVOTS = ("min_fill", "min_index")


def _root(spec: dict[str, Any], ctx: Context) -> str:
    root = (spec["execution"].get("args") or {}).get("root", "../crypto-autoresearcher")
    return root if os.path.isabs(root) else os.path.normpath(os.path.join(ctx["repo_root"], root))


def _git(root: str, *args: str) -> str | None:
    try:
        out = subprocess.run(["git", "-C", root, *args], capture_output=True, text=True, timeout=30)
    except (OSError, subprocess.SubprocessError):
        return None
    return out.stdout.strip() if out.returncode == 0 else None


class AutoresearcherIc:
    name = "autoresearcher.index_calculus"

    def __init__(self, registry: Registry | None = None):
        self.reg = registry or Registry.load()

    # ------------------------------------------------------------------ check
    def check(self, spec: dict[str, Any]) -> list[str]:
        p: list[str] = []
        curve, wl = spec["instance"]["curve"], spec["instance"]["workload"]
        if curve["regime"] != "prime" or "bits" not in curve:
            p.append("the autoresearcher solver draws prime-field curves by size: instance.curve {regime: prime, bits}")
        if "ref" in curve or "explicit" in curve:
            p.append("the solver draws its curve from bits and the seed; it does not resolve instance.curve.ref or explicit")
        if not curve.get("selected_by_seed"):
            p.append("the curve is drawn from the workload seed: set instance.curve.selected_by_seed = true")
        if wl["targets"] != 1 or wl["law"] != "known_answer":
            p.append("solve plants one target from the seed: workload {targets: 1, law: known_answer}")
        fb = spec["factor_base"]
        if fb["family"] not in FB_KINDS:
            p.append(f"factor_base.family must be one of {sorted(FB_KINDS)} (the subgroup base needs the "
                     "--subgroup-prime filter, which the spec cannot express yet)")
        if fb.get("quotient") != ["negation"]:
            p.append("one unknown per abscissa: factor_base.quotient must be [negation]")
        if fb.get("large_primes", {}).get("mode", "none") != "none" or fb.get("partition"):
            p.append("the autoresearcher solver has no large-prime or partitioned base")
        dec = spec["decomposition"]
        if dec["method"] not in ENGINES:
            p.append(f"decomposition.method must be one of {sorted(ENGINES)}; msolve is refused until its version is pinned")
        if (dec.get("params") or {}).get("accelerate"):
            p.append("the numpy scan changes no counts but is not exercised by this adapter: drop params.accelerate")
        extra = set(dec.get("params") or {}) - {"table_arity", "accelerate"}
        if extra:
            p.append(f"unknown decomposition.params {sorted(extra)}")
        if dec["method"] != "mitm" and "table_arity" in (dec.get("params") or {}):
            p.append("table_arity is a mitm parameter")
        if dec.get("encoding", "none") != "none" or dec.get("splitting", "none") != "none":
            p.append("the solver has no encoding or splitting choice: decomposition.encoding and splitting must be none")
        if spec["linear_algebra"].get("filtering"):
            p.append("the solver does no relation filtering: drop linear_algebra.filtering")
        if spec["relations"]["stop"] != "first_log" or spec["relations"]["collector"] != "random":
            p.append("solve stops at the first relation set that fixes k: relations {collector: random, stop: first_log}")
        la = spec["linear_algebra"]
        if la["method"] != "gauss" or la.get("pivot") not in PIVOTS:
            p.append(f"linear_algebra must be {{method: gauss, pivot: one of {list(PIVOTS)}}}; the pivot rule moves la_ops")
        if (spec.get("descent") or {}).get("method", "embedded") != "embedded":
            p.append("the target is written into every relation: descent.method must be embedded")
        ref = (spec.get("reference") or {}).get("rho")
        if ref not in (None, "rho.plain"):
            p.append("the solver's rho is a Teske walk without the negation map: reference.rho must be rho.plain")
        if spec["accounting"]["unit"] not in UNITS:
            p.append(f"the solver reports {list(UNITS)}")
        if spec["measurement"]["window"] not in ("cold_end_to_end", "whole_process"):
            p.append("the solver has no separate online clock: window must be cold_end_to_end or whole_process")
        if spec["execution"]["threads"] != 1:
            p.append("execution.threads must be 1")
        return p

    # ---------------------------------------------------------------- command
    def command(self, spec: dict[str, Any], workload: dict[str, Any], ctx: Context) -> dict[str, Any]:
        root = _root(spec, ctx)
        seed = str(workload["record"]["seed"])
        dec, fb = spec["decomposition"], spec["factor_base"]
        argv = ["python3", "-m", "crypto_autoresearcher.index_calculus", "solve",
                "--bits", str(spec["instance"]["curve"]["bits"]), "--m", str(dec["arity"]),
                "--fb", FB_KINDS[fb["family"]], "--fb-size", str(fb["params"]["size"]),
                "--engine", ENGINES[dec["method"]],
                "--curve-seed", seed, "--target-seed", seed, "--seed", seed,
                "--la-pivot", spec["linear_algebra"]["pivot"], "--no-accel"]
        if "table_arity" in (dec.get("params") or {}):
            argv += ["--table-arity", str(dec["params"]["table_arity"])]
        if (spec.get("reference") or {}).get("rho"):
            argv.append("--rho")
        env = dict(spec["execution"].get("env") or {})
        env["PYTHONPATH"] = os.path.join(root, "src")
        return {"argv": argv, "cwd": root, "env": env}

    # ------------------------------------------------------------- provenance
    def provenance(self, spec: dict[str, Any], ctx: Context) -> dict[str, Any]:
        root = _root(spec, ctx)
        pkg = os.path.join(root, "src", "crypto_autoresearcher", "index_calculus")
        status = _git(root, "status", "--porcelain", "--untracked-files=no")
        sources = {}
        if os.path.isdir(pkg):
            sources = {f: sha256_file(os.path.join(pkg, f)) for f in sorted(os.listdir(pkg)) if f.endswith(".py")}
        try:
            import numpy  # noqa: F401
            numpy_version = numpy.__version__
        except ImportError:
            numpy_version = None
        return {"autoresearcher": {"commit": _git(root, "rev-parse", "HEAD"),
                                   "dirty": bool(status) if status is not None else None,
                                   "remote": _git(root, "config", "--get", "remote.origin.url")},
                "index_calculus_sources_sha256": sources, "numpy_version_in_runner": numpy_version}

    # ------------------------------------------------------------------ parse
    def parse(self, spec: dict[str, Any], workload: dict[str, Any], ctx: Context, stdout_path: str) -> dict[str, Any]:
        try:
            with open(stdout_path, encoding="utf-8") as fh:
                doc = json.load(fh)
        except (OSError, ValueError) as exc:
            return {"outcome": {"status": "error", "verified": False, "reason": f"unparseable solver output: {exc}"}}
        ic, curve, rho = doc["index_calculus"], doc["curve"], doc.get("rho")
        verified = bool(ic.get("verified")) and bool(ic.get("correct"))
        order = int(curve["order"])
        fbi = ic.get("factor_base") or {}
        size = fbi.get("size")
        # The producer's s3_solves already includes the tail table's solves
        # (solver.py: stats.s3_solves + table.s3_solves); table_s3_solves is
        # the same work reported again on its own.
        s3 = ic.get("s3_solves")
        table_s3 = ic.get("table_s3_solves") or 0
        sec = {k: ic.get(f"seconds_{k}") for k in ("factor_base", "table", "relations", "linalg", "total")}

        def ns(x: float | None) -> float | None:
            return None if x is None else x * 1e9

        units = {
            "count.s3_solves": {"total": s3 if verified else None, "deterministic": True, "host_dependent_because": [],
                                "note": "the producer's s3_solves, table solves included (one modular square root each)"},
            "count.group_additions": {"total": ic.get("group_ops") if verified else None, "deterministic": True,
                                      "host_dependent_because": [],
                                      "note": "group_ops, target_ops included; S3 solves and LA field operations excluded"},
            "time.wall_ns": {"total": ns(sec["total"]), "deterministic": False,
                             "host_dependent_because": ["producer perf_counter seconds_total on this host"]},
        }
        relation_note = "queries, decompositions and relation checks on one producer clock (seconds_relations)"
        phases = {
            "factor_base": {"ops": None, "ops_unit": None, "wall_ns": ns(sec["factor_base"]), "native": {}},
            "precompute": {"ops": None, "ops_unit": None, "wall_ns": ns(sec["table"]),
                           "native": {"table_entries": ic.get("table_entries"), "table_s3_solves": ic.get("table_s3_solves")}},
            "pdp": {"ops": None, "ops_unit": None, "wall_ns": ns(sec["relations"]),
                    "native": {"includes": relation_note,
                               "s3_solves": None if s3 is None else s3 - table_s3,
                               "membership_tests": ic.get("membership_tests")}},
            "relation_la": {"ops": None, "ops_unit": None, "wall_ns": ns(sec["linalg"]),
                            "native": {"la_ops": ic.get("la_ops"), "la_pivot": ic.get("la_pivot")}},
        }
        reference = None
        if rho is not None:
            reference = {"id": "rho.plain", "method": "Teske r-adding walk, no negation map, one run",
                         "verified": rho.get("verified"), "group_ops": rho.get("group_ops"), "walk_ops": rho.get("walk_ops"),
                         "setup_ops": rho.get("setup_ops"), "walks": rho.get("walks"), "runs": 1,
                         "ratio_group_ops": (ic["group_ops"] / rho["group_ops"]) if rho.get("group_ops") else None,
                         "note": "one rho run is one sample of a random variable; its ratio is descriptive"}
        return {
            "outcome": {"status": "complete" if verified else "error", "verified": verified,
                        **({} if verified else {"reason": "solver did not verify the planted logarithm"})},
            "units": units,
            "reference": reference,
            "metrics": {
                "instance": {"p": str(curve["p"]), "a": str(curve["a"]), "b": str(curve["b"]), "r": str(order),
                             "log2_r": math.log2(order), "bits": spec["instance"]["curve"]["bits"]},
                "factor_base": {"family": spec["factor_base"]["family"], "producer_name": fbi.get("kind"),
                                "nominal_dimension": None,
                                "usable_points": 2 * size if size is not None and order % 2 else None,
                                "signed_points": 2 * size if size is not None and order % 2 else None,
                                "abscissae_with_points": size,
                                "abscissae_allowed": fbi.get("bound") if fbi.get("kind") == "small_x" else None,
                                "columns": size, "orbit_representatives": size,
                                "points_per_column": 2.0 if size is not None and order % 2 else None,
                                "native": {k: v for k, v in fbi.items() if k not in ("kind", "size")}},
                "decomposition": {"oracle": ic.get("engine"), "arity": ic.get("m"), "targets_tried": ic.get("attempts"),
                                  "relations_found": ic.get("relations"),
                                  "hit_rate": (ic["relations"] / ic["attempts"]) if ic.get("attempts") else None,
                                  "native": {k: ic.get(k) for k in ("s3_solves", "membership_tests", "table_arity",
                                                                    "table_entries", "table_s3_solves", "censored_attempts",
                                                                    "accelerated", "target_ops")}},
                "system": None,
                "solver": None,
                "linear_algebra": {"method": "gauss", "rows": ic.get("relations"), "rank": ic.get("rank"),
                                   "dependent": None, "work": ic.get("la_ops"), "work_unit": "la_ops (field operations)",
                                   "native": {"pivot": ic.get("la_pivot")}},
            },
            "phases": phases,
            "windows": {
                "cold_end_to_end": {"conformance": "derived", "ops": ic.get("group_ops") if verified else None,
                                    "ops_unit": "count.group_additions", "wall_ns": ns(sec["total"]),
                                    "composition": {"factor_base": "seconds_factor_base", "precompute": "seconds_table",
                                                    "pdp": relation_note, "relation_la": "seconds_linalg"},
                                    "why_derived": "the producer counts group_ops before its final [k]P = Q check, so the "
                                                   "recovery check's additions are not in the total, and it has no "
                                                   "separate recovery_check clock"},
                "whole_process": {"conformance": "runner", "ops": None, "ops_unit": None, "wall_ns": None},
            },
            "consistency": [],
            "producer": {"operation": "index_calculus solve", "k": ic.get("k"), "P": doc.get("P"), "Q": doc.get("Q"),
                         "msolve_calls": ic.get("msolve_calls")},
        }
