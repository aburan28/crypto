"""Adapter for the crypto repository's pluggable IC framework: ``ic bench``.

``ic bench`` runs one configuration end to end on one planted target and
prices every stage in group-addition equivalents (GAE).  Its report row
(src/cryptanalysis/ic_framework/mod.rs RunReport) carries the factor-base
sizes, the decomposition counters, the system shape and solver report for
algebraic oracles, the relation-matrix work, the verification, S and the
matched counted rho.

What this adapter has to be honest about, each recorded in the output:

* Every shipped polynomial solver is priced from host wall time
  (``priced_by: measured``), and binary-field table lookups have no pinned
  price.  Such a total is host-dependent; the record says so and the unit
  carries the env_class_id.
* ``ic bench`` stops when the target column is pinned (``first_log``), not at
  full rank.  A spec asking for another stop rule is refused.
* Relation search writes the target into every relation, so the
  ``online_one_target`` window is *derived* from coarser phases, not timed on
  the five exclusive online clocks.  It is usable for operation counts and
  refused as a primary wall-clock speedup.
* Thread count: RAYON_NUM_THREADS is fixed by the runner; F4-family solvers
  parallelise above one thread and would otherwise move wall-priced S.
"""
from __future__ import annotations

import json
import math
import os
from typing import Any

from ..canonical import sha256_file
from ..registry import Registry
from . import Context

METHOD_TO_ORACLE = {"subtract": "subtract", "mitm": "mitm", "pair_table": "mitm",
                    "mitm_frobenius": "mitm-frobenius", "algebraic": "descent-algebraic",
                    "symmetrised": "symmetrised"}
LINALG = {"incremental_gauss": "incremental-gauss", "structured_gauss": "structured-gauss"}
# The symmetrised oracle ignores --solver: it builds its own engine from its own
# parameters (plugins.rs SymmetrisedOracle::prepare), so the spec's solver is
# translated into those parameters.
SYMMETRISED_ENGINES = ("inherited-f4", "matrix-f4", "matrix-f5")
UNITS = ("crypto.S.gae_pinned", "time.wall_ns", "time.cpu_ns")
WINDOWS = ("cold_end_to_end", "online_one_target", "whole_process")


def _plugin(name: str, params: dict[str, Any]) -> str:
    if not params:
        return name
    parts = []
    for k in sorted(params):
        v = params[k]
        if isinstance(v, (list, tuple)):
            v = ";".join(str(x) for x in v)
        elif isinstance(v, bool):
            v = "1" if v else "0"
        parts.append(f"{k}={v}")
    return f"{name}:{','.join(parts)}"


class CryptoIcBench:
    name = "crypto.ic_bench"

    def __init__(self, registry: Registry | None = None):
        self.reg = registry or Registry.load()

    # ------------------------------------------------------------------ check
    def check(self, spec: dict[str, Any]) -> list[str]:
        p: list[str] = []
        curve, wl = spec["instance"]["curve"], spec["instance"]["workload"]
        if curve["regime"] not in ("prime", "char2", "koblitz"):
            p.append(f"ic bench has no {curve['regime']} instances")
        if curve["regime"] == "prime" and "bits" not in curve:
            p.append("ic bench draws prime-field curves by subgroup size: set instance.curve.bits")
        if curve["regime"] == "prime" and "family" not in curve:
            p.append("ic bench defaults a prime curve's family to generic: set instance.curve.family explicitly")
        if "ref" in curve or "explicit" in curve:
            p.append("ic bench builds its curve from regime and degree or bits; it does not resolve instance.curve.ref "
                     "or explicit, so a spec naming one would run a different curve than it says")
        if curve["regime"] in ("char2", "koblitz") and "degree" not in curve:
            p.append("ic bench needs instance.curve.degree")
        if curve["regime"] == "char2" and not curve.get("selected_by_seed"):
            p.append("a random char2 curve is drawn from the seed: set instance.curve.selected_by_seed = true")
        if wl["targets"] != 1:
            p.append("ic bench solves one planted target per run; workload.targets must be 1")
        if wl["law"] != "known_answer":
            p.append("ic bench plants its target from the seed: workload.law must be known_answer")
        fam = self.reg.family(spec["factor_base"]["family"]) or {}
        if "crypto.ic_bench" not in fam:
            p.append(f"ic bench has no factor base for family {spec['factor_base']['family']!r}")
        fbp = spec["factor_base"].get("params") or {}
        if spec["factor_base"]["family"] == "binary_subspace":
            if fbp.get("basis_law") != "prefix":
                p.append("ic bench's binary-subspace base is span{1, z, ...}: factor_base.params.basis_law must be prefix")
            if "basis_seed" in fbp:
                p.append("the prefix basis has no sampler: drop factor_base.params.basis_seed")
        if spec["factor_base"].get("large_primes", {}).get("mode", "none") != "none":
            p.append("ic bench has no large-prime variant")
        if spec["factor_base"].get("partition"):
            p.append("ic bench has no partitioned factor base")
        dec = spec["decomposition"]
        if dec["method"] not in METHOD_TO_ORACLE:
            p.append(f"ic bench has no oracle for method {dec['method']!r}")
        if dec["method"] == "algebraic" and dec.get("encoding") != "expanded_semaev":
            p.append("ic bench's descent-algebraic oracle descends the expanded summation polynomial: encoding must be expanded_semaev")
        if dec["method"] == "symmetrised" and dec.get("encoding") != "symmetrized":
            p.append("the symmetrised oracle needs encoding symmetrized")
        if dec["method"] in ("subtract",) and dec["arity"] != 2:
            p.append("subtract is the m = 2 oracle")
        if dec.get("splitting", "none") not in ("none", "mitm_table"):
            p.append("ic bench has no split decomposition")
        if "max_trials" not in spec["relations"]:
            p.append("ic bench gives up after a trial budget (default 2,000,000): set relations.max_trials")
        if dec["method"] in ("algebraic", "symmetrised") and (dec.get("limits") or {}).get("per_call_seconds") is None:
            p.append("ic bench gives each solver call a wall-clock budget (default 120 s): set "
                     "decomposition.limits.per_call_seconds (0 for none)")
        if dec["method"] == "symmetrised" and dec.get("solver"):
            sv = dec["solver"]
            if sv["name"] not in SYMMETRISED_ENGINES:
                p.append(f"the symmetrised oracle runs its own engine: solver.name must be one of {list(SYMMETRISED_ENGINES)}")
            opts = sv.get("options") or {}
            if opts.get("split", "auto") != "auto":
                p.append("the symmetrised oracle always uses the default split rule: solver.options.split must be auto")
            extra = set(opts) - {"max_degree", "node_budget", "split"}
            if extra:
                p.append(f"the symmetrised oracle reads only max_degree and node_budget, not {sorted(extra)}")
        if spec["relations"]["stop"] != "first_log":
            p.append("ic bench stops when the target column is pinned: relations.stop must be first_log")
        if spec["relations"]["collector"] not in ("walk", "random"):
            p.append("ic bench targets are walk or random")
        if spec["linear_algebra"]["method"] not in LINALG:
            p.append(f"ic bench has no relation solver {spec['linear_algebra']['method']!r}")
        if spec["accounting"]["unit"] not in UNITS:
            p.append(f"ic bench reports {UNITS}, not {spec['accounting']['unit']!r}")
        if spec["measurement"]["window"] not in WINDOWS:
            p.append(f"ic bench supports windows {WINDOWS}")
        rho = (spec.get("reference") or {}).get("rho")
        if rho and rho != "rho.measured_matched":
            p.append("ic bench prices its own matched counted rho: reference.rho must be rho.measured_matched")
        if spec["execution"]["threads"] != 1:
            p.append("the IC framework is single-threaded by design (FRAMEWORK.md section 8); set threads = 1")
        return p

    # ---------------------------------------------------------------- command
    def binary(self, spec: dict[str, Any], ctx: Context) -> str:
        rel = (spec["execution"].get("binary") or {}).get("path", "target/release/ic")
        return rel if os.path.isabs(rel) else os.path.join(ctx["repo_root"], rel)

    def command(self, spec: dict[str, Any], workload: dict[str, Any], ctx: Context) -> dict[str, Any]:
        curve, fb, dec = spec["instance"]["curve"], spec["factor_base"], spec["decomposition"]
        argv = [self.binary(spec, ctx), "bench"]
        if curve["regime"] == "prime":
            argv += ["--bits", str(curve["bits"]), "--family", curve["family"]]
        elif curve["regime"] == "char2":
            argv += ["--char2-degree", str(curve["degree"])]
            if "max_cofactor" in curve:
                argv += ["--max-cofactor", str(curve["max_cofactor"])]
        else:
            argv += ["--koblitz-degree", str(curve["degree"]), "--koblitz-a", str(curve["koblitz_a"])]
        fam = self.reg.family(fb["family"])
        fb_params = dict(fb.get("params") or {})
        fb_params.pop("basis_law", None)  # checked to be prefix, ic bench's only basis
        if fb["family"] in ("frobenius_stable_subspace", "frobenius_symmetrised", "glv_orbit", "gls_line"):
            folds_frob = "frobenius" in fb.get("quotient", []) or "automorphism" in fb.get("quotient", []) or "gls" in fb.get("quotient", [])
            if not folds_frob:
                fb_params["no_fold"] = True
        argv += ["--factor-base", _plugin(fam["crypto.ic_bench"], fb_params)]
        oparams = dict(dec.get("params") or {})
        if dec["method"] != "subtract":
            oparams.setdefault("m", dec["arity"])
        solver = dec.get("solver")
        if dec["method"] == "symmetrised" and solver:
            opts = solver.get("options") or {}
            oparams["engine"] = solver["name"]
            oparams.update({k: opts[k] for k in ("max_degree", "node_budget") if k in opts})
            solver = None
        argv += ["--oracle", _plugin(METHOD_TO_ORACLE[dec["method"]], oparams)]
        if solver:
            argv += ["--solver", _plugin(solver["name"], solver.get("options") or {})]
        argv += ["--targets", spec["relations"]["collector"],
                 "--linalg", LINALG[spec["linear_algebra"]["method"]],
                 "--seed", str(workload["record"]["seed"]),
                 "--repeats", "1",
                 "--rho-runs", str((spec.get("reference") or {}).get("rho_runs", 16)),
                 "--json"]
        if "max_trials" in spec["relations"]:
            argv += ["--max-trials", str(spec["relations"]["max_trials"])]
        budget = (dec.get("limits") or {}).get("per_call_seconds")
        if budget is not None:
            argv += ["--solver-budget-seconds", str(budget)]
        argv += [str(x) for x in (spec["execution"].get("args") or {}).get("extra_argv", [])]
        return {"argv": argv, "cwd": ctx["repo_root"], "env": dict(spec["execution"].get("env") or {})}

    # ------------------------------------------------------------- provenance
    def provenance(self, spec: dict[str, Any], ctx: Context) -> dict[str, Any]:
        b = self.binary(spec, ctx)
        cal = os.path.join(ctx["repo_root"], "docs", "ic", "calibration.json")
        lock = os.path.join(ctx["repo_root"], "Cargo.lock")
        return {"binary": {"path": os.path.relpath(b, ctx["repo_root"]), "sha256": sha256_file(b)},
                "calibration_json_sha256": sha256_file(cal),
                "cargo_lock_sha256": sha256_file(lock),
                "cargo_lock_note": "Cargo.lock is gitignored in this repository; this is the lock file present in the build tree"}

    # ------------------------------------------------------------------ parse
    def parse(self, spec: dict[str, Any], workload: dict[str, Any], ctx: Context, stdout_path: str) -> dict[str, Any]:
        try:
            with open(stdout_path, encoding="utf-8") as fh:
                rep = json.load(fh)
        except (OSError, ValueError) as exc:
            return {"outcome": {"status": "error", "verified": False, "reason": f"unparseable report: {exc}"}}
        rows = rep.get("rows") or []
        if len(rows) != 1:
            return {"outcome": {"status": "error", "verified": False,
                                "reason": f"expected exactly one row, got {len(rows)}; skipped: {rep.get('configurations_skipped')}"},
                    "producer": {"operation": rep.get("operation"), "status": rep.get("status")}}
        row = rows[0]
        fb, dec, la = row["factor_base"], row["decomposition"], row["linear_algebra"]
        solver = dec.get("solver")
        system = dec.get("system")
        cal = rep.get("calibration") or {}
        measured_units = (rep.get("calibration_pins") or {}).get("measured") or []
        solver_wall_priced = bool(solver) and solver.get("priced_by") != "pinned"
        native = (dec.get("cost") or {}).get("native") or {}
        unpriced_counters = []
        for counter, unit_key in (("lookups", "ns_per_lookup"), ("target_guard_probes", "ns_per_lookup"),
                                  ("sqrt_solves", "ns_per_sqrt"), ("inversions", "ns_per_inversion"),
                                  ("s4_pairs", "ns_per_s4_pair")):
            if native.get(counter) and not cal.get(unit_key):
                unpriced_counters.append(counter)
        deterministic = not measured_units and not solver_wall_priced
        status = "complete" if row.get("verified") else ("budget" if row.get("exhausted") else "error")

        ppc = fb.get("points_per_column")
        observed_quotient = []
        if ppc is not None:
            if ppc >= 1.5:
                observed_quotient.append("negation")
            if ppc > 2.5:
                observed_quotient.append("frobenius_or_automorphism")
        declared = set(spec["factor_base"].get("quotient") or [])
        consistency = []
        if ("negation" in declared) != ("negation" in observed_quotient):
            consistency.append({"field": "factor_base.quotient", "declared": sorted(declared),
                                "observed_points_per_column": ppc,
                                "problem": "the negation fold the spec declares is not what the report shows"})
        if bool(declared & {"frobenius", "automorphism", "gls"}) != ("frobenius_or_automorphism" in observed_quotient):
            consistency.append({"field": "factor_base.quotient", "declared": sorted(declared),
                                "observed_points_per_column": ppc,
                                "problem": "the Frobenius/automorphism fold the spec declares is not what the report shows"})
        dim = (spec["factor_base"].get("params") or {}).get("dimension")

        def phase(cost: dict[str, Any] | None) -> dict[str, Any] | None:
            if not cost:
                return None
            return {"ops": cost.get("gae"), "ops_unit": "crypto.S.gae_pinned", "wall_ns": cost.get("wall_ns"),
                    "group_ops": cost.get("group_ops"), "native": cost.get("native")}

        phases = {
            "factor_base": phase(fb.get("cost")),
            "precompute": phase(dec.get("setup")),
            "target_pdp": phase(dec.get("cost")),
            "relation_la": phase(la.get("cost")),
            "recovery_check": phase(row.get("verify")),
        }
        online_parts = ["target_pdp", "relation_la", "recovery_check"]
        cold_parts = ["factor_base", "precompute", "target_pdp", "relation_la", "recovery_check"]

        def total(parts: list[str], key: str) -> float | None:
            vals = [(phases[p] or {}).get(key) for p in parts]
            return None if any(v is None for v in vals) else float(sum(vals))

        r = row.get("r")
        sqrt_r = math.sqrt(int(r)) if r else None
        rho = rep.get("rho_reference") or {}
        windows = {
            "cold_end_to_end": {"conformance": "exact_operations",
                                "composition": {p: f"ic bench {p}" for p in cold_parts},
                                "ops": row.get("total_gae"), "ops_unit": "crypto.S.gae_pinned",
                                "wall_ns": total(cold_parts, "wall_ns")},
            "online_one_target": {"conformance": "derived",
                                  "composition": {"target_pdp": "decomposition.cost (relation search writes the target into every relation; target queries are walk steps inside it)",
                                                  "relation_la": "linear_algebra.cost (the log is read off as rows arrive)",
                                                  "recovery_check": "verify"},
                                  "ops": total(online_parts, "ops"), "ops_unit": "crypto.S.gae_pinned",
                                  "wall_ns": total(online_parts, "wall_ns"),
                                  "why_derived": "ic bench has no exclusive five-phase online clock; use ic price --single-target for a primary wall-clock claim"},
            "whole_process": {"conformance": "runner", "ops": None, "ops_unit": None, "wall_ns": None},
        }
        return {
            "outcome": {"status": status, "verified": bool(row.get("verified")), "recovered": row.get("recovered"),
                        "exhausted": row.get("exhausted"), "producer_status": rep.get("status")},
            "units": {
                "crypto.S.gae_pinned": {
                    "total": row.get("total_gae"), "S": row.get("s"), "sqrt_r": sqrt_r,
                    "deterministic": deterministic,
                    "host_dependent_because": ([f"calibration ratios measured on this host: {measured_units}"] if measured_units else [])
                                              + (["solver priced from host wall time (priced_by=%s)" % solver.get("priced_by")] if solver_wall_priced else []),
                    "unpriced_counters": unpriced_counters,
                    "calibration_in_report": cal,
                },
            },
            "reference": {"id": "rho.measured_matched", "method": rho.get("method"), "automorphisms": rho.get("automorphisms"),
                          "runs": rho.get("runs"), "mean_s": rho.get("mean_s"), "all_verified": rho.get("all_verified"),
                          "rule": rep.get("rho_reference_rule"), "s_over_rho": row.get("s_over_rho")},
            "metrics": {
                "instance": {"slug": row.get("instance") or rep.get("instance"), "r": str(r) if r is not None else None,
                             "log2_r": row.get("log2_r"), "group_order": str(row.get("group_order")) if row.get("group_order") is not None else None},
                "factor_base": {"family": spec["factor_base"]["family"], "producer_name": fb.get("name"),
                                "nominal_dimension": dim, "usable_points": fb.get("signed_points"),
                                "signed_points": fb.get("signed_points"), "abscissae_with_points": fb.get("abscissae"),
                                "abscissae_allowed": None, "columns": fb.get("columns"),
                                "points_per_column": ppc, "orbit_representatives": fb.get("columns")},
                "decomposition": {"oracle": dec.get("name"), "arity": dec.get("summands"), "targets_tried": dec.get("targets_tried"),
                                  "relations_found": dec.get("relations_found"), "hit_rate": dec.get("hit_rate"),
                                  "native": native},
                "system": None if not system else {
                    "n_vars": system.get("n_vars"), "n_equations": system.get("n_equations"),
                    "max_degree": max(system.get("degrees") or [0]) or None, "degrees": system.get("degrees"),
                    "semi_regular_degree": system.get("semi_regular_degree"),
                    "anf_monomials": None, "cnf_variables": None, "cnf_aux_variables": None,
                    "cnf_clauses": (solver or {}).get("extra", {}).get("original_clauses"),
                    "xor_rows": (solver or {}).get("extra", {}).get("native_xor_rows"), "literals": None,
                    "note": "cnf_clauses and xor_rows are summed over every solver call (ic_framework SolverTotals)"},
                "solver": None if not solver else {
                    "name": solver.get("name"), "calls": solver.get("calls"), "ops": solver.get("ops"),
                    "op_unit": solver.get("op_unit"), "priced_by": solver.get("priced_by"),
                    "wall_ns": solver.get("wall_ns"), "budget_exceeded": solver.get("budget_exceeded"),
                    "solving_degree_mean": solver.get("solving_degree_mean"), "solving_degree_max": solver.get("solving_degree_max"),
                    "peak_bytes": solver.get("peak_bytes") or None,
                    "conflicts": solver.get("ops") if solver.get("op_unit") == "conflicts" else None,
                    "decisions": solver.get("extra", {}).get("decisions"), "propagations": solver.get("extra", {}).get("propagations"),
                    "restarts": solver.get("extra", {}).get("restarts"), "extra": solver.get("extra")},
                "linear_algebra": {"method": la.get("name"), "rows": la.get("rows"), "rank": la.get("rank"),
                                   "dependent": la.get("dependent"), "work": la.get("work"), "work_unit": la.get("work_unit")},
            },
            "phases": phases,
            "windows": windows,
            "consistency": consistency,
            "producer": {"operation": rep.get("operation"), "schema_version": rep.get("schema_version"),
                         "label": row.get("label"), "spec": row.get("spec"), "software": rep.get("software"),
                         "producer_wall_seconds": row.get("wall_seconds")},
        }
