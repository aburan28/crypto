#!/usr/bin/env python3
"""Build the lab browser's data file from the repository's frozen records.

The lab browser (docs/browser/) is a searchable index of everything the
index-calculus workstream names: every curve in the ICV1 registry with its
invariants and identities, every ecbench method and factor base that a
committed session measured, every tournament candidate identity (IC1) the
repository mentions, every tournament round, every committed ecbench
session, and the vocabulary of oracles, solvers and factor-base families.

Like the leaderboard (AGENTS.md §7b) it is generated, never edited by hand,
and reads only committed files, so a number that is not in one cannot reach
the page.  `--check` fails when docs/browser/data.json is stale.

This is site tooling, not research or performance code: it joins JSON
records on the identities the repository already computes and computes no
measurement of its own.
"""

from __future__ import annotations

import argparse
import glob
import hashlib
import json
import os
import re
import sys
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
OUT = os.path.join(ROOT, "docs", "browser", "data.json")

REGISTRY = "docs/curves/registry.json"
COVERS = "docs/curves/covers.json"
LEADERBOARD = "docs/ic/leaderboard.json"
TOURNAMENT_RUNS = "research/ic_candidate_tournament_20260915/runs"
ECBENCH_SESSIONS_GLOB = "research/ecbench_*/sessions/*"
IC1_SCAN_GLOBS = (
    "docs/ic/**/*.json",
    "docs/ic/**/*.md",
    "research/ic_candidate_tournament_20260915/**/*.json",
    "research/ic_candidate_tournament_20260915/**/*.md",
    "research/f6_ic_geometric_closure_20261003/**/*.json",
    "research/f6_ic_geometric_closure_20261003/**/*.md",
    "research/notes/ecc2k130/**/*.json",
    "research/notes/ecc2k130/**/*.md",
)
IC1_RE = re.compile(r"IC1N(\d+)C([a-z][a-z0-9]*?)fb(\d+)PDP(\d+)([a-z0-9]+?)RC([a-z0-9]+?)LA([a-z0-9]+?)TD([a-z0-9]+?)ISO(\d+)h([a-f0-9]{12})")
IC1_FULL_RE = re.compile(r"IC1N\d+C[a-z][a-z0-9]*fb\d+PDP\d+[a-z0-9]+RC[a-z0-9]+LA[a-z0-9]+TD[a-z0-9]+ISO\d+h[a-f0-9]{12}")
CELL_RE = re.compile(r"^n(\d+)a([01])$")
MAX_SCAN_BYTES = 8 * 1024 * 1024

# The pipeline's named parts.  These are the strings the code accepts; the
# descriptions paraphrase the source they are implemented in.  Names are the
# only identity these have (an oracle or solver enters a hashed identity only
# through an ecbench method's params or a tournament IC1 code).
VOCABULARY = {
    "oracles": [
        {"name": "subtract", "where": "ic.pipeline (ecbench), ic_run", "what": "Decompose a point by subtracting factor-base points directly; no table."},
        {"name": "mitm", "where": "ic.pipeline, ic_run", "what": "Meet-in-the-middle over pair sums of the base; `m=2|3`; `negation_folded=1` halves the table under the negation map (prime and binary bases)."},
        {"name": "mitm-frobenius", "where": "ic.pipeline, ic_run", "what": "Meet-in-the-middle with the table folded under Frobenius; needs a Frobenius-closed base (Koblitz orbit bases)."},
        {"name": "descent-algebraic", "where": "ic_run", "what": "Summation-polynomial decomposition solved algebraically; needs `--solver`. Nondeterministic under a solver budget, so excluded from cross-method ecbench baselines."},
        {"name": "matrix-f4", "where": "ic_run (hybrid kinds matrix-f4, matrix-f5, inherited-f4)", "what": "Boolean Macaulay matrices to a fixed degree, propagation and splitting (koblitz_groebner MatrixF4)."},
        {"name": "pair_table / pair_or_triple", "where": "tournament candidate configs", "what": "Pair (or pair-then-triple) sum table decomposition in the tournament's native solver; the incumbent's oracle."},
        {"name": "semaev_s3_roots, semaev_s4_pairs_and_solve", "where": "ledger runs (docs/ic/runs)", "what": "Third and fourth summation polynomials solved by root finding; historical ladder variants."},
    ],
    "solvers": [
        {"name": "buchberger-f2", "what": "Gröbner basis over F2 by Buchberger's algorithm."},
        {"name": "f4-f2", "what": "Matrix F4 over F2."},
        {"name": "xl-f2", "what": "Extended linearisation over F2."},
        {"name": "fes-f2, fes-f2-wide", "what": "Fast exhaustive search over F2, Gray-code enumeration."},
        {"name": "crossbred-f2", "what": "Crossbred: partial linearisation then exhaustive search."},
        {"name": "sat-cdcl", "what": "Native CDCL SAT with XOR reasoning on the Boolean encoding."},
        {"name": "exhaustive", "what": "Plain enumeration; the reference for correctness."},
        {"name": "wdsat, cryptominisat, magma-f4", "what": "External solvers named in the SAT factor-base review; historical, not in the native pipeline."},
    ],
    "factor_base_families": [
        {"name": "prime-abscissa", "params": "size", "what": "The points whose abscissa is below a bound on a prime-field curve."},
        {"name": "koblitz-orbit", "params": "divisor", "what": "A Frobenius-invariant subspace of abscissae named by factors of x^n - 1; Frobenius-closed, so the table folds."},
        {"name": "binary-subspace", "params": "dimension", "what": "Abscissae in a polynomial-basis F2-subspace; needs no invariant factor, so not Frobenius-closed."},
        {"name": "compact-orbit-scan", "params": "", "what": "Orbit base built by scanning a compact representative set."},
        {"name": "glv-orbit", "params": "", "what": "Orbit base under a GLV endomorphism on a prime-field curve."},
        {"name": "gls-line", "params": "", "what": "Line base under a GLS endomorphism on an extension-field curve."},
        {"name": "symmetrised / koblitz-symmetrised", "params": "", "what": "Bases closed under negation and Frobenius so that one column stands for an orbit."},
    ],
}


def rel(path: str) -> str:
    return os.path.relpath(path, ROOT)


def load_file(path: str):
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)


def load(path_rel: str):
    return load_file(os.path.join(ROOT, path_rel))


def sha256_file(path_rel: str) -> str:
    h = hashlib.sha256()
    with open(os.path.join(ROOT, path_rel), "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def int_or_none(value):
    try:
        return int(str(value), 0)
    except (TypeError, ValueError):
        return None


def bits(value) -> float | None:
    n = int_or_none(value)
    if n is None or n <= 0:
        return None
    return round(n.bit_length() - 1 + (n / (1 << (n.bit_length() - 1)) - 1), 3)


def curve_rows(registry: dict, leaderboard: dict) -> list[dict]:
    roster = {row["slug"]: row for row in leaderboard.get("roster", [])}
    board = defaultdict(list)
    for row in leaderboard.get("board", []):
        best = row.get("best") or {}
        board[row["slug"]].append(
            {
                "regime": row.get("regime"),
                "log2_r": row.get("log2_r"),
                "reference": row.get("reference"),
                "best_recipe": best.get("variant"),
                "best_s": best.get("s"),
                "ratio_rho": best.get("ratio_rho"),
                "ratio_floor": best.get("ratio_floor"),
                "verified": best.get("verified"),
            }
        )
    rows = []
    for c in registry["curves"]:
        params = c.get("params", {})
        family = c["family"]
        if family == "prime":
            p = params.get("p")
            field = {"characteristic": "p", "p": str(p), "bits": (int_or_none(p) or 0).bit_length(), "degree": 1, "modulus": None}
        else:
            degree = params.get("n") or params.get("m")
            field = {"characteristic": 2, "p": None, "bits": degree, "degree": degree, "modulus": params.get("modulus")}
        try:
            model = json.loads(c.get("model_json") or "{}")
        except ValueError:
            model = {}
        a = params.get("a") if params.get("a") is not None else model.get("a")
        b = params.get("b") if params.get("b") is not None else model.get("b")
        rep = (c.get("representations") or [None])[0]
        subgroup = rep["curve"] if rep else {}
        r = subgroup.get("subgroup_order")
        rows.append(
            {
                "slug": c["slug"],
                "icv1": c["icv1"],
                "family": family,
                "standard_names": c.get("standard_names", []),
                "field": field,
                "a": str(a) if a is not None else None,
                "b": str(b) if b is not None else None,
                "construction": params.get("generator_call"),
                "trace": None if c.get("trace") is None else str(c.get("trace")),
                "order": c.get("order"),
                "order_bits": bits(c.get("order")),
                "r": r,
                "r_bits": bits(r),
                "cofactor": subgroup.get("cofactor"),
                "generator": subgroup.get("generator"),
                "target_group": subgroup.get("target_group"),
                "j": c.get("j"),
                "endomorphism_discriminant": None if c.get("end") in (None, "unk") else c.get("end"),
                "model_json": c.get("model_json"),
                "ec1": [x["ec1"] for x in c.get("representations") or []],
                "curve_uid": [x["curve_uid"] for x in c.get("representations") or []],
                "ec1_unresolved": c.get("ec1_unresolved"),
                "retired_names": c.get("aliases", []),
                "sources": c.get("sources", []),
                "on_leaderboard": bool(roster.get(c["slug"], {}).get("on_board")),
                "leaderboard": board.get(c["slug"], []),
                "ecbench": [],
                "factor_bases": [],
                "tournament_cells": [],
            }
        )
    return rows


def ecbench_sessions(curves_by_slug: dict) -> tuple[list[dict], list[dict], list[dict]]:
    sessions, methods, fbs, yields = [], {}, {}, []
    per_curve_arm = defaultdict(lambda: {"n": 0, "verified": 0, "s_sum": 0.0, "floor_sum": 0.0})
    for sdir in sorted(glob.glob(os.path.join(ROOT, ECBENCH_SESSIONS_GLOB))):
        session_path = os.path.join(sdir, "session.json")
        records_path = os.path.join(sdir, "records.jsonl")
        if not (os.path.exists(session_path) and os.path.exists(records_path)):
            continue
        s = load_file(session_path)
        spec = load_file(os.path.join(sdir, "spec.json")) if os.path.exists(os.path.join(sdir, "spec.json")) else {}
        slugs, arms, statuses = set(), {}, Counter()
        with open(records_path, encoding="utf-8") as fh:
            for line in fh:
                r = json.loads(line)
                statuses[r["outcome"]["status"]] += 1
                curve = r["workload"].get("curve") or {}
                slug = curve.get("slug")
                if slug:
                    slugs.add(slug)
                m = r["method"]
                arms[r["arm"]] = {"arm": r["arm"], "role": r.get("role"), "method": m["id"], "method_id": m["method_id"]}
                entry = methods.setdefault(
                    m["method_id"],
                    {"method_id": m["method_id"], "method_sha256": m.get("method_sha256"), "method": m["id"], "family": m.get("family"), "params": m.get("params", {}), "sessions": set(), "curves": set(), "runs": 0, "verified": 0},
                )
                entry["sessions"].add(s["session_id"])
                entry["runs"] += 1
                if slug:
                    entry["curves"].add(slug)
                if r["outcome"]["status"] == "verified":
                    entry["verified"] += 1
                fb = r.get("factor_base")
                if fb:
                    f = fbs.setdefault(
                        fb["fb_id"],
                        {**{k: fb.get(k) for k in ("fb_id", "fb_sha256", "family", "params", "description", "signed_points", "abscissae", "columns", "dimension", "points_sha256")}, "curve": slug, "sessions": set(), "methods": set()},
                    )
                    f["sessions"].add(s["session_id"])
                    f["methods"].add(m["method_id"])
                if r.get("warmup") or not slug:
                    continue
                if m.get("family") == "ic":
                    yields.append(yield_row(r, s["session_id"]))
                key = (slug, s["session_id"], r["arm"])
                agg = per_curve_arm[key]
                agg["n"] += 1
                agg["method"] = m["id"]
                agg["method_id"] = m["method_id"]
                if r["outcome"]["status"] == "verified":
                    agg["verified"] += 1
                    agg["s_sum"] += float(r["cost"]["s"])
                    agg["floor_sum"] += float(r["boundaries"]["floor_s"])
        sessions.append(
            {
                "session_id": s["session_id"],
                "dir": rel(sdir),
                "label": s.get("label") or spec.get("label"),
                "status": s.get("status"),
                "spec_id": s.get("spec_id"),
                "env_class_id": s.get("env_class_id"),
                "binary_sha256": s.get("binary_sha256"),
                "git_commit": s.get("git_commit"),
                "records": sum(statuses.values()),
                "status_counts": dict(statuses),
                "started_unix_ms": s.get("started_unix_ms"),
                "curves": sorted(slugs),
                "arms": sorted(arms.values(), key=lambda a: a["arm"]),
            }
        )
    for (slug, session_id, arm), agg in per_curve_arm.items():
        if slug in curves_by_slug:
            curves_by_slug[slug]["ecbench"].append(
                {
                    "session_id": session_id,
                    "arm": arm,
                    "method": agg["method"],
                    "method_id": agg["method_id"],
                    "runs": agg["n"],
                    "verified": agg["verified"],
                    "mean_s": round(agg["s_sum"] / agg["verified"], 4) if agg["verified"] else None,
                    "mean_ratio_to_floor": round(agg["s_sum"] / agg["floor_sum"], 3) if agg["floor_sum"] else None,
                }
            )
    for f in fbs.values():
        if f["curve"] in curves_by_slug:
            curves_by_slug[f["curve"]]["factor_bases"].append(f["fb_id"])
    for c in curves_by_slug.values():
        c["ecbench"].sort(key=lambda e: (e["session_id"], e["arm"]))
    method_rows = sorted(
        ({**m, "sessions": sorted(m["sessions"]), "curves": sorted(m["curves"])} for m in methods.values()),
        key=lambda m: (m["family"] or "", m["method"], m["method_id"]),
    )
    fb_rows = sorted(({**f, "sessions": sorted(f["sessions"]), "methods": sorted(f["methods"])} for f in fbs.values()), key=lambda f: (f["family"] or "", f["curve"] or "", f["fb_id"]))
    yields.sort(key=lambda y: (y["curve"], y["fb_id"] or "", y["oracle"] or "", y["session_id"], y["arm"], y["target_index"] if y["target_index"] is not None else -1, y["round"]))
    return sessions, method_rows, fb_rows, yields


def yield_row(r: dict, session_id: str) -> dict:
    """One row of the yield ledger: what one IC run yielded on one target.

    Every figure is the run's own counter; `yield` is relations / trials and
    `lookups_per_relation` is lookups / relations, both left null when the
    denominator is zero.  The `solver` block is copied from the record when
    the harness wrote one (algebraic and SAT oracles); older records have
    none and the column stays null.
    """
    phases = {p["name"]: p for p in r.get("phases", [])}
    rel = (phases.get("relations") or {}).get("native") or {}
    la = (phases.get("linear_algebra") or {}).get("native") or {}
    fb = r.get("factor_base") or {}
    detail = r.get("detail") or {}
    dec = detail.get("decomposition") or {}
    params = r["method"].get("params", {})
    trials = rel.get("trials")
    # ic.pipeline writes `relations`; ic.shared_rank writes `hits` for the
    # same count (a decomposition that became a relation-log row).
    relations = rel.get("relations", rel.get("hits"))
    lookups = rel.get("lookups")
    workload = r["workload"]
    curve = workload.get("curve") or {}
    phase_gae = {p["name"]: p.get("gae") for p in r.get("phases", [])}
    return {
        "run_id": r["run_id"],
        "record_id": r["record_id"],
        "session_id": session_id,
        "arm": r["arm"],
        "round": r.get("round"),
        "status": r["outcome"]["status"],
        "curve": curve.get("slug"),
        "log2_r": None if curve.get("r") is None else bits(curve.get("r")),
        "target_index": workload.get("target_index"),
        "target": workload.get("target"),
        "method_id": r["method"]["method_id"],
        "fb_id": fb.get("fb_id"),
        "fb_family": fb.get("family"),
        "fb_params": fb.get("params"),
        "fb_columns": fb.get("columns"),
        "fb_signed_points": fb.get("signed_points"),
        "oracle": params.get("oracle") or r["method"]["id"],
        "solver_name": params.get("solver") or dec.get("solver"),
        "summands": dec.get("summands"),
        "hit_rate": dec.get("hit_rate"),
        "trials": trials,
        "relations": relations,
        "yield": (relations / trials) if trials and relations is not None else None,
        "lookups": lookups,
        "lookups_per_relation": (lookups / relations) if relations and lookups is not None else None,
        "walk_steps": rel.get("walk_steps"),
        "lift_failures": rel.get("lift_failures"),
        "unliftable_systems": rel.get("unliftable_systems"),
        "frobfold_mismatches": rel.get("frobfold_mismatches"),
        "matrix_rows": la.get("rows"),
        "matrix_columns": la.get("columns"),
        "matrix_rank": la.get("rank"),
        "row_ops": la.get("row_ops"),
        "gae": {k: phase_gae.get(k) for k in ("factor_base", "oracle_setup", "relations", "linear_algebra", "verify")},
        "total_gae": (r.get("cost") or {}).get("total_gae"),
        "s": (r.get("cost") or {}).get("s"),
        "deterministic": (r.get("cost") or {}).get("deterministic"),
        "solver": r.get("solver"),
    }


def tournament_rounds(curves: list[dict]) -> list[dict]:
    koblitz = {(c["field"]["degree"], int(c["a"])): c["slug"] for c in curves if c["family"] == "koblitz" and c["a"] in ("0", "1")}
    by_slug = {c["slug"]: c for c in curves}
    rounds = []
    base = os.path.join(ROOT, TOURNAMENT_RUNS)
    for d in sorted(os.listdir(base)):
        path = os.path.join(base, d)
        if not os.path.isdir(path):
            continue
        row = {"round": d, "dir": rel(path), "report": rel(os.path.join(path, "REPORT.md")) if os.path.exists(os.path.join(path, "REPORT.md")) else None}
        failure = os.path.join(path, "prepare_failure.json")
        if os.path.exists(failure):
            f = load_file(failure)
            row.update({"status": "prepare failed", "reason": f.get("reason")})
            rounds.append(row)
            continue
        dec_path = os.path.join(path, "decision.json")
        if os.path.exists(dec_path):
            dec = load_file(dec_path)
            row.update(
                {
                    "status": dec.get("status"),
                    "winner": dec.get("winner"),
                    "classification": dec.get("classification"),
                    "beats_rho": dec.get("beats_rho"),
                    "beats_rho_strict": dec.get("beats_rho_strict"),
                    "rho_parity": dec.get("rho_parity"),
                    "rho_over_winner": dec.get("rho_over_winner"),
                    "unit": dec.get("unit"),
                    "objective": dec.get("objective"),
                    "scope": dec.get("scope"),
                }
            )
        else:
            row["status"] = "no decision"
        cand_path = os.path.join(path, "candidates.json")
        if os.path.exists(cand_path):
            row["candidates"] = [{"id": c.get("id"), "config": c.get("config"), "parent": c.get("parent"), "hypothesis": c.get("hypothesis"), "configuration_sha256": c.get("configuration_sha256")} for c in load_file(cand_path)]
        fix_path = os.path.join(path, "fixtures.json")
        cells = {}
        if os.path.exists(fix_path):
            fixtures = load_file(fix_path)
            for stage, items in fixtures.items():
                if not isinstance(items, list):
                    continue
                for item in items:
                    cell = item.get("cell")
                    fx = item.get("fixture") or {}
                    if not cell:
                        continue
                    entry = cells.setdefault(cell, {"cell": cell, "stages": set(), "degree": fx.get("degree"), "curve_a": fx.get("curve_a"), "subgroup_order": fx.get("subgroup_order"), "cofactor": fx.get("cofactor"), "slug": None})
                    entry["stages"].add(stage)
                    m = CELL_RE.match(cell)
                    if m:
                        entry["slug"] = koblitz.get((int(m.group(1)), int(m.group(2))))
        row["cells"] = [{**c, "stages": sorted(c["stages"])} for c in sorted(cells.values(), key=lambda c: (c["degree"] or 0, c["cell"]))]
        for c in row["cells"]:
            if c["slug"] and c["slug"] in by_slug:
                by_slug[c["slug"]]["tournament_cells"].append({"round": d, "cell": c["cell"], "stages": c["stages"]})
        rounds.append(row)
    return rounds


def candidate_identities() -> list[dict]:
    found = defaultdict(lambda: {"count": 0, "files": Counter()})
    for pattern in IC1_SCAN_GLOBS:
        for path in sorted(glob.glob(os.path.join(ROOT, pattern), recursive=True)):
            if not os.path.isfile(path) or os.path.getsize(path) > MAX_SCAN_BYTES:
                continue
            with open(path, encoding="utf-8", errors="replace") as fh:
                text = fh.read()
            for ident in IC1_FULL_RE.findall(text):
                found[ident]["count"] += 1
                found[ident]["files"][rel(path)] += 1
    rows = []
    for ident, info in found.items():
        m = IC1_RE.match(ident)
        parts = {}
        if m:
            parts = {
                "n": int(m.group(1)),
                "curve_tag": m.group(2),
                "factor_base_points": int(m.group(3)),
                "summands": int(m.group(4)),
                "solver": m.group(5),
                "collector": m.group(6),
                "linear_algebra": m.group(7),
                "descent": m.group(8),
                "isogeny": int(m.group(9)),
                "record_sha12": m.group(10),
            }
        top = sorted(info["files"].items(), key=lambda kv: (-kv[1], kv[0]))[:6]
        rows.append({"candidate_id": ident, **parts, "mentions": info["count"], "files": len(info["files"]), "where": [p for p, _ in top]})
    rows.sort(key=lambda r: (r.get("n") or 0, r.get("curve_tag") or "", r["candidate_id"]))
    return rows


def attach_covers(curves: list[dict], report: dict) -> None:
    """Join preverified metadata only; mathematical checks are native Rust."""
    if report.get("schema_version") != "curve-covers/v1" or report.get("registry_sha256") != sha256_file(REGISTRY):
        raise ValueError("cover report schema or registry digest is stale")
    findings = report.get("curves", [])
    by_slug = {r["slug"]: r for r in findings}
    if len(by_slug) != len(findings) or set(by_slug) != {c["slug"] for c in curves}:
        raise ValueError("cover report must have exactly one finding per catalog model")
    for curve in curves:
        finding = by_slug[curve["slug"]]
        if finding.get("model_sha256") != hashlib.sha256(curve["model_json"].encode()).hexdigest():
            raise ValueError(f"cover model digest mismatch: {curve['slug']}")
        curve["hyperelliptic_cover"] = finding


def build() -> dict:
    registry = load(REGISTRY)
    leaderboard = load(LEADERBOARD)
    curves = curve_rows(registry, leaderboard)
    attach_covers(curves, load(COVERS))
    by_slug = {c["slug"]: c for c in curves}
    sessions, methods, fbs, yields = ecbench_sessions(by_slug)
    rounds = tournament_rounds(curves)
    candidates = candidate_identities()
    tag_to_slug = {}
    for c in curves:
        if c["family"] == "koblitz" and c["a"] in ("0", "1"):
            tag_to_slug[(c["field"]["degree"], f"kb{c['a']}")] = c["slug"]
    for cand in candidates:
        cand["slug"] = tag_to_slug.get((cand.get("n"), cand.get("curve_tag")))
    sources = {
        REGISTRY: sha256_file(REGISTRY),
        COVERS: sha256_file(COVERS),
        LEADERBOARD: sha256_file(LEADERBOARD),
        **{os.path.join(s["dir"], "records.jsonl"): sha256_file(os.path.join(s["dir"], "records.jsonl")) for s in sessions},
    }
    return {
        "schema_version": "lab-browser/1",
        "generated_by": "scripts/build_lab_browser.py",
        "what_this_is": "A searchable index of the curves, methods, factor bases, candidate identities, tournament rounds and ecbench sessions this repository names, joined on the identities it already computes (ICV1 slugs, EC1 aliases, ECM1 method ids, FB1 factor-base ids, IC1 candidate ids). Every figure is quoted from a committed file named in `sources` or in the row; nothing is computed here beyond means over a session's own verified runs.",
        "counts": {"curves": len(curves), "methods": len(methods), "factor_bases": len(fbs), "candidates": len(candidates), "rounds": len(rounds), "sessions": len(sessions), "yields": len(yields)},
        "sources": sources,
        "vocabulary": VOCABULARY,
        "curves": curves,
        "methods": methods,
        "factor_bases": fbs,
        "candidates": candidates,
        "rounds": rounds,
        "sessions": sessions,
        "yields": yields,
    }


def render(data: dict) -> str:
    return json.dumps(data, indent=1, ensure_ascii=False, sort_keys=False) + "\n"


def comparable(data: dict) -> dict:
    """The index without its volatile parts, for the staleness check.

    How many files mention a candidate identity, and which, changes
    whenever any PR adds a report or a note that writes one, so a check on
    those counts would go red on a merge commit for reasons unrelated to the
    PR under test.  They are informational and refresh on the next
    regeneration; everything else (curves, methods, factor bases, rounds,
    sessions, the identities themselves and their decoded fields) is
    compared exactly.
    """
    out = json.loads(json.dumps(data))
    for cand in out.get("candidates", []):
        for key in ("mentions", "files", "where"):
            cand.pop(key, None)
    return out


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--check", action="store_true", help="exit 1 if docs/browser/data.json is stale")
    parser.add_argument("--out", default=OUT)
    args = parser.parse_args(argv)
    fresh = build()
    text = render(fresh)
    if args.check:
        current = ""
        if os.path.exists(args.out):
            with open(args.out, encoding="utf-8") as fh:
                current = fh.read()
        try:
            same = comparable(json.loads(current)) == comparable(fresh)
        except ValueError:
            same = False
        if not same:
            print(f"{rel(args.out)} is stale; run python3 scripts/build_lab_browser.py", file=sys.stderr)
            return 1
        print(f"{rel(args.out)} is current")
        return 0
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w", encoding="utf-8") as fh:
        fh.write(text)
    data = json.loads(text)
    print(json.dumps(data["counts"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
