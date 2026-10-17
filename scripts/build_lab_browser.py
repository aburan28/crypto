#!/usr/bin/env python3
"""Build the lab browser's data file from the repository's frozen records.

The lab browser (docs/browser/) is a searchable index of everything the
index-calculus workstream names: every curve in the ICV1 registry with its
invariants and identities, every ecbench method and factor base that a
committed session measured, every tournament candidate identity (IC1) the
repository mentions, every tournament round, every committed ecbench
session, and the vocabulary of oracles, solvers and factor-base families.

Like the leaderboard (AGENTS.md §7b) it is generated, never edited by hand.
Sparse-excluded tracked inputs are read from committed HEAD, while materialized
inputs use their worktree bytes.  `--check` fails when docs/browser/data.json
is stale.

This is site tooling, not research or performance code: it joins JSON
records on the identities the repository already computes and computes no
measurement of its own.
"""

from __future__ import annotations

import argparse
import contextlib
import fnmatch
import glob
import hashlib
import json
import os
import re
import subprocess
import sys
from collections import Counter, defaultdict
from functools import lru_cache

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
OUT = os.path.join(ROOT, "docs", "browser", "data.json")

REGISTRY = "docs/curves/registry.json"
COVERS = "docs/curves/covers.json"
COVER_LINKS = "docs/curves/cover-links.yaml"
STANDARDS_COVERAGE = "docs/curves/standards/coverage.json"
LEADERBOARD = "docs/ic/leaderboard.json"
TOURNAMENT_RUNS = "research/ic_candidate_tournament_20260915/runs"
ECBENCH_SESSIONS_GLOBS = ("research/ecbench_*/sessions/*", "research/isogeny_conductor_gap_*/sessions/*", "research/fiber_ic_scaling_*/sessions/*", "research/full_volcano_graphwide_ic_*/execution/campaigns/*/sessions/*")
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
        {"name": "prime-projected-abscissa", "params": "size,scan_cap", "what": "Cofactor-projected rational points, checked in the selected prime subgroup, with distinct negation pairs folded into columns."},
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


def _glob_match(path: str, pattern: str) -> bool:
    """Match one repository path with Python glob's path-component rules."""
    path_parts = path.split("/")
    pattern_parts = pattern.split("/")

    @lru_cache(maxsize=None)
    def match(path_i: int, pattern_i: int) -> bool:
        if pattern_i == len(pattern_parts):
            return path_i == len(path_parts)
        part = pattern_parts[pattern_i]
        if part == "**":
            return match(path_i, pattern_i + 1) or (
                path_i < len(path_parts)
                and not path_parts[path_i].startswith(".")
                and match(path_i + 1, pattern_i)
            )
        return (
            path_i < len(path_parts)
            and (
                not path_parts[path_i].startswith(".")
                or part.startswith(".")
            )
            and fnmatch.fnmatchcase(path_parts[path_i], part)
            and match(path_i + 1, pattern_i + 1)
        )

    return match(0, 0)


class RepositoryInputs:
    """Read a sparse worktree as the worktree overlaid on committed HEAD.

    Git marks sparse-excluded index entries with ``skip-worktree``.  Those
    files are absent from the filesystem, but remain inputs to this generated
    index.  Materialized files always use their worktree bytes, including
    local modifications.  Only an absent skip-worktree entry falls back to
    its blob in HEAD; an ordinary tracked deletion therefore stays deleted.
    """

    def __init__(self, root: str):
        self.root = root
        self.git_env = os.environ.copy()
        # A generator must not silently turn a sparse validation into a broad
        # promisor download.  Callers may provide the required objects through
        # the repository object store or an alternate object directory.
        self.git_env["GIT_NO_LAZY_FETCH"] = "1"
        try:
            raw = subprocess.check_output(
                ["git", "ls-files", "-t", "-z"],
                cwd=root,
                env=self.git_env,
                stderr=subprocess.PIPE,
            )
        except (FileNotFoundError, subprocess.CalledProcessError) as error:
            raise RuntimeError(
                f"cannot enumerate committed lab-browser inputs in Git checkout {root}"
            ) from error
        entries = [item for item in raw.split(b"\0") if item]
        self.tracked_files = {
            item[2:].decode("utf-8", "surrogateescape") for item in entries
        }
        self.sparse_files = {
            item[2:].decode("utf-8", "surrogateescape")
            for item in entries
            if item.startswith(b"S ")
        }

    def _relative(self, path: str) -> str:
        if os.path.isabs(path):
            path = os.path.relpath(path, self.root)
        return path.replace(os.sep, "/")

    def _absolute(self, path: str) -> str:
        return os.path.join(self.root, self._relative(path))

    def exists(self, path: str) -> bool:
        relative = self._relative(path)
        return (
            relative in self.tracked_files
            and (
                os.path.isfile(self._absolute(relative))
                or relative in self.sparse_files
            )
        )

    def is_directory(self, path: str) -> bool:
        relative = self._relative(path).rstrip("/")
        prefix = relative + "/"
        return any(
            item.startswith(prefix) and self.exists(item)
            for item in self.tracked_files
        )

    def glob_files(self, pattern: str) -> list[str]:
        pattern = self._relative(pattern)
        materialized = {
            self._relative(path)
            for path in glob.glob(os.path.join(self.root, pattern), recursive=True)
            if (
                os.path.isfile(path)
                and self._relative(path) in self.tracked_files
            )
        }
        materialized.update(path for path in self.sparse_files if _glob_match(path, pattern))
        return sorted(materialized)

    def child_directories(self, base: str) -> list[str]:
        base = self._relative(base).rstrip("/")
        children = set()
        prefix = base + "/"
        for path in self.tracked_files:
            if path.startswith(prefix) and self.exists(path):
                remainder = path[len(prefix):]
                if "/" in remainder:
                    children.add(remainder.split("/", 1)[0])
        return sorted(children)

    def size(self, path: str) -> int:
        relative = self._relative(path)
        absolute = os.path.join(self.root, relative)
        if os.path.isfile(absolute):
            return os.path.getsize(absolute)
        if relative not in self.sparse_files:
            raise FileNotFoundError(absolute)
        result = subprocess.run(
            ["git", "cat-file", "-s", f"HEAD:{relative}"],
            cwd=self.root,
            env=self.git_env,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        if result.returncode:
            detail = result.stderr.strip() or "object is unavailable"
            raise OSError(f"cannot read sparse input HEAD:{relative}: {detail}")
        return int(result.stdout.strip())

    @contextlib.contextmanager
    def open_binary(self, path: str):
        relative = self._relative(path)
        absolute = os.path.join(self.root, relative)
        if os.path.isfile(absolute):
            with open(absolute, "rb") as fh:
                yield fh
            return
        if relative not in self.sparse_files:
            raise FileNotFoundError(absolute)
        proc = subprocess.Popen(
            ["git", "cat-file", "blob", f"HEAD:{relative}"],
            cwd=self.root,
            env=self.git_env,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        assert proc.stdout is not None
        completed = False
        try:
            yield proc.stdout
            completed = True
        finally:
            proc.stdout.close()
            if not completed and proc.poll() is None:
                proc.kill()
            returncode = proc.wait()
            assert proc.stderr is not None
            stderr = proc.stderr.read()
            proc.stderr.close()
        if returncode:
            detail = (
                stderr.decode("utf-8", "replace").strip()
                or "object is unavailable"
            )
            raise OSError(f"cannot read sparse input HEAD:{relative}: {detail}")

    def read_bytes(self, path: str) -> bytes:
        with self.open_binary(path) as fh:
            return fh.read()

    def read_text(self, path: str, errors: str = "strict") -> str:
        return self.read_bytes(path).decode("utf-8", errors)


def load_file(path: str, source: RepositoryInputs | None = None):
    if source is not None:
        return json.loads(source.read_text(path))
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)


def load(path_rel: str, source: RepositoryInputs | None = None):
    return load_file(os.path.join(ROOT, path_rel), source)


def sha256_file(path_rel: str, source: RepositoryInputs | None = None) -> str:
    h = hashlib.sha256()
    opener = (
        source.open_binary(path_rel)
        if source is not None
        else open(os.path.join(ROOT, path_rel), "rb")
    )
    with opener as fh:
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
        elif family == "extension":
            # GF(p^k) = GF(p)[t]/(t^k + low(t)): bits are the field's, p^k's.
            p, k = params.get("p"), int_or_none(params.get("k"))
            field = {"characteristic": "p", "p": str(p), "bits": ((int_or_none(p) or 0) ** (k or 0)).bit_length(), "degree": k, "modulus": params.get("modulus")}
        else:
            degree = params.get("n") or params.get("m")
            field = {"characteristic": 2, "p": None, "bits": degree, "degree": degree, "modulus": params.get("modulus")}
        try:
            model = json.loads(c.get("model_json") or "{}")
        except ValueError:
            model = {}
        a = params.get("a") if params.get("a") is not None else model.get("a")
        b = params.get("b") if params.get("b") is not None else model.get("b")
        # An extension coefficient is its list over 1, t, ..., t^(k-1).
        a, b = ("[" + ", ".join(v) + "]" if isinstance(v, list) else v for v in (a, b))
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


def ecbench_sessions(
    curves_by_slug: dict, source: RepositoryInputs
) -> tuple[list[dict], list[dict], list[dict], list[dict]]:
    sessions, methods, fbs, yields = [], {}, {}, []
    per_curve_arm = defaultdict(lambda: {"n": 0, "verified": 0, "s_sum": 0.0, "floor_sum": 0.0})
    session_paths = {
        path
        for pattern in ECBENCH_SESSIONS_GLOBS
        for path in source.glob_files(pattern.rstrip("/") + "/session.json")
    }
    for session_path in sorted(session_paths):
        sdir = os.path.dirname(session_path)
        records_path = os.path.join(sdir, "records.jsonl")
        if not source.exists(records_path):
            continue
        s = load_file(session_path, source)
        spec_path = os.path.join(sdir, "spec.json")
        spec = load_file(spec_path, source) if source.exists(spec_path) else {}
        slugs, arms, statuses = set(), {}, Counter()
        with source.open_binary(records_path) as raw:
            lines = (line.decode("utf-8") for line in raw)
            for line in lines:
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
                        {**{k: fb.get(k) for k in ("fb_id", "fb_sha256", "family", "params", "description", "signed_points", "abscissae", "columns", "dimension", "points_sha256")}, "curve": slug if slug in curves_by_slug else None, "sessions": set(), "methods": set()},
                    )
                    if slug and slug not in curves_by_slug:
                        # A frozen record can name a model whose exact subgroup
                        # representation is still pending registration. Preserve
                        # that identity without emitting a broken registry link.
                        f["recorded_curve"] = slug
                        f["registration_status"] = "pending"
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
                    # A zero charged lower bound with unpriced native work is
                    # not a complete zero-cost solve or a sortable S=0 result.
                    if r["cost"].get("unpriced") and r["cost"].get("total_gae") == 0:
                        agg["unpriced_zero_runs"] += 1
        sessions.append(
            {
                "session_id": s["session_id"],
                "dir": sdir,
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
                    "mean_s": round(agg["s_sum"] / agg["verified"], 4) if agg["verified"] and not agg["unpriced_zero_runs"] else None,
                    "mean_ratio_to_floor": round(agg["s_sum"] / agg["floor_sum"], 3) if agg["floor_sum"] and not agg["unpriced_zero_runs"] else None,
                    **({"unpriced_zero_runs": agg["unpriced_zero_runs"], "score_status": "unknown complete price; zero charged lower bound"} if agg["unpriced_zero_runs"] else {}),
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
    cost = r.get("cost") or {}
    uninformative_zero = bool(cost.get("unpriced")) and cost.get("total_gae") == 0
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
        "total_gae": None if uninformative_zero else cost.get("total_gae"),
        "s": None if uninformative_zero else cost.get("s"),
        **({"charged_gae_lower_bound": 0, "charged_s_lower_bound": 0, "score_status": "unknown complete price; zero charged lower bound", "unpriced": cost["unpriced"]} if uninformative_zero else {}),
        "deterministic": (r.get("cost") or {}).get("deterministic"),
        "solver": r.get("solver"),
    }


def tournament_rounds(curves: list[dict], source: RepositoryInputs) -> list[dict]:
    koblitz = {(c["field"]["degree"], int(c["a"])): c["slug"] for c in curves if c["family"] == "koblitz" and c["a"] in ("0", "1")}
    by_slug = {c["slug"]: c for c in curves}
    rounds = []
    base = TOURNAMENT_RUNS
    for d in source.child_directories(base):
        path = os.path.join(base, d)
        row = {"round": d, "dir": path}
        failure = os.path.join(path, "prepare_failure.json")
        report_path = os.path.join(path, "REPORT.md")
        row["report"] = report_path if source.exists(report_path) else None
        if source.exists(failure):
            f = load_file(failure, source)
            row.update({"status": "prepare failed", "reason": f.get("reason")})
            rounds.append(row)
            continue
        dec_path = os.path.join(path, "decision.json")
        if source.exists(dec_path):
            dec = load_file(dec_path, source)
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
        if source.exists(cand_path):
            row["candidates"] = [
                {
                    "id": c.get("id"),
                    "config": c.get("config"),
                    "parent": c.get("parent"),
                    "hypothesis": c.get("hypothesis"),
                    "configuration_sha256": c.get("configuration_sha256"),
                }
                for c in load_file(cand_path, source)
            ]
        fix_path = os.path.join(path, "fixtures.json")
        cells = {}
        if source.exists(fix_path):
            fixtures = load_file(fix_path, source)
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


def candidate_identities(source: RepositoryInputs) -> list[dict]:
    found = defaultdict(lambda: {"count": 0, "files": Counter()})
    for pattern in IC1_SCAN_GLOBS:
        for path in source.glob_files(pattern):
            if source.size(path) > MAX_SCAN_BYTES:
                continue
            text = source.read_text(path, errors="replace")
            for ident in IC1_FULL_RE.findall(text):
                found[ident]["count"] += 1
                found[ident]["files"][path] += 1
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


def attach_covers(
    curves: list[dict],
    report: dict,
    source: RepositoryInputs | None = None,
) -> None:
    """Join preverified metadata only; mathematical checks are native Rust."""
    if (
        report.get("schema_version") != "curve-covers/v1"
        or report.get("registry_sha256") != sha256_file(REGISTRY, source)
    ):
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
    source = RepositoryInputs(ROOT)
    registry = load(REGISTRY, source)
    leaderboard = load(LEADERBOARD, source)
    curves = curve_rows(registry, leaderboard)
    attach_covers(curves, load(COVERS, source), source)
    graph = load(COVER_LINKS, source)  # JSON syntax is a subset of YAML 1.2.
    if (
        graph.get("schema_version") != "curve-cover-graph/v1"
        or graph.get("registry_sha256") != sha256_file(REGISTRY, source)
    ):
        raise ValueError("stale cover graph")
    nodes = {node["icv1_identity"]["slug"]: node for node in graph["curves"]}
    if len(nodes) != len(graph["curves"]) or set(nodes) != {c["slug"] for c in curves}:
        raise ValueError("cover graph must bind every catalog model exactly once")
    for curve in curves:
        node = nodes[curve["slug"]]
        if node["model_json"] != curve["model_json"]:
            raise ValueError("cover graph model mismatch")
        curve["cover_links"] = [dict(link, cover_id=graph["covers"][link["cover_uid"]]["cover_id"],
                                     map_id=graph["maps"][link["map_uid"]]["map_id"]) for link in node["cover_links"]]
        curve["standards_provenance"] = node.get("standards_provenance") or []
    by_slug = {c["slug"]: c for c in curves}
    sessions, methods, fbs, yields = ecbench_sessions(by_slug, source)
    rounds = tournament_rounds(curves, source)
    candidates = candidate_identities(source)
    tag_to_slug = {}
    for c in curves:
        if c["family"] == "koblitz" and c["a"] in ("0", "1"):
            tag_to_slug[(c["field"]["degree"], f"kb{c['a']}")] = c["slug"]
    for cand in candidates:
        cand["slug"] = tag_to_slug.get((cand.get("n"), cand.get("curve_tag")))
    sources = {
        REGISTRY: sha256_file(REGISTRY, source),
        COVERS: sha256_file(COVERS, source),
        COVER_LINKS: sha256_file(COVER_LINKS, source),
        STANDARDS_COVERAGE: sha256_file(STANDARDS_COVERAGE, source),
        LEADERBOARD: sha256_file(LEADERBOARD, source),
        **{
            os.path.join(s["dir"], "records.jsonl"): sha256_file(
                os.path.join(s["dir"], "records.jsonl"), source
            )
            for s in sessions
        },
    }
    return {
        "schema_version": "lab-browser/1",
        "generated_by": "scripts/build_lab_browser.py",
        "what_this_is": "A searchable index of the curves, methods, factor bases, candidate identities, tournament rounds and ecbench sessions this repository names, joined on the identities it already computes (ICV1 slugs, EC1 aliases, ECM1 method ids, FB1 factor-base ids, IC1 candidate ids). Every figure is quoted from a committed file named in `sources` or in the row; nothing is computed here beyond means over a session's own verified runs.",
        "counts": {"curves": len(curves), "methods": len(methods), "factor_bases": len(fbs), "candidates": len(candidates), "rounds": len(rounds), "sessions": len(sessions), "yields": len(yields)},
        "sources": sources,
        "vocabulary": VOCABULARY,
        "standards_coverage": load(STANDARDS_COVERAGE, source),
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
