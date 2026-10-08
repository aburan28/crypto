#!/usr/bin/env python3
"""R01's figures (PROTOCOL.md, "Figures"), from runs/ only.

    tar -xJf runs.tar.xz && python3 analyse.py > analysis.json
"""
from __future__ import annotations

import json
import math
import os
import statistics
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "harness"))
import bench  # noqa: E402
import stats  # noqa: E402

RUNS = (HERE / os.environ.get("IC_RUNS", "runs")).resolve()
ROUNDS = 5
S23_ANALYSIS = bench.ROOT / "research" / "ic_single_target_20260930" / "analysis.json"


def load(path: Path) -> dict | None:
    return json.loads(path.read_text()) if path.exists() else None


def sizes() -> list[tuple[str, list[dict]]]:
    out: dict[str, list[dict]] = {}
    for r in bench.suite_rows("S"):
        out.setdefault(f"k{r['a']}n{r['n']}", []).append(r)
    return list(out.items())


def figure(out: Path) -> tuple[dict | None, Path]:
    path = bench.figure_path(out)
    rep = load(path)
    ok = rep is not None and rep.get("status") == "complete" and bench.clean(path)
    return (rep if ok else None), path


def accounting(d: Path) -> dict:
    """Every attempt under d: clean, contended and failed processes."""
    counts = {"processes": 0, "contended": 0, "failed": 0}
    for rec in d.rglob("*.isolation.jsonl"):
        for line in rec.read_text().splitlines():
            run = json.loads(line)["run"]
            counts["processes"] += 1
            counts["contended"] += bool(run["contended"])
            counts["failed"] += run["exit_status"] != 0
    return counts


def profile_row(size: str, rows: list[dict]) -> dict:
    reps = []
    for r in rows:
        rep, _ = figure(RUNS / "profile" / "v0" / r["id"] / "r1.price.json")
        if rep is not None:
            reps.append((r, rep))
    if not reps:
        return {"size": size, "rows": 0}
    r0 = reps[0][0]
    sqrt_r = math.sqrt(r0["r"])
    unit = statistics.median(stats.median_rep(rep, "unit_ns") for _, rep in reps)
    setup = [stats.median_rep(rep, "setup_ns") for _, rep in reps]
    online = [stats.ic_online_ns(rep) for _, rep in reps]
    cold = [stats.ic_cold_ns(rep) for _, rep in reps]
    rho = [stats.rho_online_ns(rep) for _, rep in reps]
    phase_tot: dict[str, float] = {}
    for _, rep in reps:
        for k, v in stats.setup_phase_ns(rep).items():
            phase_tot[k] = phase_tot.get(k, 0.0) + v
    total = sum(phase_tot.values())
    online_tot: dict[str, float] = {}
    for _, rep in reps:
        for k, v in stats.online_phase_ns(rep).items():
            online_tot[k] = online_tot.get(k, 0.0) + v
    otot = sum(online_tot.values())
    per_summand = [stats.setup_phase_ns(rep)["collect"] / rep["counts"]["pass"]["collect"]["summands_scanned"]
                   for _, rep in reps]
    per_pair = [stats.setup_phase_ns(rep)["build"] / rep["counts"]["pass"]["build"]["stored_pairs"]
                for _, rep in reps]
    return {
        "size": size, "a": r0["a"], "n": r0["n"], "log2_r": r0["log2_r"], "rows": len(reps),
        "columns": r0["columns"], "descent_summands": r0["descent_summands"],
        "unit_ns_v0": unit,
        "setup_ms_median": statistics.median(setup) / 1e6,
        "online_ms_median": statistics.median(online) / 1e6,
        "cold_ms_median": statistics.median(cold) / 1e6,
        "rho_online_ms_median": statistics.median(rho) / 1e6,
        "s_cold_v0_unit": statistics.median(cold) / unit / sqrt_r,
        "s_setup_v0_unit": statistics.median(setup) / unit / sqrt_r,
        "s_online_v0_unit": statistics.median(online) / unit / sqrt_r,
        "s_rho_online_v0_unit": statistics.median(rho) / unit / sqrt_r,
        "setup_phase_shares": {k: v / total for k, v in sorted(phase_tot.items(), key=lambda kv: -kv[1])},
        "online_phase_shares": {k: v / otot for k, v in sorted(online_tot.items(), key=lambda kv: -kv[1])},
        "collect_ns_per_summand": statistics.median(per_summand),
        "collect_units_per_summand": statistics.median(per_summand) / unit,
        "build_ns_per_pair": statistics.median(per_pair),
        "build_units_per_pair": statistics.median(per_pair) / unit,
        "stored_pairs": [rep["counts"]["pass"]["build"]["stored_pairs"] for _, rep in reps],
        "summands_scanned": [rep["counts"]["pass"]["collect"]["summands_scanned"] for _, rep in reps],
        "all_verified": all(rep["all_verified"] and rep["ic_and_rho_agree"] for _, rep in reps),
        "max_rss_mib": max(bench.run_record(p)["max_rss_kib"] for p in
                           (figure(RUNS / "profile" / "v0" / r["id"] / "r1.price.json")[1] for r, _ in reps)) / 1024,
    }


def paired(d: Path, arms: tuple[str, str], rows: list[dict], measure) -> dict:
    ratios, missing = [], 0
    for r in rows:
        for k in range(1, ROUNDS + 1):
            a, _ = figure(d / arms[0] / r["id"] / f"r{k}.price.json")
            b, _ = figure(d / arms[1] / r["id"] / f"r{k}.price.json")
            if a is None or b is None:
                missing += 1
                continue
            ratios.append(measure(a) / measure(b))
    return {**stats.geo_ci(ratios), "missing_pairs": missing}


def aa_rows() -> list[dict]:
    out = []
    for size, rows in sizes():
        m1 = [r for r in rows if r["recipe_seed"] == 201]
        out.append({"size": size, "log2_r": m1[0]["log2_r"],
                    "cold": paired(RUNS / "aa", ("A", "A2"), m1, stats.ic_cold_ns),
                    "online": paired(RUNS / "aa", ("A", "A2"), m1, stats.ic_online_ns),
                    "rho_online": paired(RUNS / "aa", ("A", "A2"), m1, stats.rho_online_ns)})
    return out


def cpu_seconds(path: Path) -> tuple[float, float]:
    run = bench.run_record(path)
    return run["user_seconds"], run["system_seconds"]


def thp_rows() -> list[dict]:
    out = []
    for size, rows in sizes():
        m1 = [r for r in rows if r["recipe_seed"] == 201]
        if not (RUNS / "thp" / "4k" / m1[0]["id"]).exists():
            continue
        user = {"4k": [], "thp": []}
        system = {"4k": [], "thp": []}
        for r in m1:
            for k in range(1, ROUNDS + 1):
                for arm in ("4k", "thp"):
                    rep, path = figure(RUNS / "thp" / arm / r["id"] / f"r{k}.price.json")
                    if rep is not None:
                        u, s = cpu_seconds(path)
                        user[arm].append(u)
                        system[arm].append(s)
        out.append({
            "size": size, "log2_r": m1[0]["log2_r"],
            "cold_4k_over_thp": paired(RUNS / "thp", ("4k", "thp"), m1, stats.ic_cold_ns),
            "collect_4k_over_thp": paired(RUNS / "thp", ("4k", "thp"), m1,
                                          lambda rep: stats.setup_phase_ns(rep)["collect"]),
            "build_4k_over_thp": paired(RUNS / "thp", ("4k", "thp"), m1,
                                        lambda rep: stats.setup_phase_ns(rep)["build"]),
            "process_user_s_median": {k: statistics.median(v) for k, v in user.items() if v},
            "process_system_s_median": {k: statistics.median(v) for k, v in system.items() if v},
        })
    return out


def calib_rows() -> dict:
    out: dict = {}
    for mode in ("4k", "thp"):
        d = RUNS / "calib" / mode
        if not d.exists():
            continue
        lat = {}
        for f in sorted(d.glob("chase-*.json"), key=lambda p: int(p.stem.split("-")[1])):
            j = json.loads(f.read_text() or "{}")
            lat[f.stem.split("-")[1]] = j.get("ns_per_step_median")
        mlp = {}
        for f in sorted(d.glob("mlp-*.json")):
            j = json.loads(f.read_text() or "{}")
            mlp[f.stem.split("-k")[1]] = j.get("ns_per_step_median")
        out[mode] = {"latency_ns_by_bytes": lat, "ns_per_group_step_256MiB_by_k": mlp}
    return out


def callgrind_rows() -> dict:
    """The per-phase split of each callgrind run (`../../harness/callgrind_phases.py`).

    The declared step's own `*.annotate.txt` read only the base output file,
    which holds none of the phase parts; the split re-reads every part.
    """
    out = {}
    d = RUNS / "callgrind"
    for f in sorted(d.glob("*.phases.json")) if d.exists() else []:
        out[f.name.split(".")[0]] = json.loads(f.read_text())
    return out


def main() -> None:
    profile = [profile_row(s, rows) for s, rows in sizes()]
    s23 = load(S23_ANALYSIS)
    doc = {
        "what_this_is": "R01, baseline v0: every figure the round's README and the ledger quote",
        "host": load(RUNS / "host.json"),
        "pin": load(RUNS / "pin" / "pin.json"),
        "smoke": load(RUNS / "smoke" / "smoke.json"),
        "accounting": {step: accounting(RUNS / step) for step in ("profile", "aa", "thp", "calib")
                       if (RUNS / step).exists()},
        "profile": profile,
        "aa": aa_rows(),
        "thp_probe": thp_rows(),
        "calibration": calib_rows(),
        "callgrind": callgrind_rows(),
        "rule_comparison_from_s23": [
            {k: s[k] for k in ("a", "n", "log2_r", "s_ic_online_mean", "s_rho_online_mean", "s_setup")}
            | {"online_speedup": s["online_speedup"]["mean_ratio"],
               "cold_ratio_ic_over_rho": s["cold_ratio_ic_over_rho"]["value"]}
            for s in (s23 or {}).get("sizes", [])],
    }
    print(json.dumps(doc, indent=1))


if __name__ == "__main__":
    main()
