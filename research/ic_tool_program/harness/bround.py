#!/usr/bin/env python3
"""Track B's measurements (IC_TOOL_PROGRAM.md §9): one runner for B0, B1,
B2 and B3, whose protocols declare the same steps.

    IC_BASE=<base> IC_BASE_COMMIT=<sha> IC_CAND=<candidate> IC_CAND_COMMIT=<sha> \\
    IC_RUNS=<run tree> python3 bround.py --steps B0,B1,B3,B2 <step> [<step> ...]

The steps:
- `manifest`: the host and both binaries.
- `conformance`: `conformance/run.py --steps <steps>` on the candidate,
  with `--build-commit`, and on the base, whose failures record what the
  base could not do.
- `pin`: the candidate on all 90 suite rows from their v1 files,
  untimed. Every output must equal v0's from R01's profile pass.
- `translate` (B1 on): `ic check --translate` on each row must equal the
  design §8 translation written in Python (`conformance/v2/make_cases.py`,
  C009's), and `ic price` on it must give the row's v1 outputs.
- `timing`: the base against the candidate on `M1`'s 22 rows, five
  rounds ABAB, isolated (the programme's runner).
- `v2timing` (B1's measurement 6): the candidate on `M1`'s 22 rows, the
  v1 file against its v2 translation, five rounds ABAB.
- `chain` (the B steps' amendment of 2026-10-01): the timing check of
  every step at once. One interleave over the newest baseline and each
  step's arm in the queue's order (`IC_CHAIN`, a JSON object of arm name
  to binary, and `IC_CHAIN_COMMITS`, of arm name to commit), on `M1`'s 22
  rows, five rounds, the order reversed every other round. Each step's
  figure is its paired ratio against the arm before it in the chain.
- `analyse`: the verdicts, as JSON on stdout.

Every output stays on disk; every step resumes where it stopped.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import subprocess
import sys
import tarfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import bench  # noqa: E402
import stats  # noqa: E402

PROGRAMME = HERE.parent
ROOT = bench.ROOT
R01_DIR = PROGRAMME / "rounds" / "R01-baseline-v0"
ROUNDS = 5
CHAIN_ORDER = ("base", "B0", "B1", "B3", "B2", "B2b", "B7a", "B3b")


def arm(name: str) -> Path:
    path = os.environ.get(name)
    if not path or not Path(path).exists():
        raise SystemExit(f"set {name}")
    return Path(path)


def runs() -> Path:
    path = os.environ.get("IC_RUNS")
    if not path:
        raise SystemExit("set IC_RUNS: the round's run tree")
    return Path(path).resolve()


def r01_runs() -> Path:
    """R01's run tree, checked against its SHA-256 and extracted on first use."""
    tree = R01_DIR / "runs"
    if not tree.exists():
        want = (R01_DIR / "runs.tar.xz.sha256").read_text().split()[0]
        if bench.sha256(R01_DIR / "runs.tar.xz") != want:
            raise SystemExit("R01's runs.tar.xz does not match runs.tar.xz.sha256")
        with tarfile.open(R01_DIR / "runs.tar.xz") as tar:
            tar.extractall(R01_DIR, filter="data")
    return tree


def all_rows() -> list[dict]:
    return bench.slug_rows(bench.suite_rows("S") + bench.suite_rows("smoke"))


def m1() -> list[dict]:
    return [r for r in bench.slug_rows(bench.suite_rows("S")) if r["recipe_seed"] == 201]


def v0_outputs(row: dict) -> dict:
    r01 = r01_runs()
    if row["tier"] == "S":
        old = bench.load(bench.figure_path(r01 / "profile" / "v0" / row["suite_id"] / "r1.price.json"))
    else:
        old = bench.load(r01 / "smoke" / f"{row['suite_id']}.price.json")
    return bench.outputs(old)


# ── Steps ──────────────────────────────────────────────────────────


def manifest(_steps: str) -> dict:
    binaries = {"base": ("IC_BASE", os.environ.get("IC_BASE_COMMIT")),
                "candidate": ("IC_CAND", os.environ.get("IC_CAND_COMMIT"))}
    doc = bench.host_manifest(runs() / "host.json", {
        name: {"path_basename": arm(env).name, "sha256": bench.sha256(arm(env)), "built_from": commit}
        for name, (env, commit) in binaries.items()})
    r01 = json.loads((r01_runs() / "host.json").read_text())
    same = all(doc[k] == r01[k] for k in ("cpu_model", "cpu_flags_relevant", "logical_cores", "memory",
                                          "transparent_hugepage", "os"))
    out = {"host_matches_r01": same, "aa": "R01's" if same else "the round's own (not run by this harness)"}
    (runs() / "aa-source.json").write_text(json.dumps(out, indent=1) + "\n")
    return out


def conformance(steps: str) -> dict:
    out = {}
    for name, env, commit in (("candidate", "IC_CAND", os.environ.get("IC_CAND_COMMIT")),
                              ("base", "IC_BASE", None)):
        report = runs() / "conformance" / f"{name}.json"
        if not report.exists():
            report.parent.mkdir(parents=True, exist_ok=True)
            cmd = [sys.executable, str(PROGRAMME / "conformance" / "run.py"), "--ic", str(arm(env)),
                   "--steps", steps, "--out", str(report)]
            if commit:
                cmd += ["--build-commit", commit]
            subprocess.run(cmd, check=False, stdout=subprocess.DEVNULL)
        doc = json.loads(report.read_text())
        results = doc.get("results", [])
        out[name] = {"passed": doc.get("passed"), "cases": len(results),
                     "failed": [c["id"] for c in results if not c["pass"]]}
    return out


def pin(_steps: str) -> dict:
    out = runs() / "pin" / "pin.json"
    if out.exists():
        return json.loads(out.read_text())
    result = []
    for row in all_rows():
        new = bench.untimed(arm("IC_CAND"), row, runs() / "pin" / "candidate" / row["id"] / "pin.price.json")
        a, b = bench.outputs(new), v0_outputs(row)
        result.append({"row": row["id"], "suite_id": row["suite_id"], "equal": a == b,
                       "differs_in": sorted(k for k in a if a[k] != b[k])})
    doc = {"rows": result, "held": all(e["equal"] for e in result), "against": "v0's outputs, R01's profile pass"}
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(doc, indent=1) + "\n")
    return doc


def python_translation(row: dict) -> dict:
    """The design §8 translation of a suite row, as C009's generator writes it."""
    sys.path.insert(0, str(PROGRAMME / "conformance" / "v2"))
    import make_cases  # noqa: E402

    v1 = make_cases.v1_row(row["params"])
    a, n = row["a"], row["n"]
    f = make_cases.curve_id.find_irreducible_sparse(n)
    order = make_cases.koblitz_order(a, n)
    r = make_cases.largest_prime_factor(order)
    target = v1["targets"][0]
    return make_cases.document(v1["name"], make_cases.Curve(n, f, a, 1), {"form": "koblitz", "a": a}, r,
                               order // r, {"rule": "koblitz_search_v1"},
                               {k: target[k] for k in ("public_hash_seed", "known_log", "random_seed") if k in target},
                               make_cases.recipe_of(v1), row["rho_seed"])


def translation_path(row: dict) -> Path:
    return runs() / "translate" / row["id"] / "v2.json"


def translate_row(row: dict) -> dict:
    """`ic check --translate` on the row's v1 file, written beside the run."""
    path = translation_path(row)
    if not path.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
        got = subprocess.run([str(arm("IC_CAND")), "check", "--params", str(bench.SUITE / row["params"]),
                              "--translate", "--rho-seed", str(row["rho_seed"]), "--json"],
                             capture_output=True, text=True, check=False, env=bench.ENV)
        docs = json.loads(got.stdout).get("documents", [])
        path.write_text(json.dumps(docs[0] if len(docs) == 1 else docs, indent=1) + "\n")
    return json.loads(path.read_text())


def translate(_steps: str) -> dict:
    out = runs() / "translate" / "translate.json"
    if out.exists():
        return json.loads(out.read_text())
    result = []
    for row in all_rows():
        doc = translate_row(row)
        same_json = doc == python_translation(row)
        report = runs() / "translate" / row["id"] / "v2.price.json"
        if not report.exists():
            subprocess.run(["taskset", "-c", bench.CPUS, str(arm("IC_CAND")), "price", "--params",
                            str(translation_path(row)), "--json", "--out", str(report)],
                           env=bench.ENV, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=False)
        a, b = bench.outputs(bench.load(report)), v0_outputs(row)
        result.append({"row": row["id"], "same_json_as_python": same_json, "same_outputs": a == b,
                       "differs_in": sorted(k for k in a if a[k] != b[k])})
    doc = {"rows": result, "held": all(e["same_json_as_python"] and e["same_outputs"] for e in result)}
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(doc, indent=1) + "\n")
    return doc


def timing(_steps: str) -> dict:
    bench.interleave({"base": arm("IC_BASE"), "cand": arm("IC_CAND")}, m1(), ROUNDS, runs() / "timing")
    return {"done": True}


def chain_arms() -> dict[str, Path]:
    """The chain's arms in the queue's order: the baseline, then each step
    present.  A rejected step is left out, and the step after it is then
    measured against the arm before it."""
    given = json.loads(os.environ.get("IC_CHAIN") or "{}")
    if "base" not in given or set(given) - set(CHAIN_ORDER):
        raise SystemExit(f"set IC_CHAIN to a JSON object of arm to binary, with base, from {CHAIN_ORDER}")
    arms = {name: Path(given[name]) for name in CHAIN_ORDER if name in given}
    for name, path in arms.items():
        if not path.exists():
            raise SystemExit(f"IC_CHAIN's {name} binary does not exist: {path}")
    return arms


def chain(_steps: str) -> dict:
    arms = chain_arms()
    commits = json.loads(os.environ.get("IC_CHAIN_COMMITS") or "{}")
    bench.host_manifest(runs() / "host.json", {
        name: {"path_basename": path.name, "sha256": bench.sha256(path), "built_from": commits.get(name)}
        for name, path in arms.items()})
    bench.interleave(arms, m1(), ROUNDS, runs() / "chain")
    return {"done": True, "arms": list(arms)}


def v2timing(_steps: str) -> dict:
    """The candidate on each M1 row, v1 file against its translation."""
    d = runs() / "v2timing"
    names = ("v1", "v2")
    for k in range(1, ROUNDS + 1):
        for row in m1():
            translate_row(row)
            for name in (names if k % 2 else names[::-1]):
                out = d / name / row["id"] / f"r{k}.price.json"
                if out.exists():
                    continue
                out.parent.mkdir(parents=True, exist_ok=True)
                cmd = (bench.price_cmd(arm("IC_CAND"), row, out) if name == "v1" else
                       [str(arm("IC_CAND")), "price", "--params", str(translation_path(row)), "--json",
                        "--out", str(out)])
                bench.launch(cmd, out, d)
                print(f"r{k} {row['id']} {name}: {bench.load(out).get('status')}", flush=True)
    return {"done": True}


# ── The analysis ───────────────────────────────────────────────────


def r01_aa() -> dict[str, dict]:
    doc = json.loads((R01_DIR / "analysis.json").read_text())
    return {bench.curve_slug(int(row["size"][1]), int(row["size"][3:])): row["cold"] for row in doc["aa"]}


def paired(d: Path, arms: tuple[str, str], rows: list[dict]) -> dict:
    ratios, missing = [], 0
    for r in rows:
        for k in range(1, ROUNDS + 1):
            pa = bench.figure_path(d / arms[0] / r["id"] / f"r{k}.price.json")
            pb = bench.figure_path(d / arms[1] / r["id"] / f"r{k}.price.json")
            a, b = bench.load(pa), bench.load(pb)
            if not (bench.clean(pa) and bench.clean(pb) and a.get("status") == b.get("status") == "complete"):
                missing += 1
                continue
            ratios.append(stats.ic_cold_ns(a) / stats.ic_cold_ns(b))
    out = stats.geo_ci(ratios)
    out["missing_pairs"] = missing
    return out


def by_size(rows: list[dict]) -> dict[tuple[int, int], list[dict]]:
    out: dict[tuple[int, int], list[dict]] = {}
    for r in rows:
        out.setdefault((r["a"], r["n"]), []).append(r)
    return out


def accounting(d: Path) -> dict:
    counts = {"processes": 0, "contended": 0, "failed": 0}
    for rec in d.rglob("*.isolation.jsonl"):
        for line in rec.read_text().splitlines():
            run = json.loads(line)["run"]
            counts["processes"] += 1
            counts["contended"] += bool(run["contended"])
            counts["failed"] += run["exit_status"] != 0
    return counts


def sizes(d: Path, arms: tuple[str, str]) -> list[dict]:
    aa = r01_aa()
    out = []
    for (a, n), rs in sorted(by_size(m1()).items(), key=lambda kv: kv[1][0]["r"]):
        slug = bench.curve_slug(a, n)
        row = {"slug": slug, "log2_r": round(math.log2(rs[0]["r"]), 3), "cold": paired(d, arms, rs)}
        band = aa.get(slug)
        if band and "lo" in band and "hi" in row["cold"]:
            row["aa_cold"] = {k: band[k] for k in ("geomean", "lo", "hi")}
            row["regresses_beyond_aa"] = row["cold"]["hi"] < band["lo"]
        out.append(row)
    return out


def analyse(steps: str) -> dict:
    out: dict = {"steps": steps}
    for name, path in (("conformance", runs() / "conformance"), ("pin", runs() / "pin" / "pin.json"),
                       ("translate", runs() / "translate" / "translate.json")):
        if path.exists():
            out[name] = conformance(steps) if name == "conformance" else json.loads(path.read_text())
            if name != "conformance":
                out[name] = {"held": out[name]["held"],
                             "failing_rows": [e["row"] for e in out[name]["rows"]
                                              if not e.get("equal", e.get("same_outputs")) or
                                              not e.get("same_json_as_python", True)]}
    if (runs() / "timing").exists():
        rows = sizes(runs() / "timing", ("base", "cand"))
        out["timing"] = {"what": "paired cold-time ratio, base over candidate, per size; above 1 is faster",
                         "accounting": accounting(runs() / "timing"), "sizes": rows,
                         "any_regression": any(r.get("regresses_beyond_aa") for r in rows)}
    if (runs() / "chain").exists():
        names = [n for n in CHAIN_ORDER if (runs() / "chain" / n).exists()]
        steps_out = {}
        for prev, cur in zip(names, names[1:]):
            rows = sizes(runs() / "chain", (prev, cur))
            steps_out[cur] = {"against": prev, "sizes": rows,
                              "any_regression": any(r.get("regresses_beyond_aa") for r in rows)}
        out["chain"] = {"what": "per step, the paired cold-time ratio of the arm before it over the step's "
                                "arm, per size; above 1 is faster",
                        "arms": names, "accounting": accounting(runs() / "chain"), "steps": steps_out}
    if (runs() / "v2timing").exists():
        rows = sizes(runs() / "v2timing", ("v1", "v2"))
        out["v2timing"] = {"what": "paired cold-time ratio, v1 file over v2 translation, per size",
                           "accounting": accounting(runs() / "v2timing"), "sizes": rows,
                           "any_regression": any(r.get("regresses_beyond_aa") for r in rows)}
    return out


STEPS = {"manifest": manifest, "conformance": conformance, "pin": pin, "translate": translate,
         "timing": timing, "v2timing": v2timing, "chain": chain, "analyse": analyse}


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--steps", required=True, help="the conformance steps, e.g. B0,B1,B3,B2")
    p.add_argument("step", nargs="+", choices=list(STEPS))
    args = p.parse_args()
    runs().mkdir(parents=True, exist_ok=True)
    for step in args.step:
        result = STEPS[step](args.steps)
        print(json.dumps({step: result}, indent=1), flush=True)


if __name__ == "__main__":
    main()
