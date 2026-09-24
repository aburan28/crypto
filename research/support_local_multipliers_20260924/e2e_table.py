#!/usr/bin/env python3
"""Tabulate the whole-logarithm runs of RESEARCH_SUPPORT_LOCAL_MULTIPLIERS.md §4 and §6.

    python3 research/support_local_multipliers_20260924/e2e_table.py

Two sources, never merged:
  e2e/{reference,candidate}/                 the registered runs (binary of 1f0d9751,
                                             arms back to back, one run per seed);
  e2e_postfix/<binary>-<arm>/rep{1,2,3}/     the amendment's runs (both binaries,
                                             four arms interleaved, three repetitions).
Per cell and comparison: runs, verified runs, trials and relations (which must
be equal run for run), oracle word operations and their per-run ratio, and the
wall ratio (`elapsed_seconds`, the median over repetitions per seed) as a
geometric mean over seeds with a 95% paired bootstrap interval (10,000
resamples of seeds, fixed RNG seed).  Wall is a practicality note.
"""
import json, math, pathlib, random, statistics

HERE = pathlib.Path(__file__).resolve().parent
CELLS = {
    "`K_0/2^13`, seeds 201–210, holdout 301–305": ("K0_2^13", list(range(201, 211)) + list(range(301, 306))),
    "`K_0/2^9`, seeds 201–205": ("K0_2^9", list(range(201, 206))),
}


def load(paths):
    docs = [json.loads(p.read_text()) for p in paths]
    for d in docs:
        assert d["result"]["verified"] and d["result"]["expected"] == d["result"]["recovered"], d["arguments"]
    return docs


def per_seed(root, cell, seeds):
    """seed -> (trials, relations, word ops, median wall) for one arm directory."""
    out = {}
    for s in seeds:
        paths = sorted(root.glob(f"**/{cell}_seed{s}.json"))
        docs = load(paths)
        key = {(d["counts"]["trials"], d["counts"]["relations"], d["counts"]["f4_word_ops"]) for d in docs}
        assert len(key) == 1, (root, cell, s, key)
        t, r, w = key.pop()
        out[s] = (t, r, w, statistics.median(d["elapsed_seconds"] for d in docs), len(docs))
    return out


def geo_ci(ratios, rng):
    logs = [math.log(x) for x in ratios]
    boots = sorted(
        math.exp(statistics.fmean(rng.choice(logs) for _ in logs)) for _ in range(10_000)
    )
    return math.exp(statistics.fmean(logs)), boots[249], boots[9_749]


def row(label, a_dir, b_dir):
    rng = random.Random(20260924)
    lines = []
    for name, (cell, seeds) in CELLS.items():
        a, b = per_seed(a_dir, cell, seeds), per_seed(b_dir, cell, seeds)
        runs = sum(v[4] for v in a.values()), sum(v[4] for v in b.values())
        same = all(a[s][:2] == b[s][:2] for s in seeds)
        trials = sum(v[0] for v in a.values()), sum(v[0] for v in b.values())
        rels = sum(v[1] for v in a.values()), sum(v[1] for v in b.values())
        wops = sum(v[2] for v in a.values()), sum(v[2] for v in b.values())
        per_run = [a[s][2] / b[s][2] for s in seeds]
        wall = sum(v[3] for v in a.values()), sum(v[3] for v in b.values())
        g, lo, hi = geo_ci([a[s][3] / b[s][3] for s in seeds], rng)
        lines.append(
            f"| {label} | {name} | {runs[0]} + {runs[1]} | {trials[0]} / {trials[1]}"
            f"{' (every seed equal)' if same else ' (seeds differ)'} | {rels[0]} / {rels[1]} "
            f"| {wops[0]:,} → {wops[1]:,} | **{wops[0] / wops[1]:.2f}×** ({min(per_run):.2f}–{max(per_run):.2f}) "
            f"| {wall[0]:.2f} → {wall[1]:.2f} s | {g:.2f}× [{lo:.2f}, {hi:.2f}] |"
        )
    return lines


def main():
    print(
        "| comparison (A → B) | cell | runs A + B | trials A / B | relations A / B | oracle word ops A → B "
        "| ratio (per run) | wall A → B (sum of per-seed medians) | wall ratio A/B, geometric mean [95% paired bootstrap] |"
    )
    print("|:--|:--|--:|:--|:--|:--|:--|:--|:--|")
    e2e, post = HERE / "e2e", HERE / "e2e_postfix"
    for line in row("registered, back to back: reference → candidate", e2e / "reference", e2e / "candidate"):
        print(line)
    if post.exists():
        # The amendment's admissibility check for whole logarithms: every
        # amendment run of an arm has the registered run's trials, relations
        # and oracle word operations, seed by seed.
        for arm in ("reference", "candidate"):
            for cell, seeds in CELLS.values():
                registered = per_seed(e2e / arm, cell, seeds)
                for binary in ("registered", "fixed"):
                    rerun = per_seed(post / f"{binary}-{arm}", cell, seeds)
                    assert all(rerun[s][:3] == registered[s][:3] for s in seeds), (binary, arm, cell)
        for label, a, b in [
            ("registered binary, interleaved: reference → candidate", "registered-reference", "registered-candidate"),
            ("fixed binary, interleaved: reference → candidate", "fixed-reference", "fixed-candidate"),
            ("hash fix alone, reference arm: registered → fixed", "registered-reference", "fixed-reference"),
            ("hash fix alone, candidate arm: registered → fixed", "registered-candidate", "fixed-candidate"),
            ("round 2 as shipped: registered reference → fixed candidate", "registered-reference", "fixed-candidate"),
        ]:
            for line in row(label, post / a, post / b):
                print(line)


if __name__ == "__main__":
    main()
