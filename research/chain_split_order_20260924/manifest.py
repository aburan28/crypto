#!/usr/bin/env python3
"""Write manifest.json for the registered runs of RESEARCH_CHAIN_SPLIT_ORDER.md.

    python3 research/chain_split_order_20260924/manifest.py

Records the commit the runs were built from, the blake3/sha256 of every
measured source and of both binaries, the host, and the commands; refuses to
overwrite an existing manifest.
"""
import hashlib, json, os, pathlib, platform, subprocess

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SOURCES = [
    "src/cryptanalysis/koblitz_groebner.rs",
    "src/cryptanalysis/inherited_f4.rs",
    "src/cryptanalysis/koblitz_index_calculus.rs",
    "src/cryptanalysis/pq_groebner_f2.rs",
    "examples/groebner_stage_bench.rs",
]
BINARIES = ["target/release/examples/groebner_stage_bench", "target/release/ic"]


def sha256(path):
    return hashlib.sha256((ROOT / path).read_bytes()).hexdigest()


def main():
    out = HERE / "manifest.json"
    if out.exists():
        raise SystemExit(f"{out} exists; never overwritten")
    git = lambda *a: subprocess.check_output(["git", *a], cwd=ROOT, text=True).strip()
    cpu = next(
        (l.split(":", 1)[1].strip() for l in open("/proc/cpuinfo") if l.startswith("model name")),
        platform.processor(),
    )
    doc = {
        "round": "chain-split-order-20260924",
        "note": "research/notes/ecc2k130/RESEARCH_CHAIN_SPLIT_ORDER.md",
        "scope": "decomposition-oracle Groebner stage (stage diagnostic) plus whole-logarithm runs; see the note",
        "git_head": git("rev-parse", "HEAD"),
        # The commits the measured binaries were built from.  No file under
        # src/ differs between them (`git diff --stat 2809b498 0fdb05d0 -- src`
        # is empty); the harness gained the chain-holdout-2 ladder.
        "built_from": {
            "run.sh (frozen, chain, chain-holdout)": "2809b498",
            "run_t1prime.sh, oracle_ladder.sh, e2e.sh": "0fdb05d0",
        },
        "src_identical_between_builds": git("diff", "--stat", "2809b498", "0fdb05d0", "--", "src") == "",
        "working_tree_clean_for_sources": git("status", "--porcelain", "--", *SOURCES) == "",
        "measured_sources_unchanged_since_build": git("diff", "--stat", "0fdb05d0", "HEAD", "--", *SOURCES[:4]) == "",
        "sources_sha256_at_head": {s: sha256(s) for s in SOURCES},
        "binaries_sha256": {b: sha256(b) for b in BINARIES if (ROOT / b).exists()},
        "host": {
            "cpu": cpu,
            "logical_cores": os.cpu_count(),
            "os": platform.platform(),
            "note": "shared cloud host, no CPU pinning; wall time is a practicality note, not the metric",
        },
        "arms": {
            "reference": {"KIC_CHAIN_ORDER": "layout", "KIC_LINEAR_ELIM": "0", "KIC_F4_DROP": "complete"},
            "candidate": {"KIC_CHAIN_ORDER": "interleaved", "KIC_LINEAR_ELIM": "1", "KIC_F4_DROP": "complete"},
            "factorial": "run.sh: O (order) x L (linear elimination) x D (degree-drop rule), every variable explicit",
        },
        "commands": [
            "cargo build --release --example groebner_stage_bench --bin ic",
            "research/chain_split_order_20260924/run.sh",
            "research/chain_split_order_20260924/run_t1prime.sh",
            "research/chain_split_order_20260924/compare_all.sh",
            "research/chain_split_order_20260924/oracle_ladder.sh",
            "research/chain_split_order_20260924/e2e.sh",
            "python3 research/chain_split_order_20260924/table.py > research/chain_split_order_20260924/tables.md",
        ],
    }
    out.write_text(json.dumps(doc, indent=2) + "\n")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
