#!/usr/bin/env python3
"""Write manifest.json for the runs of RESEARCH_SUPPORT_LOCAL_MULTIPLIERS.md.

    OLD=<target/release built from 1f0d9751> python3 research/support_local_multipliers_20260924/manifest.py

Records the commits the runs were built from, the sha256 of every measured
source and of the binaries, whether the amendment's reruns reproduce every
registered counter, the host, and the commands; refuses to overwrite an
existing manifest.
"""
import hashlib, json, os, pathlib, platform, subprocess

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SOURCES = [
    "src/cryptanalysis/koblitz_groebner.rs",
    "src/cryptanalysis/inherited_f4.rs",
    "src/cryptanalysis/fx_hash.rs",
    "src/cryptanalysis/koblitz_index_calculus.rs",
    "src/cryptanalysis/pq_groebner_f2.rs",
    "examples/groebner_stage_bench.rs",
    "examples/chain_ladder_screen.rs",
]
BINARIES = ["examples/groebner_stage_bench", "ic"]
REGISTERED, FIXED = "1f0d9751", "d69936e1"


def sha256(path):
    return hashlib.sha256(pathlib.Path(path).read_bytes()).hexdigest()


def main():
    out = HERE / "manifest.json"
    if out.exists():
        raise SystemExit(f"{out} exists; never overwritten")
    git = lambda *a: subprocess.check_output(["git", *a], cwd=ROOT, text=True).strip()
    cpu = next(
        (l.split(":", 1)[1].strip() for l in open("/proc/cpuinfo") if l.startswith("model name")),
        platform.processor(),
    )
    old = os.environ.get("OLD")
    identity = subprocess.run(["python3", str(HERE / "check_identity.py")], capture_output=True, text=True)
    e2e = subprocess.run(["python3", str(HERE / "e2e_table.py")], capture_output=True, text=True)
    doc = {
        "round": "support-local-multipliers-20260924",
        "note": "research/notes/ecc2k130/RESEARCH_SUPPORT_LOCAL_MULTIPLIERS.md",
        "scope": "decomposition-oracle Groebner stage (stage diagnostic) plus whole-logarithm runs; see the note",
        "git_head": git("rev-parse", "HEAD"),
        "built_from": {
            "run.sh, oracle_ladder.sh, e2e.sh (the registered runs)": REGISTERED,
            "run_postfix.sh 'registered-*' arms (the same commit, rebuilt in a worktree)": REGISTERED,
            "run_postfix.sh 'fixed-*' arms (the monomial-mask hash fix)": FIXED,
        },
        "src_identical_between_builds": git("diff", "--stat", REGISTERED, FIXED, "--", "src") == "",
        "src_files_differing_between_builds": git("diff", "--name-only", REGISTERED, FIXED, "--", "src").split(),
        "registered_binary_rebuild_note": (
            "the worktree rebuild of 1f0d9751 embeds a different path and so hashes differently from the binary "
            "the registered runs recorded (software.binary_blake3 in e2e/); check_identity.py and e2e_table.py "
            "show it reproduces every registered counter"
        ),
        "amendment_counters_identical_to_registered": {
            "stage_ladders (check_identity.py, every field but *_ns, 5 suites x 4 arms x 3 reps)": identity.returncode == 0,
            "whole_logarithms (e2e_table.py: trials, relations, oracle word ops, every seed, both binaries)": e2e.returncode == 0,
        },
        "measured_sources_unchanged_since_fixed_build": git("diff", "--stat", FIXED, "HEAD", "--", *SOURCES[:5]) == "",
        "sources_sha256_at_head": {s: sha256(ROOT / s) for s in SOURCES},
        "binaries_sha256": {
            "fixed (target/release)": {b: sha256(ROOT / "target/release" / b) for b in BINARIES if (ROOT / "target/release" / b).exists()},
            "registered rebuild (OLD)": {b: sha256(pathlib.Path(old) / b) for b in BINARIES if old and (pathlib.Path(old) / b).exists()},
        },
        "host": {
            "cpu": cpu,
            "logical_cores": os.cpu_count(),
            "os": platform.platform(),
            "note": "shared cloud host, no CPU pinning; wall time is a practicality note, not the metric",
        },
        "arms": {
            "reference": {"KIC_F4_MULTIPLIERS": "occurring", "KIC_CHAIN_ORDER": "interleaved", "KIC_LINEAR_ELIM": "1", "KIC_F4_DROP": "complete"},
            "candidate": {"KIC_F4_MULTIPLIERS": "support", "KIC_CHAIN_ORDER": "interleaved", "KIC_LINEAR_ELIM": "1", "KIC_F4_DROP": "complete"},
        },
        "commands": [
            "cargo build --release --example groebner_stage_bench --example chain_ladder_screen --bin ic",
            "./target/release/examples/chain_ladder_screen --wide > research/support_local_multipliers_20260924/screen_wide.json",
            "research/support_local_multipliers_20260924/run.sh",
            "research/support_local_multipliers_20260924/oracle_ladder.sh",
            "research/support_local_multipliers_20260924/e2e.sh",
            "python3 research/support_local_multipliers_20260924/table.py > research/support_local_multipliers_20260924/tables.md",
            "OLD=<1f0d9751 build>/release research/support_local_multipliers_20260924/run_postfix.sh",
            "python3 research/support_local_multipliers_20260924/check_identity.py > research/support_local_multipliers_20260924/postfix_identity.md",
            "python3 research/support_local_multipliers_20260924/e2e_table.py > research/support_local_multipliers_20260924/e2e_tables.md",
        ],
    }
    out.write_text(json.dumps(doc, indent=2) + "\n")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
