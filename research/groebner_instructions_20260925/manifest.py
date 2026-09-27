#!/usr/bin/env python3
"""Write manifest.json for the runs of RESEARCH_GROEBNER_STAGE_INSTRUCTIONS.md.

    python3 research/groebner_instructions_20260925/manifest.py

Instructions are build-dependent, so the manifest pins the binaries, the
valgrind version and the commit they were built from; refuses to overwrite.
"""
import hashlib, json, os, pathlib, platform, subprocess

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parents[1]
BINARIES = ["target/release/examples/groebner_stage_bench", "target/release/examples/ir_calibration", "target/release/ic"]
SOURCES = [
    "src/cryptanalysis/koblitz_groebner.rs",
    "src/cryptanalysis/inherited_f4.rs",
    "src/cryptanalysis/koblitz_index_calculus.rs",
    "src/cryptanalysis/ic_boundary.rs",
    "examples/groebner_stage_bench.rs",
    "examples/ir_calibration.rs",
]
BUILT_FROM = "31bedbdf"


def sha256(path):
    return hashlib.sha256((ROOT / path).read_bytes()).hexdigest()


def main():
    out = HERE / "manifest.json"
    if out.exists():
        raise SystemExit(f"{out} exists; never overwritten")
    git = lambda *a: subprocess.check_output(["git", *a], cwd=ROOT, text=True).strip()
    cpu = next((l.split(":", 1)[1].strip() for l in open("/proc/cpuinfo") if l.startswith("model name")), "")
    gates = subprocess.run(["python3", str(HERE / "table.py")], capture_output=True, text=True)
    doc = {
        "round": "groebner-instructions-20260925",
        "note": "research/notes/ecc2k130/RESEARCH_GROEBNER_STAGE_INSTRUCTIONS.md",
        "scope": "accounting: the decomposition oracle and whole logarithms priced in instructions retired (callgrind)",
        "git_head": git("rev-parse", "HEAD"),
        "built_from": BUILT_FROM,
        "measured_sources_unchanged_since_build": git("diff", "--stat", BUILT_FROM, "HEAD", "--", *SOURCES) == "",
        "sources_sha256_at_head": {s: sha256(s) for s in SOURCES},
        "binaries_sha256": {b: sha256(b) for b in BINARIES if (ROOT / b).exists()},
        "valgrind": subprocess.check_output(["valgrind", "--version"], text=True).strip(),
        "gates_exit_code": gates.returncode,
        "gates_note": "non-zero because G2 failed on one rung (thread scheduling); see tables.md and the note §5.1",
        "host": {
            "cpu": cpu,
            "logical_cores": os.cpu_count(),
            "os": platform.platform(),
            "note": "instructions do not depend on the host; they do depend on the binary",
        },
        "arms": {
            "A0": "KIC_CHAIN_ORDER=layout KIC_LINEAR_ELIM=0 KIC_F4_MULTIPLIERS=occurring KIC_F4_DROP=complete",
            "A1": "KIC_CHAIN_ORDER=interleaved KIC_LINEAR_ELIM=1 KIC_F4_MULTIPLIERS=occurring KIC_F4_DROP=complete",
            "A2": "KIC_CHAIN_ORDER=interleaved KIC_LINEAR_ELIM=1 KIC_F4_MULTIPLIERS=support KIC_F4_DROP=complete",
        },
        "commands": [
            "cargo build --release --example groebner_stage_bench --example ir_calibration --bin ic",
            "research/groebner_instructions_20260925/run.sh",
            "python3 research/groebner_instructions_20260925/table.py > research/groebner_instructions_20260925/tables.md",
        ],
    }
    out.write_text(json.dumps(doc, indent=2) + "\n")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
