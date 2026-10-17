"""Packaging only. Native inventory/summary are the authority for all counts.

Usage: CHECKPOINT SUMMARY REPORT GRAPH_DIR PINNED_SOURCE HELPERS FRESH_STAGE
       FRESH_ARCHIVE REPAIR_EVIDENCE
"""
import json
import shutil
import sys
from pathlib import Path
from package_checkpoint import archive_tree, digest


def main():
    checkpoint, summary, report, graphs, pinned, helpers, stage, out, repair = map(Path, sys.argv[1:])
    if stage.exists() or out.exists():
        raise FileExistsError("fresh stage and archive required")
    inv = json.loads((checkpoint / "inventory.json").read_text())
    doc = json.loads((summary / "summary.json").read_text())
    manifest = json.loads((checkpoint / "campaign/campaign.json").read_text())
    if not inv["ok"] or inv["outcome_totals"] != doc["outcome_totals"]:
        raise ValueError("native inventory and summary must agree")
    if manifest["ecbench_sha256"] != doc["binary_sha256"]:
        raise ValueError("summary binary differs from frozen campaign")
    shutil.copytree(checkpoint, stage)
    shutil.copytree(summary, stage / "summary")
    shutil.copytree(report, stage / "report")
    shutil.copytree(repair, stage / "diagnostics")
    (stage / "graphs").mkdir()
    expected = {}
    for line in (stage / "campaign/queue.jsonl").read_text().splitlines():
        row = json.loads(line)
        expected[row["graph"]] = row["input"]["graph_sha256"]
    for graph, sha in expected.items():
        source = graphs / f"{graph}.sqlite"
        if digest(source) != sha:
            raise ValueError(f"graph hash mismatch: {graph}")
        shutil.copy2(source, stage / "graphs" / source.name)
    study = pinned / "research/full_volcano_graphwide_ic_20261017"
    (stage / "source").mkdir()
    for relative, name in [
        ("src/cryptanalysis/ic_framework/plugins.rs", "plugins.rs"),
        ("src/cryptanalysis/ecbench/methods.rs", "methods.rs"),
        ("src/cryptanalysis/ecbench/runner.rs", "runner-before-repair.rs"),
        ("Cargo.lock", "Cargo.lock"),
    ]:
        shutil.copy2(pinned / relative, stage / "source" / name)
    for name in ["PROTOCOL.md", "OPTIMIZED_PROTOCOL.md", "driver.rs", "summarize.rs",
                 "checkpoint.rs", "match_profiles.rs", "BUILD.sh"]:
        shutil.copy2(study / "execution" / name, stage / "source" / name)
    shutil.copy2(study / "PROTOCOL.md", stage / "source/full-matrix-protocol.md")
    shutil.copy2(Path(__file__), stage / "source/package_release_checkpoint.py")
    shutil.copy2(Path(__file__).with_name("package_checkpoint.py"), stage / "source/package_checkpoint.py")
    for name in ["volcano-ic-checkpoint", "volcano-ic-summarize-v2"]:
        shutil.copy2(helpers / name, stage / "bin" / name)
    (stage / "README.txt").write_text(f"""Stopped RELEASE survey checkpoint, run 20261017T0540Z.
{inv['vertices_audited']} audited vertices of {inv['vertices_planned']}; {inv['strata_audited']} of {inv['strata_planned']} strata sampled.
Accepted {doc['measured_records']} attempts; exact outcome totals: {json.dumps(doc['outcome_totals'], sort_keys=True)}.
Ledger states: {json.dumps(inv['states'], sort_keys=True)}.

The failed three-v10-attempt-001 capsule and its failed standard audit are
preserved under campaign/. They are NOT accepted observations and are excluded
from summary/ and report/. Record 10 combines timeout and a checked scalar.
The worker halted after 432 elapsed seconds. Do not silently resume this old
binary: the watchdog recording repair requires a separately frozen campaign.
diagnostics/ describes the failure, tests and access probes; these were not
part of this historical measured execution. No sealed record was modified.

inventory.json binds the accepted capsules' files, record seals and existing
replay receipts. It is not a new or independent solver replay. Native summary
v2 supplies per-graph/stratum/arm/status min/max/mean with explicit denominators.
Timeout counts are unknown, not zero. GAE with unpriced work is a lower bound.
All wall observations are unqualified. This release dataset is separate from
the earlier dev checkpoint and from controls/pilots; never pool repetitions.

graphs/ contains complete graphs only for the three declared base fields and
finite prime-degree sets: 24,823 vertices, 51 strata, 69,300 directed maps.
These do not complete original A--E or enumerate infinitely many compositions.
Source snapshot {manifest['source_head']}; campaign.json binds exact binaries,
queue and input hashes. Source snapshots do not imply reproducible compilation.
The original full benchmark matrix, original A--E completion, finite-field IC
and independent implementation replay remain unresolved.

Inspection: bin/volcano-ic-driver status campaign
Summary: bin/volcano-ic-summarize-v2 /fresh/summary campaign
Never edit a sealed session. This archive contains executables for Linux x86-64.
""")
    lines = [f"{digest(p)}  {p.relative_to(stage).as_posix()}"
             for p in sorted(stage.rglob("*")) if p.is_file()]
    (stage / "MANIFEST.sha256").write_text("\n".join(lines) + "\n")
    archive_tree(stage, out)
    print(json.dumps({"path": str(out), "bytes": out.stat().st_size, "sha256": digest(out)}))


if __name__ == "__main__":
    main()
