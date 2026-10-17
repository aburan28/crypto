"""Metadata/copy/archive only: never solve, run experiments or aggregate counts.

Usage: package_checkpoint.py CHECKPOINT SUMMARY REPORT GRAPH_DIR PINNED_SOURCE
       HELPER_BIN_DIR FRESH_STAGE FRESH_ARCHIVE
"""
import hashlib
import json
import shutil
import sys
import tarfile
from pathlib import Path


def digest(path):
    with path.open("rb") as src:
        return hashlib.file_digest(src, "sha256").hexdigest()


def archive_tree(stage, destination):
    # Stable order, metadata, permissions and xz stream, independent of checkout
    # timestamps and extraction ownership. No links or external filesystem data.
    with tarfile.open(destination, "x:xz", format=tarfile.GNU_FORMAT,
                      preset=6) as archive:
        for path in sorted(stage.rglob("*")):
            if path.is_symlink():
                raise ValueError(f"links are not package inputs: {path}")
            info = archive.gettarinfo(path, path.relative_to(stage).as_posix())
            info.uid = info.gid = info.mtime = 0
            info.uname = info.gname = ""
            info.mode = 0o755 if path.is_dir() or path.stat().st_mode & 0o111 else 0o644
            if path.is_file():
                with path.open("rb") as src:
                    archive.addfile(info, src)
            else:
                archive.addfile(info)


def main():
    checkpoint, summary, report, graphs, pinned, helpers, stage, out = map(Path, sys.argv[1:])
    if stage.exists() or out.exists():
        raise FileExistsError("fresh stage and archive paths required")
    shutil.copytree(checkpoint, stage)
    shutil.copytree(summary, stage / "summary")
    shutil.copytree(report, stage / "report")
    (stage / "graphs").mkdir()
    # Bind supplied graph files to hashes already certified in the native queue.
    expected = {}
    with (stage / "campaign/queue.jsonl").open() as queue:
        for line in queue:
            row = json.loads(line)
            expected[row["graph"]] = row["input"]["graph_sha256"]
    for graph, sha in expected.items():
        source = graphs / f"{graph}.sqlite"
        if digest(source) != sha:
            raise ValueError(f"graph hash mismatch: {source}")
        shutil.copy2(source, stage / "graphs" / source.name)
    study = pinned / "research/full_volcano_graphwide_ic_20261017"
    (stage / "source").mkdir()
    for relative, name in [
        ("src/cryptanalysis/ic_framework/plugins.rs", "plugins.rs"),
        ("src/cryptanalysis/ecbench/methods.rs", "methods.rs"),
        ("Cargo.lock", "Cargo.lock"),
    ]:
        shutil.copy2(pinned / relative, stage / "source" / name)
    for name in ["PROTOCOL.md", "OPTIMIZED_PROTOCOL.md", "driver.rs", "summarize.rs",
                 "checkpoint.rs", "match_profiles.rs", "BUILD.sh", "package_checkpoint.py"]:
        # The packaging script may be newer; execution sources come from pinned
        # source, not the checkout undergoing a base-branch merge.
        source = Path(__file__) if name == "package_checkpoint.py" else study / "execution" / name
        shutil.copy2(source, stage / "source" / name)
    shutil.copy2(study / "PROTOCOL.md", stage / "source/full-matrix-protocol.md")
    for name in ["volcano-ic-checkpoint", "volcano-ic-summarize-v2"]:
        shutil.copy2(helpers / name, stage / "bin" / name)
    (stage / "README.txt").write_text("""Measured stopped checkpoint: 99 / 24,823 vertices, all 51 strata.
1,584 attempts: 1,065 verified, 422 exhausted, 97 timed out. Not full execution.
Graphs are complete only for the declared base field and prime-degree sets.
They do not complete the original A--E request or enumerate all compositions.

campaign/ is the native byte-preserving stopped campaign, including its frozen
ecbench, full queue, ledger, source/protocol, sealed sessions and replay receipts.
bin/volcano-ic-driver is its exact frozen driver. The control/pilot and optimized
build are separate datasets; do not pool their timings or repeated targets here.
inventory.json binds sealed files, identities and previously issued replay
receipts. It is not a fresh solver replay or independent implementation check.

summary/ is native summary schema v2. Trials include skipped repeated queries;
actual decomposition calls subtract those skips. Both denominators and yields
are explicit. Timeouts have unknown final counters, not zero. All wall times
are descriptive and unqualified. GAE values with unpriced work are lower bounds.
report/ contains the inspected three-page PDF and all 51 stratum plot data.

Linux x86-64 resume: bin/volcano-ic-driver run campaign 24823 86400 1073741824
uses the exact frozen dev build and stops on guards/errors. The copied campaign
has no STOP file; do not run this command merely to inspect the archive.
Use bin/volcano-ic-driver status campaign to inspect, or regenerate a summary
with bin/volcano-ic-summarize-v2 /fresh/summary campaign. The lock is recreated
by the driver. For a read-only receipt inventory after extraction, create an
empty campaign/writer.lock and invoke bin/volcano-ic-checkpoint campaign
bin/volcano-ic-driver /fresh/checkpoint. Never edit a sealed session.

Source snapshot: public aburan28/crypto commit 9bca26ae2; campaign.json records
the earlier measured binary/driver provenance. The supplied source snapshot
does not claim byte-reproducible compiler output. Original inputs' NOT_RUN
fields are historical preparation status; ledger.sqlite is execution status.
Finite-field IC at these embedding degrees, full matrix, exact registration,
independent replay and original A--E completion remain open.
""")
    manifest = [f"{digest(p)}  {p.relative_to(stage).as_posix()}"
                for p in sorted(stage.rglob("*")) if p.is_file()]
    (stage / "MANIFEST.sha256").write_text("\n".join(manifest) + "\n")
    archive_tree(stage, out)
    print(json.dumps({"path": str(out), "bytes": out.stat().st_size,
                      "sha256": digest(out)}))


if __name__ == "__main__":
    main()
