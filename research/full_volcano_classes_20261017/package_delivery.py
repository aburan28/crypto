"""Artifact packaging only: copy files, hash bytes and create/check the ZIP."""
import hashlib
import json
from pathlib import Path
import shutil
import zipfile

STUDY = Path(__file__).resolve().parent
DEST = Path("/workspace/complete_volcano_delivery_20261017")
SOURCE = DEST / "source/research/full_volcano_classes_20261017"
for subdir in ["inputs", "reports", "generator"]:
    shutil.copytree(STUDY / subdir, SOURCE / subdir, dirs_exist_ok=True)
for filename in ["PROTOCOL.md", "README.md", "run.sh", "render_report.py"]:
    shutil.copy2(STUDY / filename, SOURCE / filename)
    shutil.copy2(STUDY / filename, DEST / filename)
evidence = SOURCE / "evidence"
evidence.mkdir(exist_ok=True)
for path in (STUDY / "evidence").iterdir():
    if path.name.startswith(("cover-", "curve-traits-", "browser-", "leaderboard-")):
        continue  # Auxiliary catalogue jobs can still be active while packaging.
    if path.is_file():
        shutil.copy2(path, evidence / path.name)
shutil.copytree(evidence, DEST / "evidence", dirs_exist_ok=True)
manifest = {
    "schema": "complete-volcano-delivery/v1",
    "generation_commit": "946ab7419",
    "native_fixture_export_source_commit": "77d8ed0e5",
    "cli_and_isogeny_algos_snapshot_commit": "7c2e8f46aa5aaa78e8b3ef3e0bcde93c85d3d9f2",
    "cli_source_sha256": "d3c77f50dc2e737de38561203ce3693f7899fe3db502b8e863c7ac69e5818bed",
    "binary_sha256": "bbacacca68712e32e60c6dac778ce99e43e00f8309628a9f58a994d3b0681dce",
    "binary_platform": "Linux x86_64 glibc; bundled SQLite",
    "complete_cases": json.loads((STUDY / "evidence/coverage.json").read_text()),
    "folder_url": "https://drive.google.com/drive/folders/1aOfT7QDntQZ2rxnYAauit0SYn86atA3g",
    "independent_implementation_replay": "NOT_RUN",
    "solver_measurements": "NOT_RUN",
    "original_A_C_D_E": "INCOMPLETE; these three classes are additional examples"
}
(DEST / "BUILD.json").write_text(json.dumps(manifest, indent=2) + "\n")
(DEST / "delivery.json").write_text(json.dumps({"folder_url": manifest["folder_url"], "archive_receipt": "Repository research/full_volcano_classes_20261017/delivery.json is written after upload"}, indent=2) + "\n")
(DEST / "START_HERE.txt").write_text("""Three complete ordinary isogeny graphs, 2026-10-17

Graphs: graphs/{two,three,mixed}.sqlite
Readable full exports: graphs/{two,three,mixed}/nodes.csv and edges.jsonl
Each stratum: stratum-<conductor>-representative.json (explicit rational maps)
Report: reports/run-20261017/report.pdf
Inputs: inputs/*.json; *-floor.json starts at the maximum End(E) conductor.
Evidence: evidence/; complete counts: evidence/coverage.json

Linux x86_64:
  ./bin/isogeny-class status --input graphs/mixed.sqlite
  ./bin/isogeny-class vertex --input graphs/mixed.sqlite --id 0
  cp graphs/mixed.sqlite replay.sqlite
  ./bin/isogeny-class verify --input replay.sqlite
  ./bin/isogeny-class build --curve inputs/mixed-floor.json --out fresh.sqlite --steps all

verify updates metadata; use copies to preserve byte-frozen originals.
Graph generation does not run ECDLP or index calculus.
The entire vertex class is present for each of these THREE classes.
Edges cover the declared prime degrees, not infinitely many compositions.
Original A/C/D/E are still partial; the prior full B dataset remains separate.

Rebuild the delivered binary from source/:
  cargo build --release --locked --manifest-path source/tools/volcano-cli/Cargo.toml
The source snapshot is 7c2e8f46a; BUILD.json records provenance and SHA256SUMS
records every delivered file. The native fixture/export sources are additional.
""")
files = []
for path in sorted(DEST.rglob("*")):
    if path.is_file() and path.name != "SHA256SUMS":
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        files.append((str(path.relative_to(DEST)), digest))
(DEST / "SHA256SUMS").write_text("".join(f"{digest}  {name}\n" for name, digest in files))
archive = Path("/workspace/three_complete_volcano_graphs_20261017.zip")
assert not archive.exists(), "preserve existing delivery snapshots"
with zipfile.ZipFile(archive, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=6) as z:
    for path in sorted(DEST.rglob("*")):
        if path.is_file():
            z.write(path, path.relative_to(DEST))
with zipfile.ZipFile(archive) as z:
    assert z.testzip() is None
receipt = {"path": str(archive), "bytes": archive.stat().st_size,
           "sha256": hashlib.sha256(archive.read_bytes()).hexdigest(),
           "zip_crc_verified": True, "files": len(files) + 1}
(STUDY / "evidence/archive-local.json").write_text(json.dumps(receipt, indent=2) + "\n")
print(json.dumps(receipt))
