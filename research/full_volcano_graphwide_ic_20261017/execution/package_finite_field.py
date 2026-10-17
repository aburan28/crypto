"""Metadata/hash/copy/archive only; native outputs supply every scientific value.

Usage: REPO DATASET BUILD GRAPH_DIR FRESH_STAGE FRESH_ARCHIVE
"""
import json
import shutil
import sys
from pathlib import Path
from package_checkpoint import archive_tree, digest

repo, dataset, build, graphs, stage, archive = map(Path, sys.argv[1:])
if stage.exists() or archive.exists():
    raise FileExistsError("fresh delivery paths required")
study = repo / "research/full_volcano_graphwide_ic_20261017"
execution = study / "execution"
prefix = json.loads((dataset / "runs/prefix.json").read_text())
if prefix["source_sha256"] != digest(execution / "finite_field_k1.rs"):
    raise ValueError("compiled source differs from supplied source")
for name in ["root", "prefix"]:
    receipt = json.loads((dataset / f"{name}-replay-match.json").read_text())
    if not receipt["ok"] or receipt["base_sha256"] != digest(dataset / f"runs/{name}.json") or receipt["candidate_sha256"] != digest(dataset / f"runs/{name}-reexecution.json"):
        raise ValueError("replay binding failed")
stage.mkdir()
shutil.copytree(execution / "finite-field-k1-evidence", stage / "evidence")
shutil.copytree(execution / "local-strata-evidence", stage / "local-strata")
shutil.copytree(study / "reports/run-20261017T0626Z-finite-field-C-prefix", stage / "report")
(stage / "inputs").mkdir()
shutil.copy2(repo / "research/full_volcano_classes_20261016/inputs/C.json", stage / "inputs/C.json")
source_db = Path("/workspace/isogeny_class_cli_delivery/graphs/partial-prefixes/C.sqlite")
if digest(source_db) != prefix["prefix"]["sha256"]:
    raise ValueError("prefix hash differs")
shutil.copy2(source_db, stage / "inputs/C-prefix.sqlite")
shutil.copy2(execution / "finite-field-k1-summary.json", stage / "summary.json")
(stage / "bin").mkdir()
for name in ["finite-field-k1", "local-strata"]:
    shutil.copy2(build / name, stage / "bin" / name)
for name in ["start-utc.txt", "finish-utc.txt", "bootstrap-note.txt", "missing-time-launch-stderr.log", "missing-time-launch-start-utc.txt"]:
    shutil.copy2(dataset / name, stage / name)
(stage / "source").mkdir()
for name in ["finite_field_k1.rs", "local_strata.rs", "BUILD_FINITE_FIELD.sh",
             "FINITE_FIELD_K1_PROTOCOL.md", "LOCAL_STRATA_PROTOCOL.md",
             "render_finite_field.py", "package_finite_field.py", "package_checkpoint.py"]:
    shutil.copy2(execution / name, stage / "source" / name)
shutil.copy2(repo / "Cargo.lock", stage / "source/Cargo.lock")
shutil.copy2(execution / "watchdog-repair/release-full-serial.log", stage / "source/release-full-serial.log")
shutil.copy2(study / "evidence/original-applicability.json", stage / "source/original-applicability.json")
calibration = json.loads((stage / "local-strata/calibration.json").read_text())
expected = {row["graph"]: row["db_sha256"] for row in calibration["calibrations"]}
(stage / "calibration-graphs").mkdir()
for label, sha in expected.items():
    source = graphs / f"{label}.sqlite"
    if digest(source) != sha:
        raise ValueError(f"calibration graph differs: {label}")
    shutil.copy2(source, stage / "calibration-graphs" / source.name)
(stage / "README.txt").write_text("""Original C finite-field reduction diagnostic; not the full C isogeny graph.
3738/3738 public target attempts verified at all 1869 retained prefix vertices.
The separate source-model run has 2/2 verified; never pool its repeated target
with the prefix population. Exact outcome/relation/counter re-execution matches
are in evidence/. Time fields are excluded, not asserted equal.

The genuine k=1 reduced Tate/Miller map feeds a projected-prime-subgroup
two-summand factor-base relation method, with unknown base and target logs.
This is NOT NFS or integer-smoothness index calculus. Its largest prime order
is 7, even though the field is 61 bits. No scaling or security claim follows.
All 15042 equations / 9171 destination queries are retained. More than one
equation per query is allowed; query hit rate and equations/query are distinct.

local-strata/ certifies 517 vertices across 24 observed conductor strata,
2985984 through 967458816; source conductor 967458816. 1352 vertices remain
unlabeled. Native per-stratum summaries retain the null-conductor bucket.
The local bounds were calibrated at five visibility cuts of each complete
graph. They reuse the standard volcano theorem and isogeny arithmetic, not an
independent formal or arithmetic proof. The original prefix DB and global
NOT_ESTABLISHED stratum receipt are unchanged. Full C has 1934901248 expected
vertices and 165 expected strata. Original A--E completion remains open.

Linux x86-64 reproduction (use fresh output names; no AWS needed):
  tar -xJf evidence/runs.tar.xz
  bin/finite-field-k1 inputs/C.json /fresh/root.json
  bin/finite-field-k1 inputs/C.json /fresh/prefix.json inputs/C-prefix.sqlite
  bin/finite-field-k1 compare runs/prefix.json /fresh/prefix.json /fresh/match.json
  bin/local-strata calibration-graphs inputs/C-prefix.sqlite /fresh/local-cert
  bin/local-strata summarize runs/prefix.json local-strata/C-local-strata.json /fresh/strata.json

The native summary and CSV expose phase operations, min/max/arithmetic mean,
actual n and missing observations. Native per-target elapsed intervals include
setup/projection/pairing/factor-base/relations/matrix/verification, excluding
process startup, input loading and output serialization/I/O. Bash process logs
include elapsed/user/system observations. Peak RSS is unavailable (GNU time is
absent). All times are unqualified; no external/reference speed claim is made.
Unpriced factorization, sqrt internals, allocation and matrix inversions stay
explicit. Exact generator representations are not promoted to the registry.
source/ plus toolchain/library hashes bind provenance; byte-reproducible
compilation is not asserted. Exact executed binaries are included.
""")
lines = [f"{digest(p)}  {p.relative_to(stage).as_posix()}" for p in sorted(stage.rglob("*")) if p.is_file()]
(stage / "MANIFEST.sha256").write_text("\n".join(lines) + "\n")
archive_tree(stage, archive)
print(json.dumps({"path":str(archive),"bytes":archive.stat().st_size,"sha256":digest(archive)}))
