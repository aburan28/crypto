"""Packaging/metadata only. No arithmetic, solvers, measurements or orchestration."""
import hashlib
import json
import shutil
import sys
import zipfile
from pathlib import Path

repo, out = map(Path, sys.argv[1:3])
stage = out.with_suffix("")
stage.mkdir(exist_ok=False, parents=True)
study = repo / "research/full_volcano_graphwide_ic_20261017"
shutil.copytree(study / "execution", stage / "execution")
shutil.copytree(study / "reports", stage / "reports")
shutil.copy(study / "PROTOCOL.md", stage / "full-matrix-protocol.md")
shutil.copy(study / "README.md", stage / "study-README.md")
bin_dir = stage / "bin"
bin_dir.mkdir()
for src, name in [
    ("/workspace/volcano_ic_survey_20261017/ecbench", "ecbench"),
    ("/workspace/scratch/volcano-ic-driver", "volcano-ic-driver"),
    ("/workspace/scratch/volcano-ic-driver-initial", "volcano-ic-driver-initial"),
    ("/workspace/scratch/volcano-ic-summarize", "volcano-ic-summarize"),
]:
    shutil.copy(src, bin_dir / name)
shutil.copy("/workspace/volcano_ic_survey_20261017/queue.jsonl", stage / "execution/campaigns/survey-prefix-12/queue.jsonl")
shutil.copy("/workspace/volcano_ic_survey_20261017/driver.rs", stage / "execution/campaigns/survey-prefix-12/driver.rs")
shutil.copy("/workspace/volcano_ic_survey_20261017/execution-protocol.md", stage / "execution/campaigns/survey-prefix-12/execution-protocol.md")
source = stage / "measured-library-source"
source.mkdir()
shutil.copy(repo / "src/cryptanalysis/ic_framework/plugins.rs", source / "plugins.rs")
shutil.copy(repo / "src/cryptanalysis/ecbench/methods.rs", source / "methods.rs")
shutil.copy(repo / "Cargo.lock", source / "Cargo.lock")
shutil.copy("/workspace/scratch/volcano-lab-browser-current.json", stage / "staged-browser-data.json")
(stage / "README.txt").write_text("""MEASURED IC DATA: control + pilot + first 12 survey vertices, 2026-10-17.
Survey checkpoint: 192 attempts, 119 verified, 58 exhausted, 15 timed out.
Pilot/control are separate repeated workloads; do not pool as independent targets.
This is NOT the completed full matrix, all 24,823 solved vertices, or all A--E.

Reports: reports/run-20261017T0240Z-control-pilot/report.pdf and
reports/run-20261017T0245Z-survey-prefix-12/report.pdf.
Detailed statistics: execution/{pilot-summary-v2,survey-prefix-summary-v2}.
Raw sealed sessions/audits/checkpoints: execution/campaigns/.
The survey queue contains every planned vertex, including queued/not-run rows.

Linux x86-64 binaries are supplied once in bin/ to avoid duplicate copies.
Copy bin/ecbench into each campaign directory before driver resume, then chmod
+x on bin/* and each copied ecbench if your ZIP extractor loses permissions.
Control/pilot froze bin/volcano-ic-driver-initial; the survey freezes the current
bin/volcano-ic-driver. Their manifest hashes must match. Standard ecbench audits
use bin/ecbench. Exhaustion receipts supplement its successful-only replay.
Never overwrite a sealed session. A new failed-attempt retry gets a new directory.

The survey snapshot's SQLite ledger includes 12 audited and 24,811 queued
vertices. The first 51 queue entries cover one representative per stratum.
Full matrix launch is gated as documented; the survey is bounded to 4,096
trials, two targets, one round and a 3-second watchdog per child.
All wall timings are L0/debug and unqualified. Min/max/mean use explicit n,
status and present-value denominators; timeout counts stay unknown.
Finite-field IC at the actual recorded embedding degrees is not implemented.

Measured library source commit: 4eecf6b53 (public aburan28/crypto).
Current frozen survey driver: ee54ad603fef4768817fbca2f133564e81ce84f1.
Complete repository source can be checked out at those public commits.
Exact binary/input/queue hashes are in campaign.json and MANIFEST.sha256.
Prepared browser data is staged, not the canonical conflicted page. Canonical
browser/leaderboard generated files await the existing Copilot conflict request.
Original graph and input Drive uploads remain unchanged.
""")
manifest = []
for p in sorted(stage.rglob("*")):
    if p.is_file():
        manifest.append(f"{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.relative_to(stage).as_posix()}")
(stage / "MANIFEST.sha256").write_text("\n".join(manifest) + "\n")
with zipfile.ZipFile(out, "x", compression=zipfile.ZIP_DEFLATED, compresslevel=6) as archive:
    for p in sorted(stage.rglob("*")):
        if p.is_file():
            archive.write(p, p.relative_to(stage).as_posix())
with zipfile.ZipFile(out) as archive:
    bad = archive.testzip()
    assert bad is None, bad
print(json.dumps({"path": str(out), "bytes": out.stat().st_size,
                  "sha256": hashlib.sha256(out.read_bytes()).hexdigest(),
                  "local_zip_crc": "PASS"}))
