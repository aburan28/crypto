"""Copy/hash/package frozen QA campaigns; no solver or statistical aggregation.

Usage: REPO CONTROL PILOT SUMMARY REPORT FRESH_STAGE FRESH_ARCHIVE
"""
import json
import shutil
import sys
from pathlib import Path
from package_checkpoint import archive_tree, digest

repo, control, pilot, summary, report, stage, archive = map(Path, sys.argv[1:])
if stage.exists() or archive.exists():
    raise FileExistsError("fresh delivery paths required")
execution = repo / "research/full_volcano_graphwide_ic_20261017/execution"
stage.mkdir()
(stage / "campaigns").mkdir()
for name, source in [("control",control),("pilot",pilot)]:
    wal = source / "ledger.sqlite-wal"
    if wal.exists() and wal.stat().st_size:
        raise ValueError("nonempty WAL: checkpoint first")
    for audit in source.glob("*.audit.json"):
        if not json.loads(audit.read_text())["ok"]:
            raise ValueError(f"failed QA audit: {audit}")
    shutil.copytree(source,stage / "campaigns" / name,
        ignore=shutil.ignore_patterns("writer.lock","ledger.sqlite-wal","ledger.sqlite-shm"))
shutil.copytree(summary,stage / "summary")
shutil.copytree(report,stage / "report")
(stage / "receipts").mkdir()
for name in ["control","pilot"]:
    source=Path(f"/workspace/volcano-repaired-{name}-match.json")
    if not json.loads(source.read_text())["ok"]:
        raise ValueError("failed match receipt")
    shutil.copy2(source,stage / "receipts" / source.name)
shutil.copytree(execution / "watchdog-repair",stage / "repair-evidence")
(stage / "source").mkdir()
for name in ["driver.rs","summarize.rs","match_profiles.rs","PROTOCOL.md",
             "OPTIMIZED_PROTOCOL.md","WATCHDOG_REPAIR_PROTOCOL.md",
             "package_repaired_qa.py","package_checkpoint.py","render_repaired_qa.py"]:
    shutil.copy2(execution / name,stage / "source" / name)
shutil.copy2(repo / "src/cryptanalysis/ecbench/runner.rs",stage / "source/runner.rs")
(stage / "README.txt").write_text("""Repaired release runner: control and three-root pilot, not the full survey.
60 attempts: 40 verified, 18 exhausted, 2 timeout. Four audited capsules.
Native comparison to the earlier release gates: 58 exact completed pairs,
2 censored pairs not compared. Successful and exhausted reports reproduce.
No frozen historical record, solver arithmetic or audit rule was changed.
Previous release QA is separately archived on Drive (delivery-release-qa.json
in the study). This archive does not pool its repeated targets or timing.

Each campaigns/ directory carries its exact ecbench and driver, immutable
queue, source/protocol, ledger, sealed sessions and receipts. Inspection:
  campaigns/control/volcano-ic-driver status campaigns/control
  campaigns/pilot/volcano-ic-driver status campaigns/pilot
Never edit a sealed session. Same-implementation replay is not independent
implementation validation. Wall observations are unqualified. GAE with
unpriced work is a lower bound; exhausted costs are not successful solve costs.
Timeout final costs/counters are missing, not zero. Raw process observations
are separate from completed solve timings.

These gates permit a new whole-queue bounded survey with 24823 vertices,
two public targets, one round, no warmups, 4096 trials and 3-second child
watchdogs. Actual launch/current progress require its separate receipt and
ledger; they are not established by this QA archive. Original A--E completion,
the full eight-target/five-round/million-trial matrix and large-extension
finite-field IC remain open. The C k=1 diagnostic is a separate dataset.
""")
lines=[f"{digest(p)}  {p.relative_to(stage).as_posix()}" for p in sorted(stage.rglob("*")) if p.is_file()]
(stage / "MANIFEST.sha256").write_text("\n".join(lines)+"\n")
archive_tree(stage,archive)
print(json.dumps({"path":str(archive),"bytes":archive.stat().st_size,"sha256":digest(archive)}))
