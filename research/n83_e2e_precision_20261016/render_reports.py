#!/usr/bin/env python3
"""Render the provider-bound N83 precision result and benchmark PDFs."""

import json
import shutil
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from reportlab.lib import colors
from reportlab.lib.pagesizes import letter
from reportlab.pdfgen import canvas

HERE = Path(__file__).resolve().parent
BENCH = HERE.parents[1] / "benchmarks/20261017-n83-e2e-precision"
RUN_ID = 38970960376
analysis = json.loads((HERE / f"analysis_{RUN_ID}.json").read_text())
rows = json.loads((HERE / f"ic1_rows_{RUN_ID}.json").read_text())
provider = json.loads((HERE / f"provider_{RUN_ID}.json").read_text())

assert analysis["evidence_status"] == "PASS_INDEPENDENT_REPLAY"
assert analysis["selection_admitted"] is True
assert analysis["point_count"] == 16 and rows["evidence_status"] == "PASS_INDEPENDENT_REPLAY"
assert rows["raw_artifact_sha256"] == provider["artifact_digest"].removeprefix("sha256:")
assert rows["raw_artifact_run_id"] == provider["run_id"] == RUN_ID

plt.rcParams.update({"font.size": 10, "svg.hashsalt": "n83-precision-38970960376"})
differences = analysis["paired_interval"]["differences_ms"]
interval = analysis["paired_interval"]
fig, ax = plt.subplots(figsize=(9.2, 4.5))
ax.bar(range(1, 17), differences,
       color=["#176b85" if value >= 0 else "#bd6842" for value in differences])
ax.axhline(0, color="#323c43", linewidth=1)
ax.axhline(interval["mean_ms"], color="#114f66", linestyle="--", linewidth=1.5,
           label=f"Mean: {interval['mean_ms']:.1f} ms")
ax.fill_between([0.3, 16.7], interval["lower_95_ms"], interval["upper_95_ms"],
                color="#114f66", alpha=0.11, label="95% paired mean interval")
ax.set_xlim(0.3, 16.7)
ax.set_xticks(range(1, 17))
ax.set_xlabel("Sealed public point")
ax.set_ylabel("K1104 minus K1120 cold process (ms)")
ax.set_title("N83 factor-base cold precision: sixteen fresh points")
ax.legend(loc="upper right", frameon=False)
fig.tight_layout()
for suffix in ("svg", "png", "pdf"):
    fig.savefig(HERE / f"precision.{suffix}", dpi=180,
                metadata={"Creator": "N83 precision report"})
plt.close(fig)


def report(name, heading, lines, table):
    path = HERE / name
    pdf = canvas.Canvas(str(path), pagesize=letter)
    pdf.setTitle(heading)
    pdf.setAuthor("N83 factor-base research")
    width, height = letter
    pdf.setFillColor(colors.HexColor("#163f55"))
    pdf.rect(0, height - 12, width, 12, stroke=0, fill=1)
    pdf.setFont("Helvetica-Bold", 16)
    pdf.drawString(36, height - 48, heading)
    y = height - 78
    pdf.setFillColor(colors.black)
    pdf.setFont("Helvetica", 9)
    for line in lines:
        pdf.drawString(36, y, line)
        y -= 15
    y -= 10
    pdf.setFillColor(colors.HexColor("#eef4f6"))
    pdf.rect(36, y - 8 - 23 * len(table), width - 72, 24 * len(table), fill=1, stroke=0)
    pdf.setFillColor(colors.black)
    for label, value in table:
        pdf.setFont("Helvetica-Bold", 9)
        pdf.drawString(44, y, label)
        pdf.setFont("Helvetica", 9)
        pdf.drawRightString(width - 44, y, str(value))
        y -= 23
    y -= 22
    pdf.drawImage(str(HERE / "precision.png"), 36, max(76, y - 300),
                  width=width - 72, height=264, preserveAspectRatio=True, anchor="c")
    pdf.setFont("Helvetica", 8)
    pdf.drawString(36, 43, "Source: GitHub Actions run 38970960376; provider SHA-256 559578cd...fd844d")
    pdf.drawString(36, 30, "Research record: research/n83_e2e_precision_20261016")
    pdf.save()


report(
    "RUN.pdf", "N83 sixteen-point precision run",
    ["Tier 1 cloud instance; sealed public points; one fresh process per candidate repeat.",
     "Each point: K1104 A/B, K1120 A/B, and one independently solved same-point rho.",
     "The provider ZIP, all extracted members, 64 IC arms and 16 rho solves were replayed.",
     "Cold process includes base and index setup, rank, and target work; online is separate."],
    [("Provider run", RUN_ID), ("Provider ZIP SHA-256", rows["raw_artifact_sha256"]),
     ("Bound ZIP result members", rows["bound_result_members"]),
     ("Evidence", rows["evidence_status"]),
     ("Selection status", rows["selection_status"])],
)
report(
    "RESULT.pdf", "N83 cold difference bounded on sixteen points",
    ["K1104/root11 and K1120/root11 were compared on 16 fresh public points.",
     "The paired 95% mean interval crosses zero, so the frozen selector remains null.",
     "All 64 IC and 16 same-point rho computations passed independent replay."],
    [("K1104 mean cold process", f"{rows['pipeline_summaries']['k1104']['mean_ms']['cold_process_ms']/1000:.3f} s"),
     ("K1120 mean cold process", f"{rows['pipeline_summaries']['k1120']['mean_ms']['cold_process_ms']/1000:.3f} s"),
     ("Mean K1104 minus K1120", f"{interval['mean_ms']:.3f} ms"),
     ("Paired 95% interval", f"[{interval['lower_95_ms']:.3f}, {interval['upper_95_ms']:.3f}] ms"),
     ("Selected candidate", "null")],
)

BENCH.mkdir(parents=True, exist_ok=True)
for name in ("RUN.pdf", "RESULT.pdf", "precision.pdf", "precision.png", "precision.svg"):
    shutil.copyfile(HERE / name, BENCH / name)
