"""Render frozen N83 receipts to a review PDF; never executes research code."""

from __future__ import annotations

import html
import json
from pathlib import Path

from reportlab.graphics.shapes import Drawing, Line, Rect, String
from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import mm
from reportlab.platypus import (
    KeepTogether, PageBreak, Paragraph, SimpleDocTemplate, Spacer, Table, TableStyle,
)


ROOT = Path(__file__).resolve().parent
DATA = json.loads((ROOT / "report-data.json").read_text())
PDF = ROOT / "REPORT.pdf"
PAGE_W, PAGE_H = A4
styles = getSampleStyleSheet()
styles.add(ParagraphStyle(name="T", parent=styles["Title"], fontName="Helvetica-Bold", fontSize=20, leading=24, textColor=colors.HexColor("#172b4d"), spaceAfter=10))
styles.add(ParagraphStyle(name="H", parent=styles["Heading2"], fontName="Helvetica-Bold", fontSize=13, leading=16, textColor=colors.HexColor("#172b4d"), spaceBefore=12, spaceAfter=6))
styles.add(ParagraphStyle(name="B", parent=styles["BodyText"], fontName="Helvetica", fontSize=9, leading=13, spaceAfter=6))
styles.add(ParagraphStyle(name="S", parent=styles["BodyText"], fontName="Helvetica", fontSize=7.6, leading=10.5, spaceAfter=3))
styles.add(ParagraphStyle(name="C", parent=styles["BodyText"], fontName="Helvetica-Bold", fontSize=9, leading=12, alignment=TA_CENTER))


def p(value: object, sty: str = "B") -> Paragraph:
    return Paragraph(html.escape(str(value)).replace("\n", "<br/>"), styles[sty])


def table(rows: list[list[str]], widths: list[float]) -> Table:
    body = [[p(cell, "S") for cell in row] for row in rows]
    item = Table(body, colWidths=widths, repeatRows=1, hAlign="LEFT")
    item.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#e5edf5")),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#f6f8fb")]),
        ("GRID", (0, 0), (-1, -1), 0.4, colors.HexColor("#b8c6d6")),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("LEFTPADDING", (0, 0), (-1, -1), 6),
        ("RIGHTPADDING", (0, 0), (-1, -1), 6),
        ("TOPPADDING", (0, 0), (-1, -1), 5),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 5),
    ]))
    return item


def figure() -> Drawing:
    width = 470
    drawing = Drawing(width, 240)
    navy, blue, green, gold = [colors.HexColor(v) for v in ("#172b4d", "#2463a2", "#27836c", "#b45309")]
    drawing.add(String(0, 221, "Measured cgroup peak memory - one run per model", fontName="Helvetica-Bold", fontSize=11, fillColor=navy))
    entries = DATA["capacity_rows"]
    axis_x, axis_w = 166, 264
    for tick in range(5):
        x = axis_x + axis_w * tick / 4
        drawing.add(Line(x, 188, x, 48, strokeColor=colors.HexColor("#d8e0e8"), strokeWidth=0.5))
        drawing.add(String(x - 4, 194, str(tick), fontName="Helvetica", fontSize=8, fillColor=navy))
    drawing.add(String(axis_x + axis_w - 62, 39, "4 GiB limit", fontName="Helvetica", fontSize=8, fillColor=gold))
    for ix, row in enumerate(entries):
        y = 166 - ix * 35
        peak = row["peak_bytes"]
        drawing.add(String(0, y + 6, row["label"], fontName="Helvetica", fontSize=9, fillColor=navy))
        drawing.add(Rect(axis_x, y, peak / 4294967296 * axis_w, 17, fillColor=blue if ix == 0 else green, strokeColor=None))
        drawing.add(String(axis_x + peak / 4294967296 * axis_w + 5, y + 5, f"{peak / 1073741824:.3f} GiB", fontName="Helvetica", fontSize=8, fillColor=navy))
    drawing.add(String(0, 22, "Construction admission only; shared-host timings are not compared.", fontName="Helvetica", fontSize=8, fillColor=navy))
    return drawing


def footer(canvas, doc):
    canvas.saveState()
    canvas.setStrokeColor(colors.HexColor("#c6d2de"))
    canvas.line(15 * mm, 15 * mm, PAGE_W - 15 * mm, 15 * mm)
    canvas.setFont("Helvetica", 7)
    canvas.setFillColor(colors.HexColor("#526782"))
    canvas.drawString(15 * mm, 10 * mm, "N83 factor-base study - source-linked resource and storage evidence")
    canvas.drawRightString(PAGE_W - 15 * mm, 10 * mm, str(doc.page))
    canvas.restoreState()


story = [
    p("N83 retained-domain gates and K1182 storage", "T"),
    p("2026-10-10 | Public known-answer research | Frozen source 85e14930ff0aab0d50e8fb12ee7ca382e5e2f8e2", "S"),
    Spacer(1, 5),
    p("The minimum complete cold index-calculus runtime remains unmeasured. These receipts establish exact storage and bounded model/search admission, not a factor-base winner."),
    p("Exact workload", "H"),
    p("icv1-f2m83-tm6151469093347-debefd74: y^2 + xy = x^3 + 1 over z^83 + z^45 + z^2 + z + 1. Subgroup order 2417851639230796216685689; cofactor 4. The 53-bit diagnostic arm is a separate workload."),
    p("Owner budget: 7,200 seconds total, including 3,593.092994294 seconds charged before this extension. The extension ceiling is 3,606.907005706 seconds. The source-matched cross-build charged 872 seconds. Native guarded runs use one pinned CPU, 4 GiB, zero swap, no network, read-only inputs/root and fixed wall caps. The shared arm64 VM runs the x86-64 static worker under emulation; elapsed values are budget diagnostics only."),
    p("Observed gates", "H"),
]
observed = [["Gate", "Measured observation", "Outcome"]]
for row in DATA["observation_rows"]:
    observed.append([row["gate"], row["observation"], row["outcome"]])
story += [table(observed, [113, 241, 115]), Spacer(1, 8),
          p("Each model has one memory sample, without an uncertainty estimate. Measured finite-domain models and historical algebraic source bounds have different clause scopes.", "S"),
          PageBreak(), p("Capacity and stored-object coverage", "H"), figure(),
          p("All point records verified", "H"),
          p("K1182 contains 196,212 signed-Frobenius point records in 1,182 columns. Generic multi-limb BinaryCurve replay checked every point and representative against the native Gf2_128 producer on the same host/repository. The compressed object is 5,956,114 bytes."),
          p("Point-set BLAKE3: 429cde4bc514bf1026e4993e269192f858cb7410ade777178a1c0db7dbf0001a", "S"),
          p("Compressed BLAKE3: 76e1cbe4ea5a9d26e5b24a2d27f9401271dc471846529190a0416d8e09837a92", "S"),
          p("The native uploader downloaded the S3 object and checked the compressed bytes. Full objects stay outside Git; content-addressed receipts are committed. Storage coverage is 55 constructed objects, while external-host arithmetic replay is pending."),
          p("Bounded search and retained failures", "H"),
          table([["Case", "Final observation"]] + [[v["case"], v["observation"]] for v in DATA["search_rows"]], [113, 356]),
          Spacer(1, 8),
          p("The first K1182 launch failed at container Git ownership admission; a process-local safe.directory fixed the fresh retry. A host build receipt wrapper failed after compilation on a reserved zsh variable; the clean retry is retained. The six-summand Docker EOF is an infrastructure producer failure. Its original cleanup flag is invalid; the exact exited container state was captured and removal independently verified. No solver verdict follows from that attempt."),
          PageBreak(), p("Coverage against the requested objective", "H"),
          table([["Requirement", "Status and exact gap"]] + [[v["requirement"], v["status"]] for v in DATA["requirement_rows"]], [150, 319]),
          Spacer(1, 9),
          p("What remains open", "H"),
          p("The finite v1/v2 design addresses 57,024,000 tuples, most unexecuted. WDSat/FES width and system-shape gates, a useful large-prime partial producer, full-width sparse LA, natural independent rank, individual logarithm, group verification and fully charged matched cold comparisons remain open. A K1182 construction-only chained-S3 entry point also needs extension before that larger object can enter its search gate. Censored UNKNOWN outcomes are not refutations."),
          p("Canonical graph treatment", "H"),
          p("The index-calculus scoreboard receives a source-linked capacity/storage coverage panel. docs/ic/progress-timeline.json, the cold leaderboard, docs/curves and docs/performance-gains were checked. No admitted runtime ratio, exponent or curve-map result is available for their plotted series. No new algebraic identity is claimed."),
          p("Evidence locations", "H"),
          p("Frozen protocol, commands, source/input hashes, image IDs, worker/outer JSON, logs, container states, recovery audit, tests and this PDF source are in research/koblitz_n83_factor_base_sweep_20261008/verification/two-hour-extension-20261010/. The durable object is in s3://crypto-autoresearcher/factor-bases/icv1/etc/koblitz_n83_factor_base_sweep_20261008/v2-size-frontier/a0/objects/ with filename equal to the compressed digest plus .jsonl.gz. Review PR #1599 for the same-branch source and evidence change.", "S"),
]

doc = SimpleDocTemplate(str(PDF), pagesize=A4, rightMargin=15 * mm, leftMargin=15 * mm, topMargin=15 * mm, bottomMargin=20 * mm, title="N83 retained-domain gates and K1182 storage", author="Crypto Autoresearcher")
doc.build(story, onFirstPage=footer, onLaterPages=footer)
print(PDF)
