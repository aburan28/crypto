"""Typeset the adjacent Markdown report; no research computation occurs here."""

from html import escape, unescape
from pathlib import Path
import re

from reportlab.lib import colors
from reportlab.lib.styles import ParagraphStyle
from reportlab.platypus import (
    Flowable,
    Paragraph,
    SimpleDocTemplate,
    Spacer,
    Table,
    TableStyle,
)


HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[3]
REPORT_BRANCH = "codex/ecdlp-index-calculus-lead-assessment-20261011"
PAGE_W = 612
CONTENT_W = 526

body = ParagraphStyle(
    "body", fontName="Helvetica", fontSize=8.3, leading=10.8, spaceAfter=3,
    splitLongWords=1,
)
title = ParagraphStyle(
    "title", parent=body, fontName="Helvetica-Bold", fontSize=17,
    leading=21, spaceAfter=13,
)
heading = ParagraphStyle(
    "heading", parent=body, fontName="Helvetica-Bold", fontSize=12,
    leading=15, spaceBefore=11, spaceAfter=6,
)
cell = ParagraphStyle(
    "cell", parent=body, fontSize=7.3, leading=9.3, spaceAfter=0,
    splitLongWords=1,
)
def inline(md):
    value = escape(md)

    def link(match):
        label, target = match.groups()
        target = unescape(target)
        if target.startswith("https://"):
            url = target
        else:
            path = (HERE / target).resolve().relative_to(REPO_ROOT)
            branch = REPORT_BRANCH if path.is_relative_to(HERE.relative_to(REPO_ROOT)) else "main"
            url = f"https://github.com/aburan28/crypto/blob/{branch}/{path.as_posix()}"
        return f'<link href="{escape(url, quote=True)}" color="#215a91">{label}</link>'

    value = re.sub(r"\[([^\]]+)\]\(([^)]+)\)", link, value)
    value = re.sub(r"\*\*([^*]+)\*\*", r"<b>\1</b>", value)
    value = re.sub(r"\*([^*]+)\*", r"<i>\1</i>", value)
    value = re.sub(r"\x60([^\x60]+)\x60", r'<font name="Courier">\1</font>', value)
    return value


class LeadMap(Flowable):
    def __init__(self):
        super().__init__()
        self.width = CONTENT_W
        self.height = 115

    def draw(self):
        canvas = self.canv
        canvas.scale(1, self.height / 145)
        canvas.setFillColor(colors.HexColor("#f7fafc"))
        canvas.roundRect(0, 0, self.width, self.height, 8, stroke=0, fill=1)
        canvas.setFont("Helvetica-Bold", 9)
        canvas.setFillColor(colors.HexColor("#142536"))
        canvas.drawString(12, 131, "Where each lead enters the charged pipeline")

        stages = [
            ("Public", "point"),
            ("Base /", "index"),
            ("Ordinary", "relations"),
            ("Useful", "rank"),
            ("Target /", "verify"),
        ]
        xs = [10, 116, 222, 328, 434]
        for idx, ((line1, line2), x) in enumerate(zip(stages, xs)):
            canvas.setFillColor(colors.HexColor("#eaf2fb"))
            canvas.setStrokeColor(colors.HexColor("#315d87"))
            canvas.roundRect(x, 65, 82, 43, 5, stroke=1, fill=1)
            canvas.setFont("Helvetica", 8)
            canvas.setFillColor(colors.HexColor("#142536"))
            canvas.drawCentredString(x + 41, 89, line1)
            canvas.drawCentredString(x + 41, 76, line2)
            if idx < 4:
                canvas.setStrokeColor(colors.HexColor("#315d87"))
                canvas.line(x + 83, 86, x + 105, 86)
                canvas.line(x + 101, 89, x + 105, 86)
                canvas.line(x + 101, 83, x + 105, 86)

        leads = [
            (10, "Cover / isogeny", 51),
            (116, "Compact orbit", 157),
            (222, "Higher arity", 263),
            (328, "Early pruning", 263),
            (434, "Hardware / reuse", 475),
        ]
        for x, label, dest_x in leads:
            canvas.setFillColor(colors.HexColor("#fff7e7"))
            canvas.setStrokeColor(colors.HexColor("#a66b11"))
            canvas.roundRect(x, 8, 82, 32, 5, stroke=1, fill=1)
            canvas.setFillColor(colors.HexColor("#142536"))
            canvas.setFont("Helvetica", 6.7)
            canvas.drawCentredString(x + 41, 25, label)
            canvas.setStrokeColor(colors.HexColor("#a66b11"))
            canvas.line(x + 41, 40, dest_x, 62)
        canvas.setFillColor(colors.HexColor("#2d7b49"))
        canvas.setFont("Helvetica-Bold", 7.6)
        canvas.drawString(12, 113, "Reference: same-point automorphism-aware rho")


def footer(canvas, doc):
    canvas.saveState()
    canvas.setFont("Helvetica", 8)
    canvas.setFillColor(colors.HexColor("#536879"))
    canvas.drawString(43, 28, "ECDLP index-calculus lead assessment | 2026-10-12")
    canvas.drawRightString(PAGE_W - 43, 28, str(doc.page))
    canvas.restoreState()


lines = (HERE / "REPORT.md").read_text(encoding="utf-8").splitlines()
story = []
paragraph = []
idx = 0


def flush():
    if paragraph:
        story.append(Paragraph(inline(" ".join(paragraph)), body))
        paragraph.clear()


while idx < len(lines):
    line = lines[idx].strip()
    if not line:
        flush()
        idx += 1
        continue
    if line.startswith("# "):
        flush()
        story.append(Paragraph(inline(line[2:]), title))
        idx += 1
        continue
    if line.startswith("## "):
        flush()
        story.append(Paragraph(inline(line[3:]), heading))
        idx += 1
        continue
    if line.startswith("|"):
        flush()
        rows = []
        while idx < len(lines) and lines[idx].lstrip().startswith("|"):
            parts = [part.strip() for part in lines[idx].strip().strip("|").split("|")]
            if not all(re.fullmatch(r":?-{2,}:?", part) for part in parts):
                rows.append([Paragraph(inline(part), cell) for part in parts])
            idx += 1
        widths = [142, 183, 201]
        table = Table(rows, colWidths=widths, repeatRows=1, hAlign="LEFT")
        table.setStyle(TableStyle([
            ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#dce9f5")),
            ("GRID", (0, 0), (-1, -1), 0.3, colors.HexColor("#9bb0c2")),
            ("VALIGN", (0, 0), (-1, -1), "TOP"),
            ("LEFTPADDING", (0, 0), (-1, -1), 5),
            ("RIGHTPADDING", (0, 0), (-1, -1), 5),
            ("TOPPADDING", (0, 0), (-1, -1), 3),
            ("BOTTOMPADDING", (0, 0), (-1, -1), 3),
        ]))
        story.extend([table, Spacer(1, 7)])
        continue
    if re.match(r"^[1-5]\. ", line):
        flush()
        story.append(Paragraph(inline(line), body))
        idx += 1
        continue
    paragraph.append(line)
    if "lead-map.svg" in line:
        flush()
        story.extend([
            LeadMap(),
        ])
    idx += 1
flush()

doc = SimpleDocTemplate(
    str(HERE / "REPORT.pdf"), pagesize=(612, 792),
    leftMargin=43, rightMargin=43, topMargin=43, bottomMargin=44,
    title="ECDLP index-calculus leads in the conversation export",
    author="Research assessment",
)
doc.build(story, onFirstPage=footer, onLaterPages=footer)
