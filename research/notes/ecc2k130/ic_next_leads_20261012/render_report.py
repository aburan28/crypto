"""Render the adjacent recorded search report and lead diagram. No research computation."""

from html import escape, unescape
from pathlib import Path
import re

from reportlab.lib import colors
from reportlab.lib.styles import ParagraphStyle
from reportlab.platypus import Flowable, Paragraph, SimpleDocTemplate, Spacer, Table, TableStyle


HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[3]
WIDTH = 526
BODY = ParagraphStyle("body", fontName="Helvetica", fontSize=8.4, leading=11.2, spaceAfter=3)
TITLE = ParagraphStyle("title", parent=BODY, fontName="Helvetica-Bold", fontSize=16, leading=20, spaceAfter=10)
HEAD = ParagraphStyle("head", parent=BODY, fontName="Helvetica-Bold", fontSize=11, leading=14, spaceBefore=9, spaceAfter=5)
CELL = ParagraphStyle("cell", parent=BODY, fontSize=7.2, leading=9.3, spaceAfter=0, splitLongWords=1)


def inline(value):
    value = escape(value)

    def link(match):
        label, target = match.groups()
        target = unescape(target)
        if target.startswith("https://"):
            url = target
        else:
            path = (HERE / target).resolve().relative_to(REPO_ROOT)
            url = f"https://github.com/aburan28/crypto/blob/codex/ic-next-leads-20261012/{path.as_posix()}"
        return f'<link href="{escape(url, quote=True)}" color="#205b91">{label}</link>'

    value = re.sub(r"\[([^\]]+)\]\(([^)]+)\)", link, value)
    value = re.sub(r"\*\*([^*]+)\*\*", r"<b>\1</b>", value)
    value = re.sub(r"`([^`]+)`", r'<font name="Courier">\1</font>', value)
    return value


class Diagram(Flowable):
    def __init__(self):
        super().__init__()
        self.width = WIDTH
        self.height = 132

    def draw(self):
        c = self.canv
        c.setFillColor(colors.HexColor("#f7fafc"))
        c.roundRect(0, 0, self.width, self.height, 7, stroke=0, fill=1)
        c.setFillColor(colors.HexColor("#213648"))
        c.setFont("Helvetica-Bold", 9)
        c.drawString(10, 117, "Where each proposal enters the charged pipeline")
        boxes = [(8, "Public Q"), (112, "Base / index"), (216, "Relations"), (320, "Useful rank"), (424, "Target PDP")]
        for i, (x, label) in enumerate(boxes):
            c.setFillColor(colors.HexColor("#e8f2fc"))
            c.setStrokeColor(colors.HexColor("#315e87"))
            c.roundRect(x, 54, 94, 28, 4, stroke=1, fill=1)
            c.setFillColor(colors.HexColor("#213648"))
            c.setFont("Helvetica", 7)
            c.drawCentredString(x + 47, 65, label)
            if i < 4:
                c.setStrokeColor(colors.HexColor("#315e87"))
                c.line(x + 95, 68, x + 103, 68)
        for x, label, destination in [(112, "Affine bases", 159), (216, "Joint families", 263), (320, "Large-prime", 367)]:
            c.setFillColor(colors.HexColor("#fff6df"))
            c.setStrokeColor(colors.HexColor("#a36a12"))
            c.roundRect(x, 91, 94, 17, 3, stroke=1, fill=1)
            c.setFillColor(colors.HexColor("#213648"))
            c.setFont("Helvetica", 6.5)
            c.drawCentredString(x + 47, 97, label)
            c.line(destination, 91, destination, 83)
        for x, label, destination in [(112, "N53 Bloom", 159), (216, "Prime template", 263)]:
            c.setFillColor(colors.HexColor("#fff6df"))
            c.setStrokeColor(colors.HexColor("#a36a12"))
            c.roundRect(x, 19, 94, 17, 3, stroke=1, fill=1)
            c.setFillColor(colors.HexColor("#213648"))
            c.setFont("Helvetica", 6.5)
            c.drawCentredString(x + 47, 25, label)
            c.line(destination, 36, destination, 53)
        c.setFillColor(colors.HexColor("#52687a"))
        c.setFont("Helvetica", 6.8)
        c.drawString(319, 23, "Proposals; compare cold S with same-Q rho")


def footer(canvas, doc):
    canvas.saveState()
    canvas.setFont("Helvetica", 7.5)
    canvas.setFillColor(colors.HexColor("#52687a"))
    canvas.drawString(43, 27, "IC follow-on leads | 2026-10-12 | search and design only")
    canvas.drawRightString(569, 27, str(doc.page))
    canvas.restoreState()


lines = (HERE / "REPORT.md").read_text(encoding="utf-8").splitlines()
story = []
paragraph = []


def flush():
    if paragraph:
        story.append(Paragraph(inline(" ".join(paragraph)), BODY))
        paragraph.clear()


i = 0
while i < len(lines):
    line = lines[i].strip()
    if not line:
        flush()
        i += 1
        continue
    if line.startswith("# "):
        flush()
        story.append(Paragraph(inline(line[2:]), TITLE))
        i += 1
        continue
    if line.startswith("## "):
        flush()
        story.append(Paragraph(inline(line[3:]), HEAD))
        i += 1
        continue
    if line.startswith("|"):
        flush()
        rows = []
        while i < len(lines) and lines[i].lstrip().startswith("|"):
            cells = [s.strip() for s in lines[i].strip().strip("|").split("|")]
            if not all(re.fullmatch(r":?-{2,}:?", s) for s in cells):
                rows.append([Paragraph(inline(s), CELL) for s in cells])
            i += 1
        widths = [WIDTH / len(rows[0])] * len(rows[0])
        table = Table(rows, colWidths=widths, repeatRows=1, hAlign="LEFT")
        table.setStyle(TableStyle([
            ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#dce9f5")),
            ("GRID", (0, 0), (-1, -1), 0.3, colors.HexColor("#9bb0c2")),
            ("VALIGN", (0, 0), (-1, -1), "TOP"),
            ("LEFTPADDING", (0, 0), (-1, -1), 4),
            ("RIGHTPADDING", (0, 0), (-1, -1), 4),
            ("TOPPADDING", (0, 0), (-1, -1), 3),
            ("BOTTOMPADDING", (0, 0), (-1, -1), 3),
        ]))
        story.extend([table, Spacer(1, 6)])
        continue
    paragraph.append(line)
    if "lead-gates.svg" in line:
        flush()
        story.append(Diagram())
    i += 1
flush()

doc = SimpleDocTemplate(
    str(HERE / "REPORT.pdf"), pagesize=(612, 792),
    leftMargin=43, rightMargin=43, topMargin=43, bottomMargin=44,
    title="ECDLP index-calculus follow-on leads and experiment gates",
    author="Research assessment",
)
doc.build(story, onFirstPage=footer, onLaterPages=footer)
