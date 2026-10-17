"""Render recorded prose and a control-flow diagram; no research computation."""
from pathlib import Path
from xml.sax.saxutils import escape
import re

from reportlab.graphics import renderSVG
from reportlab.graphics.shapes import Drawing, Rect, String, Line, Polygon
from reportlab.lib import colors
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont


ROOT = Path(__file__).resolve().parent
FONT_DIR = Path("/usr/share/fonts/truetype/liberation")
pdfmetrics.registerFont(TTFont("ReportSans", str(FONT_DIR / "LiberationSans-Regular.ttf")))
pdfmetrics.registerFont(TTFont("ReportSans-Bold", str(FONT_DIR / "LiberationSans-Bold.ttf")))
pdfmetrics.registerFontFamily("ReportSans", normal="ReportSans", bold="ReportSans-Bold")


def diagram():
    d = Drawing(490, 232)
    navy = colors.HexColor("#16364d")
    blue = colors.HexColor("#e8f1fa")
    green = colors.HexColor("#e4f4ed")

    def box(x, y, w, h, lines, fill):
        d.add(Rect(x, y, w, h, rx=5, ry=5, fillColor=fill, strokeColor=navy))
        for i, text in enumerate(lines):
            d.add(String(x+w/2, y+h/2+(len(lines)-1)*6-i*12, text,
                         fontName="ReportSans", fontSize=9, textAnchor="middle", fillColor=navy))

    def arrow(x1, y1, x2, y2):
        d.add(Line(x1, y1, x2, y2, strokeColor=navy))
        if y1 == y2:
            d.add(Polygon([x2, y2, x2-5, y2+3, x2-5, y2-3], fillColor=navy, strokeColor=navy))
        else:
            d.add(Polygon([x2, y2, x2-3, y2+5, x2+3, y2+5], fillColor=navy, strokeColor=navy))

    d.add(String(0, 219, "Implemented point-count changes", fontName="ReportSans-Bold", fontSize=11, fillColor=navy))
    box(0, 155, 140, 46, ["Complete walk / BSGS", "keep every matching index"], blue)
    box(174, 155, 138, 46, ["Surviving candidates", "sorted intersection"], blue)
    box(346, 155, 144, 46, ["One candidate", "public order certification"], green)
    arrow(140, 178, 174, 178)
    arrow(312, 178, 346, 178)
    box(174, 86, 220, 43, ["If ambiguous: sample another point", "choose scalar tests or full search"], blue)
    arrow(243, 155, 243, 129)
    d.add(String(0, 60, "Existing CM generation (specialized research constructor)",
                 fontName="ReportSans-Bold", fontSize=11, fillColor=navy))
    box(0, 0, 140, 43, ["Hilbert class polynomial", "fixed-point complex evaluation"], green)
    box(174, 0, 138, 43, ["Roots modulo p", "build curve / quadratic twist"], green)
    box(346, 0, 144, 43, ["Sampled order checks", "certify separately in api"], green)
    arrow(140, 21, 174, 21)
    arrow(312, 21, 346, 21)
    return d


def main():
    styles = getSampleStyleSheet()
    styles.add(ParagraphStyle(name="ReportBody", fontName="ReportSans", fontSize=9, leading=12, spaceAfter=5,
                              splitLongWords=True))
    styles["Title"].fontName = "ReportSans-Bold"
    styles["Heading2"].fontName = "ReportSans-Bold"
    styles["Title"].fontSize = 19
    styles["Title"].leading = 23
    styles["Heading2"].fontSize = 12
    styles["Heading2"].leading = 15
    flow = diagram()
    renderSVG.drawToFile(flow, str(ROOT / "flow.svg"))
    story = []
    pending = []

    def paragraph():
        if not pending:
            return
        text = escape(" ".join(pending))
        text = re.sub(r"\*\*(.*?)\*\*", r"<b>\1</b>", text)
        story.append(Paragraph(text, styles["ReportBody"]))
        pending.clear()

    for line in (ROOT / "REPORT.md").read_text().splitlines():
        if line.startswith("# "):
            paragraph()
            story.append(Paragraph(escape(line[2:]), styles["Title"]))
            story.append(flow)
            story.append(Spacer(1, 12))
        elif line.startswith("## "):
            paragraph()
            story.append(Paragraph(escape(line[3:]), styles["Heading2"]))
        elif not line.strip():
            paragraph()
        elif line.startswith("- "):
            paragraph()
            pending.append(line)
            paragraph()
        else:
            pending.append(line)
    paragraph()

    def footer(canvas, doc):
        canvas.setFont("ReportSans", 8)
        canvas.drawString(40, 24, "SEA / CM review - 2026-10-15 - correctness evidence only")
        canvas.drawRightString(572, 24, str(doc.page))

    SimpleDocTemplate(str(ROOT / "report.pdf"), pagesize=(612, 792),
                      leftMargin=50, rightMargin=50, topMargin=36, bottomMargin=40).build(
                          story, onFirstPage=footer, onLaterPages=footer)


if __name__ == "__main__":
    main()
