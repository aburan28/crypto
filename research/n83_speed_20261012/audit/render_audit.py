"""Document layout only: render written audit sources and their SVG diagram.

No algorithm, experiment, benchmark, statistical calculation or verifier runs
here. Numerical values and verification status come from cited native receipts.
"""
from pathlib import Path
import html
import re
import sys
from urllib.parse import urljoin

sys.path.insert(0, "/workspace/n83_speed/audit_layout")
import cairosvg
from reportlab.lib import colors
from reportlab.lib.enums import TA_LEFT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import inch
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import (
    SimpleDocTemplate, Paragraph, Spacer, PageBreak, Image, Table, TableStyle,
)

root = Path(__file__).resolve().parent
fonts = Path("/usr/share/fonts/truetype/dejavu")
for name, filename in [
    ("AuditSans", "DejaVuSans.ttf"),
    ("AuditBold", "DejaVuSans-Bold.ttf"),
    ("AuditMono", "DejaVuSansMono.ttf"),
]:
    pdfmetrics.registerFont(TTFont(name, str(fonts / filename)))
pdfmetrics.registerFontFamily("AuditSans", normal="AuditSans", bold="AuditBold")

svg = root / "arithmetic_proof_map.svg"
png = root / "arithmetic_proof_map.png"
cairosvg.svg2png(url=str(svg), write_to=str(png), dpi=150, scale=150 / 96)

styles = getSampleStyleSheet()
styles.add(ParagraphStyle("AuditBody", fontName="AuditSans", fontSize=9.3,
    leading=12.6, spaceAfter=6, wordWrap="CJK", alignment=TA_LEFT))
styles.add(ParagraphStyle("AuditTable", parent=styles["AuditBody"], fontSize=7.8,
    leading=10.3, spaceAfter=0))
styles.add(ParagraphStyle("AuditTitle", fontName="AuditBold", fontSize=22,
    leading=28, spaceAfter=14))
styles.add(ParagraphStyle("AuditH1", fontName="AuditBold", fontSize=16,
    leading=22, spaceAfter=11, keepWithNext=True))
styles.add(ParagraphStyle("AuditH2", fontName="AuditBold", fontSize=11,
    leading=16, spaceBefore=10, spaceAfter=7, keepWithNext=True))
base_url = "https://github.com/aburan28/crypto/blob/codex/n83-rho-speed-20261012/research/n83_speed_20261012/audit/"

def markup(text):
    tokens = []
    def link(match):
        label, target = match.group(1), match.group(2)
        url = target if target.startswith("http") else urljoin(base_url, target)
        token = f"LINKTOKEN{len(tokens)}END"
        tokens.append((token, f'<link href="{html.escape(url, quote=True)}" color="#23558c">{html.escape(label)}</link>'))
        return token
    text = re.sub(r"\[([^\]]+)\]\(([^)]+)\)", link, text)
    text = html.escape(text)
    text = re.sub(r"`([^`]+)`", r'<font name="AuditMono">\1</font>', text)
    for token, replacement in tokens:
        text = text.replace(token, replacement)
    return text

def paragraph(text, style="AuditBody"):
    return Paragraph(markup(text), styles[style])

width = A4[0] - 80
story = [Paragraph("Independent degree-83 scalar audit", styles["AuditTitle"]),
    paragraph("2026-10-12 · Exact arithmetic, formal mathematical statements and fair measurement protocol"),
    paragraph("The native campaign retained the reference scalar kernel and selected 256 lanes with DP6. Correctness of the Barrett candidate is separate from its measured performance; it did not improve arithmetic-stage process cost."),
    Image(str(png), width=width, height=width * 670 / 1100),
    Spacer(1, 10),
    paragraph("Sources and raw receipts are linked throughout. The four Lean declarations are checked mathematical statements; the compiled Rust program and word-operation bridges have the explicitly stated separate verification scope."),
    PageBreak()]

def append_markdown(path):
    lines = path.read_text().splitlines()
    pending = []
    def flush():
        if pending:
            story.append(paragraph(" ".join(pending)))
            pending.clear()
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if not line:
            flush()
        elif line.startswith("# "):
            flush()
            story.append(paragraph(line[2:], "AuditH1"))
        elif line.startswith("## "):
            flush()
            story.append(paragraph(line[3:], "AuditH2"))
        elif line.startswith("|"):
            flush()
            rows = []
            while i < len(lines) and lines[i].strip().startswith("|"):
                cells = [c.strip() for c in lines[i].strip().strip("|").split("|")]
                if not all(re.fullmatch(r":?-+:?", c) for c in cells):
                    rows.append([paragraph(c, "AuditTable") for c in cells])
                i += 1
            count = len(rows[0])
            fractions = [.19, .28, .34, .19] if count == 4 else [.22, .52, .26]
            table = Table(rows, colWidths=[width * f for f in fractions], repeatRows=1,
                hAlign="LEFT")
            table.setStyle(TableStyle([
                ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#dce7f0")),
                ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#f3f6f9")]),
                ("VALIGN", (0, 0), (-1, -1), "TOP"),
                ("LEFTPADDING", (0, 0), (-1, -1), 6),
                ("RIGHTPADDING", (0, 0), (-1, -1), 6),
                ("TOPPADDING", (0, 0), (-1, -1), 5),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 5),
                ("LINEBELOW", (0, 0), (-1, 0), .5, colors.HexColor("#8ea2b6")),
            ]))
            story.extend([table, Spacer(1, 10)])
            continue
        elif line.startswith("- "):
            flush()
            story.append(paragraph("• " + line[2:]))
        else:
            pending.append(line)
        i += 1
    flush()

for name in ["ARITHMETIC_PROOF.md", "PROOF_MAP.md", "PROTOCOL_AUDIT.md"]:
    append_markdown(root / name)
    if name != "PROTOCOL_AUDIT.md":
        story.append(PageBreak())

def footer(canvas, doc):
    canvas.setFont("AuditSans", 7)
    canvas.setFillColor(colors.HexColor("#5b6e7f"))
    canvas.drawString(40, 24, "N83 scalar audit · source-linked correctness and measurement scope")
    canvas.drawRightString(A4[0] - 40, 24, str(doc.page))

pdf = root / "independent_scalar_audit.pdf"
doc = SimpleDocTemplate(str(pdf), pagesize=A4, leftMargin=40, rightMargin=40,
    topMargin=38, bottomMargin=42,
    title="Independent degree-83 scalar arithmetic and protocol audit",
    author="Research audit")
doc.build(story, onFirstPage=footer, onLaterPages=footer)
print(str(pdf))
