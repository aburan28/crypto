"""Layout only (copied from ../04-family-sweep/make_pdf.py): render RESULTS.md + recorded figures into report.pdf (reportlab)."""
import re, sys
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import cm
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image, Table, TableStyle, Preformatted
from reportlab.lib import colors

d = sys.argv[1] if len(sys.argv) > 1 else "."
ss = getSampleStyleSheet()
body = ParagraphStyle("b", parent=ss["BodyText"], fontSize=8.5, leading=11)
cell = ParagraphStyle("c", parent=body, fontSize=6.5, leading=8)
code = ParagraphStyle("code", parent=ss["Code"], fontSize=6.5, leading=8)

def inline(s):
    s = s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
    s = re.sub(r"\*\*(.+?)\*\*", r"<b>\1</b>", s)
    s = re.sub(r"`(.+?)`", r"<font face='Courier'>\1</font>", s)
    return s

story = []
lines = open(f"{d}/RESULTS.md").read().splitlines()
i = 0
while i < len(lines):
    l = lines[i]
    if l.startswith("```"):
        j = i + 1; buf = []
        while not lines[j].startswith("```"): buf.append(lines[j]); j += 1
        story.append(Preformatted("\n".join(buf), code)); i = j + 1; continue
    if l.startswith("|"):
        tb = []
        while i < len(lines) and lines[i].startswith("|"):
            if not re.match(r"^\|[\s\-:|]+\|$", lines[i]):
                tb.append([Paragraph(inline(c.strip()), cell) for c in lines[i].strip("|").split("|")])
            i += 1
        t = Table(tb, repeatRows=1)
        t.setStyle(TableStyle([("GRID", (0, 0), (-1, -1), 0.3, colors.grey), ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#eeeeee")), ("VALIGN", (0, 0), (-1, -1), "TOP")]))
        story += [t, Spacer(1, 6)]; continue
    m = re.match(r"!\[(.*)\]\((.+\.)svg\)", l)
    if m:
        img = Image(f"{d}/{m.group(2)}png"); w = 17 * cm; img.drawHeight = img.drawHeight * w / img.drawWidth; img.drawWidth = w
        story += [img, Paragraph(inline(m.group(1)), cell), Spacer(1, 6)]; i += 1; continue
    if l.startswith("#"):
        lvl = len(l) - len(l.lstrip("#")); story.append(Paragraph(inline(l.lstrip("# ")), ss["Heading%d" % min(lvl + 0, 4)]))
    elif l.strip().startswith(("- ", "* ")):
        story.append(Paragraph("&bull; " + inline(l.strip()[2:]), body))
    elif l.strip():
        story.append(Paragraph(inline(l), body))
    i += 1
SimpleDocTemplate(f"{d}/report.pdf", pagesize=A4, leftMargin=1.6 * cm, rightMargin=1.6 * cm, topMargin=1.5 * cm, bottomMargin=1.5 * cm, title="05-qr-window report").build(story)
print("report.pdf written")
