#!/usr/bin/env python3
"""Lay out RESULTS.md + fig_yield.png into report.pdf (no computation)."""
import re
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image, Preformatted
from reportlab.lib.units import cm
ss = getSampleStyleSheet()
doc = SimpleDocTemplate("report.pdf", pagesize=A4, leftMargin=1.8*cm, rightMargin=1.8*cm,
                        title="06-pkm-rational-map report")
story, para, table, code = [], [], [], False
def esc(s):
    s = s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
    s = re.sub(r"\*\*(.+?)\*\*", r"<b>\1</b>", s)
    return re.sub(r"`(.+?)`", r"<font face='Courier'>\1</font>", s)
def flush():
    global para
    if para:
        story.append(Paragraph(esc(" ".join(para)), ss["BodyText"])); para = []
for line in open("RESULTS.md").read().splitlines():
    if line.startswith("```"):
        flush(); code = not code; continue
    if code:
        story.append(Preformatted(line, ss["Code"])); continue
    if line.startswith("|"):
        flush(); table.append(line); continue
    if table:
        story.append(Preformatted("\n".join(t for t in table if not set(t) <= set("|-")), ss["Code"], maxLineLength=200))
        story.append(Spacer(1, 6)); table = []
    if line.startswith("#"):
        flush(); lvl = len(line) - len(line.lstrip("#"))
        story.append(Paragraph(esc(line.lstrip("# ")), ss["Heading%d" % min(lvl, 3)]))
        if line.startswith("## Figures"):
            story.append(Image("fig_pkm.png", width=17.4*cm, height=5.4*cm))
            story.append(Paragraph("Figure 1. Left: mean CPU time of one m=2 relation system "
                "(pure-Python Buchberger) per arm vs the full Pollard-rho DLP on the same curve. Middle: "
                "max lcm degree of a nonzero-reducing pair vs PKM Table 4. Right: projected IC relation "
                "stage / rho with 95% bootstrap CI.", ss["Italic"]))
            story.append(Image("fig_yield.png", width=14*cm, height=8.4*cm))
            story.append(Paragraph("Figure 2. Relation yield over heuristic 2|F|^2/N, pooled 12 targets, "
                "95% Garwood CI; planted known-log base is the positive control.", ss["Italic"]))
    elif line.startswith("- ") or re.match(r"^\d+\. ", line):
        flush(); para = [line]
    elif not line.strip():
        flush()
    else:
        para.append(line.strip())
flush()
doc.build(story)
