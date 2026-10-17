#!/usr/bin/env python3
"""Lay out RESULTS.md + fig_coverage.png into report.pdf (no computation)."""
import re
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image, Preformatted
from reportlab.lib.units import cm
ss = getSampleStyleSheet()
doc = SimpleDocTemplate("report.pdf", pagesize=A4, leftMargin=1.8*cm, rightMargin=1.8*cm,
                        title="12-model-small-x report")
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
        if line.startswith("## Verdict"):
            pass
    elif line.startswith("- ") or re.match(r"^\d+\. ", line):
        flush(); para = [line]
    elif not line.strip():
        flush()
    else:
        para.append(line.strip())
flush()
story.insert(3, Image("fig_coverage.png", width=17.5*cm, height=6.5*cm))
story.insert(4, Paragraph("Figure 1. Relative m=2 coverage difference |F+F|/N of each arm against its matched random control "
    "(negation-closed arms vs R-sym, E-x and H-x vs R-asym), mean over curves with Bonferroni (21 tests) bootstrap CI; "
    "grey band is the preregistered +-1% equivalence margin. Right: positive controls PC1/PC2. Data: results.json.",
    ss["Italic"]))
doc.build(story)
