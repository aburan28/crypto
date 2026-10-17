"""Layout and plot recorded native results. No solver, benchmarking or statistics."""
from pathlib import Path
from xml.sax.saxutils import escape
import json
import re
import shutil
import subprocess

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from reportlab.graphics import renderSVG, renderPDF
from reportlab.graphics.shapes import Drawing, Rect, String, Line
from reportlab.lib import colors
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image, Table, TableStyle

STUDY = Path(__file__).resolve().parent
RUN = "run-20261016T0205Z"
DATA = STUDY / "evidence" / RUN
BROWSE = STUDY.parents[1] / "benchmarks" / "20261016-sea-point-count"
FONT_DIR = Path("/usr/share/fonts/truetype/liberation")
pdfmetrics.registerFont(TTFont("ReportSans", str(FONT_DIR / "LiberationSans-Regular.ttf")))
pdfmetrics.registerFont(TTFont("ReportSans-Bold", str(FONT_DIR / "LiberationSans-Bold.ttf")))
pdfmetrics.registerFontFamily("ReportSans", normal="ReportSans", bold="ReportSans-Bold")


def plots(rows, output):
    plt.rcParams.update({"svg.fonttype": "none", "font.size": 10})
    fig, ax = plt.subplots(figsize=(9, 4.8))
    fig.subplots_adjust(left=.20, right=.98, top=.90, bottom=.20)
    for i, row in enumerate(rows):
        noise = row["aa_max_deviation"]
        ax.fill_betweenx([i-.32, i+.32], 1-noise, 1+noise, color="#d9e2ed", alpha=.8)
        ax.plot([row["ratio_min"], row["ratio_max"]], [i, i], color="#22618c", lw=2)
        ax.plot(row["ratio_median"], i, "o", color="#15384f", ms=6)
    ax.axvline(1, color="#555", ls="--", lw=1)
    ax.set_yticks(range(len(rows)), [r["case"] for r in rows])
    ax.invert_yaxis()
    ax.set_xlabel("Baseline / candidate elapsed time per call (above 1: candidate faster)")
    ax.set_title("Matched point counts: paired median and observed range")
    ax.grid(axis="x", alpha=.2)
    ax.legend(handles=[Patch(color="#d9e2ed", label="A/A maximum paired-deviation envelope")], loc="best", fontsize=8)
    fig.text(.99, .005, RUN + "; five A/B pairs per workload; pinned VM, PSI unavailable",
             ha="right", va="bottom", fontsize=7)
    fig.savefig(output / "speed-ratios.svg")
    fig.savefig(output / "speed-ratios.png", dpi=150)
    plt.close(fig)


def proof_diagram(output):
    d = Drawing(490, 142)
    navy = colors.HexColor("#15384f")
    def box(x, y, w, h, lines, color):
        d.add(Rect(x,y,w,h,rx=4,ry=4,fillColor=colors.HexColor(color),strokeColor=navy))
        for i,t in enumerate(lines):
            d.add(String(x+w/2,y+h/2+(len(lines)-1)*6-i*12,t,
                         fontName="ReportSans",fontSize=9,textAnchor="middle",fillColor=navy))
    box(0,86,225,48,["TLA+: finite transition model", "complete hits, ordering, refinement, return"],"#e6eff8")
    box(265,86,225,48,["Lean: unbounded abstract proofs", "group search, preservation, interval uniqueness"],"#e3f2e9")
    box(105,2,280,48,["Open: concrete Rust / arithmetic refinement", "group laws, trace constraints, divisor/primality bridge"],"#fff0d6")
    d.add(Line(112,86,245,50,strokeColor=navy,dashArray=[3,3]))
    d.add(Line(377,86,245,50,strokeColor=navy,dashArray=[3,3]))
    renderSVG.drawToFile(d,str(output / "verification-scope.svg"))
    # Both the PNG browsing copy and embedded PDF use this exact vector drawing.
    renderPDF.drawToFile(d,str(output / "verification-scope.pdf"))
    subprocess.run(["pdftoppm","-singlefile","-r","150","-png",
                    str(output / "verification-scope.pdf"),str(output / "verification-scope")],check=True)
    return d


def pdf(markdown, destination, graph, diagram=None):
    styles = getSampleStyleSheet()
    styles.add(ParagraphStyle(name="Body",fontName="ReportSans",fontSize=9,leading=12,spaceAfter=5))
    styles.add(ParagraphStyle(name="Cell",fontName="ReportSans",fontSize=7,leading=9))
    for h in ["Title","Heading1","Heading2","Heading3"]:
        styles[h].fontName="ReportSans-Bold"
    styles["Title"].fontSize=19
    styles["Title"].leading=23
    story=[]; paragraph=[]; table=[]
    def flush():
        if paragraph:
            t=escape(" ".join(paragraph))
            t=re.sub(r"\*\*(.*?)\*\*",r"<b>\1</b>",t)
            t=re.sub(r"`([^`]+)`",r"<font name='ReportSans'>\1</font>",t)
            t=re.sub(r"\[([^]]+)\]\(([^)]+)\)",r"<link href='\2'>\1</link>",t)
            story.append(Paragraph(t,styles["Body"]));paragraph.clear()
        if table:
            cells=[[Paragraph(escape(c),styles["Cell"]) for c in row] for row in table]
            widths=[490/len(table[0])]*len(table[0])
            t=Table(cells,colWidths=widths,repeatRows=1,hAlign="LEFT")
            t.setStyle(TableStyle([("VALIGN",(0,0),(-1,-1),"TOP"),("BACKGROUND",(0,0),(-1,0),colors.HexColor("#e6eff8")),("LINEBELOW",(0,0),(-1,0),.5,colors.grey),("BOTTOMPADDING",(0,0),(-1,-1),5)]))
            story.extend([t,Spacer(1,10)]);table.clear()
    for line in markdown.read_text().splitlines():
        if line.startswith("| "):
            if not table: flush()
            row=[c.strip() for c in line.strip().strip("|").split("|")]
            if not all(re.fullmatch(r"[-: ]+", c) for c in row): table.append(row)
        elif line.startswith("# "):
            flush();story.append(Paragraph(escape(line[2:]),styles["Title"]))
            story.extend([Image(str(graph),width=490,height=261),Spacer(1,10)])
            if diagram is not None:story.extend([diagram,Spacer(1,10)])
        elif line.startswith("## "):
            flush();story.append(Paragraph(escape(line[3:]),styles["Heading2"]))
        elif line.startswith("### "):
            flush();story.append(Paragraph(escape(line[4:]),styles["Heading3"]))
        elif not line.strip():flush()
        else:paragraph.append(line)
    flush()
    def footer(canvas,doc):
        canvas.setFont("ReportSans",8)
        canvas.drawString(45,22,"SEA measurements and formal scope - 2026-10-16")
        canvas.drawRightString(570,22,str(doc.page))
    SimpleDocTemplate(str(destination),pagesize=(612,792),leftMargin=50,rightMargin=50,topMargin=35,bottomMargin=38).build(story,onFirstPage=footer,onLaterPages=footer)


def main():
    output=STUDY / "reports" / RUN
    output.mkdir(parents=True,exist_ok=True)
    rows=json.loads((DATA / "summary.json").read_text())
    plots(rows,output)
    drawing=proof_diagram(output)
    graph=output / "speed-ratios.png"
    pdf(output / "REPORT.md",output / "report.pdf",graph)
    pdf(STUDY / "REPORT.md",STUDY / "REPORT.pdf",graph,drawing)
    BROWSE.mkdir(parents=True,exist_ok=True)
    for src,name in [(output/"report.pdf","run-report.pdf"),(STUDY/"REPORT.pdf","results-report.pdf"),
                     (output/"speed-ratios.svg","speed-ratios.svg"),(graph,"speed-ratios.png"),
                     (output/"verification-scope.svg","verification-scope.svg")]:
        shutil.copy2(src,BROWSE/name)
    shutil.copy2(output/"verification-scope.png",BROWSE/"verification-scope.png")


if __name__ == "__main__":main()
