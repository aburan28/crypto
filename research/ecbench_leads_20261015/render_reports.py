"""Layout only: all outcomes, counts, summaries and ratios come from Rust."""
import os
os.environ.setdefault("MPLCONFIGDIR", "/tmp/ecbench-leads-mpl")
import json
import shutil
from pathlib import Path
from xml.sax.saxutils import escape
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image, PageBreak, LongTable, TableStyle

ROOT=Path(__file__).resolve().parent
DATA=json.loads((ROOT/"render-data.json").read_text())
FINDINGS=json.loads((ROOT/"findings.json").read_text())
BENCH=ROOT.parents[1]/"benchmarks/20261015-ecbench-leads"
FONT=Path("/usr/share/fonts/truetype/dejavu")
for name,file in [("DV","DejaVuSans.ttf"),("DVB","DejaVuSans-Bold.ttf")]:
    pdfmetrics.registerFont(TTFont(name,str(FONT/file)))
pdfmetrics.registerFontFamily("DV",normal="DV",bold="DVB",italic="DV")
NAVY,TEAL,GRAY="#17354a","#187f83","#586775"
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":8,"svg.fonttype":"none","axes.spines.top":False,"axes.spines.right":False})
BODY=ParagraphStyle("body",fontName="DV",fontSize=8.3,leading=11.3,spaceAfter=7,textColor=colors.HexColor(NAVY))
SMALL=ParagraphStyle("small",parent=BODY,fontSize=6.8,leading=9)
TITLE=ParagraphStyle("title",parent=BODY,fontName="DVB",fontSize=17,leading=21,spaceAfter=12)
HEAD=ParagraphStyle("head",parent=BODY,fontName="DVB",fontSize=11,leading=14,spaceBefore=8)
def p(text,style=BODY): return Paragraph(escape(str(text)),style)
def footer(c,doc):
    c.setFont("DV",6.5)
    c.drawString(36,20,"Native ecbench • 2026-10-15 • bounded counted diagnostics • no full-cost claim")
    c.drawRightString(A4[0]-36,20,str(doc.page))
def pdf(path,story):
    SimpleDocTemplate(str(path),pagesize=A4,topMargin=34,bottomMargin=36,leftMargin=36,rightMargin=36).build(story,onFirstPage=footer,onLaterPages=footer)
def figure(path,max_height=460):
    item=Image(str(path))
    item._restrictSize(510,max_height)
    return item
def figsave(fig,path):
    fig.savefig(path.with_suffix(".svg"),bbox_inches="tight")
    fig.savefig(path.with_suffix(".png"),bbox_inches="tight",dpi=150)
    svg=path.with_suffix(".svg")
    svg.write_text("\n".join(s.rstrip() for s in svg.read_text().splitlines())+"\n")
    plt.close(fig)

def shape_graph(path):
    rows=[r for r in DATA["shapes"]["rows"] if r["m"]==3]
    fig,axes=plt.subplots(1,2,figsize=(8.1,2.55))
    for ax,old,new,label in [(axes[0],"legacy_vars","constructed_vars","Boolean variables"),(axes[1],"legacy_equations","constructed_equations","Coordinate equations")]:
        ax.plot([r["n"] for r in rows],[r[old] for r in rows],"o--",color=GRAY,label="Legacy metadata")
        ax.plot([r["n"] for r in rows],[r[new] for r in rows],"o-",color=TEAL,label="Constructed / v2")
        ax.set_xlabel("Extension degree n; ell=3, m=3")
        ax.set_ylabel(label)
        ax.legend(fontsize=7)
        ax.grid(alpha=.18)
    fig.tight_layout()
    figsave(fig,path)

def plot_session(s):
    path=ROOT/s["path"]
    rows=s["rows"]
    if s["session"]=="binary-shape":
        shape_graph(path/"cost")
    elif s["session"]=="wide-n83":
        fig,ax=plt.subplots(figsize=(8.1,2.8))
        for i,r in enumerate(rows):
            c=r["counter_sums"]
            ax.text(.03,.72-i*.36,f'{r["arm"]}: {r["verified"]}/{r["records"]} recovered; '
                f'{c.get("trials",0)} total queries; {c.get("point_checked_relations",0)} checked relations; '
                f'{c.get("table_entries",0)} total stored entries',transform=ax.transAxes,fontsize=10,
                bbox={"boxstyle":"round,pad=.7","facecolor":"#e7eff5","edgecolor":"#9baebb"})
        ax.text(.03,.06,"Exact E0 GF(2^83) → unknown-log base → quotient pair table → bounded target pin search",transform=ax.transAxes)
        ax.axis("off")
        figsave(fig,path/"cost")
    elif not any(r["mean_s"] is not None for r in rows):
        fig,ax=plt.subplots(figsize=(8.1,2.6))
        ax.text(.05,.6,"Legacy counted pipeline → wide dispatch refusal",transform=ax.transAxes,fontsize=15)
        ax.text(.05,.25,"One preserved implementation error; no scalar or cost comparison",transform=ax.transAxes,fontsize=11)
        ax.axis("off")
        figsave(fig,path/"cost")
    else:
        fig,ax=plt.subplots(figsize=(8.1,max(2.6,len(rows)*.24+1)))
        for i,r in enumerate(rows):
            if r["mean_s"] is None:
                ax.text(.985,i,"censored / no verified recovery",ha="right",va="center",transform=ax.get_yaxis_transform(),fontsize=6.5,color=GRAY)
            else:
                color=NAVY if r["arm"].startswith("rho") else TEAL
                ax.scatter(r["mean_s"],i,s=20,color=color)
                ax.hlines(i,r["min_s"],r["max_s"],color=color,linewidth=1)
        ax.set_yticks(range(len(rows)),[f'{r["curve"]}\n{r["arm"]}  {r["verified"]}/{r["records"]}' for r in rows],fontsize=6.5)
        ax.set_ylim(len(rows)-.5,-.5)
        ax.set_xscale("log")
        ax.set_xlabel("Counted S = GAE / sqrt(r); native work remains unpriced")
        ax.grid(axis="x",alpha=.2)
        ax.set_title(s["ledger"]["session_id"],loc="left",fontsize=10)
        figsave(fig,path/"cost")
    ledger,meta=s["ledger"],s["metadata"]
    story=[p(s["session"]+": ecbench follow-up",TITLE),p("Session "+ledger["session_id"],SMALL),
        p("Source "+str(meta["git_commit"]),SMALL),p("Binary "+meta["binary_sha256"],SMALL),
        p(f'{ledger["records"]} attempts; {ledger["verified_records"]} verified recoveries; {ledger["replays"]} exact standard replays.'),
        p("Fresh public single-target workloads, seed20261015. Cold setup repeated; no planted scalar. Host "+meta["env_class_id"]+", L0; wall times inadmissible. Exhausted/error outcomes have a separate native replay."),
        p("Points and ranges use verified runs only: native mean and observed min/max, not confidence intervals. Failed work remains in the resource vectors. Count quotients do not bound full-cost ratios. Shape panels plot independently constructed coordinate counts; wide panels show recorded bounded outcomes.",SMALL),
        figure(path/"cost.png"),
        p("Evidence: research/ecbench_leads_20261015/sessions/"+s["session"]+"; "+s["session"]+".audit.json; failed-replays-v2.json; "+s["session"]+".resources.json.",SMALL)]
    pdf(path/"REPORT.pdf",story)
    dest=BENCH/s["session"]
    dest.mkdir(parents=True,exist_ok=True)
    for name in ["cost.svg","cost.png","REPORT.pdf"]: shutil.copyfile(path/name,dest/name)

BENCH.mkdir(parents=True,exist_ok=True)
for s in DATA["sessions"]: plot_session(s)
shape_graph(ROOT/"shape-correction")
res=DATA["results"]
story=[p("Pursuing the ecbench leads",TITLE),
    p(f'{res["attempts"]} fresh public single-target attempts; {res["verified"]} verified scalar recoveries. Outcomes: '+str(res["statuses"])+"."),
    p(FINDINGS["paragraphs"][1]),
    p(FINDINGS["paragraphs"][2]),
    p("All session audits passed; unsuccessful deterministic outcomes replayed separately. The new method uses unknown factor-base logs and exact point witnesses at both integer widths. Versioned metadata includes all intermediate coordinates and equation blocks."),
    p("The exact E0 degree83 subgroup is retained. Bounded execution does not satisfy the N83 recovery gate. All ten mathematical questions remain partial. L0; no admitted wall time, fully priced speedup, asymptotic exponent or N131 transfer."),
    p("System shape: independently constructed inputs",HEAD),figure(ROOT/"shape-correction.png",175),
    p("ell=3, m=3. The plot uses native constructed counts and preserved legacy metadata; it is not a solving-degree or solver-cost claim.",SMALL),
    p("Sources: PROTOCOL-v1.md; FINDINGS.md/findings.json; results.json; shape-controls.json; per-target.json; execution-manifest.json; parent literature/REPORT.md (29 references). The audit-reader amendment preserves the initial N83 reader failure and validates its correction against a pinned historical reader. See REPORT.md for source-linked run reports and exact remaining obligations.",SMALL),
    Paragraph('<link href="https://github.com/aburan28/crypto/blob/research/ecbench-ten-questions-20261014/research/ecbench_leads_20261015/REPORT.md">Source-linked report and immutable evidence</link>',SMALL),PageBreak(),p("Measured per-curve diagnostics",HEAD)]
cells=[[p(x,SMALL) for x in ["Curve","Arm","Verified","Mean S","Count / rho"]]]
for r in DATA["rows"]:
    cells.append([p(r["curve"],SMALL),p(r["arm"],SMALL),p(f'{r["verified"]}/{r["records"]}',SMALL),
        p("—" if r["mean_s"] is None else f'{r["mean_s"]:.5f}',SMALL),
        p("—" if r["paired_complete_counted_ratio_candidate_over_rho"] is None else f'{r["paired_complete_counted_ratio_candidate_over_rho"]:.5f}',SMALL)])
table=LongTable(cells,colWidths=[172,130,54,74,80],repeatRows=1,hAlign="LEFT")
table.setStyle(TableStyle([("VALIGN",(0,0),(-1,-1),"TOP"),("BACKGROUND",(0,0),(-1,0),colors.HexColor("#e7eff5")),("LINEBELOW",(0,0),(-1,0),.5,colors.HexColor(GRAY)),("TOPPADDING",(0,0),(-1,-1),2),("BOTTOMPADDING",(0,0),(-1,-1),2)]))
story.extend([p("Verified-only means; ratios are missing on censored cells. Native work unpriced: these are count quotients, not full-cost ratios. Every target and failure remains in per-target.json.",SMALL),table,PageBreak(),p("Original questions and remaining gates",HEAD)])
for q,text in [(1,"Prime-field base-size balance tested on fresh targets; scalable cheap decomposition/full cost open."),(2,"Actual chained metadata corrected and checked; solving-degree growth open."),(3,"Enumerated quotient witnesses lift exactly; compact algebraic correspondence open."),(4,"Parent nonlinear-base evidence retained; scalable membership/full cost open."),(5,"Parent quasi-subfield census retained; useful general family/full IC open."),(6,"Intermediate-variable accounting corrected; higher-arity solver/lifting advantage open."),(7,"Width-generic modular rank implemented; general useful-rank and descent law open."),(8,"Paid certified isogeny/Hecke transport open."),(9,"General norm and cheaper elliptic-component solver open."),(10,"Successful N83 recovery, full costs, isolated independent hardware and growing-size crossover open.")]:
    story.append(p(f"Q{q}. {text}"))
story.append(p("No new mathematical theorem. Native Rust sources and sealed records derive all values; this Python file only draws them and lays out the report. PR1853 remains subject to review and applicable CI. No research jobs continue after this completed local campaign.",SMALL))
pdf(ROOT/"REPORT.pdf",story)
for name in ["REPORT.pdf","shape-correction.svg","shape-correction.png"]: shutil.copyfile(ROOT/name,BENCH/name)
