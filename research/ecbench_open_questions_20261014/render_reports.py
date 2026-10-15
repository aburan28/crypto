"""Rendering only: draw and lay out values already recorded by native Rust.

No experiments, cost derivations, statistics, counts, ratios or verdicts here.
"""
import os
os.environ.setdefault("MPLCONFIGDIR", "/tmp/ecbench-ten-question-mpl")
from pathlib import Path
import json
import re
import shutil
from xml.sax.saxutils import escape
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib import colors
from reportlab.platypus import (
    SimpleDocTemplate, Paragraph, Spacer, Image, PageBreak, Table, TableStyle,
)
from reportlab.lib.pagesizes import A4

ROOT = Path(__file__).resolve().parent
DATA = json.loads((ROOT / "render-data.json").read_text())
NAVY, TEAL, GRAY = "#17354a", "#187f83", "#586775"
FONTS = Path("/usr/share/fonts/truetype/dejavu")
pdfmetrics.registerFont(TTFont("DV", str(FONTS / "DejaVuSans.ttf")))
pdfmetrics.registerFont(TTFont("DVB", str(FONTS / "DejaVuSans-Bold.ttf")))
pdfmetrics.registerFontFamily("DV", normal="DV", bold="DVB", italic="DV")
plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 9, "svg.fonttype": "none",
    "axes.spines.top": False, "axes.spines.right": False,
    "text.color": NAVY, "axes.labelcolor": NAVY,
})
BODY = ParagraphStyle("body", fontName="DV", fontSize=8.7, leading=12,
                      textColor=colors.HexColor(NAVY), spaceAfter=7,
                      splitLongWords=True)
TITLE = ParagraphStyle("title", parent=BODY, fontName="DVB", fontSize=17,
                       leading=21, spaceAfter=12)
HEAD = ParagraphStyle("heading", parent=BODY, fontName="DVB", fontSize=11.5,
                      leading=15, spaceBefore=9, spaceAfter=7)
SMALL = ParagraphStyle("small", parent=BODY, fontSize=7.4, leading=10)
BENCH = ROOT.parents[1] / "benchmarks/20261014-ecbench-ten-questions"
BENCH.mkdir(parents=True, exist_ok=True)

def para(text, style=BODY):
    return Paragraph(escape(str(text)), style)

def figsave(fig, directory, name):
    directory.mkdir(parents=True, exist_ok=True)
    for suffix in ["svg", "png"]:
        fig.savefig(directory / f"{name}.{suffix}", bbox_inches="tight",
                    **({"dpi": 150} if suffix == "png" else {}))
        if suffix == "svg":
            path = directory / f"{name}.{suffix}"
            path.write_text("\n".join(line.rstrip() for line in path.read_text().splitlines())+"\n")
    plt.close(fig)

def footer(canvas, doc):
    canvas.setFont("DV", 7)
    canvas.setFillColor(colors.HexColor(GRAY))
    canvas.drawString(36, 20, "Native ecbench • 2026-10-14 • counted diagnostics • no full-cost claim")
    canvas.drawRightString(A4[0]-36, 20, str(doc.page))

def pdf(path, story):
    SimpleDocTemplate(str(path), pagesize=A4, rightMargin=36, leftMargin=36,
                      topMargin=34, bottomMargin=36).build(
                          story, onFirstPage=footer, onLaterPages=footer)

def plot_session(session):
    out = ROOT / session["path"]
    rows = session["rows"]
    fig, ax = plt.subplots(figsize=(8.4, max(2.5, len(rows)*0.29 + 1.2)))
    numeric = [r for r in rows if r["mean_s"] is not None]
    if not numeric:
        ax.axis("off")
        boxes = [
            (0.03, "Registered N83 E0\npublic target", "#e7eff5"),
            (0.37, "IC width boundary\nimplementation missing", "#f5e8df"),
            (0.71, "No scalar\nno cost result", "#edf0f3"),
        ]
        for x, label, color in boxes:
            ax.text(x, 0.52, label, transform=ax.transAxes, va="center",
                    ha="left", fontsize=11,
                    bbox={"boxstyle":"round,pad=0.7", "facecolor":color,
                          "edgecolor":"#9baebb"})
        ax.text(0.03, 0.08, "One recorded applicability error; no comparable S values.",
                transform=ax.transAxes, fontsize=10)
    else:
        for i, row in enumerate(rows):
            if row["mean_s"] is not None:
                color = NAVY if row["arm"].startswith("rho") else TEAL
                ax.scatter(row["mean_s"], i, color=color, s=24, zorder=3)
                ax.hlines(i, row["min_s"], row["max_s"], color=color, linewidth=1.1)
            else:
                ax.text(0.985, i, "no verified scalar",
                        transform=ax.get_yaxis_transform(), ha="right", va="center",
                        fontsize=7.5, color=GRAY)
        ax.set_xscale("log")
        ax.set_yticks(range(len(rows)), [
            f'{r["curve"]}\n{r["arm"]}   {r["verified"]}/{r["records"]} verified'
            for r in rows
        ], fontsize=7.5)
        ax.set_ylim(len(rows)-0.5, -0.5)
        ax.set_xlabel("Cold counted S = GAE / √r; native work remains unpriced")
        ax.grid(axis="x", alpha=0.2)
    ax.set_title(f'{session["session"]} • {session["ledger"]["session_id"]}',
                 loc="left", fontweight="bold", fontsize=11, pad=12)
    figsave(fig, out, "cost")
    story = [para(session["session"] + ": native ecbench run", TITLE),
             para("Session " + session["ledger"]["session_id"]),
             para("Source " + session["metadata"]["git_commit"], SMALL),
             para("Binary " + session["metadata"]["binary_sha256"], SMALL),
             para("Host " + session["metadata"]["env_class_id"] +
                  "; earned L0. Public points, no planted scalar, paired rho/A/A, "
                  "seed 20261014, cold setup repeated. Wall clocks are inadmissible."),
             para(f'{session["ledger"]["records"]} attempts; '
                  f'{session["ledger"]["verified_records"]} verified recoveries; '
                  f'{session["ledger"]["replays"]} exact standard replays.'),
             para("Exhausted outcomes and the N83 error were replayed separately. "
                  "Successful-only means/ranges do not include failed work; failed "
                  "attempts retain their counters and resource vectors. Ratios are "
                  "missing on censored curve/arm cells."),
             para("Graph: points are recorded mean S; lines are observed min/max "
                  "over verified recoveries, not confidence intervals. Sample counts "
                  "appear beside every arm. A quotient of partial counts does not "
                  "bound full cost. This is no speedup or scaling claim.", SMALL)]
    graph_height = min(490, max(180, len(rows)*17 + 70))
    story += [Image(str(out/"cost.png"), width=510, height=graph_height),
              para("Raw records: " + STUDY_PATH + "/sessions/" + session["session"] +
                   "/records.jsonl", SMALL),
              para("Receipts: " + session["session"] +
                   ".audit.json; failed-replays-v2.json; resource vector " +
                   session["session"] + ".resources.json.", SMALL)]
    pdf(out/"REPORT.pdf", story)
    target = BENCH / session["session"]
    target.mkdir(exist_ok=True)
    for name in ["REPORT.pdf", "cost.svg", "cost.png"]:
        shutil.copyfile(out/name, target/name)

STUDY_PATH = "research/ecbench_open_questions_20261014"
for session in DATA["sessions"]:
    plot_session(session)

# The quantities below are quoted from native summaries, never derived here.
fig, ax = plt.subplots(figsize=(9, 4.9))
names = [s["session"] for s in DATA["sessions"]]
ax.barh(names, [s["ledger"]["records"] for s in DATA["sessions"]],
        color="#c6d9e5", label="Recorded attempts")
ax.barh(names, [s["ledger"]["verified_records"] for s in DATA["sessions"]],
        color=TEAL, label="Verified recoveries")
ax.invert_yaxis()
ax.set_xlabel("Executions • failures are retained")
ax.set_title("Nine ecbench sessions covering all ten research tracks",
             fontweight="bold", loc="left")
ax.legend(frameon=False)
figsave(fig, ROOT, "execution-ledger")

control = DATA["results"]["controls"]
control_dir = ROOT / "reports/run-20261014T1954Z-extended-controls"
fig, axes = plt.subplots(1, 2, figsize=(9, 3.8))
qrows = control["quasi_subfield"]["data"]
labels = [f'n={r["n"]}, n₀={r["n0"]}' for r in qrows]
axes[0].bar(labels, [r["coefficient_vectors_examined"] for r in qrows], color=NAVY)
axes[0].set_yscale("log")
axes[0].set_ylabel("Coefficient vectors examined (exact finite domain)")
axes[1].bar(labels, [r["split_count"] for r in qrows], color=TEAL)
axes[1].set_ylabel("Splitting polynomials found (recorded count)")
for ax in axes:
    ax.tick_params(axis="x", rotation=25)
fig.suptitle("Quasi-subfield census • native controls-v3.json",
             fontweight="bold", color=NAVY)
figsave(fig, control_dir, "quasi-census")
pdf(control_dir/"REPORT.pdf", [
    para("Extended structural controls", TITLE),
    para("Frozen protocol v1 plus structural/replay amendments v2/v3. "
         "Receipt: controls-v3.json. These are construction, algebraic and "
         "functional checks; no timing or full-cost comparison was made."),
    para("Quasi-subfield census", HEAD),
    Image(str(control_dir/"quasi-census.png"), width=510, height=230),
    para("Every retained witness was checked by independent root enumeration. "
         "Zero results apply only to the specified coefficient domain."),
    para("Nonlinear bases, degree profiles and QUBO", HEAD),
    para("Six gamma-weight bases passed membership/order/orbit coefficient checks. "
         "Macaulay rank profiles are recorded through degree 5 at small dimensions; "
         "they are not a solving-degree theorem. All 8,192 original QUBO assignments "
         "passed faithful completion checks. No annealer or quantum hardware ran."),
    para("Actual covers and functional recovery", HEAD),
    para("Degree 7 and degree 83 E0/E1 divisor pullbacks passed checked arithmetic "
         "controls. The degree 7 public point recovered scalar " +
         control["covers"]["data"][0]["public_target_destination_recovery"]["recovered"] +
         " through generic Jacobian BSGS and verified on both groups. No measured "
         "speed claim follows. General norm and a cheaper destination IC remain open."),
    para("The E1 cofactor-two constructor failure and first replay-checker identity "
         "failure are preserved, with explicit corrections. F101 genus-2 certificate "
         "replay passed; its two-infinity Jacobian backend is unsupported.", SMALL),
])
for file in ["REPORT.pdf", "quasi-census.svg", "quasi-census.png"]:
    (BENCH/"extended-controls").mkdir(exist_ok=True)
    shutil.copyfile(control_dir/file, BENCH/"extended-controls"/file)

results = DATA["results"]
story = [
    para("All ten index-calculus questions: ecbench evidence", TITLE),
    para("Prime fields • ordinary binary fields • binary Koblitz curves"),
    para(f'{results["attempts"]} attempts; {results["verified_recoveries"]} verified '
         f'recoveries; {results["statuses"]["exhausted"]} exhausted searches; '
         f'{results["statuses"]["error"]} N83 applicability error.'),
    para(f'Nine standard audits passed. {results["standard_verified_replays"]} '
         f'verified runs and {results["supplemental_failed_replays"]} failure outcomes '
         'reproduced. All ten original mathematical questions remain open.'),
    para("The native harness is aburan28/crypto. Source c716e462934128c75e99dd458748e08623a2708c; "
         "binary cc5bba39d33465ddc1b1d31af8ac25c6680d41193fc76012d9a1019b6046f02d; "
         "host ECBENV2h687984d4e73d, virtualized AMD EPYC 9V45, L0.", SMALL),
    Image(str(ROOT/"execution-ledger.png"), width=510, height=280),
    para("Counts and failures only. No full-cost speedup, asymptotic advance, "
         "independent host, admitted wall comparison, GPU or quantum run is claimed. "
         "A quotient of partial counts does not bound full cost. The rho controls "
         "use the same public points and algorithm seeds; setup is charged cold."),
    para("Scope and sources", HEAD),
    para("The ranking and original resolutions come from the source-linked "
         "literature/questions.json and literature/REPORT.md, reviewed 2026-10-14. "
         "Recent sources include Caminata–Gorla (2023), Bernstein–Lange (2024), "
         "the September 2026 Hamming-ideal and annealing papers, plus the "
         "2020/2021 subfield and quasi-subfield work. Gamma-weight point bases "
         "are distinct from syndrome-decoding Hamming ideals."),
    PageBreak(),
]
for requirement in results["requirements"]:
    story += [para(str(requirement["question"]) + ". " + requirement["title"], HEAD),
              para("Status: partial. " + requirement["remaining"])]
    citations = []
    for ident in requirement["sources"]:
        source = next(s for s in DATA["literature"]["sources"] if s["id"] == ident)
        citations.append(f'<link href="{escape(source["url"])}">{escape(ident)}</link>')
    story.append(Paragraph("Repo/paper basis: " + "; ".join(citations), SMALL))
    if requirement["question"] == 5:
        story.append(PageBreak())
story += [PageBreak(), para("Structural observations and limits", TITLE),
          Image(str(control_dir/"quasi-census.png"), width=510, height=230),
          para("Independent root checks verify finite witnesses, not a general "
               "family or natural-query complexity law. Actual subgroup support "
               "is recorded per witness."),
          para("The actual chained system has intermediate field variables beyond "
               "the factor-base coordinates. First fall degree, Boolean roots "
               "and full point-decomposition success are separate observations."),
          para("Degree-83 actual genus-3 image classes passed annihilation and sampled "
               "homomorphism/push-pull checks. The prime-field genus-2 certificate "
               "is valid, while its arithmetic and transport remain unsupported. "
               "Generic toy Jacobian recovery supplies a functional check only."),
          para("N83 IC was attempted and returned the missing-wide-path error. "
               "Degree 31 has ord₃₁(2)=5 and does not replace N83/N131 module "
               "structure. E1's smaller degree-83 subgroup does not replace E0's "
               "required subgroup."),
          para("The exact ten-summand toy experiment ran without downscaling. "
               "The raw sessions preserve every failed attempt and native "
               "resource vector. Primary online windows are recorded, but their "
               "L0 wall clocks are inadmissible."),
          para("Reproduce and review", HEAD),
          para("Protocols, sealed records, standard audits, failure replay inputs/"
               "outputs, locked dependencies, source hashes, native tools and "
               "per-run reports are beside this report. No new mathematical theorem "
               "or independent source-order certificate is claimed."),
          Paragraph('<link href="https://github.com/aburan28/crypto/tree/research/'
                    'ecbench-ten-questions-20261014/research/'
                    'ecbench_open_questions_20261014">Repository evidence directory</link>', SMALL)]
pdf(ROOT/"REPORT.pdf", story)
for file in ["REPORT.pdf", "execution-ledger.svg", "execution-ledger.png"]:
    shutil.copyfile(ROOT/file, BENCH/file)
print("Rendered nine session PDFs, structural-control PDF, result PDF and 150dpi figures.")
