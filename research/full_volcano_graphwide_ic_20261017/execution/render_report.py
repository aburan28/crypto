"""Presentation only: draw statistics already recorded by the Rust summarizer."""
import html
import json
import sys
from pathlib import Path

import typst

data_dir, out = map(Path, sys.argv[1:3])
out.mkdir(exist_ok=False, parents=True)
doc = json.loads((data_dir / "summary.json").read_text())
(out / "plot-data.json").write_text(json.dumps(doc, indent=2) + "\n")
groups = {(r["graph"], r["conductor_0_means_all"], r["arm"], r["status"]): r
          for r in doc["groups"]}
classes = [("two", "icv1-fp57-t536870914-6cabb0c9", "r = 61,751"),
           ("three", "icv1-fp57-tm536870941-e6d0c766", "r = 19,768,889,568,671"),
           ("mixed", "icv1-fp59-t1073741954-4dd0b1ef", "r = 4,253")]
arms = ["projected-16-mitm", "projected-16-subtract", "projected-64-mitm", "projected-64-subtract"]
colors = ["#256f9a", "#7faec7", "#ba661c", "#dfab7e"]
svg = ['<svg xmlns="http://www.w3.org/2000/svg" width="1080" height="410" viewBox="0 0 1080 410">',
       '<rect width="1080" height="410" fill="white"/>',
       '<g font-family="sans-serif" fill="#18364c">',
       '<text x="20" y="28" font-size="21" font-weight="bold">Direct EC relation yield: audited root-vertex pilot</text>',
       '<text x="20" y="51" font-size="13">Two public targets per curve and arm. Mean with observed min–max; no confidence interval.</text>']
table = []
for ci, (graph, slug, subgroup) in enumerate(classes):
    left = 65 + ci * 350
    svg += [f'<text x="{left-20}" y="80" font-size="12">{html.escape(slug)}</text>',
            f'<text x="{left-20}" y="99" font-size="12">{subgroup}</text>']
    for tick in [0, .25, .5, .75, 1]:
        y = 300 - 180 * tick
        svg += [f'<line x1="{left}" y1="{y}" x2="{left+255}" y2="{y}" stroke="#dce5ea"/>',
                f'<text x="{left-35}" y="{y+4}" font-size="11">{tick:.2f}</text>']
    for ai, arm in enumerate(arms):
        row = groups.get((graph, 0, arm, "verified")) or groups.get((graph, 0, arm, "exhausted"))
        metric = row["metrics"]["relation_yield"]
        x = left + 15 + 60 * ai
        mean, low, high = (metric[k] for k in ["mean", "min", "max"])
        svg += [f'<rect x="{x}" y="{300-180*mean}" width="34" height="{180*mean}" fill="{colors[ai]}"/>',
                f'<line x1="{x+17}" x2="{x+17}" y1="{300-180*high}" y2="{300-180*low}" stroke="#132c3a" stroke-width="2"/>',
                f'<text x="{x+17}" y="{max(116,291-180*high)}" text-anchor="middle" font-size="11">{mean:.3f}</text>',
                f'<text transform="translate({x+19},315) rotate(35)" font-size="11">{arm.removeprefix("projected-")}</text>']
        gae = row["metrics"]["cold_gae_lower_bound"]
        trial = row["metrics"]["relation_trials"]
        table.append(f'[{slug}], [{arm.removeprefix("projected-")}], [{row["status"]}], [{row["records"]}], [{gae["min"]:,.0f}], [{gae["mean"]:,.1f}], [{gae["max"]:,.0f}], [{trial["sum"]:,.0f}], [{row["metrics"]["relations_found"]["sum"]:,.0f}],')
svg += ['<text x="20" y="394" font-size="11">Measured relation counters; denominator is attempted decompositions. Zero is observed, not a claim of impossibility.</text>', '</g></svg>']
(out / "relation-yield.svg").write_text("\n".join(svg) + "\n")
source = r'''#set page(paper: "a4", margin: 13mm)
#set text(size: 8.3pt)
#set heading(numbering: none)
= Volcano index calculus: measured control and pilot
Run 20261017-control-pilot. Counts and correctness diagnostic; no runtime or speedup claim.

#image("relation-yield.svg", width: 100%)

*Actual scope:* 3 root vertices (conductor 1), 2 public targets each, 8 arms,
1 round, 0 warmups: 48 pilot attempts. Additional small-subgroup control:
12 attempts. This is not graph-wide completion and does not replace the
original A--E requirement. A separate resumable survey targets all 24,823
vertices and 51 strata; consult its current ledger rather than this frozen report.

*Observed:* pilot: 30 verified recoveries, 14 exhausted budgets, 4 watchdog
timeouts. Control: 8 verified, 4 exhausted. Both projected-base methods solved
both control targets. The raw-abscissa control exhausted its budget in every
pilot target. On the actual 45-bit subgroup the projected arms found zero
relations in 4,096 trials per target; no smaller subgroup was substituted.

*Validation:* all four sealed sessions pass the ecbench audit. All 38 successful
deterministic records replay identically. Supplemental native replay also
reproduces all 18 exhausted records, including every phase/counter, field
operation and base identity. Four timeouts have no final counted report;
their unknown relation counts are not treated as zero. Every recovered scalar
is checked by multiplication in the runner and auditor. The replays use the
same implementation and host; independent implementation replay is still open.

*Method:* cofactor-project each raw abscissa point, reject infinity, check
prime-subgroup membership and deduplicate negation pairs. Factor-base sizes
16 and 64; scan cap 20,000; subtraction and two-summand MITM; checked full-rank
incremental elimination. Public-target seed 20261017; order seed 1, alternate.
Cold base/table setup and verification are charged on every one-target run.
References are rho.negation and textbook BSGS, with matched rho A/A. A 3-second
child watchdog censors unfinished runs. Both methods use actual largest r.

#pagebreak()
= Counted cold costs and yields
GAE values below are *lower bounds with unpriced native work*. Each row has
two measured target attempts; min/mean/max describe those observations only.
Exhausted rows are kept separate and do not claim recovered DLP cost.

#set text(size: 6.6pt)
#table(columns: (2.6fr, 1.3fr, .9fr, .3fr, .8fr, .8fr, .8fr, .8fr, .7fr), inset: 2.5pt,
  [Source model], [Projected base/oracle], [Status], [n], [GAE min], [GAE mean], [GAE max], [Trials sum], [Relations sum],
''' + "\n".join(table) + r''')
#set text(size: 8.3pt)

*Complete accounting:* native summary.json and statistics.csv include min,
max, arithmetic mean and explicit n by graph, conductor, arm, status and
phase; observations.jsonl retains record/workload ids and all phase-native
counters. Raw sessions retain every field operation, phase, setup, matrix
rank/dependency, target window and failure. pricing-band.json applies common
0/0.1/1 GAE native-work weights to IC and reference arms alike. No calibrated
versus uncalibrated speed ranking is inferred here.

*Host and limits:* Intel Xeon Platinum 8573C, Linux x86-64 container, four-CPU
quota and 16 GiB RAM. ecbench reserves CPU 4; every measured record earns L0:
unoptimized build, movable-thread/PSI and host-isolation limitations. Raw
in-process and end-to-end elapsed observations are retained. None is an
admitted optimized runtime distribution. Two targets/one round also do not
qualify wall-clock comparisons. PMU instructions/cycles were unavailable.

*Provenance:* projected-base source/tests and first runner commit 4eecf6b53;
exhaustion replay and survey ordering commit ee54ad603. Exact ecbench SHA-256:
#text(size: 6.7pt)[1869ca2e7be271794df8128c08a6dca0201038e2c81a3e0ce6434f7ec0d8379d].
Per-campaign manifests bind source, binary, queue and input hashes. Exact
control/pilot protocols, frozen specs, plans, records and replay receipts
accompany this report. Compiler/test exit logs are retained. New subgroup
base tests pass 2/2, including every nonzero control logarithm without a
planted-answer input. Conductor reservation was unavailable (exit 127).

*Open requirements:* exact per-vertex representation registry, optimized and
qualified full-matrix execution, all targets/rounds/base sizes, transported
edge experiments, composite-order PH, independent implementation replay and
finite-field IC. Actual embedding degrees are 30,875; 3,953,777,913,734; and
2,126. The existing supersingular degree-2 teaching MOV backend cannot
execute finite-field IC in those destinations. This is an implementation
gap, not proof that mathematical reductions do not exist. Original A/C/D/E
full generation remains storage-blocked; no original case is replaced here.

No new mathematical theorem or fastest-method claim is made. The established
IC/rho frontier remains unchanged. These results establish the measured
control/pilot behavior, not unrun strata or vertices.
'''
(out / "report.typ").write_text(source)
typst.compile(str(out / "report.typ"), output=str(out / "report.pdf"))
(out / "figure.typ").write_text('#set page(width: 1080pt, height: 410pt, margin: 0pt)\n#image("relation-yield.svg", width: 100%)\n')
typst.compile(str(out / "figure.typ"), output=str(out / "figure.pdf"))
print(out / "report.pdf")
