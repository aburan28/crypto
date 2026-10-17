"""Presentation only: native checkpoint/summary values, no new measurements."""
import html
import json
import sys
from pathlib import Path
import typst

data, out = map(Path, sys.argv[1:3])
doc = json.loads((data / "summary.json").read_text())
out.mkdir(parents=True, exist_ok=False)
(out / "plot-data.json").write_text(json.dumps(doc, indent=2) + "\n")
svg = ['<svg xmlns="http://www.w3.org/2000/svg" width="1050" height="360" viewBox="0 0 1050 360">',
       '<rect width="1050" height="360" fill="white"/><g font-family="sans-serif" fill="#18364c">',
       '<text x="20" y="28" font-size="21" font-weight="bold">First 12 strata: projected size-64 MITM relation yield</text>',
       '<text x="20" y="51" font-size="13">Two public targets per vertex; mean with observed min–max. All strata are separate graph/stratum pairs.</text>']
table = []
for ci, (graph, slug) in enumerate([
    ("two", "icv1-fp57-t536870914-6cabb0c9"),
    ("three", "icv1-fp57-tm536870941-e6d0c766"),
    ("mixed", "icv1-fp59-t1073741954-4dd0b1ef")]):
    left = 60 + ci * 345
    svg.append(f'<text x="{left-25}" y="85" font-size="12">{html.escape(slug)}</text>')
    for t in [0, .25, .5, .75, 1]:
        y = 280 - 170 * t
        svg += [f'<line x1="{left}" x2="{left+250}" y1="{y}" y2="{y}" stroke="#dce5ea"/>',
                f'<text x="{left-35}" y="{y+4}" font-size="11">{t:.2f}</text>']
    rows = [r for r in doc["groups"] if r["graph"] == graph and r["conductor_0_means_all"] > 0 and r["arm"] == "projected-64-mitm"]
    for index, row in enumerate(rows):
        m = row["metrics"]["relation_yield"]
        x = left + 15 + index * 60
        svg += [f'<rect x="{x}" y="{280-170*m["mean"]}" width="32" height="{170*m["mean"]}" fill="#256f9a"/>',
                f'<line x1="{x+16}" x2="{x+16}" y1="{280-170*m["max"]}" y2="{280-170*m["min"]}" stroke="#132c3a" stroke-width="2"/>',
                f'<text x="{x+16}" y="{max(106,270-170*m["max"])}" text-anchor="middle" font-size="11">{m["mean"]:.3f}</text>',
                f'<text x="{x+16}" y="300" text-anchor="middle" font-size="11">f={row["conductor_0_means_all"]}</text>']
        cost = row["metrics"]["cold_gae_lower_bound"]
        table.append(f'[{graph}], [{row["conductor_0_means_all"]}], [{row["status"]}], [{row["records"]}], [{cost["min"]:,.0f}], [{cost["mean"]:,.1f}], [{cost["max"]:,.0f}], [{m["mean"]:.6f}],')
svg += ['<text x="20" y="336" font-size="12">f = conductor of End(E). One vertex per displayed stratum; variation is not a graph-wide distribution.</text>', '</g></svg>']
(out / "strata-yield.svg").write_text("\n".join(svg) + "\n")
source = r'''#set page(paper: "a4", margin: 13mm)
#set text(size: 8.5pt)
= All-vertex survey: first 12 audited strata
Run 20261017-survey-prefix-12. Frozen checkpoint, not completion.

#image("strata-yield.svg", width: 100%)

*Native recorded coverage:* 12 audited vertices in 12 distinct graph/stratum
pairs. Four vertices in each graph. The immutable queue contains all 24,823
vertices and 51 strata; the first phase visits one vertex per stratum before
all remaining vertices. Pilot/root reruns use the same targets and seeds;
they are repeated observations, not additional independent targets.

*Observed outcomes:* 192 measured attempts: 119 verified, 58 exhausted,
15 watchdog timeouts. All 12 standard audits pass. Each successful deterministic
record replays, and each exhausted deterministic record passes the additional
native replay receipt. No wrong answers were recorded. Timeouts have unknown
final relation and operation counters. Their elapsed observations are retained.

*Recorded direct-IC counters:* 301,715 attempted decompositions, 6,368 relations
found across the 120 direct-IC attempts with final reports. Numerator and
denominator include exhausted reports and the raw-base control. They are not
120 successful DLPs or independent identically distributed trials. Phase,
rank, dependency, failure, memory and unpriced-work records remain available.

*Subset below:* size-64 projected MITM only. Cold GAE are lower bounds with
unpriced native work; no fully priced cost or runtime ranking follows.
Each row is two measured one-target attempts on one exact vertex.

#set text(size: 7.5pt)
#table(columns: (.8fr,.5fr,1fr,.3fr,1fr,1fr,1fr,1fr), inset: 3pt,
  [Graph key], [f], [Status], [n], [GAE min], [GAE mean], [GAE max], [Mean yield],
''' + "\n".join(table) + r''')
#set text(size: 8.5pt)

*Conditions:* actual largest-prime subgroups; projected sizes 16/64 with
MITM/subtraction, unchanged raw-base MITM, rho.negation, textbook BSGS and rho
A/A. Two public targets, one round, no warmup, 4,096 trials, 3-second watchdog.
Prime-subgroup membership and all recovered scalars are checked. Every run
charges its cold solver setup; no reusable logs are assumed free. Ordinary
finite-field IC remains unsupported at the actual recorded embedding degrees.

*Host:* Intel Xeon Platinum 8573C, Linux x86-64, four-CPU quota, 16 GiB RAM.
Unoptimized ecbench SHA-256 1869ca2e7be271794df8128c08a6dca0201038e2c81a3e0ce6434f7ec0d8379d;
native driver source commit ee54ad603. All records earn L0; raw solve/process
times are descriptive, not admitted performance. Short driver rejection checks
ran during this counts-only block; no build or test suite ran during it.
The observed invocation took 272 elapsed seconds, including runner overhead
and audits. This is a local scheduling observation, not an all-vertex forecast.

*Evidence:* execution/campaigns/survey-prefix-12, immutable full queue in the
delivery archive, campaign.json hashes, SQLite checkpoint, standard audits
and exhaustion-replay receipts. The native summary v2 includes queued/audited
coverage by conductor and explicit sample counts for every min/max/mean.
Common 0/0.1/1 GAE pricing-band output is separate. Same-implementation replay
is not independent implementation validation. A second writer and changed
driver were rejected; zero-work resume preserved the 12 audited sessions.

The original full matrix (million trials, eight targets, five rounds, warmup,
four base sizes), original A--E full graphs, exact representation registration,
transported targets, full composite PH, qualified runtimes, finite-field IC
and independent implementation replay remain open. A real queue is not a
completed solve. Canonical best-method/frontier claims are unchanged.
'''
(out / "report.typ").write_text(source)
typst.compile(str(out / "report.typ"), output=str(out / "report.pdf"))
(out / "figure.typ").write_text('#set page(width: 1050pt, height: 360pt, margin: 0pt)\n#image("strata-yield.svg", width: 100%)\n')
typst.compile(str(out / "figure.typ"), output=str(out / "figure.pdf"))
print(out / "report.pdf")
