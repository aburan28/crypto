"""Presentation only: lay out frozen native counts and native stratum edges."""
import html
import json
from pathlib import Path
import typst

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "reports/run-20261017"
OUT.mkdir(parents=True, exist_ok=True)
rows = json.loads((ROOT / "evidence/coverage.json").read_text())
(OUT / "plot-data.json").write_text(json.dumps(rows, indent=2) + "\n")
svg = ['<svg xmlns="http://www.w3.org/2000/svg" width="1120" height="810" viewBox="0 0 1120 810">',
       '<rect width="1120" height="810" fill="white"/>',
       '<g font-family="sans-serif" fill="#173047">',
       '<text x="24" y="27" font-size="19" font-weight="bold">Complete conductor-stratum coverage</text>',
       '<text x="24" y="49" font-size="12">Each box aggregates one endomorphism conductor; labels are conductor f and exact vertex population n.</text>']
for index, row in enumerate(rows):
    label = row["case"]
    strata = json.loads((ROOT / f"evidence/{label}-strata.json").read_text())
    edges = json.loads((ROOT / f"evidence/{label}-stratum-edges.json").read_text())
    left = [30, 265, 520][index]
    svg += [f'<text x="{left}" y="80" font-size="13" font-weight="bold">{html.escape(row["source_identity"]["slug"])}</text>',
            f'<text x="{left}" y="100" font-size="12">{row["vertices"]:,} vertices / {row["directed_maps"]:,} directed maps</text>']
    positions = {}
    for s in strata:
        if label == "mixed":
            x, y = 587 + 140 * s["level3"], 144 + 86 * s["level2"]
        else:
            x, y = left + 95, 144 + 43 * s["level2" if label == "two" else "level3"]
        positions[s["conductor"]] = (x, y)
    for edge in edges:
        sf, tf = edge["source_conductor"], edge["target_conductor"]
        if sf >= tf:
            continue
        x1, y1 = positions[sf]
        x2, y2 = positions[tf]
        color = "#2b6a93" if edge["degree"] == 2 else "#ad5b20"
        svg += [f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="{color}" stroke-width="2"/>']
    for s in strata:
        x, y = positions[s["conductor"]]
        svg += [f'<rect x="{x-57}" y="{y-17}" width="114" height="34" rx="5" fill="#edf3f8" stroke="#91a9bb"/>',
                f'<text x="{x}" y="{y-3}" text-anchor="middle" font-size="11">f={s["conductor"]}</text>',
                f'<text x="{x}" y="{y+11}" text-anchor="middle" font-size="11">n={s["vertices"]:,}</text>']
svg += ['<text x="24" y="755" font-size="12">Blue: degree 2. Orange: degree 3. Lines collapse reciprocal maps; horizontal loops and multiplicities are in the data exports.</text>',
        '<text x="24" y="778" font-size="11">Run 20261017. Sources: evidence/*-strata.json and *-stratum-edges.json. Exact counts; no statistical uncertainty.</text>', '</g></svg>']
(OUT / "strata.svg").write_text("\n".join(svg) + "\n")
table_md = ["| Source model | p | Maximum conductor | Vertices | Directed maps | Strata |", "| --- | ---: | ---: | ---: | ---: | ---: |"]
table_typ = []
for row in rows:
    table_md.append(f'| {row["source_identity"]["slug"]} | {row["field_p"]} | {row["maximum_conductor"]} | {row["vertices"]} | {row["directed_maps"]} | {row["strata"]} |')
    table_typ.append(f'[{row["source_identity"]["slug"]}], [{row["vertices"]:,}], [{row["directed_maps"]:,}], [{row["strata"]}],')
report = """# Three additional complete ordinary isogeny classes

All three full vertex classes and their selected-degree kernel graphs were
generated and replayed. These are additional tractable 57/59-bit examples;
the original A/C/D/E full datasets remain unfinished. Maximum conductors
are 8192, 6561 and 1728, smaller than the original large-case conductors.

""" + "\n".join(table_md) + """

![Exact conductor-stratum populations](strata.svg)

The 24,823 vertices, 69,300 directed maps and 51 strata are exact recorded
counts. Degree sets are respectively {2}, {3}, {2,3}. Every vertex is
expanded, the frontier is empty, every kernel set is regenerated, exact
rational-map identities replay, and every conductor population agrees with
the native preflight. A second replay used database copies and left the
original graph hashes unchanged. Replay shares the construction arithmetic;
an independent implementation audit is NOT_RUN. The examples are ordinary
curves with D_K = -7, -11 and -19; no new structural theorem is claimed.

`inputs/*-floor.json` gives a verified model with End(E) conductor equal to
the listed maximum. `inputs/{two,three,mixed}.json` starts at conductor 1.
Each stratum has an explicit representative and outgoing rational maps in
the downloadable archive. Model coefficients, j-invariants and conductors
for all vertices are in nodes.csv. Every directed kernel, target and target
isomorphism is in edges.jsonl and the SQLite graph. Parallel kernels and
loops are retained. All isogeny degrees and their compositions form an
infinite family; this archive enumerates the stated finite degree sets.

Reproduction: see ../../run.sh, ../../generator/, ../../PROTOCOL.md and
../../evidence/. Generation protocol commit 3b7f4b39c; frozen input and
execution commit 946ab7419. The delivered executable is the already tested
7c2e8f46a implementation with source SHA d3c77f50dc2e737de38561203ce3693f7899fe3db502b8e863c7ac69e5818bed.
Host: Linux x86_64, 4 CPU quota, 16 GiB RAM. This is a correctness/census
run; no solver runtime or performance comparison is admitted.

ECDLP min/max/mean, Pollard rho/PH/BSGS campaigns, direct elliptic-curve IC,
factor-base relation yield, and finite-field reductions remain NOT_RUN for
these three classes. No IC/rho ratios or timing data changed, so the
canonical IC progress timeline and scoreboard retain their prior values.
The curve roster, cover index and browser are refreshed separately.

The Conductor executable is unavailable (127); no task reservation was
obtained. See the prior study for the separate remaining storage blocker
for A/C/D/E. No process continues generation after this completed run.
"""
(OUT / "REPORT.md").write_text(report)
source = r'''#set page(paper: "a4", margin: 12mm)
#set text(size: 8.5pt)
#set heading(numbering: none)
= Three complete ordinary isogeny classes
Additional 57/59-bit examples; run 20261017. Original A/C/D/E remain incomplete.

#table(columns: (3fr, 1fr, 1fr, .6fr), inset: 4pt,
  [Source model], [Vertices], [Directed maps], [Strata],
''' + "\n".join(table_typ) + r''')

#image("strata.svg", width: 100%)

*Verified coverage:* 24,823 vertices, 69,300 directed maps, 51 conductor strata.
All selected-degree outgoing kernels and exact map identities replayed;
all conductor populations agree. Original database hashes survived replay
on copies. Independent implementation replay: NOT_RUN.

Maximum endomorphism conductors are 8192, 6561 and 1728. Floor-model inputs
are included. Sources start at conductor 1. SQLite and CSV/JSONL exports
contain the entire vertex class and all degree-2/3 maps in the stated sets.
Every stratum has a representative with explicit outgoing rational maps.
Arbitrary-degree maps are an infinite family, outside this finite export.

*Evidence:* evidence/coverage.json, *-build.json, *-verify.json,
*-strata.json and *-stratum-edges.json. Exact integer census; uncertainty
intervals do not apply. Protocol: PROTOCOL.md, commit 3b7f4b39c; frozen
inputs/run commands: 946ab7419. Linux x86_64, 4 CPU quota, 16 GiB RAM.
Construction/replay share native arithmetic. Conductor executable missing (127).

*Still open:* full A/C/D/E generation, independent implementation audit,
graph-wide ECDLP min/max/mean and operation counts, direct EC index calculus,
factor-base yield and finite-field reductions. No solver timings or speedup
claims are made; IC scoreboards retain their prior performance values.
Full report, exact inputs, commands, hashes and provenance accompany this PDF.
'''
(OUT / "report.typ").write_text(source)
typst.compile(str(OUT / "report.typ"), output=str(OUT / "report.pdf"))
print(OUT / "report.pdf")
