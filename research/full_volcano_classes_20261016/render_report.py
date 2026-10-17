"""Document layout only. Reads frozen native outputs; no curve arithmetic or research execution."""
import json
from pathlib import Path
import typst

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "reports" / "run-20261016T1700Z"
OUT.mkdir(parents=True, exist_ok=True)
rows = []
for case in "ABCDE":
    plan_path = ROOT / "evidence" / f"{case}-plan.json"
    plan = json.loads(plan_path.read_text())
    prefix = json.loads((ROOT / "evidence" / "prefixes" / f"{case}-verify.json").read_text())
    rows.append({"historical_case": case, "slug": plan["source_identity"]["slug"],
                 "expected_vertices": plan["expected_vertices"], "strata": len(plan["strata"]),
                 "prefix_vertices": prefix["vertices"], "prefix_maps": prefix["directed_kernel_edges"],
                 "source": str(plan_path.relative_to(ROOT))})
(OUT / "plot-data.json").write_text(json.dumps(rows, indent=2) + "\n")

# Positions and scales below lay out recorded values only.
svg = ['<svg xmlns="http://www.w3.org/2000/svg" width="940" height="365" viewBox="0 0 940 365">',
       '<rect width="940" height="365" fill="white"/>',
       '<g font-family="sans-serif" fill="#152b40">',
       '<text x="18" y="27" font-size="18" font-weight="bold">Expected full-class populations — native exact preflights</text>',
       '<text x="18" y="48" font-size="12">Vertices, not generated counts. Cases preserve historical labels; exact source model slugs identify each row.</text>']
maximum = max(row["expected_vertices"] for row in rows)
for i, row in enumerate(rows):
    y = 84 + i * 51
    width = 360 * row["expected_vertices"] / maximum
    svg += [f'<text x="18" y="{y}" font-size="12">{row["slug"]}</text>',
            f'<rect x="325" y="{y-14}" width="{width:.3f}" height="21" rx="3" fill="#2b6a93"/>',
            f'<text x="695" y="{y}" font-size="12">{row["expected_vertices"]:,} vertices</text>',
            f'<text x="325" y="{y+22}" font-size="11">{row["strata"]} expected strata; current prefix: {row["prefix_vertices"]:,} vertices / {row["prefix_maps"]:,} maps</text>']
svg += ['<text x="18" y="351" font-size="11">Run 20261016T1700Z. Source: evidence/{A..E}-plan.json and prefixes/*-verify.json. No uncertainty interval: exact integer counts.</text>', '</g></svg>']
(OUT / "coverage.svg").write_text("\n".join(svg) + "\n")

source = r'''#set page(paper: "a4", margin: 15mm)
#set text(size: 9pt)
#set heading(numbering: none)
= Complete isogeny classes: implementation and coverage
Native enumeration controls and fixed-input preflights, 16 October 2026.

The CLI is implemented and the bounded runs below passed. The full five-case
dataset request remains open. A complete vertex class and a complete graph
for selected prime degrees are separate claims. Arbitrary isogeny degrees
form an infinite family.

#image("coverage.svg", width: 100%)

== Observed correctness evidence
#table(columns: (2.5fr, 1fr, 1fr, 2fr), inset: 4pt,
  [Source model (ICV1)], [Vertices], [Directed maps], [Coverage],
  [icv1-fp7-t2-4910f2d0], [8], [16], [Complete degree-2 class; 4 strata],
  [icv1-fp8-tm1-1ed3feb7], [9], [18], [Complete degree-3 class; 3 strata],
  [icv1-fp7-t4-0c09e92f], [10], [36], [Complete degree-2/3 class; 4 strata])

All three vertex sets agree with an independent exhaustive census of every
nonsingular model in the corresponding small field. Six native release tests
also cover a nontrivial surface class group, checkpoint/resume equivalence,
corrupt-map rejection, missing-degree components, BSGS and subgroup/budget checks.
Two binary controls also check a killed writer's recovery/resume against an
uninterrupted graph and preservation of existing output files.
Raw test output: evidence/native-tests.txt. No method timing or speed comparison
is admitted by these correctness tests.

Each exact large input was expanded 1,024 times and replayed. Those outputs are
partial; unexpanded frontier vertices remain. Replay regenerated every outgoing
kernel on expanded vertices and checked every exact rational-map identity and
target scaling. It shares native arithmetic with construction. An independent
implementation audit of the new generic backend has not been performed.

== What the command does
`crypto isogeny class build --curve curve.json --out graph.sqlite --steps all`
uses certified root order, exact class populations, complete rational-kernel
enumeration and atomic SQLite vertex expansions. The equivalent standalone
commands are `isogeny class ...` and `isogeny-class ...`.
Automatic completion triggers replay and conductor-population checks. Resume
keeps committed expansions. Vertex inspection returns explicit rational maps.

Backend limits: ordinary prime fields 3 < p < 2^61, -10,000,000 <= D_K < -4,
prime degrees <=31 different from p. Other cases fail explicitly. The source
endomorphism conductor is not inferred from the Frobenius conductor alone.

== Full datasets, storage and DLP obligations
The historical degree-3 case has all 129,140,163 vertices, 18 strata and
258,280,326 directed kernel maps in its separately distributed lossless archive.
Its compressed bytes and hashes are preserved in evidence/historical-B/.
Historical timings were produced outside this repository's measurement protocol.

The five full classes contain 7,833,658,967 vertices and 515 strata. A 24-byte
layout for a,b,j alone occupies 188,007,815,208 bytes, excluding edges, indexes,
journals and benchmark evidence. This is a fixed-layout baseline, not a measured
SQLite size estimate or an information-theoretic lower bound. Compressed or
generative encodings may be smaller. The workspace has 32 GiB total.
Full remaining class generation needs an appropriate durable execution volume
or another validated storage backend; it is not declared impossible.

The solve command verifies one explicit cyclic-subgroup BSGS logarithm. It does
not complete the requested graph-wide rho/BSGS/PH and direct/finite-field IC
campaign. All-vertex end-to-end timing, total operations, relation yields and
min/max/mean remain pending. No sampled result is labeled exhaustive.

Evidence and reproduction: PROTOCOL.md, run-controls.sh, run-prefixes.sh,
tools/volcano-cli/README.md, native plans, build/replay JSON and the named storage
problem. The existing class-number and ordinary-volcano theorems are reused;
no new mathematical theorem or Lean proof is claimed.
'''
(OUT / "report.typ").write_text(source)
typst.compile(str(OUT / "report.typ"), output=str(OUT / "report.pdf"))
