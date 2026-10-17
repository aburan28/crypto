"""Presentation only: plot the certified counts already emitted by Rust."""
import html
import json
from pathlib import Path
import typst

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "reports/run-20261017-inputs"
OUT.mkdir(parents=True, exist_ok=True)
doc = json.loads((ROOT / "evidence/inventory.json").read_text())
rows = doc["cases"]
svg = ['<svg xmlns="http://www.w3.org/2000/svg" width="980" height="300" viewBox="0 0 980 300">',
       '<rect width="980" height="300" fill="white"/>',
       '<g font-family="sans-serif" fill="#18364c">',
       '<text x="20" y="30" font-size="20" font-weight="bold">Certified workload inputs across every stored vertex</text>',
       '<text x="20" y="55" font-size="13">Inputs only. Solver execution: NOT_RUN. No timing or relation-yield measurements.</text>']
for index, row in enumerate(rows):
    y = 100 + 63 * index
    width = row["vertices_certified"] / 25
    svg += [f'<text x="20" y="{y}" font-size="12">{html.escape(row["source_model"])}</text>',
            f'<rect x="390" y="{y-17}" width="{width}" height="25" fill="#3378a0"/>',
            f'<text x="815" y="{y}" font-size="12">{row["vertices_certified"]:,} vertices</text>',
            f'<text x="390" y="{y+27}" font-size="11">{row["strata_certified"]} strata; r={row["largest_prime_subgroup"]}; k={row["embedding_degree"]}</text>']
svg += ['<text x="20" y="287" font-size="11">Source: native evidence/inventory.json, run 20261017-inputs. Exact integer populations; no statistical uncertainty.</text>', '</g></svg>']
(OUT / "coverage.svg").write_text("\n".join(svg) + "\n")
(OUT / "plot-data.json").write_text(json.dumps(doc, indent=2) + "\n")
table = "\n".join(f'[{r["source_model"]}], [{r["vertices_certified"]}], [{r["largest_prime_subgroup"]}], [{r["embedding_degree"]}],' for r in rows)
source = r'''#set page(paper: "a4", margin: 15mm)
#set text(size: 9pt)
#set heading(numbering: none)
= Full graph-wide IC: certified workload preparation
Run 20261017-inputs. This report contains *zero solver executions*.

#image("coverage.svg", width: 100%)

#table(columns: (2.8fr, .8fr, 1.3fr, 1.3fr), inset: 4pt,
  [Source model], [Vertices], [Subgroup r], [Embedding k],
''' + table + r''')

*Observed:* native preparation covered all 24,823 vertices and 51 conductor
strata of the three additional complete graphs. Each input binds the exact
model, conductor, largest-prime subgroup and deterministic nonzero generator.
Curve membership and [r]G=O pass for every vertex. Original database SHA-256
values and every conductor population match the prior graph delivery.

*Derived by native arithmetic:* complete factorizations certify the selected
prime subgroup and the exact multiplicative order of p modulo r. The finite-field
destination has size p^k. No cheaper subgroup replaced the degree-3 class's r.

*Scope:* all-vertex direct EC relations and applicable finite-field IC, base
grids, independent one-target workloads, matched references, full cold/online
accounting, relation yields, rank, min/max/mean and explicit failures are required
by the committed protocol. Input certification does not execute that matrix.
Largest-prime-subgroup runs alone also do not measure full composite-order PH.

*Still open:* exact model/subgroup/generator representation registration,
subgroup/quotient base validation, audited native control/pilot, a resumable
full driver and provisioned execution/storage. Existing supersingular k=2 MOV
uses BSGS/PH, not finite-field IC for these ordinary large-k destinations.
This is a backend gap, not nonexistence of a mathematical reduction.

*Evidence:* evidence/inventory.json, inventory.stderr, implementation.sha256;
per-vertex JSONL records and hashes accompany the inputs-only archive.
Protocol commit 4253df8f2; native source commit e1b0e302e preceded execution.
Linux x86_64, 4 CPU quota, 16 GiB RAM. Initial link exhausted temporary disk;
stripped workspace relink succeeded. Failed outputs remain preserved.

*No performance claim:* solver runtime, operation distributions, relations and
yield remain unknown. Independent implementation replay is not performed.
IC/rho scoreboards retain their previous performance figures. Original A--E
requirements, including full A/C/D/E generation, remain open. The draft backlog
card is not a scheduled task; no solver job continues in the background.
'''
(OUT / "report.typ").write_text(source)
typst.compile(str(OUT / "report.typ"), output=str(OUT / "report.pdf"))
print(OUT / "report.pdf")
