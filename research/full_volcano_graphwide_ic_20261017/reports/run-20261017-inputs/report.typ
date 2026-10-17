#set page(paper: "a4", margin: 15mm)
#set text(size: 9pt)
#set heading(numbering: none)
= Full graph-wide IC: certified workload preparation
Run 20261017-inputs. This report contains *zero solver executions*.

#image("coverage.svg", width: 100%)

#table(columns: (2.8fr, .8fr, 1.3fr, 1.3fr), inset: 4pt,
  [Source model], [Vertices], [Subgroup r], [Embedding k],
[icv1-fp57-t536870914-6cabb0c9], [8192], [61751], [30875],
[icv1-fp57-tm536870941-e6d0c766], [6561], [19768889568671], [3953777913734],
[icv1-fp59-t1073741954-4dd0b1ef], [10070], [4253], [2126],)

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
