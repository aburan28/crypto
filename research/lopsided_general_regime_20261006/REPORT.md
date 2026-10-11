# General-regime watch: resolved quantitatively (2026-10-06)

Questions and method in the [protocol](README.md). Rounds 1–2 built the
thin-product interface and screened sixteen targets against the concrete
Corollary-26 line (merged PR #1485; open PR #1493); round 3 audited the
refutation's foundations (open PR #1494). This round closes round 2's open
check 3 — the general regime of paper Section 4 — and works the remaining
Section 6 directions. Status labels: **[proposed]** / **[derived]** /
**[verified]** / **[measured]** (none in this round).

## Headline

**The envelope widens; nothing moves.** The proven regime runs from the
concrete line (`D <= N^{1/18}`, gamma `0.063`) out to `D <= N^{0.114}`
with gamma up to `0.2412` for very thin, very sparse wants — frozen with
provenance in `table2.json` **[verified]** against a rendered image of
paper p.40. All sixteen targets stay blocked at shape, ring, or scope,
gates the envelope cannot touch: **size never binds for any repository
workload, so widening size changes no verdict.** Sparse matrix
multiplication is closed twice over (the paper's own negative plus our
ring), and the reduction-loss remarks confirm rather than disturb our
end-to-end screening rule.

## Q1: the general envelope re-screen [derived]

Theorem 25 (paper p.40) says: for every eps < eps* = 0.1204 and every
kappa > 0 there is gamma > 0 with `O(N^2/D^gamma)` whenever `D <= N^eps`
and `|W| <= N^2/D^kappa`. The six proven rows (Figure 1; `table2.json`):

| c | eps (D <= N^eps) | gamma at kappa=1/2 | gamma at kappa=1 |
|---|---|---|---|
| 40 | 0.029 | 0.1060 | 0.2412 |
| 21 | 0.056 | 0.0653 | 0.1546 |
| 18-line (Cor. 26) | 1/18 | 0.063 | — |
| 19 | 0.062 | 0.0578 | 0.1379 |
| 15 | 0.079 | 0.0389 | 0.0941 |
| 12 | 0.099 | 0.0186 | 0.0459 |
| 10.5 | 0.114 | 0.0052 | 0.0129 |

Two structural facts matter more than any cell. First, larger `c` trades a
*smaller* eps for a *larger* gamma: the biggest savings sit at the
thinnest products, and the widest reach (eps `0.114`) buys only gamma
`0.0129`. Second, the envelope moves exactly one gate — size — while
every recorded block in rounds 1–2 sits at shape, ring, or scope:

| Target group (rounds 1–2 verdict) | Blocking gate | Envelope effect |
|---|---|---|
| Relation-verification batch; factor-base queries (closest fit); sparse-LA block; target-descent batch (4 insertion points) | Shape/ring/model: polynomial systems not inner products; algebraic yield not set-theoretic; matvecs over `Z/rZ` not small-integer products; per-target recursion | **None**: no (D, W) mapping exists to resize |
| Pair-table / MITM / BSGS / rho probes; Macaulay batch (#1–3, #6) | Shape: hash probes, no shared middle | **None**: no coordinates in the plane at all |
| Batch MSM checks (#4); isogeny-route search (#10) | Ring then scale-free size (eps ≈ 1); tiny unweighted graphs | **None**: eps ≈ 1 is off every row including the widest; graphs have no n to grow |
| Semaev/MITM vs 3SUM; SAT/F4 vs 3XOR (#7–8) | Missing reduction across domains | **None**: the envelope is not a reduction |
| Convolution, real-RAM, SETH/OV et al. (#9, #11–12) | Out of scope (absent path / wrong model / paper-excluded) | **None** |

Falsification check **[verified by inspection]**: no target's *only*
block is size under the concrete line — the cell that would overturn the
headline does not exist. The pre-registered expectation holds.

## Q2: sparse matrix multiplication is closed twice over [derived]

The paper asks the question itself (Sec. 6): much work exists on sparse
matrix multiplication, yet "our techniques do not seem to give new
algorithms for sparse matrix multiplication" — the balanced case its
method does not speed up. That is the authors' negative, not ours. Ours
stacks on top: our linear algebra (block Wiedemann in `koblitz_sparse_la`)
is not sparse-times-sparse multiplication at all but sequences of sparse
matvecs over `Z/rZ` — wrong operation *and* wrong ring (round-2 row #5).
Confirmatory greps **[verified]** (2026-10-06): `girth` returns nothing;
sparse-MM hardness language returns nothing; balanced/all-edges-sparse-
triangle vocabulary outside our own study files returns nothing. There is
no sparse-MM claim in the tree for the technique to disturb.

## Q3: reduction losses and gray boxes confirm the screening rule [derived]

Two Section 6 remarks bear directly on method. First, the hardness
reductions this paper reuses as algorithms were designed to keep *some*
polynomial saving (Exact Triangle keeps half the exponent saving; 3SUM
and APSP keep a half and a third through it), and "now that they are
being used as algorithms, how much they keep is more important." That is
our `S` rule in other words: only kept, end-to-end, measured savings
count, never headline exponents — a stage-count reduction that leaves
total `S` flat is relabelling (round-2 promotion protocol). Second,
Sheffield–Vassilevska Williams–Xi show the all-edges-to-detection third
is optimal for black-box reductions, so tighter transfers "must use the
structure of the particular problem" — exactly the shape of every
round-1/2 block ("needs a new reduction, not substitution"). The gray
boxes (problems merely 3SUM/APSP-*hard*, reductions going the wrong way,
no speedup) mirror round 3's audit: hardness-flavoured language that is
not a load-bearing equivalence is unaffected.

The technique family could reach our workloads only through *two*
independent breakthroughs: (a) a new base identity lifting the ceiling
past eps* = 0.1204 — the paper's open computer-search problem, since its
analysis leans on Schoenhage-specific sparsity and leaf-sharing
(Sec. 2.4, Sec. 6); *and* (b) a reduction reshaping an algebraic IC
workload into thin small-integer products with known sparse W. Either
alone changes nothing here. That conjunction is the standing watch
condition; it replaces the open-ended "watch the general regime."

## Method implications [proposed]

1. **Parameterize the gate on the frozen envelope — after PR #1493
   merges.** Exact next action: extend `screen_workload` with the six
   `table2.json` rows so a future candidate screens against the widest
   proven line, not just the 1/18 line. Acceptance: all existing gate
   tests unchanged, plus one test that the widest row still rejects the
   MSM operating point and the synthetic fitter still passes. Owned by
   whoever merges #1493; specified here so it is a task, not a wish.
2. **Screen transfers on kept saving, end to end.** The SVX27 lesson and
   the reduction-loss ledger belong in future protocols as one line: a
   transfer claim states the fraction of the exponent saving its
   reduction keeps, priced inside total `S`. No template change beyond
   that sentence.

## Boundary, table, ratio (AGENTS.md Secs. 1–8)

- Boundaries unchanged: generic floor + matched rho; unit `S`.
- No table row produced: envelope cells are the paper's numbers,
  transcribed — not DLP answers — and re-screen verdicts are gate
  evaluations, kept out of every result panel.
- Falsification target: a repo workload whose only block is size under
  the concrete line but which passes a general-regime row. The table in
  Q1 is the search that found none.
- All phases unpriced (`null`); wall time appears nowhere.

## Graphs checked — no change [verified by inspection]

- `docs/index-calculus-scoreboard.html` + `docs/ic/progress-timeline.json`:
  no new ratio, no point added.
- `docs/ic-leaderboard.html` / `LEADERBOARD.md` / `leaderboard.json`: no
  whole-pipeline measurement.
- `docs/browser/data.json`: no new identities of any kind.
- `docs/curves/registry.json`: no curve named; family notation only.
- Figure 1 below is a new visual for this search round (envelope data
  plotted from `table2.json`), not a canonical-graph update.

## Open checks [proposed]

1. Round-2 PR #1493 and round-3 PR #1494 are still open; neither affects
   this analysis (different files), but implication 1 above is sequenced
   on #1493's merge — re-verify the citation path then.
2. The gate-parameterization patch (implication 1) and the identity-claim
   certificate format (round-3 implication 1) are specified but unowned;
   owning either means writing the code/section, not proposing it.
3. If a new base identity lifts the ceiling past eps* = 0.1204, only
   condition (a) of the watch conjunction is met — re-screening is owed
   only together with a candidate reduction (b).

## Sources

- Paper: `https://arxiv.org/pdf/2610.06783` (v1); Table 2 and Theorems
  24–25 (p.40, Sec. 4.1); Sec. 5 (reductions survey depth); Sec. 6
  (future directions, quoted); SVX27 via the paper's citation.
- Frozen envelope: `table2.json` (this directory).
- Rounds 1–3: `research/lopsided_thin_product_20261006/` (merged PR
  #1485), `research/lopsided_other_speedups_20261006/` (open PR #1493),
  `research/lopsided_implications_20261006/` (open PR #1494).
- Audit receipts: the greps quoted in Q2, run 2026-10-06 on branch
  `cursor/lopsided-general-regime-35ee`.
