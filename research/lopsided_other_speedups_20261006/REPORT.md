# What else can the lopsided technique speed up? Survey report (2026-10-06)

Question and method in the [protocol](README.md). Round 1 (interface +
four IC insertion points, all proposals) is in
[../lopsided_thin_product_20261006/](../lopsided_thin_product_20261006/REPORT.md).
Statement status is labelled **[proposed]** / **[derived]** / **[verified]**
/ **[measured]** (none in this round).

## Bottom line

Twelve candidate transfers were screened. **Zero repository workloads pass
all three gates; the reusable gate itself is the positive deliverable.**
Five fail the shape gate (hash probes, not inner products), two fail the
ring gate, none from the repo reaches the size gate with a valid shape, and
five are out of scope (paper-excluded or absent from the ECDLP path). The
one passing shape is a synthetic thin-integer workload proving the gate is
not vacuous **[verified]** by unit test. This is a negative result stated in
gate order, not in vibes: each row names the exact gate a future candidate
must clear.

## Ranked table

| # | Candidate (paper source → repo target) | Verdict [derived] | Failing gate + reason |
|---|---|---|---|
| 1 | Hinted-OMv workflow → online pair-table probes (`orbit_pair_table` m=2/m=3; `native_signed_mitm` batched lookups) | Negative | **Shape**: frozen table + online queries matches the *workflow*, but each probe is a hash lookup over an unstructured key space (`QUERY_BATCH = 4096` batched lookups in `native_signed_mitm.rs`; `table.probe` in `orbit_pair_table.rs`), not an inner product over a shared middle — nothing to prune |
| 2 | — → BSGS giant-step lookups vs baby table | Negative | **Shape**: hash probes; the table is also large (√n), the opposite of thin |
| 3 | — → rho distinguished-point / collision tables (`pollard_rho`, `koblitz_strong_rho`, `preprocessing_rho`) | Negative | **Shape**: hash probes on walk states |
| 4 | Thin batch verification → batch MSM relation checks | Negative | **Ring** (then size): the shared middle exists (factor-base points) but the work is group operations, not small-integer products; and at the standard operating point R≈F relations over F points, epsilon≈1 scale-free, far above 1/18 |
| 5 | — → block-Wiedemann sparse LA (`koblitz_sparse_la`; round-1 re-screen) | Negative | **Ring**: large sparse matvecs over Z/rZ, not dense small-integer products |
| 6 | Rectangular MM → F4/F5 Macaulay batch reduction | Negative | **Shape/Ring**: polynomial matrices over fields; no shared thin integer middle |
| 7 | 3SUM O(n^1.9992) → Semaev/MITM point decomposition | Negative | **Reduction missing**: 3SUM is over integers; our MITM is over curve points with algebraic constraints — different domains, no transfer without a new reduction |
| 8 | 3XOR truly subquadratic → Boolean SAT/F4 systems | Negative | **Reduction missing**: 3XOR instances are not our polynomial systems |
| 9 | MonoConvolution O(n^1.5−ε) → convolution stages | Out of scope | No convolution stage exists in the ECDLP path (only NTT-adjacent comment in `semaev_decomp.rs`) |
| 10 | (min,+)/APSP O(n^2.9995) → isogeny-route search (`isogeny_walk`, volcano experiments) | Negative | **Scale + model**: our isogeny graphs are tiny (hundreds of vertices) and unweighted; the paper's polynomial gain needs large n with polynomial integer weights |
| 11 | Real-valued Las Vegas 3SUM/APSP/Exact Triangle (§5.2) | Out of scope | Finite-field repository; the real-RAM model does not occur here |
| 12 | SETH, OV, k-SUM/XOR k≥4, unhinted OMv, 3SUM-indexing | Out of scope | Explicitly unaffected per the paper (§1.2, Fig. 1); recorded so no follow-up re-screens them without new evidence |

The gate evaluations are exact closed-form checks **[verified]**:
`screening_passes_a_thin_integer_shape` (N=2^18, D=2, |W|=N²/2 →
`FitsConcrete`, epsilon=1/18, kappa=1),
`screening_rejects_hash_probe_work` (`NoSharedMiddle`),
`screening_rejects_wrong_ring` (`WrongRing`),
`screening_rejects_dense_regime` (R≈F MSM shape → epsilon=1 at scales 64
and 1024), `screening_rejects_oversized_wanted_set` (kappa=0<1/2).

## Screening map

<figure>
<object data="screen_map.svg" type="image/svg+xml" style="width:100%">screening map (see screen_map.svg)</object>
<figcaption>Figure 2 — the (epsilon, kappa) plane. Green: paper's concrete
pass region. Plotted points are exact gate evaluations of workload shapes
(not measurements): the synthetic fitter, the MSM operating point
(epsilon=1 at every scale), and the oversized-wanted demo. Structural
rejects (#1–3, #6) have no coordinates and are listed, not plotted.</figcaption>
</figure>

## Promotion protocol for any future passer [proposed]

A candidate returning `FitsConcrete` earns a follow-up measurement
protocol, not a claim: freeze curve/subgroup/base/targets/seeds, run
matched baseline/candidate full-DLP suites with every phase priced in `S`
against the generic floor and matched rho, classify per AGENTS.md §3, and
update scoreboard/leaderboard/browser in that PR. A stage-count reduction
that leaves total `S` flat or worse is relabelling.

## Boundary, table, ratio (AGENTS.md §§1–8)

- Boundaries stated before measuring: generic floor + matched rho,
  unchanged; unit `S`.
- No table row produced: gate evaluations are not verified answers to a
  DLP and never enter a result panel.
- Falsification target: a repo workload returning `FitsConcrete` with a
  concrete (N, D, W) mapping would overturn this round's headline and open
  the promotion protocol above. Absent that, no thin-product speedup exists
  in the pipeline.
- All phases unpriced (`null`); wall time appears nowhere.

## Graphs checked — no change [verified by inspection]

- `docs/index-calculus-scoreboard.html` + `docs/ic/progress-timeline.json`:
  no new ratio, no point added.
- `docs/ic-leaderboard.html` / `LEADERBOARD.md` / `leaderboard.json`: no
  whole-pipeline measurement.
- `docs/browser/data.json`: no new curve, session, round, or IC1 identity.
- `docs/curves/registry.json`: no curve named; no retired spellings
  introduced (checked: family notation only).
- Figures 1–2 below are new visuals for this search round, not
  canonical-graph updates.

## Open checks [proposed]

1. Any future batch workload with N small-integer inner products over one
   shared D-dimensional middle should be run through `screen_workload`
   first; the call and verdict belong in its protocol.
2. Candidate #4's second failure (epsilon≈1) is scale-free only at R≈F;
   a workload with R≫F² (many relations over a tiny base) would move
   epsilon — no such workload exists in the pipeline today.
3. The paper's general regime (epsilon<0.1204, size-dependent gamma) is not
   encoded in the gate; a candidate failing the concrete 1/18 line but
   plausibly inside the general regime needs a Section-4 analysis first.

## Sources

- Paper: `https://arxiv.org/pdf/2610.06783` (v1); consequences in Fig. 1
  (p. 9) and §§5.1–5.4.
- Round 1: `../lopsided_thin_product_20261006/` (interface, four insertion
  points, merged PR #1485).
- Code: `screen_workload` in
  `src/cryptanalysis/lopsided_thin_product.rs` (this PR).
- Visuals: `figure.mmd` / `figure.svg`, `screen_map.svg`, embedded in
  `report.html` / `report.pdf`.
