# Pair-46 factor-base campaign (2026-10-15)

This campaign asked whether any factor-base family proposed for pair 46 gives
an index-calculus advantage on its prime-field CM(-19) curves. Pair 46 has
`p = 5t^2 + 19t + 19` with `t = 34489701`, a 53-bit prime `N = p + 1 - t`, a
`j = -96^3` crater curve and a floor curve. The campaign tested 14 families.
Each item has its own preregistered `PROTOCOL.md`, data, `RESULTS.md`, graph,
`report.pdf` and an independent red-team review (`REDTEAM.md`). Item 11 is the
exception: it did not run.

Every verdict below is the red team's corrected verdict. The red team upheld
12 of the 13 items that ran and refuted one: item 14's `exact_obstruction`
becomes `empirical_negative`.

- Summary graph: [`fig_campaign_summary.svg`](fig_campaign_summary.svg)
  ([PNG](fig_campaign_summary.png)). Source:
  [`plot_summary.py`](plot_summary.py). Data:
  [`summary_data.json`](summary_data.json).
- Campaign PDF: [`report.pdf`](report.pdf), built from
  [`report.typ`](report.typ).

## Summary table

| # | item | verdict (red-team corrected) | headline | red team | link |
|:--|:--|:--|:--|:--|:--|
| 01 | endo-chain | negative | Best arm (single ψ-chain) costs S = 1.47 [1.28, 1.69] per √N against ±BSGS at 0.955 [0.86, 1.05] on the 22-bit CM(-19) toy (ratio 1.54 [1.30, 1.82]). One ψ x-evaluation costs 0.98× an affine addition on pair 46. Exact: a known-head chain is BSGS with table entries costing c_E ≈ 1, so it cannot beat BSGS. | upheld. The rho-neg reference is inflated by restart overhead that grows with log N. The swap of the preregistered ratio_to_rho_neg should have been flagged as a deviation. | [01-endo-chain](01-endo-chain/) |
| 02 | deuring-lift | exact_obstruction | The lift is Cremona 361a1, and E(Q) = Z·(0,9) has rank 1. Reduction is onto E(F_p), so lifting a target is the DLP. Twist generators (115 twists, maximum rank 3) give 0 Q-relations. Their mod-p lattice is indistinguishable from random (λ1/GH +0.14 [-0.18, +0.26]). | upheld. The minimal-height lift is about 5.3e29 nats, not 2.1e30. Power is not shown; the verdict rests on the exact argument. | [02-deuring-lift](02-deuring-lift/) |
| 03 | floor-pullback | exact_obstruction | φ is an isomorphism when N is prime, so a pushed base has exactly the crater yield. On 8 toys the yield ratio against random-x was 0.9978 to 1.0023 (\|z\| ≤ 1.40); the positive control had \|z\| ≥ 562. Membership costs one dual-isogeny evaluation. | upheld. X0 is a consistency check, not independent evidence. The ~6e3 field-op cost is an order-of-magnitude lower bound. | [03-floor-pullback](03-floor-pullback/) |
| 04 | family-sweep | null | The family p = 5t²+19t+19 is within the control distribution at 20 to 30 bits. M1 ratio 0.984 [0.964, 1.005], M2 Δ -0.006 [-0.016, +0.004], M3 1.0011 [0.998, 1.004]. All positive controls were detected (p ≤ 0.0005). | upheld, scoped to random models of j = -96^3 crater curves, the small-x base and m = 2. Resolution is about ±2% in M1. | [04-family-sweep](04-family-sweep/) |
| 05 | qr-window | constant_only | Exactly the predicted 2× density gain: 2.00004 [1.99967, 2.00042] on the pair-46 crater window. Relations are identical to the plain window on 8000/8000 toy targets. Against any base that already tests liftability, the gain is exactly 1. | upheld. The 2× is bookkeeping. The positive control agrees with itself by construction; the independent checks are the E1 counts and the A2/A3 halving controls. | [05-qr-window](05-qr-window/) |
| 06 | pkm-rational-map | negative | PKM chained bases reproduce the low lcm degree (4 to 5), but one system costs about p^2.0 in our Buchberger (PKM's Magma times give p^1.51; rho gives p^0.52). The projected relation stage is about 4e6× rho at n1 = 5. Yield equals the flat set exactly. | upheld. The slope is a local toy-range growth rate, not an exponent. The X3 "only subgroup" wording needs fixing. | [06-pkm-rational-map](06-pkm-rational-map/) |
| 07 | torsion-division | exact_obstruction | No ψ_m (m ∈ {2,3,5,7,11,13}) has an F_p-root on E. ψ_5's 2 roots are on the twist. Preimage bases change coverage by exactly 0 and multiply the degree by m². Every arm has size/degree ρ ≤ 0.25. | upheld, scoped to the tested arms. The LVL mean size is 0.47 to 0.56, not 0.6. | [07-torsion-division](07-torsion-division/) |
| 08 | learned | exact_obstruction (counting and success only; empirical part formally inconclusive) | Learned selectors match random on held-out targets within about ±0.5%. Exact identity: Σ_R c(R) = \|S\|^m. Held-out m = 2 success is at most L(2L+1)/(N-1-T), so any advantage is a constant of about μ/(1-e^-μ). | upheld with corrections. Greedy is below random in 5 of 6 cells, not all 6. The random2 control failed at 16 bits. Bases chosen to make solving cheaper are not covered. | [08-learned](08-learned/) |
| 09 | cm-norm | exact_obstruction | A norm-form base is a 2-D generalized arithmetic progression with F_X + F_X ⊆ F_{4X}. On toys its coverage is 0.0159 [0.0156, 0.0161]× random at B = 256 and 0.0079 [0.0064, 0.0089]× at B = 1024. The floor reduces to Z, and class-group ideals never decompose a point. | upheld. The toy λ comes from the roots of x²-x+5. The Θ(√N) decomposition cost is heuristic; the sumset bound is the exact part. | [09-cm-norm](09-cm-norm/) |
| 10 | bounded-degree | exact_obstruction | Every bounded-degree (d ≤ 8) algebraic base has \|F\|/p ≥ 0.167 (168 rows). The extension-field control has exponent about 0.5. On pair 46, ρ̂ ≈ 0.63 and membership costs about 160d² multiplications. X3: \|F\| ≤ √p needs d ≥ √p/5. | upheld. The "exponent ≥ 0.839" is a constant fraction at toy p, not a measured exponent. The result covers only the tested families. | [10-bounded-degree](10-bounded-degree/) |
| 11 | ff-xedni | blocked | Did not run. Candidate B2 and function-field xedni remain open (EV-XEDN-001, GOAL-ICLIFT-001). | confirmed blocked, with no artifacts | none (no directory) |
| 12 | model-small-x | exact_obstruction | Pair 46 and its twists have no rational 2- or 3-torsion, so no Montgomery, Edwards or Hessian model exists over F_p. On 240 toy curves all 21 tests are indistinguishable from random (Bonferroni CIs within ±1%, within ±0.10% at b = 22). | upheld. The offset of negation-closed sets is about 0.77/K by a first-order heuristic, not exactly 1/K. | [12-model-small-x](12-model-small-x/) |
| 13 | xy-windows | null | Sign and parity y-windows are exactly the x-only window (8000/8000). Matched-size \|y\|-box and y-bit windows give D/D(X) of 0.986 to 0.995 with CIs in [0.956, 1.026], at 8 to 128 x-scans per element against 2. | upheld. No planted control tests for a yield gain. The E3 pass uses a post-hoc exact binomial null. | [13-xy-windows](13-xy-windows/) |
| 14 | dor-verdict | **empirical_negative** (red team refuted `exact_obstruction`) | Binary toy model (m = 2, n ≤ 23). The dimension abandon rule fires for generic bases in 14 of 16 groups (0.86 to 1.36 degrees per doubling) and never for the standard window. For prime n, r ≥ min(n, 2n'-1) is exact, and the standard window attains it. An (l,q)-matched null reproduces 40 of 45 cells. | **refuted.** The degree-2 Macaulay rank exceeds 1+l+q on 937 of 4638 targets (E1 holds only for the equation span). The headline range was wrong (the true range is 0.86 to 1.36). The 5 null failures are one-sided. Does not transfer to pair 46. | [14-dor-verdict](14-dor-verdict/) |

## Campaign-level conclusion

The conclusion is stated only as strongly as the evidence allows.

- **No item found an index-calculus advantage on pair 46, or on its toy
  analogues, from the choice of factor base.** No arm beat its matched
  reference in a direction that would help an attack.
- **Exact obstructions (02, 03, 05, 07, 08, 09, 10, 12) are scoped.** Each
  closes the family as tested: the stated metric (relation count, m = 2
  coverage or yield, membership degree) and the stated arms. None of them
  shows that "no factor base can help pair 46".
  - Item 08's identity Σ_R c(R) = \|S\|^m caps relation counts for uniform
    targets under any selector. It does not cover bases chosen to make the
    decomposition solve cheaper.
  - Item 10 does not exclude compact high-genus descriptions.
- **Empirical results are toy-scale (12 to 30 bits, m = 2) unless they say
  otherwise.** Items 04 and 13 are nulls at a resolution of about ±2 to 3%.
  Item 14 is a binary-field toy model and does not transfer to the prime-field
  pair 46.
- **Negative cost results (01, 06):**
  - Item 01's ψ-chain is a constant factor (about 1.5×) worse than ±BSGS.
    This is explained exactly by c_E ≈ 1 ≥ 1/2.
  - Item 06's chained systems grow about p^2.0 per system in a pure-Python
    Buchberger. Only the exponent gap with rho transfers; the 4e6× figure
    carries Python constants.
- **Not tested:** m ≥ 3 decompositions (except where an exact identity covers
  them), relation collection with linear algebra at pair-46 size, F4/msolve
  solvers, and item 11's function-field xedni.

## Positive-control status (item 05)

Item 05 is the campaign's instrument-level positive control for density
effects. It detected exactly the predicted effect, with a 2× density gain of
2.00004 [1.99967, 2.00042] on the pair-46 crater window. Its negative
controls (the χ(x) filter and a random half) give 1. The red team notes that
the G(A1) = 1/M1(A0) and 8000/8000 agreements hold by construction, because
they share the jacobi(f(x)) value. The independent validation of the
instrument is:

- the E1 full-field counts on 36 toy curves;
- the A2/A3 halving controls, which keep about 1/4 of the m = 2
  decompositions.

Later items should cite the A2/A3 controls as their sensitivity evidence, not
the ±0.02% Wilson width. Every other item that ran also had its own positive
control detected. The weak ones are:

- 08: the random2 null-vs-null control failed at 16 bits, which makes the
  empirical part formally inconclusive.
- 06: the preregistered sugar-degree match failed; the Amendment-1
  lcm-degree match passed.

## Blocked and cut

- **Blocked:** item 11 (ff-xedni). It produced no protocol, data or report,
  and candidate B2 and function-field xedni remain open.
- **Cut or not run, per item:**
  - 01: no cuts. Two post-hoc supplementary walks were added and are labelled
    as such.
  - 03: no isogeny evaluation at pair-46 size. The preregistered population
    held only 3 instances, so the run was extended to 8 under the frozen
    `PROTOCOL_EXT.md`.
  - 04: no floor-curve family members. No M1/M2 at 53 bits; only M3 was run
    at pair-46 scale.
  - 06: the chained arms stopped at n1 = 5 by the preregistered stop rule. No
    F4, msolve or Magma was available.
  - 07: no Gröbner invariants were computed, because no
    Sage/Singular/M2/Magma was on the host.
  - 08: m = 2 toy only.
  - 09: no membership-test cost was measured.
  - 12: the pair-46 curves themselves cannot be run because of the exact
    obstruction. Jacobi-quartic and Huff models and m ≥ 3 were not tested.
  - 13: m ≥ 3 was not tested.
  - 14: Sage was absent, so the m = 2 engine was ported to Rust and gated by
    an exact 4638/4638 replay. The m = 3, sym-S4, stacking and n = 131 checks
    were not ported. n = 23 and some n = 19 cells were not replayed.
- **Process:** no item used ecbench sessions. The standalone instruments are
  under each item's `src/` or scripts. Timings come from a shared host, at
  isolation level L0.

## Performance copies

Items 01 and 06 report timings: the ns cost of a ψ evaluation against an
affine addition, and solver time against p. Their PDFs and graphs are copied
into
[`benchmarks/20261015-pair46-factor-bases/`](../../benchmarks/20261015-pair46-factor-bases/),
together with the campaign summary.
