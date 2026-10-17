#set page(paper: "a4", margin: 1.8cm)
#set text(size: 9.5pt)
#show table: set text(size: 8.5pt)
#let md = read("RESULTS.md")
= 14-dor-verdict: does factor-base dimension or structure change D\*?
_Pair-46 factor-base campaign, item 14. Verdict on research/factor\_base\_dor\_campaign\_20261013. Protocol frozen 2026-10-15T17:00Z (PROTOCOL.md)._

*Verdict.* Dimension raises D\*. For generic bases it rises by at least one degree per doubling of n', so the campaign's own abandon rule fires. At matched n', structure acts on D\* almost entirely through two exact invariants: l = dim(V1+V2) and r = dim(V1·V2). The invariant-matched null reproduces 40/45 cells within ±0.25. For prime n, r ≥ min(n, 2n'−1) holds for every subspace (linear Cauchy–Davenport), and the plain polynomial window already attains it. So in this model no structure beats the standard window, and a subfield-closed base cannot exist when n is prime.

*Replay gate.* A dependency-free Rust port of the sage m=2 engine reproduced the recorded per-target d\_star and rank\_profile on 4638/4638 targets in 62 cells, with 239/239 arms identical in bases, targets, fractions and points. The full numbers, intervals, deviations and cuts are in RESULTS.md, reproduced below.

#figure(image("fig_dstar_vs_dim.svg", width: 100%), caption: [Recorded D\* (mean, 95% cluster bootstrap) against n', by structure and field.])
#figure(image("fig_dstar_vs_invariant.svg", width: 72%), caption: [Cell-mean D\* against the exact invariant (l + r)/n. All structures collapse onto one curve.])
#figure(image("fig_null_vs_semaev.svg", width: 60%), caption: [Semaev D\* against the N2 invariant-matched null (200 systems per cell). The band is ±0.25.])
#pagebreak()
#raw(md, lang: "markdown", block: true)
