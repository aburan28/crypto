#set page(paper: "a4", margin: 1.8cm)
#set text(size: 9.5pt)
#set heading(numbering: "1.")
#show link: set text(fill: rgb("#2a6fdb"))

= Literature search, 2025–26: ECDLP, index calculus, Gröbner solving
*Date* 2026-10-11 · *Type* literature search round, no experiment run · *Status* every statement is _proposed_, taken from abstracts that were opened and checked; nothing measured.

*Headline (negative).* None of the 61 classified ePrint/arXiv papers (2024-10 to 2026-10) claims or plausibly delivers an exponent-level improvement for ECDLP. Nothing in the window addresses summation polynomials, point decomposition, Weil descent for ECDLP, or GHS. No new break or trapdoor was found for P-curves, secp256k1, Brainpool or the sect curves.

#figure(image("landscape.svg", width: 100%), caption: [Papers by repo component and expected impact, N = 61, counts from `papers.csv`. The classification is a judgement made from the abstracts.])

== What matters most
- *Structured generic-group model* (ePrint 2026/384, EUROCRYPT 2026). Algorithms that exploit structure on a δ fraction of elements need time Ω(min(√q, 1/δ)), with analogous bounds stated for EC points. This formalises the repo's "Koblitz IC is generic" verdict.
- *Non-adaptive preprocessing bound* (arXiv 2505.00894, STOC 2026). Non-adaptive advice cannot beat √N online time unless S = Ω(√N).
- *Solving-degree diagnostics* for descended S#sub[m+1]:
  - Kudo–Yokoyama, saturation exponent and first-fall synchronisation (ePrint 2024/528, 2025/822).
  - Multi-homogeneous XL (ePrint 2025/2060).
  - Symmetric-function solving-degree proofs (arXiv 2609.14052).
  - Generic leading monomials (arXiv 2505.10246).
- *Constant-factor engineering:*
  - lane-batched tracing (arXiv 2607.06372);
  - early abort on inconsistency (arXiv 2603.24890);
  - Crossbred trimming (ePrint 2026/360, 2026/408, 2024/758);
  - Wiedemann-XL (ePrint 2026/1787);
  - TNFS Galois k#super[2] linear algebra (ePrint 2026/560), likely already captured by GGMP n#super[2];
  - binary-curve x-only isogenies (ePrint 2026/704);
  - co-Z table-free BSGS (ePrint 2026/2102).
- *S5 small-root oracle feasibility on paper:* cuso (ePrint 2024/1577) and the LP asymptotics (ePrint 2026/1027).
- *Flagged as unsupported:* ePrint 2025/015, arXiv 2607.09814, a 16-bit D-Wave ECDLP result, arXiv 2607.09163.

== Proposed falsifiers (ranked)
+ Saturation exponent and first degree fall against the measured solving degree, descended S3/S4, n = 11–23.
+ Multi-graded Hilbert function of descended S3 with point blocks, n ≈ 15.
+ δ table: measured IC cost against min(√r, 1/δ) for every configuration. A cost below the bound would indicate an accounting bug.
+ Early abort at the F4 refutation degree.
+ Lane-batched inherited-F4 replay, k = 4–8.
+ E/G check: FGHR-symmetrised S#sub[m] on E against x-only S#sub[m] on E/⟨T⟩.
+ LP/cuso asymptotic bound for S5. A bound below p#super[1/5] closes the route.

Canonical graphs checked and unchanged, because there is no new measured value: the IC scoreboard, `progress-timeline.json`, `docs/curves/`, `docs/performance-gains/` and `figures/`. Coverage caveats: ePrint search matches titles and keywords and was rate-limited (HTTP 429); ePrint 2026/1198 and 2025/271 were not opened. The full discussion is in `README.md`.

== Classified papers
#let rows = csv("papers.csv")
#set text(size: 7pt)
#table(columns: (auto, 1fr, auto, auto, auto), stroke: 0.3pt + gray, inset: 3pt,
  table.header([*ID*], [*Paper*], [*Date*], [*Component*], [*Impact*]),
  ..rows.slice(1).map(r => (link(r.at(1))[#r.at(0)], r.at(2), r.at(3), r.at(4), r.at(5))).flatten())
