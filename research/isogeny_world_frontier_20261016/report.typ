#set page(paper: "a4", margin: 19mm)
#set text(size: 10pt)
#align(center)[#text(size: 19pt, weight: "bold")[Prime-degree isogeny frontier]]
#align(center)[Active campaign · initial record survey · 15 October 2026 PDT]

The goal is to push independently verified prime-degree isogenies on all registered standard prime-field curves, separately ranked. The first proposed expanded-map rung is P192/50341, above the recorded P192 degree 10453. No new map is claimed in this survey.

#image("frontier.svg", width: 100%)

= What the primary sources establish

Jao–Soukharev (2010), sections 5.1–5.2, report ideal-based evaluations at degree 5000000029 and on Certicom ECCp-109 at the least Elkies prime above p/2. They give codomains and public-point images with conjugate consistency checks. Section 5.3 records an ECCp-239 ideal relation but leaves evaluation unfinished. These observations were read from the paper; they were not independently reproduced here.

Square-root Vélu (2020) distinguishes evaluating an isogeny from emitting its expanded rational-map coefficients. The parallel implementation (2024) benchmarks degrees through 321193 on a chosen 1792-bit prime field. These are different representation and curve categories from our standard-curve expanded-map board. PRISM supplies a further compact-representation comparison to audit.

= Frozen execution and claim gates

The native trial will run only after its Coordinator-approved protocol and execution handoff are committed and pushed. The first batch specifies a P192/13 positive, a coefficient-mutation negative, and P192/50341. Construction is bounded at 1800 seconds, replay at 3600 seconds, sampled RSS at 8 GiB and map size at 256 MiB. A timed-out attempt has no score and establishes no nonexistence result. A later result gets a separate write-once run report.

The 54 standard source rows stay separate; 47 adapters exist and seven conversions remain unavailable. No world record, performance gain or cryptanalytic conclusion is established. IC/rho and performance charts are unchanged because this survey contributes no such measurement.

= Sources and artifacts

#link("https://arxiv.org/pdf/1002.4228v2")[Jao–Soukharev, version 2 (2010), sections 5–6].
#link("https://arxiv.org/pdf/2003.10118")[Bernstein et al., square-root Vélu (2020), introduction].
#link("https://www.preprints.org/manuscript/202401.1366")[Parallel square-root Vélu (2024), section 4].
#link("https://r98inver.github.io/papers/2501_prism.pdf")[PRISM author paper].

Editable report and diagram: research/isogeny_world_frontier_20261016. The README retains the request, exact comparison scope, unresolved checks and execution gates. Runtime goal anchor: 01a11ebd-9436-73b3-b24e-63b10aad80fa.

#text(size: 9pt)[PRISM comparison caveat: the basic 2025 response has composite degree q(2^a-q). The #link("https://link.springer.com/article/10.1007/s00145-026-09585-6")[2026 journal extension, section 5] separately defines prime degree q=2^(2e-1)-s². A prime factor is not a prime total degree. Numeric-instance extraction and independent evaluation remain outstanding.]
