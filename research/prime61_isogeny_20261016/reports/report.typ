#set page(paper: "a4", margin: 17mm, numbering: "1")
#set text(font: "DejaVu Sans", size: 9.5pt)
#set par(justify: false)
= Prime61 horizontal-isogeny panel
2026-10-16 • Partial requested experiment; completed finite outgoing-edge panel.

All 15 outgoing kernel-distinct rational isogenies of prime degree at most 61 from the fixed source were constructed, independently kernel-verified and tested with the transported target. All 192 relation streams reached their cap with zero relations and zero solves. The full horizontal class remains unexhausted. There are no vertical separable isogenies defined over this base field in the class.

#image("coverage.svg", width: 100%)

== Exact input and finite coverage
Prime field and curve: #raw("p=2305843009213693951; a=2305843009213693948; b=111"). Group and subgroup order: #raw("r=2305843010818082053; h=1"). Source G: #raw("(1,1796364306280342617)"); Q: #raw("(992726311114892094,558260110083040833)"). This synthetic target satisfies Q = [539365399]G. The secret is used only for checks, never supplied to the collector. This is a new large-field fixture, not a transported screenshot instance.

The panel covers one ramified degree-3 edge and two split edges at each of degrees 7, 11, 13, 19, 23, 53 and 61. Degree 2 has no rational kernel. The inert primes 5, 17, 29, 31, 37, 41, 43, 47 and 59 have no rational edges. Degree 5 was also tried in the earlier degree-13 construction control. Distinct kernels and 15 distinct destination j-invariants are retained. Quotient maps are identified up to codomain field isomorphism; larger degrees, further hops and arbitrary postcompositions are outside this panel.

#pagebreak()
== Structure and transport
The trace is −1604388101 and the Frobenius discriminant is −6649310858224389603. Its absolute value factors as 3 × 97 × 113 × 35831 × 5643470711. The native driver reconstructs the product and deterministically checks all factors prime. The discriminant is squarefree and 1 modulo 4, hence fundamental. Since every ordinary base-field endomorphism order contains Z[pi], its conductor divides the Frobenius conductor. Both conductors are therefore one throughout the class. A vertical base-field edge would change this conductor, so none exist here. Extension-field and inseparable maps are outside this statement.

The native walker emits explicit kernel polynomials and exact models. Its separate verifier passed 16 curves and 15 edges with two audit points per edge. The new driver rechecks squarefreeness, division-polynomial divisibility and cyclic closure, evaluates the existing Vélu formula through traces in F_p[T]/h(T), and applies the target-model isomorphism. Both G and Q traverse that same map. The destination generator is the image of source G, not the independent generator chosen by the walker.

Each map passed 64 scalar checks, including infinity, plus the frozen target relation and subgroup checks. The 36 signed pair sums from the transported source base also match as multisets on each map. These are transport controls, not a new general proof of Vélu's theorem. The collected workloads use fresh native low-x three-pair bases on each destination. The example unit test compares degree-3 evaluation against the separate existing small-degree implementation, including a kernel pole.

== Recorded relation-stage result
Seeds 2026101600 through 2026101611 are paired across the source and every destination. Each stream requests 20 relations and stops at 100,000 target trials. The native summary checks all curve IDs, stream IDs, seeds, caps and the complete edge-degree histogram.

#table(columns: (3fr, 2fr), inset: 5pt, [Quantity], [Recorded value], [Workloads / streams], [16 / 192], [Capped streams], [192], [Target trials], [19,200,000], [Relations / verified solves], [0 / 0], [Scalar / pair-sum checks], [960 / 540], [Full cost, S, rho reference, speedup], [Unknown])

Every neighbour ties at zero relations within this cap. A second identified execution is byte-identical, and collection rows match the initial execution after adding EC1 metadata. This is same-host reproducibility. The earlier 96 degree-13 pilot streams and repeated executions are not pooled into these counts. The earlier source-only study provides a separate solved 20-bit control.

No complete large-field solve, candidate ranking or performance improvement is established. Target trials are stage counters; map construction, transport, setup and primitive accounting are not priced. Diagnostic unisolated timings in raw walker files are not performance results. Canonical cost dashboards and leaderboard graphs remain unchanged.

#pagebreak()
== Full class: certified lower bound and resource blocker
PARI/GP 2.15.2 returned the agreeing class number 545296208 from three routines, with conditional group structure [68162026, 2, 2, 2]. Full BNF certification hit the frozen 300-second outer cap (exit 124). The full count remains GRH-conditional; the Shanks output alone is not a proof at this discriminant size.

A separate exact reduced-form certificate establishes a lower bound without assuming that full count. For A = Qfb(7,5,237475387793728201), n = 68162026 = 2 × 34081013, exact composition gives A^n principal while A^(n/2) and A^(n/34081013) are nonprincipal. Both factors are verified prime. Any proper divisor of n would divide at least one of those exponents, so A has exact order n. The standard ordinary maximal-order class action gives at least n base-field isomorphism classes in this trace class.

The current walker retains Node and Edge vectors, whose native header sizes are 240 and 152 bytes. Even n nodes with n−1 connecting edges need 26,719,514,040 bytes, excluding all owned heap allocations, kernels, maps and hash tables. That exceeds the container limit of 17,179,869,184 bytes. The initial node-only estimate did not exceed the limit; the additional edge-header bound is retained separately. A streaming or external-memory implementation could change this blocker. No exhaustive traversal or remote fleet was launched.

The exact-order script's first interactive attempt failed because I is a reserved GP variable and printed a misleading symbolic PASS. That transcript is explicitly invalid and preserved with stderr. The replacement script checks types and values and fails closed on error; its successful stderr is empty. Full class coverage remains open in problems/prime61-horizontal-class-memory.md and the linked blocked backlog card.

== Evidence, build and remaining obligations
Source library revision: #raw("de0c9391480155ebdd5a8495d2c79e4508d4795d"). Protocol revisions: #raw("596cae7d,14298786,da85e554,95ec4bec"). Initial driver: #raw("764d53c0"); measured EC1 metadata: #raw("95ec4bec"). The exact measured EC1 and full UID of every workload are in the curve rows; model registration in ICV1 is pending. Walker UIDs with independent generators are not measured workload IDs.

Native walk verification, transport checks, example test and deterministic replay passed. Earlier unchanged-library release tests passed 4,154 with 117 ignored; this is not a fresh whole-library run. Strict Clippy found an existing unrelated mismatched-width expression; the rerun allowing that one lint passed. Failed compilation and tool attempts are retained. Conductor and Lean are unavailable. No Lean proof of the CM correspondence or Vélu theorem is claimed; the arithmetic certificates do not formalize those bridges.

Reproduction and raw records: README.md; reproduce.sh; PROTOCOL.json and addenda; one-hop-l61/isogeny_routes.json; results/walk-l61-verify.json; results/factor-base-l61-identified.jsonl; results/factor-base-l61-replay.jsonl; results/summary.json; class_order_certificate.gp; results/class-order-certificate.txt; results/min-class-node-edge-memory.json; SHA-256 lists in results. Paths are relative to this study. The source-only baseline is a separate sibling study. Host: AMD EPYC 9V74, five visible CPUs, Linux x86_64, 16 GiB limit; no wall-time admission.

References: Sutherland, Isogeny volcanoes (arXiv:1208.5370), for standard CM and volcano classification; PARI/GP arithmetic and general-number-field manuals for form composition and conditional versus certified class groups. This report records known formulas and exact finite computations, not a new polynomial identity or a proof that the full experiment is mathematically impossible.
