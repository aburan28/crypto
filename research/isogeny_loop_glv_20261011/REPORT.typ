// Rendered with Typst 0.15 (Python wheel, scratch venv); see README.md §6.
#set document(title: "Isogeny-loop GLV, round 2", date: datetime(year: 2026, month: 10, day: 11))
#set page(paper: "a4", margin: (x: 1.8cm, y: 1.7cm), numbering: "1 / 1",
  footer: context [#set text(7.5pt, fill: rgb("#52514e")); research/isogeny_loop_glv_20261011/REPORT.pdf #h(1fr) #counter(page).display("1 / 1", both: true)])
#set text(font: "DejaVu Sans", size: 9pt, fill: rgb("#0b0b0b"))
#set par(justify: true, leading: 0.55em)
#set heading(numbering: "1.")
#show heading.where(level: 1): set text(11pt)
#set table(stroke: 0.4pt + rgb("#c3c2b7"), inset: 4pt)
#show table: set text(7.6pt)
#show table: set par(justify: false)

#align(left)[
  #text(15pt, weight: "bold")[Isogeny-loop GLV, round 2: open leads, an erratum, and 4-dimensional GLV-GLS]
  #v(2pt)
  #text(8.5pt, fill: rgb("#52514e"))[2026-10-11 · study directory `research/isogeny_loop_glv_20261011/` · follows `RESEARCH_ISOGENY_LOOP_GLV.md` (PR \#1644) · status words: *measured* (native run recorded it), *derived* (follows by the stated argument), *proposed* (untested)]
]

= Corrections to round 1

*Language rule.* Round 1's sweeper is Python (cypari2), in breach of AGENTS.md "Implementation language: no Python" (in force since 2026-10-03). Its outputs stay as frozen evidence with Python provenance; every new number here comes from native Rust (`examples/isogeny_loop_degree_bound.rs`, `examples/isogeny_loop_glv_gls_lattice.rs`). Porting the sweeper needs a native subexponential class-group algorithm the crate lacks.

*Erratum (derived).* Round 1 charged both directions of a characteristic-2 two-step 2 M. Only the Frobenius $F$ is that cheap; the conjugate step is the Verschiebung $V = [2] compose F^(-1)$. Original statement (§4.3): the folklore split $lambda = t - 1$ via $V^n$ has ratio 0.49, "not established". Correction: $V^n = [2^n] compose pi^(-1)$ is $n$ doublings plus square roots, spent to save $n slash 2$, so the ratio is $>= 2(c_"dbl" + 2) slash c_"dbl" = 2.5$: a provable loss. Every other round-1 conclusion stands a fortiori, because raising a weight only raises loop costs.

= A class-group-free lower bound (derived; native table)

A loop of prime-degree steps is an endomorphism $alpha$ with $N(alpha) = product l_i$ and cost $sum c(l_i) >= kappa log_2 N(alpha)$, $kappa = min_l c(l) slash log_2 l$. A non-scalar $alpha = (A + B sqrt(D)) slash 2$ has $B != 0$, hence $N(alpha) >= d_min (D) = floor((|D| + 1) slash 4)$, and $d_min$ only grows in sub-orders. So every loop endomorphism costs at least $kappa log_2 d_min (D_K)$ — no class group needed. $kappa$ = 14.196 (paper's projective $7.5 l$), 5.047 (affine $2l+2$), 3.786 (optimistic x-only, 3-isogeny eval $4M + 2S$).

#figure(image("figures/degree_bound.svg", width: 78%), caption: [Bound (three per-step models) and the cheapest loop round 1 found, as a fraction of the doublings a 2-dimensional split saves ($8 dot ceil(log_2 r slash 2)$ M). Run `run-20261011T2258Z-degree-bound`; every $t^2 - 4q = f^2 D_K$ re-verified natively. Binary "found" values are round 1's under-priced costs, i.e. lower bounds on the corrected costs.])

For a random curve $|D_K| approx q approx r$, so the ratio is $approx 2 kappa slash c_"dbl"$: 3.27–3.55 (projective) and 1.16–1.26 (affine) on every random prime and binary challenge curve, *including CryptoPro-A*, whose class group round 1 could not compute. A profitable loop on a random curve would need a step formula below $c_"dbl" slash 2 = 4$ M per bit of degree; the optimistic x-only model (0.87–0.95) is not excluded by this bound, only by round 1's searches (cheapest loops 20–156× the saving, projective). Scope: separable steps; binary loops using $F$ are covered by round 1's exhaustive search only.

= Rho on Koblitz curves (derived; heuristic labelled)

An extra rho equivalence needs $alpha in O_K$ whose eigenvalue has finite order $k$ outside $chevron.l -1, lambda_tau chevron.r$. As $O_K^* = {plus.minus 1}$, $r <= N(alpha^k - 1) <= (N(alpha)^(k slash 2) + 1)^2$, so $N(alpha) gt.tilde r^(1 slash k)$; canonicalising an orbit costs $k - 1$ evaluations, about $kappa log_2 r approx 1800$ M per rho step (projective) to gain $sqrt(k)$, needing $k gt.tilde 5 dot 10^4$ against an 8 M addition. Heuristically the minimum norm over the $k$ cosets is $approx r slash k$, which closes the lead for every $k$. Open only for an unknown orbit invariant (as the normal-basis rotation is for $tau$).

= 4-dimensional GLV-GLS with a moderate-$D$ loop (measured lattice; derived model)

On the quadratic twist $E' slash bb(F)_(p^2)$ of a CM curve $E slash bb(F)_p$ ($4p = t^2 + |D| v^2$), $psi$ has eigenvalue $mu = (p-1) slash t$ ($mu^2 = -1$), and the loop endomorphism $alpha' = phi alpha phi^(-1)$ has $lambda = (A + B s) slash 2$, $s = (2 mu - t) slash v$, both derived from the ring map $pi |-> mu$ and asserted per instance. Native LLL (`lll_reduce_hp`) and Babai rounding over 10,000 seeded scalars per lattice, $p$ of 127 bits.

#figure(image("figures/construction.svg", width: 92%), caption: [What is computed natively, what is derived, and what is not validated (no curve is constructed; $lambda$ is not checked on points).])

#table(columns: 8, align: (right, left, right, right, right, right, right, right),
  [*D*], [*cheapest loop*], [*cost M*], [*h*], [*log#sub[2]r*], [*GLS-2 max/mean*], [*GLV-GLS-4 max/mean*], [*model ratio*],
  [−3], [aut. $omega$ (control)], [1], [1], [254], [126 / 125.5], [63 / 62.6], [1.357],
  [−4], [aut. $i$], [1], [—], [—], [no instance (twist order splits)], [—], [—],
  [−7], [2], [15], [4], [251], [125 / 124.1], [63 / 61.9], [1.340],
  [−8], [2], [15], [8], [250], [125 / 123.5], [63 / 61.7], [1.339],
  [−11], [3], [22.5], [20], [249], [124 / 123.0], [62 / 61.4], [1.332],
  [−23], [2·3], [37.5], [52], [248], [124 / 122.8], [62–63 / 61.2–61.3], [1.316],
  [−71], [2·3²], [60], [20], [249], [124 / 123.0], [62–63 / 61.4], [1.292],
  [−431], [2²·3³], [97.5], [20], [249], [124 / 122.9], [63 / 61.7–62.1], [1.248 / 1.241],
  [−619], [5²·7], [127.5], [25], [250], [125 / 123.5], [62–63 / 61.6–62.1], [1.220 / 1.228],
  [−3019], [5³·13], [210], [5], [251], [125 / 124.2], [63 / 62.0], [1.158],
)

#figure(image("figures/glvgls_bits.svg", width: 88%), caption: [Mean over 10,000 scalars of $max_i "bits"(|k_i|)$; both ring embeddings ($s$, $-s$) shown. Run `run-20261011T2258Z-glvgls-lattice`.])

#figure(image("figures/glvgls_model.svg", width: 88%), caption: [Derived operation-count model (not measured): interleaved $w = 4$ NAF, $c_"dbl" = 8$, $c_"add" = 12$, $c_psi = 4$ $bb(F)_(p^2)$ multiplications, loop charged once.])

*Measured:* cheapest loops are certified exhaustively up to $N <= 2^(C slash kappa)$ and reproduce round 1 and the paper ($D = -619$: 127.5 M vs 128; $D = -3019$: $5^3 dot 13$, 210 M, Table 10). The 4-dimensional lattice reaches 62–63 bits ($log_2 r slash 4 = 62$–$63.5$) for every $D$, as for the control. *Derived:* $D = -4$ has no quadratic-GLS instance, since with $x = t slash 2$, $(p - 1)^2 + t^2 = (x^2 + (1+v)^2)(x^2 + (1-v)^2)$ (a first run stalled 13.5 CPU-minutes there and was stopped; kept as a failed attempt). *Not established:* $lambda$ on points of a constructed $E'$, and any timing.

= Crate issue, and open leads

`lattice::lll_reduce` (f64 Gram–Schmidt) returned the 254-bit GLS-2 lattice unreduced without error (189-bit vectors vs 126 with `lll_reduce_hp`); reported, not fixed here.

#table(columns: 3, align: (left, left, left),
  [*lead*], [*status*], [*what settles it*],
  [native sweeper port (130-bit class groups)], [required by AGENTS.md before more sweeps; not attempted], [native subexponential class-group algorithm],
  [point-level validation of $lambda$ on a GLS curve], [open], [native Hilbert class polynomials for $h > 1$, twisted loop over $bb(F)_(p^2)$],
  [GLV-GLS-4 ($D = -619$) implementation timing], [open], [Rust benchmark as in `cryptopro_b_glv_chain_20261005`],
  [binary GLS analogue], [proposed], [same lattice run, binary cost constants],
  [random curves under x-only costs], [bound inconclusive (≈ 0.93×)], [native x-only re-run of the searches],
  [extra rho equivalences on Koblitz curves], [closed (evaluate-the-orbit model)], [a new orbit invariant],
  [folklore split $lambda = t - 1$ on binary curves], [closed: provable loss], [—],
)

#text(8pt, fill: rgb("#52514e"))[Sources: Koshelev, Sanso, ePrint 2024/1985 (CIFRIS25 ACTA, doi:10.69091/koine/vol-7-W09); Galbraith, Lin, Scott, ePrint 2008/194; Longa, Sica, ASIACRYPT 2012. Records: `degree_bound.{md,json}`, `glvgls_lattice.{md,json}`, `revision.txt` (base commit db0393b4, rustc 1.97.0, source SHA-256).]
