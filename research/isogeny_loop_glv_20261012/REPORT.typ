// Rendered with Typst (Python wheel, scratch venv); figures from figures/plot_figures.py (recorded values only).
#set document(title: "Isogeny-loop GLV, round 3", date: datetime(year: 2026, month: 10, day: 12))
#set page(paper: "a4", margin: (x: 1.8cm, y: 1.7cm), numbering: "1 / 1",
  footer: context [#set text(7.5pt, fill: rgb("#52514e")); research/isogeny_loop_glv_20261012/REPORT.pdf #h(1fr) #counter(page).display("1 / 1", both: true)])
#set text(font: "DejaVu Sans", size: 9pt, fill: rgb("#0b0b0b"))
#set par(justify: true, leading: 0.55em)
#set heading(numbering: "1.")
#show heading.where(level: 1): set text(11pt)
#set table(stroke: 0.4pt + rgb("#c3c2b7"), inset: 4pt)
#show table: set text(7.6pt)
#show table: set par(justify: false)

#align(left)[
  #text(15pt, weight: "bold")[Isogeny-loop GLV, round 3: constructed curves, point-level validation, timings, binary GLS]
  #v(2pt)
  #text(8.5pt, fill: rgb("#52514e"))[2026-10-12 · study directory `research/isogeny_loop_glv_20261012/` · follows round 2 (`research/isogeny_loop_glv_20261011/`, PR \#1698) and round 1 (`RESEARCH_ISOGENY_LOOP_GLV.md`, PR \#1644) · status words: *measured*, *derived*, *proposed* · all computation native Rust (`isogeny_algos/`)]
]

= Summary

Round 2 left the 4-dimensional GLV-GLS split as a lead backed only by lattices and an operation-count model. On a quadratic GLS twist $E' slash bb(F)_(p^2)$, the split combines $psi$ ($psi^2 = -1$) with a horizontal-isogeny loop endomorphism of a moderate-$|D|$ CM curve. This round builds the curves, checks the eigenvalues on points, times the split, and carries it to characteristic 2.

#table(columns: (auto, 1fr),
  [*point-level validation*], [*measured, all pass.* Curves at $D in {-619, -23, -71, -431, -3019}$, with $p$ of 127 bits, are built from a natively computed $H_D$. $psi$ acts as $+mu$, and each closing loop acts as $plus.minus lambda(plus.minus s)$, on 24 points of order $r$ (248–253 bits). The two closing walks land on conjugate embeddings.],
  [*benchmark, $D = -619$*], [*measured.* GLV-GLS-4 with a projective loop beats GLS-2 by 1.147× (run 1) and 1.091× (run 2, a different host), median/median. In run 1 it runs and 0.498× the time of wNAF; the A/A control reads 0.987. Across five $D$, the ratio is 1.06–1.21, which is 0.05–0.17 below round 2's model.],
  [*binary GLS*], [*constructed; validated at toy size.* Frobenius-assisted loops $frak(l) dot F^k$ exist on prime-order GLS twists at $m = 127$ and $131$, but not at $m = 83$ for $l <= 61$. An x-only point check at $m = 31$ over GF($2^62$) passes on 2 discriminants. *Not timed.*],
  [*native sweeper, x-only*], [§5.],
)

#figure(image("figures/status.svg", width: 92%), caption: [What round 3 computed natively (blue), validated or measured (green), and built in characteristic 2 (amber).])

= Point-level validation (measured)

`H_D` comes from `src/cm.rs`. It evaluates $j = E_4^3 slash Delta$ on reduced forms in fixed point, with $1 slash q = e^(pi sqrt(|D|) slash a)$ and the $eta$-product computed by the pentagonal series. Precision grows until the worst coefficient rounding error is below $2^(-16)$, and that bound is asserted. At this size the recorded bound is $2^(-86)$ to $2^(-118)$.

A curve over $bb(F)_p$ comes from a root of $H_D$, twisted if its trace is $-t$. A depth-first crater walk (`api::isogenies`, crater $j$ only, no backtracking) follows the certified degree sequence back to $j_0$, and a closing isomorphism finishes the loop. Checks:
- on $E(bb(F)_p)[r_E]$, each loop acts as one of $plus.minus lambda(plus.minus s)$;
- on $E'(bb(F)_(p^2))[r]$, $psi = [mu]$ with $mu = (p - 1) slash t$;
- on $E'(bb(F)_(p^2))[r]$, $phi alpha phi^(-1) = [plus.minus lambda(plus.minus s)]$, with $lambda(s) = (A + B s) slash 2$ and $s = (2 mu - t) slash v$.

#table(columns: 7, align: (right, left, left, right, left, left, left),
  [*D*], [*steps*], [*$H_D$ deg / bits / err*], [*log#sub[2]r*], [*$psi$*], [*twisted loops*], [*verdict*],
  [−619], [5·5·7], [5 / 303 / $2^(-90)$], [251], [$+mu$ 24/24], [$-lambda(s)$; $-lambda(-s)$], [PASS],
  [−23], [2·3], [3 / 150 / $2^(-86)$], [251], [$+mu$ 24/24], [$-lambda(-s)$; $-lambda(s)$], [PASS],
  [−71], [2·3·3], [7 / 238 / $2^(-92)$], [248], [$+mu$ 24/24], [$-lambda(-s)$; $-lambda(s)$], [PASS],
  [−431], [2·2·3·3·3], [21 / 631 / $2^(-118)$], [252], [$+mu$ 24/24], [$+lambda(-s)$; $+lambda(s)$], [PASS],
  [−3019], [5·5·5·13], [7 / 523 / $2^(-88)$], [253], [$+mu$ 24/24], [$-lambda(-s)$; $+lambda(s)$], [PASS],
)
Run `run-20261011T2330Z-gls-validate`. A replay on the final source reproduced the table byte for byte.

= Timing (measured)

All arms run in one binary on the same arithmetic: `Fp2<FpM<2>>` (Montgomery, Karatsuba), a generic 127-bit $p$, Jacobian doubling and mixed addition, and $w = 5$ tables normalised by one batched inversion. Each GLV-GLS-4 call computes $alpha'(P)$ and both $psi$ images. Every pair is checked before timing: all arms agree, $psi(P) = [mu]P$, $alpha'(P) = [lambda]P$, and each decomposition is $equiv k$. The runs use 256 pairs and 30 alternating rounds, pinned to one core of a cloud Xeon at 2.10 GHz, not isolated.

#table(columns: 6, align: (left, right, right, right, right, right),
  [*arm, D = −619*], [*median ns*], [*min ns*], [*wNAF / arm*], [*GLS-2 / arm*], [*note*],
  [wNAF], [142,773], [124,949], [1.000], [0.572], [],
  [GLS-2], [81,665], [71,901], [1.748], [1.000], [125-bit coefficients],
  [GLV-GLS-4, affine loop], [73,066], [67,012], [1.954], [1.118], [63-bit coefficients],
  [GLV-GLS-4, projective loop], [71,172], [64,788], [2.006], [1.147], [],
  [wNAF (A/A control)], [144,657], [125,440], [0.987], [—], [noise floor],
)
Stage medians: decomposition 2,850 ns; $alpha'(P)$ 7,776 ns affine and 5,602 ns projective.

*Confirmation run 2* (`run-20261012T0216Z-bench-confirm`): same command and seed, but a different host after a container restart (Xeon at 2.80 GHz, load average 1.23 at start). Medians: wNAF 211,895, GLS-2 114,823, GLV-GLS-4 105,273 (projective) and 106,268 (affine); A/A 0.989. That makes GLS-2 / GLV-GLS-4 = 1.091, against 1.147 in run 1. The 0.056 difference exceeds either run's A/A deviation; its cause is not isolated (host, load, and decomposition 2.9 vs 5.3 µs all changed).

#figure(image("figures/bench_arms.svg", width: 100%), caption: [Per-round values (µs per scalar multiplication, mean over 256 scalars), 30 rounds, D = −619: run 1 (left, `run-20261011T2333Z-bench`) and run 2 (right, `run-20261012T0216Z-bench-confirm`, a different host).])

#figure(image("figures/bench_trend.svg", width: 82%), caption: [GLS-2 median / GLV-GLS-4 median at five $D$, against the certified loop's model cost, beside round 2's derived operation-count model. A/A controls 0.987–1.019.])

*Observations.*
- The split was faster than GLS-2 at every $D$ tested: 1.214, 1.120, 1.115, 1.147 (1.091 in run 2) and 1.062 with the projective loop.
- Every measurement is below the model. Two costs the model left out are visible in the stage timings:
  - the decomposition, about 3.5% of a GLS-2 call;
  - $alpha'(P)$, 2–10.5% of a GLS-2 call. Its stage time rises with model cost: 1.7, 2.7, 4.2, 5.6 and 8.5 µs.
- The end-to-end ordering between $D$ is not resolved below about 3%, because GLS-2 medians alone range 79.1–82.2 µs across runs.
- Scope: variable time, a generic $p$, one host, one build, no optimised-library comparison.

= Binary-field GLS with Frobenius-assisted loops

*Construction (derived).*
- Take $E slash bb(F)_(2^m)$ ordinary, with $pi = (t + v sqrt(D)) slash 2$ and $t$ odd. Then $D equiv 1 space (mod 8)$ and $pi$ generates $frak(p)_2^m$.
- An odd step $frak(l)$ followed by $k$ squaring steps $F$ closes when $[frak(l)][frak(p)_2]^k = 1$. Its degree is $l dot 2^(k k)$ with $k k = min(k, m - k)$, and its cost is $c(l)$ plus squarings.
- Only $F$ is cheap (round 2's erratum: $V$ costs a doubling), so $alpha$ must be F-oriented. A 2-adic root $s_2$ of $D$ with $v_2(t + v s_2) > m$ fixes that direction, and $alpha$ is F-oriented iff $v_2(A + B s_2) > k k$.
- The quadratic twist over $bb(F)_(2^(2m))$ carries $psi$ as in the prime case, with $\#E' = (2^m - 1)^2 + t^2 = 2r$.

*Scans (measured).*
#table(columns: 7, align: (right, left, right, right, right, right, left),
  [*m*], [*|D| range*], [*scanned*], [*ord = m*], [*prime r*], [*loop, l ≤ 61*], [*cheapest (D: l, F steps, model M)*],
  [83], [≤ 8·10⁶], [810,576], [440], [4], [0], [—],
  [127], [≤ 4·10⁶], [405,292], [319], [3], [2], [−1532327: 3, 53, 128.5 · −2854663: 43, 21, 364.5],
  [131], [≤ 8·10⁶], [810,576], [459], [4], [3], [−283583: 13, 16, 129.5 · −359383: 7, 60, 172.5 · −1148567: 29, 59, 335.5],
)
GLV-GLS-4 lattices give 64-bit coefficients ($m = 127$) and 66-bit ($m = 131$), against 126 and 130 for GLS-2 (5,000 scalars each). $m = 127$ is the field of GLS254 (Oliveira, López, Aranha, Rodríguez-Henríquez). $m = 131$ is ECC2K-130's base field, but these are moderate-$D$ CM curves, not the $D = -7$ challenge curve.

#figure(image("figures/binary_scan.svg", width: 92%), caption: [Funnel counts per $m$ (log scale) and the cheapest loop's model cost per kept discriminant. Run `run-20261012T0000Z-binary-scan`.])

*Point-level validation at $m = 31$ (measured).* At $m = 127$, $h >= 127$, and the hit has $h = 1397$, so $H_D$ is out of reach for this code. The check therefore runs at toy size: $H_D$ mod 2, roots in GF($2^62$), Kohel kernel polynomials for odd $l$, then $F^(k k)$ back to $j_0$, and the twist $a_2 + delta$ with $"Tr"(delta) = 1$. Twist and isomorphisms fix $x$ in characteristic 2, so the checks are x-only:
- $x(psi P) = x^q$ is compared with $x([mu] P)$;
- $R(x)^(2^(k k))$ is compared with $x([lambda(plus.minus s)] P)$.

Results:
- $D = -39583$ ($h = 62$; $l = 31$, 11 F steps; $H_D$ at 4,162 bits): PASS on 24 points.
- $D = -482623$ ($h = 186$; $l = 31$, 12 F steps; $H_D$ at 14,081 bits): PASS on 24 points.

Run `run-20261012T0002Z-binary-validate`.

*Cost model (derived, untimed).*
- Over $bb(F)_(2^254)$, 4-dimensional decomposition cuts the shared doublings from about 126 to about 64, the measured coefficient sizes. At about 4 M per doubling (López–Dahab or λ-coordinates, squarings free), that saves about 250 M.
- The $D = -1532327$ loop is one 3-isogeny (22.5 M) plus 53 F steps. The F steps compose to $x |-> x^(2^53)$, which a multi-squaring table applies in about 1–2 M.
- So the binary lead looks stronger than the prime one. In the prime case, though, the model overstated the measured gain by 0.05–0.17.
- Not examined: security beyond prime order, including GHS for $bb(F)_(2^(2m))$ with $m$ prime and the moderate class number.

= Native sweeper port and x-only re-run

*Pending.* The native sweeper port and x-only re-run are in progress; see `README.md` §5.

= Open after round 3

- A constant-time GLV-GLS-4 over a special prime. The condition is on $D$ ($4p = t^2 + |D| v^2$), not on $p$, so Mersenne $p$ is allowed. Does the shortfall against the model persist once the field is fast?
- A binary GLV-GLS-4 benchmark over $bb(F)_(2^254)$ with carry-less multiplication.
- Folding the x-only loop costs into the GLV-GLS model.

Evidence: `README.md` §7 indexes every record and run report.
