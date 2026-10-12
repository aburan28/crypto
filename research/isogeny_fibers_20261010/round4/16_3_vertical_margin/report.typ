#set page(paper: "a4", margin: 2cm, numbering: "1")
#set text(size: 9.5pt)
#set par(justify: true)
#show heading: set text(size: 11pt)

= Lead 16.3: the vertical security margin
Round 4, 2026-10-11. Code `examples/r4_16_3_vertical_margin.rs` (native Rust), base commit `440c0e163`,
branch `research/isogeny-fibers-round4`. Raw records: `rational.json`, `generic.json`, `elkies.json`;
fits and derived table: `analysis.json` / `analysis.txt`. Full write-up: `RESULTS.md`.
Toy, synthetic instances only; no sub-rho or deployed-curve claim.

== Question
Above which prime $ell$ in a conductor gap does a vertical transfer (find the kernel + evaluate the
isogeny) cost more than rho on the target, $sqrt(pi r slash 4)$ group operations? Prop. 3.3 predicts
finding $tilde(O)(ell^2)$, evaluating $tilde(O)(sqrt(ell))$ (√élu) or $O(ell)$ (Vélu).

== Method (measured)
Instances: the $j = -3375$ curve $y^2 = x^3 - 35x - 98$ (CM by the maximal order of $QQ(sqrt(-7))$)
over primes $p = (t^2 + 7v^2) slash 4 in [2^60, 2^61)$ with $ell | v$ and $t equiv 2 lambda$ (mod $ell$):
the crater of an $ell$-volcano, Frobenius scalar $lambda$ on $E[ell]$, so the $ell - (-7 slash ell)$
non-horizontal lines are descending (vertical) isogenies. A codomain with $j != -3375$ certifies a vertical edge.
- *Rational kernel* ($lambda = 1$, $E[ell] subset E(FF_p)$, $ell^2 || \#E$): find = cofactor ladder; evaluate by Vélu and by √élu (BDLS 2020), $ell$ = 11 … $2^24$, 3 instances each.
- *Generic kernel* ($lambda$ primitive, kernel $x$-coordinates generate $FF_(p^k)$, $k = (ell-1) slash 2$): find = cofactor ladder on the twist over $FF_(p^k)$ + Berlekamp–Massey minimal polynomial; evaluate from the kernel polynomial. $ell$ = 5 … 433 (420 s cap).
- *Library Elkies route* (`isogeny_walk`: $Phi_ell$ mod $p$, roots, Elkies kernel): contended wall time, $ell <= 79$ (300 s cap).
Unit: M = $FF_p$ multiplications and squarings (inside polynomial and $FF_(p^k)$ arithmetic too; Karatsuba, threshold 24); I = inversions, $"Meq" = M + 91 I$.
Checks (all passed): √élu = Vélu on codomain and images (126 rows); images have order dividing $\#E$ on the codomain;
library `verify_kernel` + equal $j$ (24 rational rows, 44 generic rows); Elkies succeeded on every non-crater root of $Phi_ell (-3375, Y)$.

#figure(image("fig1_cost_vs_ell.svg", width: 100%),
  caption: [Counted cost per phase against $ell$ (measured, 61-bit $p$). Dashed: OLS fits. Dotted: rho on $r = 2^64$ at 6 M per step (derived).])

== Fitted exponents (measured; 95 % t-intervals measure power-law misfit, counts are deterministic)
#table(columns: 4, stroke: 0.4pt, inset: 4pt,
  [*phase*], [*range*], [*exponent [95 % CI]*], [*log₂ const*],
  [rational find], [11 … $2^24$], [−0.025 [−0.036, −0.015] (flat, ≤ 1,663 Meq)], [10.01],
  [Vélu, codomain + 1 image], [$>= 2^10$], [0.998 [0.997, 0.999]], [2.62],
  [√élu, codomain + 1 image], [$>= 2^10$], [0.856 [0.845, 0.868]], [4.22],
  [√élu, codomain + 1 image], [$>= 2^16$], [0.822 [0.805, 0.838]], [4.93],
  [generic find, expected], [31 … 433], [2.605 [2.578, 2.631]], [9.38],
  [generic find, expected], [$>= 150$], [2.552 [2.477, 2.628]], [9.81],
  [generic $FF_(p^k)$ setup], [31 … 433], [2.831 [2.547, 3.115]], [4.94],
  [generic eval / point (M)], [31 … 433], [0.960 [0.956, 0.963]], [0.35],
  [library $Phi_ell$, seconds (contended)], [11 … 79], [4.688 [4.633, 4.743]], [−23.25],
)
At $ell approx 2^24$: Vélu $1.01 dot 10^8$ Meq, √élu $2.45 dot 10^7$ (4.1× fewer); √élu passes Vélu near $ell approx 1000$.
√élu's 0.82–0.86 instead of 0.5 is Karatsuba inside the trees ($sqrt(ell)^(1.585)$). Generic finding 2.6 = $k log p$ ladder steps × Karatsuba $FF_(p^k)$ products.

#figure(image("fig2_crossover.svg", width: 100%),
  caption: [Derived crossover $ell^*$ — EXTRAPOLATED from the 61-bit counts (6 M per rho step, $p approx r$, ladder phases scaled by $log_2 r slash 61$). Bars: $ell_max = 2 sqrt(p)$, the largest prime a conductor gap can contain.])

== Derived crossover table (extrapolation, g = 6 M per rho step)
#table(columns: 4, stroke: 0.4pt, inset: 4pt,
  [*model*], [*$r = 2^64$*], [*$r = 2^128$*], [*$r = 2^256$*],
  [rho (F_p mults)], [$2^34.4$], [$2^66.4$], [$2^130.4$],
  [$ell_max = 2 sqrt(p)$], [$2^33$], [$2^65$], [$2^129$],
  [rational, Vélu], [$2^31.9$], [$2^63.9$], [$2^128.0$],
  [rational, √élu (fit)], [$2^35.9 > ell_max$], [$2^74.8 > ell_max$], [$2^152.7 > ell_max$],
  [rational, √élu ($sqrt(ell) log^2 ell$ model)], [$2^40.7 > ell_max$], [$2^99.5 > ell_max$], [$2^222.9 > ell_max$],
  [*generic, fit*], [*$2^9.6$*], [*$2^21.5$*], [*$2^45.7$*],
  [generic, $ell^2 log ell$ model], [$2^9.8$], [$2^24.7$], [$2^55.6$],
)
Kernel over a small extension $FF_(p^k)$ ($r = 2^128$, derived): $ell^* = 2^99.5, 2^95.1, 2^91.5, 2^87.9, 2^84.3, 2^81.4, 2^78.5$ for $k = 1, 2, 4, …, 64$ — all above $ell_max = 2^65$.

#figure(image("fig3_wall_generic.svg", width: 100%),
  caption: [Generic kernel finding, contended wall time: library $Phi_ell$ route ($ell <= 79$) and the $FF_(p^k)$ torsion route of this example ($ell <= 433$).])

== Verdict
*Confirmed, sharpened.* With a generic kernel ($k approx ell slash 2$) the transfer costs more than rho from
$ell approx 2^10$ ($r = 2^64$), $2^22$–$2^25$ ($r = 2^128$), $2^46$–$2^56$ ($r = 2^256$); the cost is driven by the
kernel's field degree $k$, not by $ell$ itself. With a kernel point over $FF_p$ or a small extension there is *no*
margin: with √élu no $ell <= 2 sqrt(p)$ reaches rho, and Vélu reaches it only at $ell approx sqrt(r)$.
Prop. 3.3's "$ell gt.tilde 2^60$ for a 128-bit group" is loose for the generic case and wrong for the rational one.
For a cofactor-1 target the rational case needs $ell | \#E$ or $\#E^"twist"$ over a small extension, which is unusual.

== Limits
61-bit toy field; larger $r$ is extrapolation (generic: measured $ell <= 433$). Karatsuba, no FFT: the
$tilde(O)$ rows are calibrated models. Not implemented: $psi_ell$ factoring, Sutherland's $Phi_ell (j, Y)$ evaluation.
Only descending edges from the crater; one discriminant ($-7$); instances restricted to $v_ell = 2$.
Additions and memory not counted. Conductor `check` was not available in this environment.
