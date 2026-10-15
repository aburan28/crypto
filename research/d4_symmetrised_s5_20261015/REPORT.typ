#set page(paper: "a4", margin: 2cm)
#set text(size: 10pt)
#set par(justify: true)

= Lead 5b: the $D_4$-symmetrised $S_5$ as four weighted unknowns

_2026-10-15, lane d4-s5. Study `research/d4_symmetrised_s5_20261015`. Protocol committed in `87fe78b63` before any measured run. Runs 06:38Z–07:43Z, one thread each, under a 2.5 GB cap, on a shared 4-vCPU container. Class: engineering, a constant of one phase. No end-to-end $S$._

== The lead was already half-run

The encyclopedia lists the $D_4$-symmetrised $S_5$ as "never run". That entry is stale. Candidate C8 (PR \#1785, `research/gaudry_quartic_d4_20261013`) ran it on 2026-10-13:

- form: five *unweighted* unknowns $(sigma_1, sigma_2, sigma_3, sigma_4, pi)$ with $pi^2 = sigma_4$, Macaulay degree 14;
- matrix: $13412 times 11628$;
- cost: $C_4 = 2.54 dot 10^10$, a factor of 48 below the plain $S_4$ form's $1.213 dot 10^12$. The projection was $approx 630$.

This study builds the *weighted* Macaulay solve that the GLV note §9.4 named and nobody built.

== Construction (derived)

*Curve and base.* The curve is over $FF_(p^4)$ with $T = (x_0, 0)$ and $f'(x_0) = c^2$. The coordinate is $Y = (x - x_0 - c)/(x - x_0 + c)$, with $Y(P+T) = -Y(P)$. The base is ${P : Y(P) in FF_p}$.

*Group and invariants.* $D_4 = (ZZ\/2)^3 ⋊ S_4$, of order 192, acts by permutations and by even sign changes, i.e. by even $T$-translates. The invariants are $sigma_k = e_k (Y^2)$ and $pi = product Y_i$.

*System.* Four Weil components in $(sigma_1, sigma_2, sigma_3, pi)$, each of weighted degree 8 under the weights $(2,2,2,1)$.

*Bézout.* Two derivations agree:
- weighted: $8^4 \/ (2 dot 2 dot 2 dot 1) = 512$;
- orbits: $4096 \/ (192\/24) = 512$.

*Macaulay bound.* $32 - 7 + 1 = 26$, at which the matrix has *4,200 columns*. C8 has 11,628 and the plain form 40,920.

== Controls (measured)

- *Invariance detector.* A wrong $c$ makes the invariant interpolation fail; the true $c$ passes. Checked at $p = 13$ and $p = 269$.
- *Generic planted weighted systems.* 9 of 9 have quotient 512 and recover the planted root, at $1.39$–$1.45 dot 10^9$.
- *Planted targets.* The planted quadruple is recovered by arm w in 41 of 42 and by C8 in 23 of 24.

== Results (measured)

#figure(image("figures/c_per_solve.svg", width: 100%), caption: [$F_p$ multiplications per solve, every target. C8 (red) and the weighted form (blue) run on the same targets. The plain level is quoted from §11.17; it was not run, because its 3.42 GB matrix is over the cap.])

#table(columns: 6, align: (right, left, left, left, left, left),
  [$p$], [$C_w$ median], [$C_"c8"$ median], [paired $C_"c8"\/C_w$], [$1.213 dot 10^12 \/ C_w$], [ground truth],
  [13–37], [$1.21$–$1.34 dot 10^9$], [$2.05$–$2.42 dot 10^10$], [17.9–18.5], [905–1005], [brute 23/24, MITM 23/24],
  [269], [$1.411 dot 10^9$], [$2.542 dot 10^10$], [18.0 (16.3–18.7)], [860 (772–893)], [MITM 20/20],
  [1033], [$1.406 dot 10^9$], [$2.554 dot 10^10$], [18.2 (17.4–18.7)], [863 (823–888)], [MITM 20/20],
  [4093], [$1.437 dot 10^9$], [$2.561 dot 10^10$], [17.8 (16.5–18.6)], [844 (778–885)], [MITM 20/20],
)

- *Lift-back.* The split, $F_p$ square roots, lifts and group checks of all 8 translates cost at most $2.8 dot 10^5$, which is at most $2 dot 10^(-4)$ of $C_w$.
- *Degree-27 reruns.* 5 of 101 solved targets had a degenerate eigenspace and reran at degree 27, charged at $2.5$–$3.3 dot 10^9$ each.
- *One target unsolved.* At $p = 17$, planted target 4 has quotient 514 and is unsolved by both arms; brute force finds 4 decompositions there.
- *Exactness.* On every other target the two arms and the ground truth agree exactly.

#figure(image("figures/phases.svg", width: 95%), caption: [Phases of one solve (target 1, seed 1). Elimination falls by a factor of 40 and normal forms by 8. The characteristic polynomial is unchanged, since $D = 512$ in both arms.])

== Verdict under the registered rule

- *W1 and W2.*
  - W1 held at $p >= 29$. At $p = 17$ the dimensions were 510–514.
  - W2 is *not met as registered*: one of 102 targets was unsolved. It was a flagged failure shared with C8, not a wrong set. The decision below therefore carries that caveat.
- *W3* held: the medians are flat in $p$ to 2.2 %.
- *Lead 5b: "projection met".*
  - The median ratio to the plain $C_4$ is 844–863 at all three sizes, inside $[315, 1260]$.
  - It is not "accounting", because the lift-back is negligible.
- *Against C8: "advances on C8".* The paired ratio is at least 7.7 on all 47 paired targets.

The gap between 48 and 630 was the presentation, not the symmetry.

#figure(image("figures/handover.svg", width: 100%), caption: [The $k = 4$ handover in E18's $tau_T$ frame. Inputs: $h = 4$, 106 residuals per relation, $r = 0.341$ against the negation walk. The formula is $n^* = 157.9 - 4 log_2 (1.213 dot 10^12 \/ C)$. Derived and extrapolated, not measured.])

== Handover (derived, extrapolation)

On 2-torsion curves the weighted solve moves the $tau_T$ handover from $2^157.9$ to $2^118.9$. C8 alone moves it to $2^135.6$.

Past the handover the method tends to $0.341 times$ the negation walk, and $r_infinity$ is unchanged. The figure rests on the exponents $n^(1\/4)$ and $n^(1\/2)$, and on a $C$ that is flat at $p <= 4093$. No $k = 4$ run below rho exists.

== Limits

- The plain arm was not paired on the same targets; it is over the cap.
- No relation phase or linear algebra was run.
- The degenerate targets' cause is unknown.
- Unit tests: 11 of 11 pass after a test-only fix. The old test planted a zero sum at $p = 13$; the test now uses $p = 29$.
- Recorded $C_w$ values were reproduced exactly from a private-target rebuild after merging main.
- No identity certificate was produced.
