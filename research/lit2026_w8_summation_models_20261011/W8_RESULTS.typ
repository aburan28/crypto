#set page(paper: "a4", margin: 1.5cm)
#set text(size: 9pt)
#set heading(numbering: "1.")

#align(center)[#text(size: 14pt, weight: "bold")[lit2026 W8: summation-polynomial model checks]
#v(2pt) 2026-10-11 · branch `lit2026/w8-summation-models` (stacked on aburan28/crypto\#1681) · program baseline `4b6e09b78`]

*Question.* Do four 2026 papers on curve models and system modelling (Dryło ePrint 2026/1299; Wu arXiv
2608.01675; Bariant ePrint 2026/1927; La Scala–Marchesin–Tiwari arXiv 2609.18866) change the cost of the
repository's Weil-descended decomposition systems over GF(2^n)? All code is Rust
(`examples/lit2026_w8_summation_models.rs`) on the unchanged library; systems are solved by the repo's
Boolean matrix-F4 with splitting, all roots enumerated; arms interleaved, 8 targets x 5 rounds per cell;
deterministic gate = F4 word XORs; wall time L0 on a host at load 80–360. Run reports:
`reports/run-20261011T2017Z/`, `reports/run-20261011T2131Z/`. Full tables: `README.md`.

= E/G check (Dryło): symmetrisation = E/⟨T⟩, plus a cover at m = 3

On a Koblitz curve the FGHR invariant is $w = u^2 + u = 1 slash X compose Phi$ with $X = x + sqrt(b) slash x$ the Vélu
abscissa of $Phi: E -> E slash "⟨T⟩"$, and $b'' = 1$ (checked on 33 relations per cell). The repo's
symmetrised polynomial is $G(w) + s H(w)$ with $s = sum u$, $s^2 + s = sum w$, and

$ "Res"_s (G + s H, s^2 + s + sum w) = G^2 + G H + H^2 sum w = S^(E slash "⟨T⟩")_(m+1)(1 slash w) dot.op product w^(2^(m-1)) quad "over" bb(F)_2, $

certified as *IDC1hc7b594938b50f160* ($m = 2$) and *IDC1hcf83176ae9a93715* ($m = 3$). At $m = 2$, $H = 1$:
the repo's system is the $E slash "⟨T⟩"$ system plus a linearly eliminable parity bit. At $m = 3$,
$s$ multiplies $H$: the repo's system is the $E slash "⟨T⟩"$ $S_4$ (Boolean degree 6) plus an
auxiliary cover that lowers the degree to 4. Root sets (as $w$-tuples) were identical across `sym`, `eg` and
`u-plain` on all 72 targets.

#figure(image("reports/run-20261011T2017Z/eg_arms.svg", width: 100%),
  caption: [Median word XORs per solve and Boolean unknowns, per arm and cell (run-20261011T2017Z).])

At K1/F2^15, m = 3 (median ms; word ops): x-direct 166.7; 9.01M → u-plain 64.5; 5.77M → eg 20.7; 0.54M →
sym 13.1; 1.02M. Interpretation: most of the repo's "symmetry" gain is the change of curve to
$E slash "⟨T⟩"$ that the invariant coordinate is; the non-relabelling part is the $s$-cover at
$m = 3$ (1.6x wall time, 1.9x more word ops, fewer splits).

= The faster variant: solve E/⟨T⟩ directly at m = 2

#figure(image("reports/run-20261011T2131Z/eg_speed_session.svg", width: 100%),
  caption: [Interleaved original (4b6e09b78 library, `sym`) / previous (`sym`) / candidate (`eg`), 5 rounds ABAB, 8 targets per cell (run-20261011T2131Z).])

#table(columns: 5, inset: 3pt,
  [*cell (ICV1)*], [*word ops orig = prev*], [*cand*], [*ratio*], [*wall prev/cand, median [min–max]*],
  [K1/F2^15 `icv1-f2m15-t275-b7f03703`], [6 750 384], [4 218 414], [1.600], [1.60 [0.59–8.62]],
  [K0/F2^15 `icv1-f2m15-tm275-2d22ff5d`], [6 235 447], [3 797 654], [1.642], [1.43 [0.80–4.45]],
  [K1/F2^17 `icv1-f2m17-tm101-00378d4e`], [169 687 396], [117 819 615], [1.440], [2.27 [1.46–2.86]],
  [K0/F2^23 `icv1-f2m23-t5197-69e76b73`], [11 054 138 493], [8 205 050 470], [1.347], [1.25 [1.21–1.26]],
  [K1/F2^23 `icv1-f2m23-tm5197-1f85e9e1`], [10 866 191 751], [8 107 487 123], [1.340], [1.22 [1.12–1.40]],
)
Word-op ratios are deterministic (identical in all rounds); root digests identical across the three arms.
Original/candidate wall ratios are in the run report; original and previous run the same code and their
cross-process spread (up to 5.4x) is the noise floor. Scope: the m = 2 solve on toy cells; lifting an `eg`
root tries both parities.

= Wu's binary state model, Frobenius auxiliary, convolution splitting

#figure(image("reports/run-20261011T2017Z/wu_models.svg", width: 100%),
  caption: [S3 per coordinate: Weierstrass x, Wu Kummer ratio k = x/d, Wu u = d/(x+d) (the FGHR frame), binary Edwards w (an E/⟨T⟩ coordinate).])

Wu's $u$ is a Möbius frame of $x$ and is the FGHR frame ($u -> u + 1$ is translation by $T$; $d = b^(1 slash 4)$;
the model needs $"Tr"(a) = 0$). Its $S_3$ has 11 terms and total degree 6 against 5 and 4; descended
monomials 66–83 against 81–100; word ops 1.2–1.4x fewer at n = 15 and 3.1x more at n = 17. Binary Edwards
$w$ is invariant under $T$ (a Dryło compression) and costs 5.7–7.5x more word ops on these bases. No model is
uniformly cheaper.

#figure(image("reports/run-20261011T2017Z/frob_lucas.svg", width: 100%),
  caption: [Frobenius auxiliary blocks (MatrixF4 and InheritedF4) and e_i splittings.])

*Frobenius auxiliary (negative).* Over the F2 descent $x -> x^(2^k)$ is linear, so the auxiliary is in the
span of the unknowns: same Boolean degree, same ideal after linear elimination. MatrixF4: 57–58x (m = 2) and
701x (m = 3) more word ops; InheritedF4 (eliminates linear generators): 1.03x (m = 2), 8.6x (m = 3).

*Convolution splitting (negative).* Degree 3 → 2 and 4–14% fewer monomials, the same split tree, 1.21–2.07x
more word ops and 2.4–2.9x the wall time at n = 9, 12, 15. Lucas' factorisation applies to symmetric functions
of F2-valued variables and has no counterpart here; m = 4 is measured only at correspondence level (no binary
symmetrised S5 in the repo).

= Status and limits

E/G check (2-torsion, m = 2, 3): complete except the K1/F2^17 m = 3 cell, stopped unfinished after 2 h.
4-torsion: no 4-torsion-symmetrised binary system exists in the repo. Wu: complete. Frobenius auxiliary:
complete (negative). Convolution splitting: partial (S4 full; m = 4 correspondence only). Speedup:
`eg` at m = 2, 1.34–1.64x word ops, interleaved with original and previous. Toy fields n ≤ 23 only; no
IC-versus-rho, exponent or deployed-curve statement. Found: the F4 engine's wide-system path switches at 24
variables, leaving the 22–23-variable `sym`/`eg` systems at n = 23 on a path roughly 4x slower in a
cross-session observation (README).
