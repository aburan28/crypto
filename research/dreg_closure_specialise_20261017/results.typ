#set page(paper: "a4", margin: 2cm)
#set text(size: 10pt)
#set par(justify: true)

= One fixed variable usually brings the ℓ = 7 draws back to degree 4

These are the results of `research/dreg_closure_specialise_20261017/`; the
full text is in `RESULTS.md`. *Class:* stage diagnostic.

#figure(image("reports/run-20261017T0211Z/specialise.svg", width: 100%),
  caption: [8 draws × 4 arms, run 2026-10-17T02:11Z.])

- *Registered verdict: uninformative for H_sum vs H_max.* The null control
  C (one intermediate coordinate fixed) refutes on 6 of 8, as often as the
  discriminating arm A.
- *Descriptive:*
  - Any single fixed variable lets degree 4 refute on 6 of 8 draws. The
    ℓ = 7 failure at degree 4 is marginal.
  - On 5 of 8 draws, both values of one guessed variable refute. A
    one-guess hybrid therefore refutes at degree 4, with
    `C(N−1, ≤4)` ≈ 92k–112k columns per branch against
    `C(N, ≤5)` ≈ 760k–1M.
- *ECC2K-130 (derived):*
  - Guess-then-degree-4 is a new candidate oracle,
    `2^k · C(N−k, ≤4)^w`.
  - How k grows with ℓ toward the priced ℓ = 33 is unmeasured. If it is
    about m·(ℓ−6), it is worse than degree 5.
  - No row moves, and no row beats rho.
- *Next:* the number of guesses needed at ℓ = 8, which sits at this host's
  memory limit; and m = 4, ℓ = 6 on a larger host.
