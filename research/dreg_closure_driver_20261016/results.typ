#set page(paper: "a4", margin: 2cm)
#set text(size: 10pt)
#set par(justify: true)

= What sets the MutantXL closure degree: not N; ℓ = 7 or m·ℓ ≥ 21 remain

These are the results of `research/dreg_closure_driver_20261016/`. The full
text is in `RESULTS.md`. *Class:* stage diagnostic. The `n = 131` statements
are *derived*.

#figure(image("reports/run-20261016T2049Z/closure_driver.svg", width: 100%),
  caption: [56 unsatisfiable draws over two studies, placed by (ℓ, N). This
  study's four cells (16 draws) are outlined in black.])

#table(columns: 3,
  [*Prediction*], [*outcome*], [*confidence*],
  [Q7: m = 3, ℓ = 6 refutes at 4 at N = 40, 42], [*holds*: 8 of 8], [about 70%],
  [Q8: m = 4, m·ℓ = 20 rises past 4], [*refuted*: 8 of 8 refute at 4], [about 50%],
)

- *N:* rejected.
- *m·ℓ with a threshold ≤ 20:* rejected.
- *Surviving explanations:* the degree rises at ℓ = 7, or once
  m·ℓ ≥ 21. Both fit all 56 draws.
- *ECC2K-130 (derived):* both surviving explanations put the m = 4 rise at
  ℓ = 6 or 7, far below the priced ℓ = 33. The m = 4 "4, constant" row
  (2^68.43) is therefore contradicted, and the degree-5 rows (m = 4 at
  2^71.13) become the optimistic case.
- *Next:* m = 4, ℓ = 6 separates the two survivors, and degree 5 past
  ℓ = 7 decides the rows. Both need more memory or a sparse closure.
