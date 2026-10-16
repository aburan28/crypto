#set page(paper: "a4", margin: 2cm)
#set text(size: 10pt)
#set par(justify: true)

= MutantXL closure at scale: degree 4 holds to N = 42 at m = 4, and rises past 4 at m = 3, ℓ = 7

These are the results of `research/dreg_closure_scale_20261016/`. The full
text is `RESULTS.md`, and the run report is `reports/run-20261016T1705Z/`.
*Class:* stage diagnostic. The `n = 131` rows are *derived extrapolations*.

#figure(image("reports/run-20261016T1705Z/closure_scale.svg", width: 100%),
  caption: [40 unsatisfiable draws, four per cell, run 2026-10-16T17:05Z.])

#table(columns: 4,
  [*Prediction*], [*outcome*], [*confidence*], [*evidence*],
  [Q4: the closure refutes at ≤ 4 in every draw], [*refuted*: 0 of 8 at m = 3, ℓ = 7; 32 of 32 elsewhere], [about 55%], [`score-output.txt`],
  [Q5: XL never refutes at 4], [*holds*: 0 of 40], [high], [same],
  [Q6: log₂ ops ÷ log₂ cols in 2.0–2.4], [*holds*: 2.18–2.33], [—], [same],
)

- *Registered verdict:* "rises at N = 40".
- *Descriptive reading (not registered):*
  - the rise tracks ℓ rather than N, because m = 4 at N = 40 and 42 stays
    at 4;
  - the refuting cells have m·ℓ ≤ 18 and the rising cells have 21.
- *Soundness control:* the one registered satisfiable draw was not
  refuted. That is thin, and the report says so.

== ECC2K-130 (derived)

#table(columns: 6,
  [*m*], [*degree*], [*w*], [*ℓ*], [*log₂ total*], [*vs rho*],
  [3], [4, constant — refuted at ℓ = 7], [1], [40], [82.27], [+21.46],
  [3], [5, constant], [1], [41], [85.18], [+24.37],
  [4], [4, constant — holds to ℓ = 4], [1], [33], [68.43], [+7.62],
  [4], [5, constant], [1], [34], [71.13], [+10.32],
)

- No row beats rho, and the units are unconverted.
- *The next cell* is m = 4, ℓ = 5 at degree 4 (about 179k columns). It is
  at this container's memory limit.
- *Degree 5 at N ≥ 40* needs a sparse closure.
