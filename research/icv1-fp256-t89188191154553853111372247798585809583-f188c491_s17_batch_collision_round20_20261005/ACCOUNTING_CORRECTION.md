# Round 20 finite-domain accounting correction

Date frozen: 2026-10-05, after the first width-ladder run and before merge.

The preregistered symmetric birthday projection

```text
N_left = N_right = sqrt(R*n/p_disjoint)
```

is not attainable for 138,031 P-256 rows.  It requires `2^136.538`
distinct left images, while the complete signed-eight domain contains only

```text
L = C(131458,8) * 2^8 = 2^128.734 images.
```

Sampling beyond `L` repeats representations and cannot create new independent
relation rows.  The initial `2^144.625` direct projection and `2^141.625`
one-addition-per-image boundary are therefore withdrawn.  Their measured
width-ladder inputs remain valid and their first frozen JSON remains preserved
in commit history; they may not be used as attack costs.

The corrected projection is asymmetric and finite-domain bounded:

1. materialise every unique signed-eight image once;
2. use the complete left domain `L` as the collision table;
3. stream target-adjusted signed-nine images; and
4. require

   ```text
   R_right = collision_events * n / (L * p_disjoint).
   ```

Use round 19's whole-factor-base success probability
`0.9640001820745058`, so the deterministic target allowance is
`138031 / success`.  Compute the exact per-target signed-17 relation mean as

```text
lambda = C(131458,17) * 2^17 / n.
```

The target pool contains `lambda * target_allowance` expected distinct
relations.  To obtain 138,031 distinct rows under uniform sampling, charge the
occupancy inversion

```text
collision_events = -population * ln(1 - 138031/population).
```

This explicitly charges duplicate relations.  The 138,031-row requirement
already includes the frozen rank allowance over 131,458 columns.  Continue to
charge seven group additions per unique left image, nine per streamed right
image, measured batch-normalisation FME per image, relation replay, target
generation, and sparse linear algebra.  The optimistic generic boundary
charges one group addition per left/right image and gives key extraction for
free.

Direct materialisation stores the complete left table only; the right side is
streamed.  A projected full-depth run remains prohibited unless the corrected
cost, per-row, storage, degree, exactness, and rho gates all pass.

No measured cell, target, seed, width, or observed collision is changed by
this correction.  It is an accounting correction discovered by checking the
projection against the exact finite representation domains.
