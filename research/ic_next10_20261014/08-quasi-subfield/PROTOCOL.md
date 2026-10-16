# Item 8 protocol: quasi-subfield polynomials at prime extension degree

Frozen 2026-10-14, before any data in this directory was collected.

## Check-first record

The item has been partly run before. What exists:

- `src/cryptanalysis/quasi_subfield.rs` and the examples `quasi_subfield_census`
  and `quasi_subfield_reach`. They hold an exhaustive census of linearised
  candidates at `n <= 13` (search capped at 2^34) and a Frobenius-stable
  divisor scan for odd `n <= 600`.
- `research/notes/index-calculus/RESEARCH_QUASI_SUBFIELD.md` reports the
  results of that code:
  - the `n mod n0` term of HKPYY Lemma 4.1 cannot be dropped;
  - at prime `n` only 13 non-subfield stable cells exist, the best being
    `n = 73`, `n0 = 9`, `j = 1`;
  - Euler–Petit Theorem 1 gives `beta >= 3/4` for every linearised
    polynomial that splits completely;
  - the HKPYY generic-beating threshold is `beta < 1/(2 kappa) = 0.1025`
    with `kappa = 4.876`.
- `DECOMPOSITION-SURVEY.md` §3.4 repeats this.

Nothing in the repository covers three things:

1. The general-coefficient linearised class at prime `n` from 17 to 23.
2. The **multiplicative** class `X^{2^{n'}} - X^a`. Euler–Petit Prop. 6 shows
   that Theorem 1 does not cover it: their `X^8 - X^3` over `F_{2^4}` has
   `beta = 0.70 < 3/4`. Nobody here has searched it at prime `n`.
3. A relation-yield and refutation-degree measurement of a QSP factor base
   against random subspaces in the binary decomposition system.

This item adds those three and re-derives the stable frontier with an
independent implementation in Python.

## Hypothesis

H1. Some prime `n <= 61` admits a quasi-subfield polynomial (additive or
multiplicative) with `beta = l*n/n'^2 < 0.1025`, where `l = log2 deg lambda`.

H0, the expected outcome, is that no such polynomial exists. The best `beta`
at prime `n` then stays at or above 0.7, and any QSP subspace base decomposes
like a random subspace of the same dimension.

## Arms

| Arm | Class | Range | Method |
|---|---|---|---|
| A1 | linearised, `F_2` coefficients (Frobenius-stable) | prime `3 <= n <= 61`, all `1 <= n0 < n`; context to `n <= 127` | enumerate every divisor `g` of `t^n - 1` over `F_2` (products of irreducible factors); `l` = second-highest exponent of `g` |
| A2 | linearised, coefficients in `F_{2^n}` | exhaustive over `j in {0,1}` for prime `n <= 23`, every `n0` in `2..n-1`; `j = 2` for prime `n <= 13` | Rust rank test (`quasi_subfield::splits_completely`), leading `lambda` coefficient scaled to 1, valid when `n` does not divide `n0 - j` (always true here) |
| A3 | multiplicative `X^{2^{n'}} - X^a` | prime `3 <= n <= 61`, every divisor `D` of `2^n - 1` | `a = 2^{n'} - D` with `n' = ceil(log2 D)`, plus the next `n'` up; root count `D + 1`; keep rows whose root fraction `(D+1)/2^{n'}` is at least 1/2 |
| D | decomposition, `m = 2`, curve `y^2 + xy = x^3 + 1` | `n = 31`, `k = 15, 16` (QSP stable subspaces); `n = 17`, `k = 9` (best stable, not QSP) | same instrument on QSP base and on 4 random subspaces of equal dimension |

Positive controls (the detector must recover each one):

- A1: `n = 7`, `n0 = 4`, `g = t^4 + t^2 + t + 1`, so `l = 2`.
- A1: `n = 31`, `n0 = 15` and `16` (Euler–Petit Type 2). `n = 73`, `n0 = 9`,
  `l = 1` (context, outside `n <= 61`).
- A2: the `n = 7` witnesses `X^2 + X` (`n0 = 3`) and `X^4 + X^2 + X`
  (`n0 = 4`) must appear. Every A1 hit at `n <= 23` with `j <= 1` must also
  appear in A2.
- A3: Euler–Petit Prop. 6 family 1 at composite `n`:
  - `X^8 - X^3` at `n = 4` (`beta = 0.704`);
  - `n = 6`, `n' = 5`, `a = 11`;
  - `n = 8`, `n' = 7`, `a = 43`;
  - `n = 8`, `n' = 6`, `a = 13`.
- D: a subfield base `F_{2^9}` in `F_{2^18}` against random 9-dimensional
  subspaces. Subfield bases are the known case where Weil descent behaves
  differently, so it checks that the instrument can see a structural
  difference.

Null control for D: 4 independent random `F_2`-subspaces of the same
dimension, built with the same code.

## Metrics

- `beta = l * n / n'^2` (Euler–Petit), and `l / n'`.
- HKPYY threshold for beating generic algorithms: `beta < 0.1025`.
  Definition 3.1 (QSP at all) requires `beta < 1`.
- D, yield: the fraction of 2000 uniform random target abscissae `x_R`
  (1000 at `n = 31`) that decompose as `x(P1) + x(P2)` with
  `x(P_i) in V`. The check is exact per target: for each `x1 in V`, solve
  the `S3` quadratic for `x2` and test `x2 in V`. Reported with a Wilson
  95% interval, together with the mean number of decompositions per target.
- D, refutation degree: on 6 targets per base that are verified
  non-decomposable, take the smallest `d` in `{3, 4}` (`{3, 4, 5}` at
  `n <= 18`) at which the degree-`d` F4 closure (`matrix_f4_f2` to a fixed
  point) contains `1`. If no such `d` exists, report "> d_max".

## Kill criterion (from the report, item 8)

Either of these closes the route:

- No family at any prime `n <= 47` (we search to 61) has
  `log2(deg lambda)/k` below the HKPYY threshold, i.e. `beta < 0.1025`.
- A family passes but its refutation degree is no lower than a random
  base of equal size.

## Budget

Everything runs single-threaded under `nice`, with `prlimit --as=700MB`.
No job runs longer than 45 min.
