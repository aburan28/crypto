# Protocol: Frobenius-orbit index calculus on Koblitz curves, fully priced, against the best rho (2026-10-14)

Frozen before any claim data. Conductor T-1. Class: **engineering**. A
constant can move; the expected outcome is that nothing else does (§6).

## 1. Question

On a Koblitz curve, Frobenius `φ` acts on the prime-order subgroup as
multiplication by `λ`. A Frobenius-stable factor base therefore needs one
column per signed orbit (`|F| / 2n` columns), and its pair table can be
stored modulo Frobenius. Rho gets the same `√(2n)` from the same automorphism
(`rho.signed_frobenius`). Expected costs for `m = 2`:

- IC with `c` orbit columns: about `n c²/2` table additions plus
  `r / (2 n² c)` trials, minimised at `c ∝ (r/n³)^{1/3}`, giving `∝ r^{2/3} / n`;
- Frobenius rho: `√(π r / 4n)`.

So **IC ÷ rho-frob ∝ r^{1/6} n^{-1/2}**: a constant-level lever, not an
exponent one. This protocol measures it instead of assuming it.

## 2. Changes under test

1. **Full pricing** (`src/cryptanalysis/ecbench/full_price.rs`). New opt-in
   method ids `ic.pipeline_full`, `rho.signed_frobenius_full` and
   `rho.signed_frobenius_strong_full` run the unchanged solvers, then price
   every native counter from `docs/ecbench/full_prices.json`, a table of
   ratios to one affine addition measured by `ecbench price-units` (median of
   15 passes, one host, 2026-10-14) for every unit and every panel curve:
   Artin–Schreier solve, hash probe (a `2^20`-entry `FastMap`, so priced at or
   above any table in the panel), row operation, Frobenius map,
   canonicalisation, field multiplication, squaring and inversion, word XOR.
   The orbit factor base's construction counters that the old methods left
   `_uncharged` (raw-x lifts, field squarings and multiplications of the
   lift and the orbit key, curve-equation checks, inversions, hash probes and
   inserts, negations, the base's bytes) are priced; so are rho-frob's
   canonicalisations and distinguished-point probes and rho-strong's
   canonicalisations, partition hashes and table probes. A counter the rules
   do not know stays unpriced and keeps the lower-bound flag. Additions and
   doublings are one unit each, as everywhere in ecbench. `icv1-f2m19-t797-b6cf2467`,
   which has no row in `docs/ic/calibration.json`, has one in the new table,
   so its IC rows are full costs. The pinned table and every existing method
   are unchanged, so committed sessions replay as before.
2. **stride-renew targets** (PR #1842, this branch's base) on the orbit factor
   base `compact-orbit-scan` with oracle `mitm-frobenius-counted:m=2`.

## 3. Exploration (not evidence)

Scratch sessions on target seeds 777101 and 777201 (2026-10-14; tables in
`exploration/`, sessions not kept) chose the column count `c` per curve
(the lowest mean S among the counts tried) and the negfold dimension `d`:

| a | n | slug | r | log2 r | cofactor | c tried → chosen | d tried → chosen |
|---|---:|---|---:|---:|---:|---|---|
| 0 | 13 | icv1-f2m13-t181-515ee569 | 2003 | 10.97 | 4 | 1, 2, 3 → **1** | 4 → **4** (d = 4 from PR #1842) |
| 1 | 17 | icv1-f2m17-tm101-00378d4e | 65587 | 16.00 | 2 | 1, 2, 3 → **1** | 4, 5, 6 → **6** |
| 0 | 19 | icv1-f2m19-t797-b6cf2467 | 130873 | 17.00 | 4 | 1–4 → **2** | 5, 6, 7 → **6** |
| 0 | 23 | icv1-f2m23-t5197-69e76b73 | 2095853 | 21.00 | 4 | 1, 2, 3, 4, 6, 8 → **3** | 7, 8, 9 → **8** |
| 1 | 29 | icv1-f2m29-tm40309-30c52b96 | 42457 | 15.37 | 12646 | 1, 2, 3, 4, 6, 9, 13 → **1** | 10, 11: S ≈ 1400–5600; excluded |
| 0 | 31 | icv1-f2m31-tm90707-c95f16f5 | 1439393 | 20.46 | 1492 | 2, 3, 5, 8, 12, 16 → **2** | 11, 12: S ≈ 1000–4000; excluded |
| 0 | 37 | icv1-f2m37-tm534059-32aad96b | 230603167 | 27.78 | 596 | 5, 8, 12 → **8** | not run; excluded |
| 0 | 41 | icv1-f2m41-tm2308219-7f48b14a | 549756390943 | 39.00 | 4 | 60, 90, 130 → **90** | not run; excluded |

The subspace base `koblitz-standard-subspace` is not projected by the
cofactor, so on the three large-cofactor curves its hit rate is set by the
whole group and it is hopeless (exploration S ≈ 10³); at n = 41 it would cost
about `r^{2/3} ≈ 2^{26}` additions a solve and a table of `~2^{25}` entries.
It runs only where the cofactor is ≤ 4 and n ≤ 23. These choices were made
on the exploratory seeds above and are disclosed as such; the claim uses fresh
seeds.

Exploratory ratios (orbit-renew ÷ rho-frob, mean S): 3.4, 2.2, 3.6, 4.4, 2.7,
3.5, 7.0, 13.9 in the row order above. Rho-strong was below rho-frob at
n = 41 only.

## 4. Sessions

One session per curve, specs from `make_specs.sh`:
`docs/ecbench/specs/claim-koblitz-frob-orbit-full-a<a>-n<n>.json`. Each has 16
public targets × 3 rounds (+1 warm-up), alternate order, fresh target seeds
202610142001–202610142015, L0 (`--cpus none --allow-busy`; this container has
no CPU reservation, so wall times are descriptive). Arms, all `*_full`:

- `rho-frob` (reference), `rho-frob-aa` (A/A control), `rho-strong`;
- `ic-orbit-c<c>-renew` (candidate): orbit base, `mitm-frobenius-counted:m=2`,
  `targets=stride-renew`;
- `ic-orbit-c<c>-walk`: the same with the incumbent `targets=walk`;
- `ic-negfold-d<d>-renew` (n ≤ 23, cofactor ≤ 4): the subspace base with
  `pdp3-koblitz:m=2,engine=rootfind-negfold`, `targets=stride-renew`.

`run.sh` runs each session, `ecbench verify --replay-all`, `ecbench table`,
and `ecbench compare --save` for the pairs below; ratios are B ÷ A paired by
target with cluster-bootstrap 95 % intervals. `ecbench trend` fits the
ratio across curves; `ecbench bound fit` fits each arm's own `gae = C · r^α`.

## 5. Hypotheses and decision rules

Every record must verify, every replay must match, and every record of every
arm must have `lower_bound = false` (no unpriced work). `rho-frob-aa ÷
rho-frob` must contain 1 on each curve, or that curve is reported unreliable.

- **H1 (full cost).** No record is a lower bound. Fails if any is.
- **H2 (stride-renew on orbits, per curve).** `orbit-renew ÷ orbit-walk`:
  helps iff the upper end < 1. Prediction: < 1 on the small curves, closer to
  1 as r grows (it removes a fixed set-up of ~45 scalar multiplications).
- **H3 (against rho-frob, per curve).** `orbit-renew ÷ rho-frob`. Prediction:
  > 1 (lower end > 1) on every curve. IC is below rho-frob on a curve iff the
  upper end < 1.
- **H4 (against rho-strong, per curve).** `orbit-renew ÷ rho-strong`,
  reported with the same rule; prediction: < 1 on the small curves, > 1 by
  log2 r ≈ 21 and above.
- **H5 (orbits against the subspace base).** `orbit-renew ÷ negfold-renew` on
  n = 13, 17, 19, 23: orbits are cheaper iff the upper end < 1.
- **H6 (trend; the verdict).** `ecbench trend` on the eight H3 comparisons:
  `free_r` (`ln ratio = c + β ln r`), `free_r_n` (`+ γ ln n`) and the model
  (`β = 1/6`, `γ = −1/2`, only `c` fitted). **Constant-only** if `β`'s 95 %
  interval in `free_r` lies above 0 (the ratio grows with r: the orbit base
  does not change IC's exponent relative to rho-frob). **Exponent-level** only
  if `β`'s interval lies at or below 0. Otherwise inconclusive. The model is
  consistent if its `χ²/dof ≤ 2`; the fitted `c` (the constant `e^c`) is
  reported either way. As a cross-check, `ecbench bound fit` of
  `ic-orbit-*-renew` and `rho-frob` across the panel: expected `α ≈ 2/3` and
  `≈ 1/2` (the per-curve `c` and `n` differ, so these are descriptive).

## 6. What cannot be claimed

`m = 2` decomposition is at best `r^{2/3}`; Frobenius divides its table and its
trials by `n`, and rho's walk by `√n`. No statement about larger fields,
ECC2K-130, deployed curves or security follows. The prices in
`full_prices.json` are one host's microbenchmark ratios; a different host
reprices the native share of each total (material for IC at small r,
where lookups and Frobenius maps are a large part of a short run) but not the
operation counts, which the records keep.
