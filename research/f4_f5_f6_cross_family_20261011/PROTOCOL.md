# F4, F5 and F6-IC across curve families: protocol

Status: **preregistered**, written before any measured run. The smoke run
that checked the plumbing is disclosed below.

## Question

How do the repository's Boolean Gröbner decomposers compare on two binary
curve families, in the same counted unit and on the same public one-target
workloads?

- **Decomposers:**
  - matrix F4, built from scratch at every node;
  - matrix F5, which is F4's matrices less the rows the F5 criterion
    predicts reduce to zero;
  - inherited F4;
  - F6-IC, which is inherited F4 with exact factor-base geometry inside the
    splitting tree.
- **Families:**
  - **Koblitz** curves, `E_0: y² + xy = x³ + 1`;
  - **generic binary** curves, `y² + xy = x³ + ax² + b` with random `b ∉ F_2`.

Prime-field curves have no algebraic decomposition oracle in `ecbench`. They
enter the report only through existing evidence, labelled as not comparable
(see *Prime curves*).

The user's focus is F6-IC. Until this round it had only ever run on
Koblitz curves.

## What changed to make it measurable

Only the `ecbench` binary's PDP3 plug-in changes
(`src/bin/ecbench/pdp3_koblitz.rs`). The library does not.

1. **Engines.** `engine=matrix-f4|matrix-f5` join `inherited-f4|f6-ic`.
   - Every engine runs under its own default splitting rule
     (`SplitRule::Auto`): the inherited engines split on the highest free
     variable, the matrix engines on the lowest. An arm is therefore a whole
     solver choice, as in #1333's three-way comparison, not one kernel on a
     shared tree.
2. **Generic binary curves.** The plug-in now runs on any binary curve with
   `a ∈ {0, 1}`.
   - **How:** a non-Koblitz instance is seen through a Frobenius-free
     `KoblitzCurve` view: `frobenius_is_endomorphism = false`, `λ = 1`, and
     the curve's own order, subgroup and generator. This is the same shape
     `KoblitzCurve::isogenous_model` gives a floor curve.
   - **Why that is enough:** nothing the decomposers read depends on
     Frobenius. The Semaev system takes only `b`, the base is folded by
     negation alone, and F6-IC's gate needs only the subspace, base
     membership and group addition.
   - **Names:** the base and the oracle get the names
     `binary-standard-subspace` and `pdp3-binary`. The old names still
     resolve, so the F6-IC ladder's method identities are unchanged.

## Curves (all in `docs/curves/registry.json`)

| family | m | slug | r | cofactor | base dimension d |
| --- | --: | --- | --: | --: | --: |
| Koblitz `E_0` | 13 | `icv1-f2m13-t181-515ee569` | ≈ 2^10.97 | 4 | 5 |
| Koblitz `E_0` | 19 | `icv1-f2m19-t797-b6cf2467` | ≈ 2^17.00 | 4 | 7 |
| random binary (registry `char2:13:59297:8`) | 13 | `icv1-f2m13-tm9-6863a838` | ≈ 2^10.42 | ≤ 8 | 5 |
| random binary (registry `char2:19:59297:8`) | 19 | `icv1-f2m19-tm649-a9133fe6` | ≈ 2^18.00 | ≤ 8 | 7 |

Within each size the two families' subgroup orders differ by up to a factor
of 2 (`r` is not matched). The factor base is the same construction on both:
every point whose abscissa lies in `span(1, z, …, z^{d−1})`. The bases'
point counts differ by curve and are reported.

## Specs (frozen in this directory)

| spec | curves | targets | arms | per-run cap |
| --- | --- | --: | --- | --: |
| `SPEC-K13` | Koblitz m = 13 | 8 | strong rho, rho-negation, matrix F4, matrix F5, inherited F4, F6-IC | 3600 s |
| `SPEC-R13` | random m = 13 | 8 | rho-negation, matrix F4, matrix F5, inherited F4, F6-IC | 3600 s |
| `SPEC-K19` | Koblitz m = 19 | 8 | strong rho, rho-negation, inherited F4, F6-IC | 3600 s |
| `SPEC-R19` | random m = 19 | 8 | rho-negation, inherited F4, F6-IC | 3600 s |
| `SPEC-X19` | both m = 19 curves | 1 (the first target of K19 and R19) | matrix F4, matrix F5 | 1800 s |

Common settings:

- **Targets:** all specs share target seed 202610111901 and measurement
  seed 202610111902. There is one round, and targets are public.
- **IC arms:** `pdp3-binary:m=3,degree=3,node_budget=8192`,
  `max_trials=200000`.
- **Strong rho:** the F6-IC ladder's reference (32 lanes, `dp_bits` 4).
- **Rho-negation:** the generic reference. It is also run on Koblitz curves,
  so both families share one reference.
- **Matrix arms at m = 19:** at m = 13 they cost 100–250 times inherited F4
  (smoke run below), so at m = 19 they get one target per family and a
  1800 s cap. A timeout there is recorded as a timeout, not as a cost.

## Units and accounting

- **`S`** = charged group-addition equivalents / √r, the `ecbench` unit.
  IC `S` is a **lower bound**: the Boolean solver's word XORs are counted
  but unpriced. F6-IC's geometric point additions are charged at one
  addition each.
- **Primary solver measure:** counted word XORs per run (`solver_ops`), which
  is deterministic.
- **Secondary:** relation-phase wall time, a practicality note only. This
  is an L1 cloud VM; `ecbench` records the isolation actually earned.
- **Boundaries:**
  - the floor `√(π/4m)` for Koblitz (automorphism group of order 2m);
  - `√(π/4)` for generic binary curves (negation only);
  - each family's rho measured on the same targets.
- **Correctness:** every IC record must recover and verify `[k]G = Q`.
  `ecbench verify --replay-all` must replay every session identically.

## Questions and decision rules

Ratios are per target over the 8 workloads; medians and ranges are
reported.

1. **F5 against F4 (matrix), m = 13, each family.**
   - Statistic: `W(matrix-F5) / W(matrix-F4)`.
   - **"F5 cuts F4's work"** on a family if the ratio is below 1 on all 8
     targets, and the relations found and solver call counts are identical
     (same tree).
   - Otherwise report the ratio.
2. **Inherited against matrix F4, m = 13, each family.**
   - Statistic: `W(matrix-F4) / W(inherited F4)`. Reported as is.
3. **F6-IC against inherited F4, each family and size.**
   - Statistic: `log₂(W(inherited) / W(F6-IC))`.
   - **"F6-IC saves solver work"** on a family and size if the median is at
     least 0.1 and all 8 targets are above 0. **"No saving"** if the median
     lies within ±0.05. **"Mixed"** otherwise.
   - Also reported: lower-bound `S`, `S/rho`, the geometric additions
     charged, and relation-phase wall time.
4. **Any arm against rho.**
   - **"IC below rho"** for an arm and cell only if its lower-bound `S` is
     below the same-target rho reference on all 8 targets.
   - A lower bound above rho proves IC slower; a lower bound below rho
     proves nothing, because the word XORs are unpriced.
5. **Matrix arms at m = 19 (X19).**
   - Each run either completes or reaches the 1800 s cap.
   - Report counts for completed runs and the cap for the others.
   - No decision rests on X19.

**Stop conditions.**

- A session in which any IC run fails verification is reported as
  failed and enters no ratio.
- A cell whose runs time out on more than half its targets gets no median.

## Prime curves

These are not measured here. The report collects the existing evidence:

- the PKM tower F4 against signature-F5 study
  (`research/notes/index-calculus/RESEARCH_PKM_TOWER_ORACLE.md` §13);
- the `E(F_{p³})` residual-walk ledger (`RESEARCH_RESIDUAL_WALKS.md` §11);
- the GLV orbit-system study (`RESEARCH_GLV_INDEX_CALCULUS.md`).

All three are labelled with their own units. There is no F6 for prime
curves: an `F_p` factor base is not an `F_2`-linear subspace, so the gate's
encoder does not apply as it stands.

## Smoke run (before this protocol; plumbing only)

One target per curve on `icv1-f2m13-t181-515ee569` and an unregistered
random m = 13 curve (seed 1). All 10 records verified. Word XORs, with
relation-phase wall time in brackets:

| m = 13 | matrix F4 | matrix F5 | inherited F4 | F6-IC |
| --- | --: | --: | --: | --: |
| Koblitz | 2.78 × 10⁹ (36.7 s) | 1.78 × 10⁹ (30.7 s) | 1.15 × 10⁷ (0.21 s) | 9.46 × 10⁶ (0.13 s) |
| random (seed 1) | 5.23 × 10⁸ (7.3 s) | 3.11 × 10⁸ (5.8 s) | 4.24 × 10⁶ (0.07 s) | 4.10 × 10⁶ (0.06 s) |

The smoke curve is not one of this round's curves, and none of its figures
enters a decision.

## Commands

```sh
cargo build --release --bin ecbench
for s in K13 R13 K19 R19 X19; do
  target/release/ecbench run --spec research/f4_f5_f6_cross_family_20261011/SPEC-$s.json \
    --out research/f4_f5_f6_cross_family_20261011/sessions/$s --wait
  target/release/ecbench verify --dir research/f4_f5_f6_cross_family_20261011/sessions/$s --replay-all
done
```
