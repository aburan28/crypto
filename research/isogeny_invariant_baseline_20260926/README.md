# Eligible isogeny derived invariant sets: exact small field benchmark

This benchmark migrates a completed local small parameter experiment into the crypto repository. The run predates this protocol file; it is historical baseline evidence, **not a preregistered performance comparison**. The frozen commands and an audit script make the figures reproducible. The two primary constructions are the torus isogeny fraction flags of [Couveignes and Lercier, §7](https://arxiv.org/abs/0802.0282), also studied for Koblitz curves by [Galbraith, Granger, Merz and Petit, §§4.2, 5.2](https://sacworkshop.org/SAC20/files/preproceedings/18-IndexCalculus.pdf).

## Question, boundary and frozen parameters

Hypothesis under test: at eligible small parameters, an isogeny derived nonlinear Frobenius invariant subset has a compact fraction description, correct summation equation lifting, and enough *actual curve points* to warrant a matched solver experiment. The comparison reference is **every** three dimensional F3 linear subspace: 40 total, including the two Frobenius invariant subspaces. All sets have 27 **field** elements. There is no measured rho or full ECDLP reference in this round; the normalized `S`, ratio to rho, total operation cost and speedup remain null.

- Small field: F3[w]/(w⁴-w³+w+1), nonsquare D=2, torus multiplication by 4, fiber above a=1, Frobenius `w³=(2w+2)/(w+2)`.
- Fraction base: V₁={(a+bw)/(c+dw): coefficients in F3, (c,d)≠(0,0)}. Distinct point keys come from the **evaluated** field set, not from repeated parameter tuples. The denominator is nonzero for every admissible pair.
- Curve controls: all six short models `y²=x³+ax+b`, a∈{1,2}, b∈{0,1,2}, over F81. These are supersingular, and do not represent all curves over F3. Each target is every affine point on the particular curve. Group identity is counted in the point addition table but excluded from the coverage denominator. Repeated points and opposite pairs are allowed; the raw file also retains a distinct-x diagnostic. No subgroup restriction.
- Construction control: the published q=13, degree=7, D=2, fiber a=8 polynomial, with Frobenius `w¹³=(4w+2)/(w+4)`; verify its k=1 cardinality and orbits. This control does not include curve point decomposition.
- Workload: m=2, all **ordered** pairs of affine base points. For every pair of liftable x coordinates and every liftable target x coordinate, test Semaev f3 against the existence of compatible rational y lifts obtained independently by group addition.

The finite field irreducibility criterion, torus fiber, full set Frobenius closure, curve table associativity and lifting equivalence are checked in the benchmark itself. Acceptance for this bounded benchmark requires six complete curves, all 40 subspaces, 325,998 lifting checks without disagreement, and the published 2,197 element fraction set with 13 fixed elements plus 312 full length seven orbits. Any disagreement, missing curve, or representation counted as a distinct point is a failure. Wall seconds are recorded only as diagnostics and are not an algorithm comparison. Point counts and target coverage are the shared units in the table below. The underlying theoretical floor for the complete relation collection and the same group rho reference have not been established for these six cases, so no attack efficiency ratio is reported.

## Exact result (run 001)

Each cell is `affine factor base points / covered affine targets` using all field values of the named 27 element set. W1 and W2 are the two invariant hyperplanes, kernels of coordinate linear forms (0,1,1,2) and (1,1,1,1) in (1,w,w²,w³).

| Curve `y²=x³+ax+b` | Affine targets | Torus V1 | Invariant W1 | Invariant W2 |
|---|---:|---:|---:|---:|
| a=1, b=1 | 63 | 5 / 7 | 21 / 63 | 30 / 63 |
| a=1, b=0 | 63 | 5 / 7 | 21 / 63 | 3 / 3 |
| a=1, b=2 | 63 | 5 / 7 | 21 / 63 | 30 / 63 |
| a=2, b=0 | 63 | 51 / 63 | 3 / 3 | 21 / 63 |
| a=2, b=1 | 90 | 30 / 90 | 30 / 90 | 30 / 90 |
| a=2, b=2 | 90 | 30 / 90 | 30 / 90 | 30 / 90 |

V1 is nonlinear (the raw result includes an additive closure counterexample). Its 27 field values have three fixed elements and six four-element Frobenius orbits. Its 36 normalized numerator/denominator tuples represent each constant four times and each nonconstant once. At q=13, degree 7, the 2,366 tuples represent each constant 14 times and each nonconstant once. Thus factor base size must come from unique evaluated field values and liftable curve points, not tuple count. The tiny curves show that point admissibility is curve dependent; the torus base ranges from 5 to 51 affine points. This is a **stage diagnostic and accounting result**, not a solver speed comparison or proof that either construction dominates.

## Reproduce and audit

```sh
python3 research/isogeny_invariant_baseline_20260926/benchmark.py \
  --output /tmp/isogeny-run-002.json
python3 research/isogeny_invariant_baseline_20260926/summarize.py \
  /tmp/isogeny-run-002.json --output /tmp/isogeny-summary-002.json
```

The programs use the Python standard library only. Outputs are new files (`x` mode); reruns cannot overwrite evidence. `results/run-001.json.gz` contains gzip-compressed complete curve point order, all target counts including failed coverage, all 40 linear comparator rows on all six curves, checks and construction receipts. `results/summary-001.json` holds 18 normalized rows with full DLP costs null and the SHA-256 of the benchmark source. Retain any future run in a separately numbered file; timing may differ, but exact counts must match the summary.

Current checked artifact hashes:

```text
86683cf33aec362accaef9ff75c6b4aeda64953e2bff165c1cd46059e4005374  benchmark.py
bc74d09adafb8b4c067a58edf9b7aa0f0398314f51dfbd2ead644e3178343c2a  summarize.py
7653ec502a446fd0bdc46044a3c81359fa0d81f5db91bdf675ee745c83fba557  results/run-001.json.gz
f889db266f0e9de1fb8f93cc8d1d7f2f23d405e82d0fcc17d7b1bc292df2241c  results/summary-001.json
```

## Next experimental gate

A future candidate needs *matched* m=3 Semaev polynomial systems at eligible parameters, with equivalent vector space bases, fixed random and planted targets, explicit grouping by Frobenius orbit, separately priced parameterization, failed solves, point lifting and verification, and an independent enumerative control. Freeze target seeds and solver configurations **before** running it. A meaningful claimed speedup additionally needs the repo's equivalent full DLP suite, matched rho, precomputation, matrix work and scalar verification; this m=2 enumeration is outside the frozen binary WDSat corpus and cannot be spliced into its `S` axis. Do not infer performance at the binary 131 degree challenge from this odd characteristic degree four control.
