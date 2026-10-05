# GLV-2 with an isogeny-chain endomorphism on GOST CryptoPro-B, measured against width-w NAF

A constructive scalar-multiplication measurement on the GOST R 34.10-2001
**CryptoPro-B** parameter set (RFC 4357 §11.4.4): `y² = x³ - 3x + b` over
`F_p`, `p = 2^255 + 3225`, prime group order `n`, cofactor 1.  The curve has
CM by the order of discriminant `-619` (class number 5), so its endomorphism
ring contains no cheap degree-small map; the endomorphism sweep in the
autoresearcher (`research/endosweep_20261005` of `crypto-autoresearcher`)
found that the cheapest usable endomorphism is the element `4 + ω`
(`ω = (1 + √-619)/2`) of norm `175 = 5·5·7`, evaluated as a chain of two
5-isogenies and one 7-isogeny through neighbouring curves and back to a
curve isomorphic to `E`, followed by the isomorphism `(x, y) → (u²x, u³y)`.
On the group it acts as multiplication by a 256-bit scalar `λ`, which gives a
2-dimensional GLV decomposition `k = k1 + k2·λ (mod n)` with
`|k1|, |k2| < 2^128`.  A second chain, the element `ω` itself of norm
`155 = 5·31` (a 5-isogeny then a 31-isogeny), is measured as a further arm.

This is **not an attack**: it computes `k·P` faster on the same curve, it does
not change any problem's difficulty.  The `S = ops/√n`, matched-rho and
scoreboard rules of `AGENTS.md` §§1–8 therefore do not apply; the one table /
one unit / correctness-column discipline (§2), the class label (§3), the
clean-baseline and isolation rules (§10) and the frozen-input rules do.

## Hypothesis, stated before the measurement

> On identical Jacobian `a = -3` arithmetic, GLV-2 with the degree-175 chain
> endomorphism of CryptoPro-B is faster than width-w NAF; a Python
> operation-count model predicts 1.45×.

The model is the autoresearcher's `costmodel.py` (endosweep Table 1):
1681 M for GLV-2 with `4 + ω` against 2439 M generic, i.e. 1.45×, with the
chain costed at `2·(8M+2S) + (12M+2S) + 20M ≈ 53 M` in affine Kohel form.
The same model ranks `4 + ω` (53 M) ahead of `ω` itself (`5·31`, 91 M as a
chain), so the secondary expectation is that the degree-175 arm beats the
degree-155 arm.  The endosweep README is explicit that its 1.45× is modelled
and that "a projective evaluation of the 5-5-7 chain in a real library is
the measurement that would settle the number"; this directory is that
measurement.

## Frozen inputs

| file | sha256 | content |
|:--|:--|:--|
| `frozen_inputs/cryptopro_b_chain_5_5_7.constants.json` | `3ad2ca805b7181cc8c183852d42b3bbbf0b4fa7bfe9a098f1173b2f3eeb59ade` | the degree-175 chain `5·5·7` (element `4 + ω`): curve, kernel polynomials `ψ`, numerators `N`, `M`, codomain coefficients per step, `iso_u`, `λ`, the LLL-reduced GLV basis, the Babai bound, four test vectors |
| `frozen_inputs/cryptopro_b_chain_5_31.constants.json` | `de3ea720d37cdb3c46e4560d848b10be49496de715bf983ea5cc3a9a932a5272` | the degree-155 chain `5·31` (element `ω`), same fields |
| `frozen_inputs/SHA256SUMS` | — | the two hashes above, `sha256sum -c` form |

The constants were generated and verified outside this repository (the
autoresearcher's `export_chain_constants.py` / `emit_rust_constants.py`;
Python provenance, kept as frozen input per `AGENTS.md`).  The Rust module
`src/ecc/cryptopro_b_chain_consts.rs` was emitted from those two JSON files
and copied here **unchanged**; its header records both hashes.  Every field
element in it is a canonical little-endian `[u64; 4]` (not Montgomery form)
and polynomials are listed low-degree first.  The curve parameters come from
the existing `CurveParams::gost_cryptopro_b()` in `src/ecc/curve_zoo.rs`;
the tests check the frozen `p, a, b, n` against it.

Nothing in this directory was computed by Python; the benchmark, the
arithmetic, the decomposition and every check are Rust in the crate.

## Reference and arms

All arms run in one binary on the same variable-time Jacobian `a = -3`
arithmetic (`src/ecc/cryptopro_b_point.rs`: EFD `dbl-2001-b` 3M+5S,
`add-2007-bl` 11M+5S, mixed `madd-2007-bl` 7M+4S) over one Montgomery-form
field (`src/ecc/cryptopro_b_field.rs`, built on the crate's `U256`
primitives).  The comparison is therefore between scalar-multiplication
*strategies*, not between arithmetic styles.

* **Reference / baseline**: width-`w` NAF, `w = 5`.  Table of the odd
  multiples `P, 3P, …, 15P` (one doubling, seven full additions), converted
  to affine with **one** batched inversion, then one doubling per digit and
  one mixed addition per non-zero digit.
* **GLV-2 total, chain `5·5·7`** (the hypothesis' arm): Babai rounding of
  `k` against the frozen basis with `num-bigint`, `φ(P)` by the chain, the
  two odd-multiple tables (two doublings, fourteen full additions) converted
  with **one** batched inversion over all sixteen points, wNAF digits of
  `|k1|` and `|k2|`, and the interleaved loop with shared doublings
  (negative digits and negative halves use the negated table point).
* **GLV-2 total, chain `5·31`**: the same with the degree-155 chain, its own
  `λ` and basis.
* **Stage diagnostics**, labelled as such and never a speed: the
  decomposition alone and `φ(P)` alone, per chain.
* **A/A control**: the baseline routine timed a second time as a separate
  arm, so its ratio to the baseline is the noise floor of the run
  (`AGENTS.md` §10).

Chain evaluation, per step, in projective `P²` coordinates `(X : Y : Z)`,
`x = X/Z`, `y = Y/Z`, with the homogenised polynomials evaluated by Horner
over precomputed powers of `Z`:
`X' = N_h·ψ_h`, `Y' = Y·M_h`, `Z' = Z·ψ_h³`; after the last step the
isomorphism `X'' = u²X`, `Y'' = u³Y`; then to Jacobian without an inversion
as `(X·Z, Y·Z², Z)`.  Degrees: `ψ` has degree `(ℓ-1)/2`, `N` degree `ℓ`, `M`
degree `(3ℓ-3)/2` (for `ℓ = 31`: 15, 31, 45).

## Cost accounting

Everything input-dependent is inside each arm's timed call.  The GLV arms'
totals include the decomposition, `φ(P)`, both tables, the shared batched
inversion, both digit expansions and the main loop.  Outside every arm, by
design and identically for all of them:

* a `GlvContext` per chain, built once per process: the chain's
  coefficients converted to Montgomery form, the basis parsed, `λ` and `n`
  as big integers — curve-constant setup, the analogue of the hardcoded
  Montgomery-form coefficients of the other point modules;
* the final affine conversion of the result (one inversion), which would be
  the same fixed cost in every arm; the timed output is a Jacobian point and
  the correctness checks compare points by cross-multiplication.

The decomposition is deliberately the plain exact-rational Babai rounding
with `num-bigint` (heap-allocating big integers, two 384-bit divisions); the
usual precomputed-rounding-constant trick was not implemented, so the
decomposition stage is an upper bound on what GLV has to pay there, and it
is reported separately.

## Protocol

`M = 256` seeded `(scalar, point)` pairs (`StdRng`, seed `20261005`; points
are `r·G` for random `r`, scalars uniform in `[1, n)`), `R = 20` rounds, in
each round every arm is timed over all `M` pairs, the arm order alternates
between rounds (forward on even rounds, reversed on odd rounds), and the
table reports the median and the minimum over rounds of the per-operation
time.  Before any timing the binary verifies, on every pair, that each GLV
arm equals the baseline, that the baseline equals the crate's textbook
`BigUint` implementation (`Point::scalar_mul_vartime`), that `φ(P) = λ·P` for
each chain, and that every decomposition satisfies `k1 + k2·λ ≡ k (mod n)`
within the Babai bound; the exit status is non-zero if any check fails.

**Success condition** (declared in advance): measured ratio baseline/GLV
> 1 for the `5·5·7` arm, with the median and the minimum of both arms on
the same side (both ratio columns above 1), and every correctness check
passing.  **Stop condition**: one benchmark run of the stated size.

## Results

**Pending** — the measured run has not been executed yet; this section is filled by a follow-up commit from the run itself.

## Class

**Pending** — filled by the follow-up commit that lands the run.

## Honesty notes

* **Variable time, everywhere.**  Identity checks, the `H = 0` branch in the
  addition law, wNAF digit scans, the `num-bigint` decomposition and the
  textbook reference all branch on their inputs.  The module docs say so.
  This is a measurement harness; nothing here is for secret scalars.
* **Unoptimised pieces.**  The field multiplication is the crate's generic
  const-generic schoolbook `U256::mont_mul` (no dedicated squaring, no
  assembly); the decomposition is `num-bigint`; the chain uses no special
  case for `Z = 1` in the first step.  Both arms share the first point, so
  the ratio is what it is on this arithmetic; a faster field would change
  both arms' absolute numbers and could move the ratio either way, since the
  decomposition is a fixed cost that does not scale with the field.
* **One machine, one hardware class.**  See the host line of the run.  A
  shared virtual CPU on a cloud container; expect ±5–10 % residual wall
  noise, which is what the A/A control row measures.  Not measured on Arm64
  or on any other x86-64 host.
* **Wall time, not operation counts.**  The unit is nanoseconds per scalar
  multiplication because that is what the hypothesis is about; the field
  operation count of each arm is fixed by the formulas above and does not
  depend on the host.

## Commands

**Pending** — filled by the follow-up commit that lands the run.

## Files

| file | what |
|:--|:--|
| `README.md` | this record |
| `frozen_inputs/*.json`, `frozen_inputs/SHA256SUMS` | the frozen chain constants |
| `bench_w5.json` | the full JSON record of the run (times per round, host, build, correctness) |
| `bench_w5.stdout.md` | the bench's standard output, verbatim |
| `isolated_bench_record.jsonl` | the `isolated_bench` conditions record of the run |
| `src/ecc/cryptopro_b_chain_consts.rs` | the constants module (copied unchanged) |
| `src/ecc/cryptopro_b_field.rs` | the Montgomery-form field |
| `src/ecc/cryptopro_b_point.rs` | Jacobian arithmetic, wNAF, the chain, the decomposition, both arms |
| `src/bin/cryptopro_b_glv_bench.rs` | the benchmark |
