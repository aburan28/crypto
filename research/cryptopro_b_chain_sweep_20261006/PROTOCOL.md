# CryptoPro-B chain sweep — protocol, stated before measuring

This file is committed before the counting run and the timing run it
describes.  The results go into `README.md` in a later commit; this file is
not edited after the runs.

## Question

`research/cryptopro_b_glv_chain_20261005` (PR #1408) measured GLV-2 on GOST
R 34.10-2001 CryptoPro-B with one chain realisation of one endomorphism
(`4 + ω`, two 5-isogenies and a 7-isogeny, projective evaluation) at one
setting of the scalar multiplication (width-5 NAF, affine tables, in both
arms): **1.252× (median) / 1.337× (min)** against a modelled 1.45×.  Its
reconciliation named three costs: the chain (≈ 128 M), the batched
inversion both arms pay (a Fermat inversion is 264 M here), and the second
table.  This experiment sweeps those parameters:

* **which chain** — every element and every order of its isogeny steps,
  within 2× of the cheapest (the frozen set below);
* **which chain evaluator** — the projective evaluator of PR #1408
  (`generic`) and a Jacobian one (`optimised`);
* **which scalar-multiplication settings** — wNAF width 3–7 and the table
  kept affine (one batched inversion, mixed additions) or Jacobian (no
  inversion, full additions), for the baseline as well as for GLV.

It is constructive (`k·P` on the same curve, faster), so the class is
**engineering** (§3); the `S = ops/√n`, matched-rho and scoreboard rules of
§§1–8 do not apply.  The one-table / one-unit / correctness-column rule
(§2), the clean-baseline and isolation rule (§10), the frozen-input rules
and the no-Python rule apply.

## Frozen inputs

| file | sha256 | content |
|:--|:--|:--|
| `frozen_inputs/chains.constants.json` | `683f92d0e56f5606f78b6c685f747a7b67106a7b5e5883e81869c8d6199d7336` | the chain set: 22 orderings of 11 elements, format `endosweep-chainsweep/1` |
| `frozen_inputs/chainsweep_model.json` | `9c93becb5a3421413d86c5577fb2d0e29c299fd87217583760a00efdad20ca0e` | the prediction: expected field-operation counts of every arm below |
| `frozen_inputs/SHA256SUMS` | — | the two hashes, `sha256sum -c` form |

Both were generated outside this repository by the autoresearcher's
`harness/endosweep/chainsweep.py` (aburan28/crypto-autoresearcher, branch
`endosweep-chainsweep-20261006`, commit `6d93b47a6a178a7070e58f417796b5ecc683012a`),
which is Python and stays there: frozen input with Python provenance.
Nothing in this directory is computed by Python.

How the set was chosen, from that generator's record: the catalogue is
every primitive non-scalar element `a + bω` of the maximal order of
discriminant −619 with norm ≤ 20 000 and every prime factor ≤ 61 (76
elements), and every distinct ordering of each element's prime-degree steps
(381 orderings), each built on the curve and verified on points.  Under the
`optimised` evaluator's operation count the cheapest is 80 M_eq, and every
chain outside the catalogue costs at least 173 M_eq.  The set is every
ordering costing at most 2 × 80 = 160 M_eq, so it contains every chain that
could come within 2× of the cheapest.

The loader (`src/ecc/cryptopro_b_chain_set.rs`) refuses a set whose curve
constants differ from this crate's, whose field elements are not canonical,
whose steps do not compose from `E` back to `E`, whose polynomials are not
monic of the right degrees, or whose GLV basis is not in `λ`'s lattice.  The
benchmark refuses to run unless the file's sha256 is the one above
(`--expect-sha256`).

## Arms

All arms run in one binary on the same variable-time Jacobian `a = −3`
arithmetic of `src/ecc/cryptopro_b_point.rs` (EFD `dbl-2001-b` 3M+5S,
`add-2007-bl` 11M+5S, `madd-2007-bl` 7M+4S) over one Montgomery-form field.

| group | arms | what |
|:--|--:|:--|
| baseline | 10 | width-`w` NAF, `w` ∈ {3, 4, 5, 6, 7}, table `affine` or `jacobian` |
| GLV grid | 20 | chain `4+1w/7.5.5`, evaluator `generic`/`optimised` × table × `w` |
| GLV per element | 10 | every other element's cheapest ordering at `optimised/jacobian/w5` |
| continuity | 1 | `4+1w/5.5.7` at `generic/affine/w5`: PR #1408's configuration (the `generic` evaluator does not depend on the order, so this is PR #1408's arm on the same walk up to conjugation) |
| stage diagnostic: `φ(P)` | 44 | every chain in the set × both evaluators |
| stage diagnostic: decomposition | 11 | the reduced basis of every element |
| A/A control | 1 | the reference baseline timed again |

The **reference** for every ratio column is the baseline `jacobian/w5`, the
cheapest baseline in the model.  The evaluators:

* `generic` — PR #1408's: projective steps, homogeneous Horner over powers
  of `Z`, `X' = N·ψ, Y' = Y·M, Z' = Z·ψ³`, isomorphism `u²X, u³Y`, then to
  Jacobian.  `11s + 2ℓ + 3` M + 1 S per `ℓ`-step (`s = (ℓ−1)/2`) and 4 M +
  1 S after the last.
* `optimised` — Jacobian steps: with the polynomials homogenised in
  `(X, Z²)` a step is `X' = N, Y' = Y·M, Z' = Z·ψ`; `ψ`, `N`, `M` are monic;
  the first step takes the affine input (`3ℓ − 4` M); later steps
  `11s + 2ℓ − 2` M + 1 S; the isomorphism is `Z·u⁻¹` (1 M).

## Units and accounting

* **Primary: counted field operations** per scalar multiplication
  (per call for stage diagnostics), mean over the `M` pairs, from a counting
  build of the same binary (`--features cryptopro-b-opcount`, a thread-local
  counter in `CryptoProBFieldElement`) on the same inputs.  `M_eq = M + S`:
  `mont_sqr` is `mont_mul` in this field.  An inversion is counted through
  its 256 squarings and 8 multiplications.  Counts are deterministic and
  host-independent.
* **Secondary: nanoseconds** per scalar multiplication, median and minimum
  over `R` rounds, from the default (non-counting) build, isolated.
* Inside every GLV arm: the decomposition (`num-bigint` Babai rounding, no
  field operations, so it appears only in the timings), `φ(P)`, both tables,
  the batched inversion when tables are affine, both digit expansions and
  the main loop.  Outside every arm, identically: the per-chain
  `GlvContext` (Montgomery-form constants, `u⁻¹`, `λ`, basis), and the final
  affine conversion of the result.

## Predictions (from `chainsweep_model.json`; derived, not measured)

* Every `φ(P)` arm's counted `M` and `S` equal the model exactly; for
  `4+1w/7.5.5`: `generic` 124 M + 4 S, `optimised` 78 M + 2 S.  The
  cheapest stage is `4+1w/7.5.5` `optimised` (80 M_eq); under `optimised`,
  `0+1w/31.5` (121) costs about half of what `5·31` would (238, not in the
  set); under `generic` every ordering of an element costs the same.
* Expected counted M_eq per scalar multiplication:

| arm | M_eq | | arm | M_eq |
|:--|--:|:--|:--|--:|
| baseline jacobian/w5 (reference) | 2 807 | | GLV 7.5.5 optimised/jacobian/w5 | 1 991 |
| baseline jacobian/w6 | 2 837 | | GLV 7.5.5 optimised/jacobian/w4 | 1 999 |
| baseline jacobian/w4 | 2 884 | | GLV 7.5.5 generic/jacobian/w5 | 2 039 |
| baseline affine/w5 (PR #1408's) | 2 916 | | GLV 7.5.5 optimised/affine/w4 | 2 066 |
| baseline affine/w4 | 2 922 | | GLV continuity 5.5.7 generic/affine/w5 | 2 203 |
| baseline jacobian/w7 | 3 016 | | GLV 9+1w/7.5.7 optimised/jacobian/w5 | 2 005 |
| baseline affine/w3 | 3 019 | | GLV 0+1w/31.5 optimised/jacobian/w5 | 2 030 |
| baseline affine/w6 | 3 031 | | GLV 37+3w/23.5.5.5 optimised/jacobian/w5 | 2 073 |
| baseline jacobian/w3 | 3 058 | | GLV 7.5.5 generic/affine/w7 | 3 187 |
| baseline affine/w7 | 3 345 | | | |

* Counted ratio, cheapest baseline over cheapest GLV: **1.41**
  (2 807 / 1 991), against 1.32 for PR #1408's configuration on the same
  counts (2 916 / 2 203).

## Success conditions (declared now)

1. **Correctness.**  Both runs pass every check: every baseline arm equals
   the textbook `BigUint` `k·P` on every pair; every GLV arm equals it too;
   every `φ(P)` arm equals `λ·P` on every point and reproduces the frozen
   test vectors; every decomposition satisfies `k1 + k2·λ ≡ k (mod n)`
   within its Babai bound.  Without this no number below is evidence.
2. **The count model is exact where it claims to be.**  Every `φ(P)` arm's
   counted `M` and `S` equal the model's; every total arm's mean counted
   M_eq is within 1 % of the model's expectation.
3. **Counted headline.**  The cheapest GLV arm's counted M_eq is at most the
   cheapest baseline arm's divided by 1.35.
4. **Wall-time headline.**  The GLV configuration declared here,
   `4+1w/7.5.5/optimised/jacobian/w5`, beats the fastest baseline arm (the
   baseline with the lowest median, picked after the run, which favours the
   baseline) by a median ratio above 1.252 and a minimum ratio above 1.337,
   PR #1408's values; and it beats the continuity arm by more than the A/A
   noise floor (`1 − A/A median ratio`) on medians.
5. **Chain ranking.**  The lowest `φ(P)` median of the 44 stage arms is
   `4+1w/7.5.5` `optimised`.

A condition that fails is reported as failed.  Inadmissible: changing the
arms, inputs, widths, seed or accounting after seeing a result; rerunning
to improve a number; dropping an arm.

## Stop conditions

One counting run (`M = 256`, seed `20261006`) and one timing run
(`M = 256`, `R = 20`, same seed) through `isolated_bench run` on CPU 3 with
`--max-other-cpu 0.5` (PR #1408 needed that budget because the agent runtime
itself exceeds the default 0.1 CPUs; the record reports what was actually
used).  A run the isolation tool refuses may be retried; a run it admits is
final.  Every attempt is recorded.

## Commands

```sh
cargo test --release --lib cryptopro_b
cargo test --release --features cryptopro-b-opcount --lib cryptopro_b
cargo build --release --bin cryptopro-b-chain-sweep --bin isolated_bench
CARGO_TARGET_DIR=target/opcount cargo build --release \
  --features cryptopro-b-opcount --bin cryptopro-b-chain-sweep

D=research/cryptopro_b_chain_sweep_20261006
SHA=683f92d0e56f5606f78b6c685f747a7b67106a7b5e5883e81869c8d6199d7336

# counting run (deterministic; no isolation needed)
target/opcount/release/cryptopro-b-chain-sweep --count-only \
  --chains $D/frozen_inputs/chains.constants.json --expect-sha256 $SHA \
  --pairs 256 --seed 20261006 --json $D/counts.json > $D/counts.stdout.md

# timing run
target/release/isolated_bench run --wait --cpus 3 --settle 5 --max-other-cpu 0.5 \
  --label cryptopro-b-chain-sweep --out $D/isolated_bench_record.jsonl -- \
  target/release/cryptopro-b-chain-sweep \
  --chains $D/frozen_inputs/chains.constants.json --expect-sha256 $SHA \
  --pairs 256 --rounds 20 --seed 20261006 --counts $D/counts.json \
  --json $D/sweep.json > $D/sweep.stdout.md
```

The defaults of `cryptopro-b-chain-sweep` are the arms above (`--widths
3,4,5,6,7`, grid chain = the least `optimised` model count, `--element-config
optimised/jacobian/5`, `--reference jacobian/5`, `--continuity-chain
4+1w/5.5.7`, `--continuity-config generic/affine/5`).

## Scope

One curve, one hardware class (the Linux x86-64 VM the run happens on), one
field implementation (the crate's generic `U256::mont_mul`, no assembly), a
variable-time harness.  Not measured: Arm64, other x86-64 hosts, other field
arithmetic, constant-time implementations.
