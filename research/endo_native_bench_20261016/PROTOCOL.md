# Native GLV-chain benchmark on Bandersnatch, Tom-384 and CP6-782 — protocol, stated before measuring

This file is committed before the counting runs and the timing runs it
describes.  The results go into `README.md` in a later commit; this file is
not edited after the runs.

## Question

`research/cryptopro_b_chain_sweep_20261006` measured GLV-2 through a cheap
isogeny-chain endomorphism natively on one curve, GOST CryptoPro-B: **1.408×
counted, 1.417× median wall time** against the best width-`w` NAF, within
0.26 % of the chain-sweep model's counts.  The autoresearcher's cross-curve
sweep (aburan28/crypto-autoresearcher, `research/endosweep_curves_20261006`,
finding 5) modelled the same construction on every curve with a small CM
discriminant, and named as the next native measurements the largest modelled
gains on a non-unit curve: **Tom-384 (modelled 1.525), Bandersnatch (1.522, on
its short Weierstrass model) and CP6-782 G1 (1.494)**.  This experiment
measures those three natively, with the same discipline, the same arms and
the same formulas the model counts, so that the model is tested on curves it
was not calibrated on: other field sizes (255, 384 and 782 bits), other
doubling formulas (general `a`, small `a`), a degree-2 step (Bandersnatch's
`√−2`), and a large cofactor (CP6-782, `n` of 377 bits in a 782-bit field).

It is constructive (`k·P` on the same curve, faster): class **engineering**.
It is not an attack and makes no security claim.  The harness is variable
time.

## The harness

Committed in `85793758e` (the commit before this file's).  New modules beside the CryptoPro-B ones, which are not changed:

* `src/ecc/endo_generic_field.rs` — `GenericField<L>`: a Montgomery field of
  `L` 64-bit limbs whose modulus is a run-time value (4 limbs for the 255-bit
  field, 6 for 384 bits, 13 for 782 bits) over the crate's `Uint<L>`
  `mont_mul`/`mont_sqr` (no dedicated squaring, no assembly); Fermat
  inversion `bits(p−2)` S + `popcount(p−2)` M; the `endo-chain-opcount`
  counter (same accounting as `cryptopro-b-opcount`).
* `src/ecc/endo_generic_point.rs` — Jacobian arithmetic with the curve's own
  doubling, chosen as the model chooses it (`curvesweep.arithmetic`): `a = 0`
  dbl-2009-l 2M + 5S; `a = −3` dbl-2001-b 3M + 5S; small `|a| ≤ 16`
  dbl-2007-bl with `a·ZZ²` by additions 1M + 8S; general `a` dbl-2007-bl
  1M + 8S + 1M; add-2007-bl 11M + 5S, madd-2007-bl 7M + 4S; the two chain
  evaluators of the CryptoPro-B harness generalised, plus Vélu 2-steps;
  width-`w` NAF, affine or Jacobian tables; Babai decomposition.
* `src/ecc/endo_generic_chain_set.rs` — the loader of the
  `endosweep-chainsweep/1` files with the curve's own constants (checks listed
  in the module docs: probable-prime `p` and `n`, cofactor inside Hasse's
  interval, `ω`'s root, norms, eigenvalues `= a' + b'·ω_root`, canonical
  field elements, composing steps, monic polynomials, Vélu 2-steps, basis of
  determinant `±n` in `λ`'s lattice, generator on the curve and of order `n`
  by the textbook arithmetic).
* `src/bin/endo_chain_bench.rs` — the benchmark (`endo-chain-bench`): the arm
  set, pair generation, correctness gate, counting protocol and timing loop of
  `cryptopro-b-chain-sweep`, instantiated at 4, 6 or 13 limbs.

### Positive control (already run, before this protocol)

The generic harness was run in counting mode on the frozen CryptoPro-B chain
set of `research/cryptopro_b_chain_sweep_20261006` with that experiment's
inputs (seed 20261006, M = 256, the crate's CryptoPro-B generator, its
continuity arm), and its record
(`positive_control/cryptopro_b.counts.json`, `.stdout.md`) was compared with
that experiment's committed `counts.json`: **all 97 arms, same keys in the
same order, identical `mul`, `sqr`, `inv`, `add`, `m_eq`, `m_eq_min` and
`m_eq_max`** (0 differences), every correctness check passing.  The
generic field, formulas, evaluators, recoding, tables and decomposition
therefore reproduce the specialised harness's counts exactly on the curve
where those were confirmed.  The unit tests also check the generic path
against the crate's textbook arithmetic on CryptoPro-B, on a general-`a`
384-bit curve, on secp256k1 (`a = 0`), on a small-`a` 384-bit curve, the
2-step formulas against affine Vélu, and every test vector and `φ = [λ]` of
the three frozen sets below at their own limb counts.

## Frozen inputs

Produced by the autoresearcher (aburan28/crypto-autoresearcher, branch
`endosweep-native-chains-20261016`, commit `b2f0ac3095`,
`research/endosweep_native_chains_20261016/`, generator
`harness/endosweep/chainsweep.py` / `curvesweep.py`; Python provenance, which
stays there), handed over through the session scratchpad, and byte-identical
to that commit's `SHA256SUMS`:

| file | sha256 | content |
|:--|:--|:--|
| `frozen_inputs/bandersnatch.chains.constants.json` | `12100c87b6ffe5a3469f3c8ab625f2e8f23411df025a01bf3c06136cd646114b` | 2 chains (2 elements): `0+1w/2` (`√−2`, one 2-isogeny), `1+1w/3` |
| `frozen_inputs/bandersnatch.model.json` | `2ca638388eff39906d29cc217f956759e3e3504ecc8d50cadfbcc2cea2b54a52` | the model's expected counts of every arm |
| `frozen_inputs/tom384.chains.constants.json` | `c93b4b4752a3bad0331d6bdc6a7056c7a8608df16db626cb968c4e013a7ed00f` | 22 chains (11 elements), `D_K = −619` |
| `frozen_inputs/tom384.model.json` | `5f8a019c35ee22371b3bfa70720eb034495f0298d560ceb1ecb43ea7236fb253` | the model |
| `frozen_inputs/cp6_782.chains.constants.json` | `9e08fe6a45643b75c513c75164ecb51d39393d65c078dc07dadd9461bfa8ac90` | 57 chains (18 elements), `D_K = −339` |
| `frozen_inputs/cp6_782.model.json` | `eebecc832b78edbafde5e9485841d61168f9c652106df27b8f2117b0cca828c3` | the model |
| `frozen_inputs/SHA256SUMS` | — | the six hashes, `sha256sum -c` form |

The constants the loader uses are certified in Lean by the autoresearcher
(`formal/CryptoResearch/Endosweep/NativeChains.lean`:
`bandersnatch_native_chain_certificate`, `tom384_native_chain_certificate`,
`cp6_782_native_chain_certificate`; primality of `p` and `n` is a stated
hypothesis there, not a Lean proof).  This harness re-checks them natively at
load time and refuses to run unless each file's sha256 is the one above
(`--expect-sha256`).

The curves, as the files state them and the loader confirms:

| curve | p | n | cofactor | `a` and doubling | `D_K` | model |
|:--|:--|:--|:--|:--|:--|:--|
| Bandersnatch | 255 bits (BLS12-381 scalar field) | 253 bits | 4 | general `a`: dbl-2007-bl 1M + 8S + 1M | −8 | **its short Weierstrass model**; the deployed twisted-Edwards arithmetic is **not** measured |
| Tom-384 | 384 bits | 384 bits | 1 | general `a`: 1M + 8S + 1M | −619 | short Weierstrass |
| CP6-782 G1 | 782 bits | 377 bits | 405 bits | `a = 5`: dbl-2007-bl 1M + 8S (`5·ZZ²` by additions) | −339 | short Weierstrass; every input point is `r·G` with `G` of order `n` (checked), so in the order-`n` subgroup |

## Arms (per curve)

All arms of one curve run in one binary on the same arithmetic.  As in the
CryptoPro-B sweep:

| group | arms | what |
|:--|--:|:--|
| baseline | 10 | width-`w` NAF, `w` ∈ {3, …, 7}, table `affine` or `jacobian` |
| GLV grid | 20 | the grid chain × evaluator `generic`/`optimised` × table × `w` |
| GLV per element | 1 per element | every element's cheapest ordering (by the model's optimised count, ties by id) at `optimised/jacobian/w5` (the grid chain's entry is shared with the grid) |
| stage diagnostic: `φ(P)` | 2 per chain | every chain in the set × both evaluators |
| stage diagnostic: decomposition | 1 per element | the reduced basis of every element |
| A/A control | 1 | the reference baseline run again |

No continuity arm (it was PR #1408's configuration on CryptoPro-B).

| curve | reference baseline (the model's cheapest) | grid chain (the model's cheapest optimised chain, ties by id) | per-element arms | `φ(P)` arms | total arms |
|:--|:--|:--|--:|--:|--:|
| Bandersnatch | **`jacobian/w5`** | `0+1w/2` (`√−2`, 2-isogeny, 6 M_eq; ties `1+1w/3` at 6) | 2 | 4 | 10 + 20 + 1 + 4 + 2 + 1 = 38 |
| Tom-384 | **`jacobian/w6`** | `4+1w/7.5.5` (80 M_eq) | 11 | 44 | 10 + 20 + 10 + 44 + 11 + 1 = 96 |
| CP6-782 G1 | **`jacobian/w6`** | `4+1w/7.3.5` (65 M_eq; ties `4+1w/7.5.3`) | 18 | 114 | 10 + 20 + 17 + 114 + 18 + 1 = 180 |

The **declared GLV configuration** of each curve is its grid chain at
`optimised/jacobian/w5` (the model's cheapest GLV arm on Tom-384 and CP6-782;
on Bandersnatch the model's cheapest is `1+1w/3` at the same configuration,
0.4 M_eq below `0+1w/2`, and both are measured).

The evaluators are those of the CryptoPro-B sweep (`generic`: projective
steps, crypto#1408's; `optimised`: Jacobian steps, affine first step,
isomorphism `Z·u⁻¹`), plus, for the 2-step: `optimised` Jacobian Vélu
`X' = (X·d + v·Z⁴)·d`, `Y' = Y·(d² − v·Z⁴)·d`, `Z' = Z·d`, `d = X − x₀Z²`
(4M + 1S on the affine input), and `generic` projective Vélu
`X' = (X·d + v·Z²)·d`, `Y' = Y·(d² − v·Z²)`, `Z' = Z·d²`, `d = X − x₀Z`
(3M + 1S on the affine input).

## Units and accounting

* **Primary: counted field operations** per scalar multiplication (per call
  for stage diagnostics), mean over the `M` pairs, from a counting build of
  the same binary (`--features endo-chain-opcount`) on the same inputs.
  `M_eq = M + S` (`mont_sqr` is `mont_mul` here); an inversion is counted
  through its squarings and multiplications (Bandersnatch 255 S + 164 M,
  Tom-384 384 S + 300 M, CP6-782 782 S + 391 M: the model's `inv`).
  Deterministic and host-independent.
* **Secondary: nanoseconds** per scalar multiplication, median and minimum
  over `R` rounds, interleaved arms (order reversed every other round), one
  pinned core, from the default build.
* Inside every GLV arm: decomposition, `φ(P)`, both tables, the batched
  inversion when affine, both recodings, the main loop.  Outside every arm:
  the per-chain context (Montgomery constants, `u⁻¹`, `λ`, basis) and the
  final affine conversion.
* Wall time per counted operation differs between curves (limb count); the
  ratios are within one curve only.  No claim about CPU cycles per field
  operation is made.

## Predictions (from the frozen model files; derived, not measured)

The model's expectations are over 400 random scalars (seed 7); the runs use
256 pairs (seed 20261016), so totals are compared within a tolerance.

### Chain counts (`φ(P)` stage arms)

Every `φ(P)` arm's counted `M` and `S` equal the chain set's `model_ops`
exactly, **except one arm declared here in advance**: Bandersnatch
`phi/0+1w/2/generic`.  The model charges a 2-step 4M + 1S under both
evaluators and then the generic tail (4M + 1S), 8M + 2S = 10 M_eq; but the
generic evaluator's state is projective, so its correct 2-step on the affine
input is the projective Vélu step above, 3M + 1S, and the harness counts
**7M + 2S = 9 M_eq** (it does not run a Jacobian step into a projective tail,
which would compute a wrong point).  Consequently every Bandersnatch GLV arm
on `0+1w/2` with the `generic` evaluator is predicted at the model's value
minus exactly 1 M_eq.  The `optimised` 2-step (5M + 1S = 6) is predicted
exact.

### Totals (expected counted M_eq per `k·P`; model)

Bandersnatch (reference `jacobian/w5`):

| baseline | model M_eq | GLV on `0+1w/2` | model M_eq |
|:--|--:|:--|--:|
| jacobian/w5 (reference) | 3 285.3 | optimised/jacobian/w5 (declared) | 2 157.6 |
| jacobian/w6 | 3 317.3 | generic/jacobian/w5 | 2 161.6 (harness: 2 160.6) |
| jacobian/w4 | 3 363.1 | optimised/jacobian/w4 | 2 167.3 |
| jacobian/w7 | 3 495.1 | optimised/jacobian/w3 | 2 305.1 |
| jacobian/w3 | 3 536.7 | optimised/jacobian/w6 | 2 314.5 |
| affine/w5 | 3 551.4 | optimised/affine/w4 | 2 389.4 |
| affine/w4 | 3 558.4 | optimised/affine/w3 | 2 437.2 |
| affine/w3 | 3 655.3 | optimised/affine/w5 | 2 477.7 |
| affine/w6 | 3 668.2 | optimised/jacobian/w7 | 2 755.0 |
| affine/w7 | 3 980.7 | optimised/affine/w7 | 3 463.0 |

Per element at optimised/jacobian/w5: `0+1w/2` 2 157.6, `1+1w/3` 2 157.2
(chains 5M + 1S and 6M + 0S).  **Model ratio, cheapest baseline over cheapest
GLV: 3 285.3 / 2 157.2 = 1.523** (sweep README: 1.522).

Tom-384 (reference `jacobian/w6`):

| baseline | model M_eq | GLV on `4+1w/7.5.5` | model M_eq |
|:--|--:|:--|--:|
| jacobian/w6 (reference) | 4 930.5 | optimised/jacobian/w5 (declared) | 3 234.6 |
| jacobian/w5 | 4 946.0 | generic/jacobian/w5 | 3 282.6 |
| jacobian/w7 | 5 069.3 | optimised/jacobian/w4 | 3 313.7 |
| jacobian/w4 | 5 091.8 | optimised/jacobian/w6 | 3 341.3 |
| affine/w5 | 5 367.9 | optimised/jacobian/w3 | 3 557.2 |
| jacobian/w3 | 5 371.3 | optimised/affine/w4 | 3 670.8 |
| affine/w4 | 5 421.7 | optimised/affine/w5 | 3 711.3 |
| affine/w6 | 5 452.7 | optimised/jacobian/w7 | 3 743.5 |
| affine/w3 | 5 591.5 | optimised/affine/w3 | 3 791.2 |
| affine/w7 | 5 738.3 | optimised/affine/w7 | 4 635.4 |

Per element at optimised/jacobian/w5 (chain M + S): `4+1w/7.5.5` 78 + 2 →
3 234.6; `9+1w/7.5.7` 93 + 2 → 3 251.6; `15+2w/7.5.5.5` 108 + 3 → 3 266.4;
`2+1w/23.7` 111 + 1 → 3 270.0; `0+1w/31.5` 120 + 1 → 3 278.7;
`20+1w/23.5.5` 126 + 2 → 3 284.7; `54+1w/5.5.5.5.5` 132 + 4 → 3 289.6;
`5+1w/37.5` 138 + 1 → 3 292.7; `25+1w/23.5.7` 141 + 2 → 3 300.8;
`39+1w/7.5.7.7` 138 + 3 → 3 312.6; `37+3w/23.5.5.5` 156 + 3 → 3 317.1.
The 22 chain counts are those of CryptoPro-B (same `D_K`, same set).
**Model ratio: 4 930.5 / 3 234.6 = 1.524** (sweep README: 1.525).

CP6-782 G1 (reference `jacobian/w6`):

| baseline | model M_eq | GLV on `4+1w/7.3.5` | model M_eq |
|:--|--:|:--|--:|
| jacobian/w6 (reference) | 4 463.9 | optimised/jacobian/w5 (declared) | 2 987.6 |
| jacobian/w5 | 4 479.5 | generic/jacobian/w5 | 3 035.6 |
| jacobian/w7 | 4 607.4 | optimised/jacobian/w4 | 3 068.2 |
| jacobian/w4 | 4 620.5 | optimised/jacobian/w6 | 3 098.2 |
| jacobian/w3 | 4 891.6 | optimised/jacobian/w3 | 3 305.3 |
| affine/w5 | 5 396.2 | optimised/jacobian/w7 | 3 498.7 |
| affine/w4 | 5 446.6 | optimised/affine/w4 | 3 917.8 |
| affine/w6 | 5 480.7 | optimised/affine/w5 | 3 957.4 |
| affine/w3 | 5 610.2 | optimised/affine/w3 | 4 033.8 |
| affine/w7 | 5 769.9 | optimised/affine/w7 | 4 883.2 |

Per element at optimised/jacobian/w5, the 18 totals range 2 987.6
(`4+1w/7.3.5`, chain 63M + 2S) to 3 051.5 (`1+2w/7.7.7`); the full table is in
the model file.  **Model ratio: 4 463.9 / 2 987.6 = 1.494** (sweep README:
1.494).

The `generic` grid rows are the `optimised` ones plus the chain difference
(48 M_eq on Tom-384 and CP6-782, 4 on Bandersnatch's model / 3 in the
harness); the omitted generic rows are in the model files and are predictions
too.

## Success conditions (declared now)

1. **Correctness.**  Both runs of every curve pass every check: every
   baseline arm equals the textbook `BigUint` `k·P` on every pair, every GLV
   arm too, every `φ(P)` equals `λ·P` on every point and reproduces the frozen
   test vectors, every decomposition satisfies `k1 + k2·λ ≡ k (mod n)` within
   its Babai bound.  Without this no number is evidence.
2. **The count model is exact where it claims to be.**  Every `φ(P)` arm's
   counted `M` and `S` equal the model's, except the one arm declared above,
   which equals 7M + 2S; every total arm's mean counted M_eq is within 1 % of
   the model's expectation (the Bandersnatch generic `0+1w/2` arms against the
   model minus 1).
3. **Counted headline, per curve.**  The counted ratio (cheapest baseline arm
   over cheapest GLV arm, both picked from the counts) is at least 99 % of the
   model's: **Bandersnatch ≥ 1.508, Tom-384 ≥ 1.509, CP6-782 ≥ 1.479**.  It
   is reported against the sweep's named 1.522 / 1.525 / 1.494 as well.
4. **Wall-time headline, per curve.**  The declared GLV configuration beats
   the fastest baseline arm by median (picked after the run, which favours the
   baseline) by a **median ratio above 1.40** on each curve, and its median
   ratio lies within 5 % of its own counted ratio against that baseline arm
   (`|timed / counted − 1| ≤ 0.05`).  The A/A floor (`1 − A/A median ratio`)
   is reported; a difference smaller than twice the floor is not called a
   difference.
5. **Chain ranking, per curve.**  The lowest `φ(P)` median of the stage arms
   is an `optimised` arm of a chain at the model's least optimised count:
   Bandersnatch `0+1w/2` or `1+1w/3`; Tom-384 `4+1w/7.5.5`; CP6-782
   `4+1w/7.3.5` or `4+1w/7.5.3`.

A condition that fails is reported as failed.  Inadmissible: changing the
arms, inputs, widths, seed, references or accounting after seeing a result;
rerunning to improve a number; dropping an arm.

## Stop conditions

Per curve, one counting run (`M = 256`, seed `20261016`) and one timing run
(`M = 256`, `R = 20`, same seed) through `isolated_bench run` on CPU 3 with
`--max-other-cpu 0.5`, the curves in the order Bandersnatch, Tom-384,
CP6-782, all counting runs first.  Before each timing run a low-priority
background job on the host is paused (the orchestrator's pause file, about
90 s of settling, `uptime`/`top` recorded), and no build or other heavy work
of this session runs during timing.  A run the isolation tool refuses may be
retried; a run it admits is final.  Every attempt is recorded.

## Commands

```sh
cargo test --release --lib endo_generic
cargo test --release --features endo-chain-opcount --lib endo_generic
cargo build --release --bin endo-chain-bench --bin isolated_bench
CARGO_TARGET_DIR=target/opcount cargo build --release \
  --features endo-chain-opcount --bin endo-chain-bench

D=research/endo_native_bench_20261016
F=$D/frozen_inputs
# per curve: C, SHA, REF
#   bandersnatch 12100c87b6ffe5a3469f3c8ab625f2e8f23411df025a01bf3c06136cd646114b jacobian/5
#   tom384       c93b4b4752a3bad0331d6bdc6a7056c7a8608df16db626cb968c4e013a7ed00f jacobian/6
#   cp6_782      9e08fe6a45643b75c513c75164ecb51d39393d65c078dc07dadd9461bfa8ac90 jacobian/6

# counting run (deterministic; no isolation needed)
target/opcount/release/endo-chain-bench --count-only \
  --chains $F/$C.chains.constants.json --expect-sha256 $SHA --model $F/$C.model.json \
  --reference $REF --pairs 256 --seed 20261016 \
  --json $D/$C.counts.json > $D/$C.counts.stdout.md

# timing run
target/release/isolated_bench run --wait --cpus 3 --settle 5 --max-other-cpu 0.5 \
  --label endo-chain-bench-$C --out $D/$C.isolated_bench_record.jsonl -- \
  target/release/endo-chain-bench \
  --chains $F/$C.chains.constants.json --expect-sha256 $SHA --model $F/$C.model.json \
  --reference $REF --pairs 256 --rounds 20 --seed 20261016 --counts $D/$C.counts.json \
  --json $D/$C.timing.json > $D/$C.timing.stdout.md
```

Defaults used: `--widths 3,4,5,6,7`, grid chain = the least optimised model
count (ties by id), `--element-config optimised/jacobian/5`, no continuity
arm.

## Scope

Three curves, one hardware class (the Linux x86-64 VM the runs happen on),
one field implementation (the crate's generic `Uint<L>::mont_mul` with a
run-time modulus, no assembly, no dedicated squaring), short Weierstrass
models only, a variable-time harness, variable-base single scalar
multiplication only.  Not measured: Bandersnatch's twisted-Edwards
arithmetic, constant-time implementations, other hosts, faster inversions
(the table-mode conclusions depend on the Fermat inversion), multi-scalar
multiplication.
