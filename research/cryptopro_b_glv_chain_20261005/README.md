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

**Success condition** (declared in advance): measured ratio baseline/GLV
> 1 for the `5·5·7` arm, with the median and the minimum of both arms on
the same side (both ratio columns above 1), and every correctness check
passing.  **Stop condition**: one benchmark run of the stated size.

## Frozen inputs

| file | sha256 | content |
|:--|:--|:--|
| `frozen_inputs/cryptopro_b_chain_5_5_7.constants.json` | `3ad2ca805b7181cc8c183852d42b3bbbf0b4fa7bfe9a098f1173b2f3eeb59ade` | the degree-175 chain `5·5·7` (element `4 + ω`, matched on points as `-5 + ω`): curve, kernel polynomials `ψ`, numerators `N`, `M`, codomain coefficients per step, `iso_u`, `λ`, the LLL-reduced GLV basis, the Babai bound, four test vectors |
| `frozen_inputs/cryptopro_b_chain_5_31.constants.json` | `de3ea720d37cdb3c46e4560d848b10be49496de715bf983ea5cc3a9a932a5272` | the degree-155 chain `5·31` (element `ω`, matched as `-ω`), same fields |
| `frozen_inputs/SHA256SUMS` | — | the two hashes above, `sha256sum -c` form |

The constants were generated and verified outside this repository (the
autoresearcher's `export_chain_constants.py` / `emit_rust_constants.py`;
Python provenance, kept as frozen input per `AGENTS.md`).  The Rust module
`src/ecc/cryptopro_b_chain_consts.rs` was emitted from those two JSON files
and copied here **unchanged except for one added line**,
`#![cfg_attr(rustfmt, rustfmt::skip)]`, so that the generated layout passes
the PR-scoped rustfmt check without being reformatted; its header records
both hashes.  Every field element in it is a canonical little-endian
`[u64; 4]` (not Montgomery form) and polynomials are listed low-degree
first.  The curve parameters come from the existing
`CurveParams::gost_cryptopro_b()` in `src/ecc/curve_zoo.rs`; the tests check
the frozen `p, a, b, n` against it.

Each JSON names its chain by an element but also records the element the
composite was matched to on points (`matched_element`: `(-5, 1)` and
`(0, -1)`).  The tests confirm that `λ` is the eigenvalue of the *matched*
element — `λ² + 9λ + 175 ≡ 0` and `λ² + λ + 155 ≡ 0 (mod n)` respectively,
and `λ₁₇₅ + λ₁₅₅ + 5 ≡ 0` — which is `4 + ω` up to sign and conjugation,
same norm, same cost.

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

## Results

One run, as declared.  Unit: **nanoseconds per scalar multiplication** (per
operation for the stage-diagnostic rows), `w = 5`, `M = 256`, `R = 20`,
seed `20261005`.  The table below is the bench's own output
([`bench_w5.stdout.md`](bench_w5.stdout.md), verbatim rows; the full record
with every per-round value is [`bench_w5.json`](bench_w5.json)).

| arm | class | ns per scalar multiplication, median of 20 rounds | min of 20 rounds | ratio baseline/arm (median) | ratio baseline/arm (min) | correctness |
|:--|:--|--:|--:|--:|--:|:--|
| baseline: width-5 NAF (affine odd-multiple table, mixed additions) | reference | 111158 | 107220 | 1.000 | 1.000 | baseline == textbook BigUint reference on all 256 pairs: yes |
| GLV-2 total via chain 5·5·7 (element 4 + 1 w, degree 175): decomposition + phi(P) + interleaved width-5 NAF | engineering | 88759 | 80199 | 1.252 | 1.337 | GLV == baseline on all 256 pairs: yes |
| stage diagnostic: decomposition alone, basis of chain 5·5·7 (num-bigint Babai rounding) | stage diagnostic | 1006 | 875 | — (stage) | — (stage) | k1 + k2·λ ≡ k (mod n) and |k1|,|k2| ≤ 128 bits on all 256 scalars: yes |
| stage diagnostic: phi(P) alone, chain 5·5·7 (element 4 + 1 w, degree 175; projective Horner, no inversion) | stage diagnostic | 4269 | 3841 | — (stage) | — (stage) | phi(P) == λ·P on all 256 points: yes |
| GLV-2 total via chain 5·31 (element 0 + 1 w, degree 155): decomposition + phi(P) + interleaved width-5 NAF | engineering | 91799 | 84544 | 1.211 | 1.268 | GLV == baseline on all 256 pairs: yes |
| stage diagnostic: decomposition alone, basis of chain 5·31 (num-bigint Babai rounding) | stage diagnostic | 1007 | 863 | — (stage) | — (stage) | k1 + k2·λ ≡ k (mod n) and |k1|,|k2| ≤ 128 bits on all 256 scalars: yes |
| stage diagnostic: phi(P) alone, chain 5·31 (element 0 + 1 w, degree 155; projective Horner, no inversion) | stage diagnostic | 9510 | 8479 | — (stage) | — (stage) | phi(P) == λ·P on all 256 points: yes |
| A/A control: the baseline arm timed again as a separate arm (identical code and inputs; its ratio to the baseline is the noise floor) | A/A control | 118769 | 108890 | 0.936 | 0.985 | same routine as the baseline, so baseline == textbook BigUint reference on all 256 pairs: yes |

All correctness checks passed (the bench prints `all correctness checks
passed: yes` and exits 0).

**Verdict against the declared success condition: met.**  Baseline/GLV for
the `5·5·7` chain is 1.252 on the medians and 1.337 on the minima, both
above 1, every check passing.  The A/A control puts the noise floor of a
median ratio at about 6 % (0.936) and of a minimum ratio at about 1.5 %
(0.985); the GLV ratios clear both.  Per-round spread (max/min − 1 over the
20 rounds, from the JSON): baseline 14.7 %, GLV `5·5·7` 18.7 %, GLV `5·31`
19.9 %, A/A 15.0 %, `φ(P)` rows about 20 %, the sub-microsecond
decomposition rows 41–55 %; the medians are the robust statistic, the minima
the least-contended one.

**Against the modelled 1.45×: short of it.**  Measured 1.25–1.34 against a
modelled 1.45.  A formula-derived count for *this* implementation (derived
from the formulas and the expected `256/(w+1) ≈ 43` non-zero wNAF digits,
counting `S` as `M` because `mont_sqr` is `mont_mul` here; not a
measurement) puts the baseline at roughly 2 940 M-equivalents (table 120,
batch conversion with its 264-M Fermat inversion 318, 255 doublings 2 040,
≈ 42 mixed additions 462) and the `5·5·7` GLV arm at roughly 2 240 plus the
measured 1.0 µs decomposition (chain ≈ 128, two tables 240, batch of
sixteen 374, 128 doublings 1 024, ≈ 43 mixed additions 473), a ratio near
1.30, which is where the measurement sits.  The implied 37.8 ns per
M-equivalent also predicts `φ(P)` at ≈ 4.8 µs for `5·5·7` (128
M-equivalents) and ≈ 10.3 µs for `5·31` (≈ 272, the degree-31 step's
polynomials have degrees 15, 31 and 45), against 4.27 µs and 9.51 µs
measured.  What compresses the modelled 1.45 is therefore accounted for by
costs the model prices differently or not at all: the inversion both arms
pay (≈ 9 % of the baseline), the second table and larger batch the GLV arm
pays, and a projective chain at ≈ 128 M against the model's 53 M in affine
Kohel form.  None of this is a measurement of the model; it is the
reconciliation a reader would otherwise have to do.

**Secondary expectation (`5·5·7` beats `5·31`): consistent, not
separately established.**  The `φ(P)` stage rows measure the chains directly:
4.27 µs against 9.51 µs (medians), the degree-155 chain costing 2.2× more to
evaluate, as the model's ordering says.  The two *total* arms differ by
3.0 µs on the medians and 4.3 µs on the minima, inside the ≈ 6 % A/A band
on 90 µs, so the totals alone do not separate the two chains; the stage
rows do.

The decomposition stage is about 1.0 µs, ≈ 1 % of the GLV total, even with
`num-bigint`; the precomputed-constant rounding trick would not change the
verdict.

## Class

**engineering** (`AGENTS.md` §3): same problem, same curve, same
arithmetic; a constant factor on `k·P`.  No boundary moves, nothing is
exponent-moving, and this is a constructive measurement, not an attack, so
there is no `S` row and no scoreboard entry to update.

## Host, build and isolation

* Host: `Intel(R) Xeon(R) Processor @ 2.10GHz`, a cloud virtual machine with
  4 logical CPUs (`cpu_flags_of_interest`: `pclmulqdq sse4_2 popcnt aes avx
  bmi1 avx2 bmi2 avx512f adx`), kernel `6.18.44-fc-v70`, Linux x86-64.  The
  bench's own host line says `logical cpus: 1` because it ran pinned to
  CPU 3 and `available_parallelism` reports the affinity mask; the
  isolation record says 4.
* Build: `rustc 1.97.0 (2d8144b78 2026-07-07)`, `release` profile, no
  `target-cpu` flags; the crate's generic const-generic `U256::mont_mul`.
* Commit at run time: `e1df564168943f7d03f9462503feb63b218fa392`.  The
  record says `git_dirty: true` only because the run's own three output
  files were untracked at that moment; `git diff HEAD` was empty (no tracked
  file differed from that commit), and the outputs were committed next.
* Isolation (`AGENTS.md` §10): run through the native `isolated_bench run`
  (`src/bin/isolated_bench.rs`), which took the benchmark lock, moved every
  movable thread off CPU 3, pinned the bench there and recorded the
  conditions in [`isolated_bench_record.jsonl`](isolated_bench_record.jsonl):
  `contended: false`, other processes used 0.12 CPU-seconds during the
  2.9 s run (the agent runtime `claude[96]` 0.11 of it), 33 involuntary
  context switches, no major faults, PSI `some avg10` 0.03 before and 0.02
  after.  The first attempt, with the tool's default budget of 0.10 CPUs for
  other processes, was refused (`other processes used 0.80 CPU s in 5.0 s
  (limit 0.50)`), the excess being the agent runtime itself, which cannot be
  stopped from inside the session; the run was admitted with
  `--max-other-cpu 0.5`.  Both the refusal and the admitted run's numbers
  are recorded above.  The VM's host neighbours and CPU frequency are
  outside the tool's reach, which is what the A/A row measures.
* One hardware class only: Linux x86-64 on this VM.  Not measured on Arm64
  (Apple silicon, Graviton), on any other x86-64 host, or with any
  vectorised or assembly field arithmetic.

## Honesty notes

* **Variable time, everywhere.**  Identity checks, the `H = 0` branch in the
  addition law, wNAF digit scans, the `num-bigint` decomposition and the
  textbook reference all branch on their inputs.  The module docs say so.
  This is a measurement harness; nothing here is for secret scalars.
* **Unoptimised pieces.**  The field multiplication is the crate's generic
  const-generic schoolbook `U256::mont_mul` (no dedicated squaring, no
  assembly); the field inversion is Fermat with the sparse exponent `p - 2`
  (256 squarings, 8 multiplications); the decomposition is `num-bigint`; the
  chain has no special case for `Z = 1` in the first step.  Both arms share
  the field, so the ratio is what it is on this arithmetic; a faster field
  would change both arms' absolute numbers and could move the ratio either
  way, since the decomposition is a fixed cost that does not scale with the
  field.
* **Wall time, not operation counts.**  The unit is nanoseconds per scalar
  multiplication because that is what the hypothesis is about; the field
  operation count of each arm is fixed by the formulas above and does not
  depend on the host.  The counts quoted under "Results" are derived, and
  labelled so.
* **Not an ICV1-registered curve.**  The curve is named by its published
  standard name (RFC 4357 CryptoPro-B), which `AGENTS.md` §11 allows; it is
  not in `docs/curves/registry.json`, and registering it is a separate
  change outside this PR's scope.

## Commands

Exactly what was run, in this order, on branch
`cryptopro-b-glv-chain-20261005` at commit `e1df5641`:

```sh
cargo test --release --lib cryptopro_b            # 32 passed, 0 failed
cargo build --release --bin cryptopro-b-glv-bench --bin isolated_bench

# refused by the isolation tool (agent runtime above the 0.10-CPU budget):
target/release/isolated_bench run --wait --cpus 3 --settle 5 \
  --label cryptopro-b-glv-bench-w5 \
  --out research/cryptopro_b_glv_chain_20261005/isolated_bench_record.jsonl -- \
  target/release/cryptopro-b-glv-bench --pairs 256 --rounds 20 --w 5 --seed 20261005 \
  --json research/cryptopro_b_glv_chain_20261005/bench_w5.json

# the run reported above:
target/release/isolated_bench run --wait --cpus 3 --settle 5 --max-other-cpu 0.5 \
  --label cryptopro-b-glv-bench-w5 \
  --out research/cryptopro_b_glv_chain_20261005/isolated_bench_record.jsonl -- \
  target/release/cryptopro-b-glv-bench --pairs 256 --rounds 20 --w 5 --seed 20261005 \
  --json research/cryptopro_b_glv_chain_20261005/bench_w5.json \
  > research/cryptopro_b_glv_chain_20261005/bench_w5.stdout.md
```

`cargo run --release --bin cryptopro-b-glv-bench -- --pairs 256 --rounds 20
--w 5 --seed 20261005 --json out.json` reproduces the same inputs and
checks on any host; the timings are the host's.

## Files

| file | what |
|:--|:--|
| `README.md` | this record |
| `frozen_inputs/*.json`, `frozen_inputs/SHA256SUMS` | the frozen chain constants |
| `bench_w5.json` | the full JSON record of the run (times per round, host, build, correctness) |
| `bench_w5.stdout.md` | the bench's standard output, verbatim |
| `isolated_bench_record.jsonl` | the `isolated_bench` conditions record of the run |
| `src/ecc/cryptopro_b_chain_consts.rs` | the constants module (copied unchanged but for the rustfmt-skip line) |
| `src/ecc/cryptopro_b_field.rs` | the Montgomery-form field |
| `src/ecc/cryptopro_b_point.rs` | Jacobian arithmetic, wNAF, the chains, the decomposition, both arms |
| `src/bin/cryptopro_b_glv_bench.rs` | the benchmark |
