# B5a: extension fields `GF(p^k)`, at F0

**Declared 2026-10-06, before any B5a measurement.** This is Track B in
`research/notes/index-calculus/IC_TOOL_PROGRAM.md` (§9, B5).
- **Its design:** [`../../design/extension-fields.md`](../../design/extension-fields.md).
- **Its cases:** [`../../conformance/v2-b5a/`](../../conformance/v2-b5a/cases.json).
- **Its instances:** [`instances.json`](instances.json).

**This declaration replaces #1178's,** which was closed before anything
ran: its generators and its harness were Python, which AGENTS.md no
longer admits in the programme's tooling. What it declared is unchanged
except where this text says so:
- **the tools are native** (`icprog`; below);
- **B5a is measured on main's head,** as every Track B step now is
  ([`../../track-b/README.md`](../../track-b/README.md));
- **measurement 6's reductions** are stated, which #1178 left open;
- **C118's row** names `G1`, the instance its frozen case uses; #1178's
  table said `G2`.

B5a's code was written after #1178 and checked on its arm (below). No
B5a measurement has run. Nothing below changes after the first
measurement, except by a dated amendment appended at the end.

## What B5a adds

**Every valid instance over `GF(p^k)` gets a route**, where today every one
is refused as `no-pipeline-for-field` (design §2):
- **`rho-negation` on `GF(p^k)`, for `q ≤ 2^62` and `r < 2^63`:** the
  matched rho's walk over the document's own field.
- **`rho-bignum` on `GF(p^k)`:** any `p` and any `k`.
- **`ic-gaudry-cubic`:** Gaudry's index calculus on `E(GF(p³))`, imported
  from the residual-walk thread's `gaudry_cubic`. It admits:
  - the modulus `t³ − c`;
  - `h = 1`;
  - `r < 2^63`.

  Elsewhere it is refused with a gate code that names the condition.
- **ICV1's extension part** (`docs/curves/ICV1.md`):
  - the specification, and the reference implementation
    (`scripts/curve_id.py extension`), in this declaration, with the
    nine curves it names registered (`docs/curves/registry.json`);
  - the Rust port (`curve_id::extension`), with the code. The report's
    `curve_id` names every extension instance.

**B5a changes no binary or prime path.** The suite's rows, and every
earlier step's instances, run exactly what they ran.

**The plan's B5 is split.** B5a is extension fields. Prime fields past one
word are B5b, declared on their own later. B5's row in the plan is done
when both are accepted.

## Class

**Robustness.** B5a adds instances the tool can solve. It claims no
speedup. Its figures are new rows, not gains.

## Arms

- **The base:** main's head, `995ea207` (v3, if R07 accepts it), with
  Track B's chain through B4: arm `tbarm-B4`, `d5761c13`. No B5a run
  starts before R07's decision and B4's results.
- **B5a's arm:** `tbarm-B5a`, `97ade3d0`. It is one merge of B5a's code
  (`b5a-int`, `538afdce`) into arm B4, by the chain's rules
  ([`../../track-b/README.md`](../../track-b/README.md)), with two
  resolutions by hand:
  - `src/cryptanalysis/mod.rs`: the arm's `koblitz_two_word` kept beside
    B5a's modules;
  - `src/cryptanalysis/gaudry_cubic.rs`: main had added
    `Fp3::with_cube_nonresidue`, the field B5a's `Fp3::with_c` builds,
    under the same conditions, so `with_c` returns its result.

  The arm is on record in
  [`../../track-b/stack-20261006-b5a.bundle`](../../track-b/stack-20261006-b5a.bundle).
- **The arm's development check** (not a measurement;
  [`../../track-b/arms-check-20261006-b5a/`](../../track-b/arms-check-20261006-b5a/)):
  - its `ic` built with its commit embedded;
  - 259 library tests, 50 of `ic`'s and 8 of `tests/curve_id.rs` pass;
  - the conformance suite through B5a passes 110 of 110 cases.

## Measurements

1. **Tests.** `cargo test --release --bin ic`, `--test curve_id`, and the
   library tests of every module B5a touches. They must include design
   §5's tests:
   - the field against `Fpk`;
   - the group law and the keys;
   - both rho pipelines on small instances;
   - sameness with `gaudry_cubic` on its own instances;
   - `curve_id::extension` against the reference's vectors.
2. **Conformance.** `icprog conformance --steps <accepted>,B5a`, with
   B5a's cases (below). The same runner also runs on the base, and the
   base's failures are recorded as what the base cannot do.
3. **The pin, untimed.** All 90 suite rows from their v1 files
   (`icprog bround pin`). Every output must equal the base's.
4. **No slowdown, timed.** By the chain's rule (`icprog bround chain`),
   with B5a the next arm after B4; this declaration appends B5a to the
   chain's order.
5. **F0 on extension fields.** `icprog f0 run --set b5a`, then
   `icprog f0 analyse --set b5a`. Each run is one process, isolated, on
   one core, with an hour (`timeout 3600`). A contended or failed attempt
   is kept and the run made again, at most twice. Each instance below
   (design §4) is run on two targets: its known-answer target and its
   public point `T001`. `H1` is C108's refusal, and is not run here.
   - `G1`–`G3`: `ic price`, paired, `ic-gaudry-cubic` against
     `rho-negation`.
   - `E2`, `E5` and `E11`: `solve: rho` with `rho-negation`.
   - C050's document, on its own known-answer target only:
     `rho-negation`.
   - `B2`: `solve: rho` with `rho-bignum`.

   Only the paired documents pass `--repeats 1 --repeats-fast 1`, as
   #1178's runner did. Recorded for each run:
   - the logarithm, checked in the run and replayed outside it, in
     `icprog`'s own arithmetic for `GF(p^k)`, every arm's certificate;
   - each arm's `S`;
   - every phase's cost and counts;
   - the host manifest and the isolation record.

   A run passes when its report is `complete`, every arm's scalar
   replays, the arms agree, and a known logarithm equals them.
6. **Calibration of the estimates.** `icprog b5a calibrate`, then
   `icprog b5a estimates`. Five isolated runs each, every instance once
   a round, on its known-answer document:
   - `rho-negation` on `E2`, `E5` and `E11`;
   - `rho-bignum` on `B2`;
   - `ic-gaudry-cubic` (paired with `rho-negation`) on `G1`–`G3`.

   Each instance's figures are its medians over its runs: rho's
   nanoseconds a step (online time over steps), and the index calculus's
   operations (`counts.total_ops`) and seconds (set-up plus online). From
   them:
   - **`rho-negation`'s step cost** is the median of `E2`'s, `E5`'s and
     `E11`'s;
   - **`rho-bignum`'s** is `B2`'s, at the width of `B2`'s field (71 bits);
   - **`ic-gaudry-cubic`'s exponent** is the least-squares slope of the
     log of the operations on the log of `r` over `G1`–`G3`, and its
     anchor is `G3`'s, the largest rung's: `p`, `log₂ r`, the operations,
     the seconds and the nanoseconds an operation.

   They replace the provisional constants in
   [`estimates.json`](estimates.json) in B5a's results pull request. No
   case's expectation depends on them.

## B5a's cases

`conformance/v2-b5a/` is written by `icprog b5a cases` from
`instances.json`, which `icprog b5a instances` writes. Their arithmetic
is `icprog`'s own (`src/bin/icprog/b5a.rs`) and shares nothing with the
tool:
- **Group orders:** found exactly, by baby steps and giant steps over the
  Hasse interval, narrowed by labelled points' orders until one is left,
  as B2's generator does.
- **Subgroups:** factored by trial division and Brent's method, each
  prime certified by B1's exact test.
- **Instances:** every property re-checked before a document is written
  (design §4.5's method 3).
- **ICV1 slugs:** computed as the reference computes them; the registry,
  which the reference builds, holds every one.

Both commands refuse to overwrite their files. With `--check`, each
derives its files again and compares the bytes, and an `icprog` test
writes the cases again on every run. They reproduce #1178's Python
generators: every instance record and all 18 parameter files byte for
byte, and the same 16 cases; only the generator each file names differs.

| case | what it checks |
|:--|:--|
| C103 | C050's successor (it `supersedes` C050, whose `until` is B5): C050's document, `GF(1009³)` with a general modulus and `h = 1524`, recovers its known logarithm with `rho-negation`, verified |
| C104–C106 | `ic-gaudry-cubic` paired with `rho-negation` at F0 on `G1`–`G3`. Both arms recover the known logarithm and are verified, and `curve_id.slug` is the registry's |
| C107 | C050's document under `paired`: `no-ic-route`, with `ic-gaudry-cubic` refused as `modulus-not-binomial` |
| C108 | `H1` under `paired`: `no-ic-route`, with `cofactor-not-one` |
| C109 | `E2` under `paired`: `no-ic-route`, with `extension-degree-not-three` |
| C110 | `E2` under `solve: rho`: `rho-negation` at `q ≈ 2^62`, verified |
| C111 | `E5` under `solve: rho`: verified |
| C112 | `E11` under `solve: rho`: verified, with 11 coefficients |
| C113 | `B2` under `solve: rho`: `rho-negation` refused as `field-wider-than-one-word` and `rho-bignum` admitted, verified |
| C114 | `G1`'s curve in `general_weierstrass` form: converted, the conversion recorded, and the logarithm verified by rho |
| C115 | `G1` under `check`: exit 0, with `ic-gaudry-cubic` and `rho-negation` admitted |
| C116 | `G1` at `fidelity: F2`: `status: estimated`, with an estimate for each arm |
| C117 | `ic-gaudry-cubic` named on a binary document: refused as `no-pipeline-for-field` |
| C118 | `G1` with `recipe: {"oracle": "groebner"}`: the logarithm recovered and verified with Gaudry's `S₄` solve |

C050 names B5 as its `until`, and B5 is now two steps, so C103 retires it
by the `supersedes` rule. The runner's list of steps gains `B5a` and
`B5b`; this declaration makes that change to `icprog conformance`.

## Acceptance

B5a is **accepted** when:
- every test passes, sameness included;
- every case the steps select passes;
- the pin holds on every row;
- no size regresses beyond its A/A band;
- every run in measurement 5 recovers its logarithm, verified in the run
  and by the replay.

It is **stopped** if a binary or prime output differs. It is **rejected**
if a test or a case fails, if a size regresses, or if a run in
measurement 5 gives a wrong answer, or none within its hour.

**What an accepted B5a may claim:** the tool solves the discrete
logarithm on curves over `GF(p^k)` at F0 where the subgroup allows:
- with rho, for any `p` and any `k`;
- with Gaudry's index calculus, on `E(GF(p³))` in the module's basis and
  of prime order;

with every answer verified. It claims nothing about the index calculus's
cost against rho beyond the rows it measures, which the residual-walk
thread priced (design §2.4).

## What B5a does not do

- Prime fields past one word (B5b).
- A basis change for other cubic moduli, `p ≡ 2 (mod 3)`, or cofactors
  in the index calculus (design §6).
- `k = 4` (`gaudry_quartic.rs`), or any other degree, for the index
  calculus.
- Characteristic 3.

## Cost, as an estimate

- **Measurement 5:** fifteen runs, each under its hour; the arm's
  development check ran the same documents' cases in minutes.
- **Measurement 6:** thirty-five runs, about half an hour, most of it
  `G3`'s index calculus and `B2`'s `rho-bignum`.
- **The chain step and the pin:** the chain's, as for every arm.

## Commands

    icprog b5a instances --check
    icprog b5a cases --check
    icprog conformance --ic <B5a's ic> --steps B0,B1,B3,B2,B2b,B7a,B3b,B4,B5a --build-commit <commit> --out <file>
    icprog f0 run --set b5a --runs <tree> --ic <B5a's ic>
    icprog f0 analyse --set b5a --runs <tree>
    icprog b5a calibrate --runs <tree> --ic <B5a's ic>
    icprog b5a estimates --runs <tree>
