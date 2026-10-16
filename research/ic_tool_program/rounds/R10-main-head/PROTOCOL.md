# R10: the programme re-based on main's head, again

**Declared 2026-10-16, before any R10 timed run.** This is Track A in
`research/notes/index-calculus/IC_TOOL_PROGRAM.md`, under the plan's
standing rule (§8, 2026-10-05): when main has moved under the
programme's baseline, a drift round like R07 runs before the next
declared round. A one-process diagnostic ran first and is disclosed
below. Nothing below changes after the first timed run, except by a
dated amendment appended at the end.

## Why a round, and why now

**The baseline is v3, main at `995ea207` (R07).** Main is now at
`92246530`, about 1,450 commits later, and 44 of them touch the three
files every timed phase runs (`koblitz_fast.rs`,
`koblitz_index_calculus.rs`, `semaev_decomp.rs`). Among them:
- `a8f4bb60`: the batched inversion stops zero-filling its scratch on
  every call;
- `a98ba95a`: the portable path keys Frobenius orbits in lockstep lanes;
- `0cfba411` (#1652): batched wide inversions and a two-word fold
  reduction;
- `8ab924b9`: isogenous (`b ≠ 1`) curves in the Koblitz pipeline;
- the SAT, F4, F6 and root-finding oracles (#1646, #1729 and others),
  which `ic price`'s pipeline does not call.

**The next declared rounds wait on this one.**
- **R06 and R08** are declared on v3 for the reference class
  (AVX-512 with VBMI2, GFNI and VPCLMULQDQ), and have waited for a host
  of it (R06 amendment 1). This container is one: an Emerald Rapids
  (`Intel(R) Xeon(R) Processor @ 2.10GHz`, family 6 model 207, kernel
  build `6.18.44-fc-v139`), where `icprog host-class r06` and `r08` find
  nothing missing.
- R06's amendment 1 says that if the baseline moves first, R06 does not
  run on its declaration: an amendment ports the candidate onto the new
  baseline before any run. The standing rule says the baseline moves
  now. So R10 runs first, and R06 is ported onto its outcome.
- **R09** is declared for the AVX-512 class without VPCLMULQDQ. This host
  is outside it, so R09 cannot run here (`icprog run r09` refuses). Its
  runs on the Cascade Lake host were lost with that container (below).

**What this round is.** It asks whether main's head computes v0's
answers, and what it costs against v3, size by size. It claims no lever
of its own: any gain or loss is main's.

## Hypothesis

- **Main's head (`92246530`) gives v0's outputs on every pinned row:**
  the counts, relations and logarithms; rho's walks and both arms'
  scalars; the verification flags.
- **It costs about what v3 costs**, within a few percent at every size:
  the drift on the timed path is housekeeping (scratch, portable path,
  wide fields), and the suite's fields are one word.

## Class

**Engineering or accounting.** The pin makes the counts identical, so the
ratio to the floor is flat. The ledger credits any change to main, not to
a lever of the programme's.

## Arms

Both arms are built on this host with its toolchain, `rustc 1.97.0`, so
the compiler is not a difference between them.
- **The base: v3, rebuilt.** Commit
  `995ea2071cc7453877a503d30eb561ae82cddab9`, built with
  `cargo build --release --bin ic` and the first `Cargo.lock` main
  tracked (`ba345bbc`, SHA-256 `4ddc5ef76119493b…`).
  - **It is not v3's binary** (`9c332039…`). That binary, and the lock it
    was built with (`4f17b356…`), lived in a scratch directory the
    container replacement erased; no commit carries that lock. A
    rebuild cannot be byte-identical anyway, because the build path is
    embedded in the binary.
  - The pin below shows the rebuild's outputs are v0's, as v3's were.
  - Binary SHA-256 `00cff0b04c3edb9a46076af20857fb27ba1b80f177b277b65abb4c4cf2c3f562`.
- **The candidate: main at `92246530c`.** Built with
  `cargo build --release --locked --bin ic` and main's own committed
  `Cargo.lock`. Binary SHA-256
  `d8daac9aa75ab26e4bac70810f52d8925d518f2d2645ac2f98306d935f23a653`.
- **The manifest** (`icprog run r10 manifest`) records both binaries'
  SHA-256 and the commits they were built from.

## Before any timed step

1. **The tests.** On the candidate's tree:
   `cargo test --release --lib -- koblitz_fast koblitz_index_calculus semaev_decomp`.
2. **The pin, untimed.** Both arms on all 90 suite rows
   (`icprog run r10 pin`, then the same with the base as the candidate).
   Every output must equal v0's from R01's profile pass, and the names
   must agree.
3. **The A/A.** The base against a byte-identical copy, on `M1`'s 22
   rows, five rounds, ABAB, isolated (`icprog run r10 aa`): 220
   processes. The bands are the round's own: the host is new.

## Rows

R07's rows, with fresh holdouts.
- **The suite rows:** `M1`'s 22 rows, and the other six suite rows at
  each of the three target sizes: 40 rows, five rounds, ABAB, isolated,
  400 processes.
- **The fresh holdouts:** eight rows at each target size, drawn by
  `icprog holdouts` (which first reproduced R07's holdouts byte for byte
  on this host):
  - recipe seeds 230, 231, 232 and 233, two targets per seed, `T151` to
    `T158`;
  - `public_hash_seed` 23151 to 23158, and rho seeds `0x230000 + 151` to
    `+ 158`.

  They are committed in [`holdouts/`](holdouts/) (`SHA256SUMS`). No round,
  exploration or diagnostic has run them. Five rounds, ABAB, isolated:
  240 processes.
- **The extension.** A set whose interval half-width exceeds 3% after
  five rounds gets rounds 6–10, pooled.
- **The callgrind cross-check:** as R07's, `M1`'s first target at
  `icv1-f2m61-t158598901-ab42b6c5` and `icv1-f2m53-tm56619371-dac20a85`,
  both arms, untimed. It decides one thing: both arms recover the same
  logarithm.

## Prediction

**From the diagnostic below**, cold time, base over candidate, one
process each:

| curve | `log₂ r` | v3 (rebuilt), ms | main, ms | ratio |
|:--|--:|--:|--:|--:|
| `icv1-f2m53-tm56619371-dac20a85` | 44.3 | 757 | 735 | 1.03 |
| `icv1-f2m59-tm943548413-98844ecc` | 44.5 | 769 | 756 | 1.02 |
| `icv1-f2m61-t158598901-ab42b6c5` | 47.2 | 3089 | 3239 | 0.95 |

- **Cold time, both sets:** within about 5% of 1 at every size. One
  process per arm cannot tell `2^47.2`'s 0.95 from noise; the A/A rule
  decides it.

## Success and stop

**Accepted, as baseline v4** (main's head), if all of the following hold:
1. The tests pass, and every pinned output is identical, for both arms.
2. No size regresses beyond its A/A band: no size whose cold interval's
   upper end lies below the lower end of the round's own A/A interval.
3. Under callgrind, both arms recover the same logarithm at both sizes.

There is **no minimum gain**: the programme moves to the code its users
build, and the figures report what that move did.

**Rejected** if a size regresses beyond its A/A band. Then v3 stays the
baseline, the regression is located and fixed on main by its own pull
request, and the round is re-declared on a head that carries the fix.

**Stopped** if a test fails, an output differs, or a verification fails.
A pin failure means main changed an output: it must be declared as an
algorithm change with its own measurement, or fixed, before the
programme re-bases.

**If accepted:**
- v4 enters `baselines.json` with its binary hash, host and `S`, beside
  v3's rebuilt `S` from the same runs.
- **R06 is ported onto v4** by an amendment, before any R06 run, and runs
  here, on its class. R08 follows R06's decision.

## Inadmissible

- Pooling the diagnostic with R10's runs, or choosing between them.
- Changing either binary, commit or lock after a run.
- Changing the suite, the recipes, the holdouts or the targets, or
  drawing new holdout seeds after a run.
- Pooling contended runs.
- Quoting collection, build or rho as the speedup: the speedup is cold
  time.
- Crediting a gain to the programme: it is main's.

## Cost

- The tests and the pins: a few minutes, and 180 untimed processes.
- The A/A: 220 timed processes. The suite rows and holdouts: 640 timed
  processes, and up to 640 more if extended.
- The callgrind cross-check: four runs under valgrind.

On the last host, each timed process waited for the CPU-pressure gate to
settle after the one before, and the queue ran one to four processes a
minute. At that pace the round takes several hours.

## Order and restarts

- **R10 runs first,** then R06 (ported), then R08.
- **Every step resumes what its run tree holds.** Each resume runs
  `icprog run r10 manifest-resumed`. A resume onto a host whose CPU
  model, flags, cores or memory differ from `host.json`'s stops the
  round; the rest re-runs from the A/A on the new host, as a dated
  amendment.

## What the last container took with it (disclosed)

Between 2026-10-08 and this declaration the container was replaced, not
restarted, and its scratch directory went with it:
- **R09's runs** (declared in #1414): its drift pin had held on all 90
  rows and its own pin had held; its A/A had started. None of it was
  committed, so R09 has no runs on record. It can resume only on a host
  of its class.
- **Track B's campaign** (the chain over B0–B5a, about 906 of 1,100
  processes, and the step-specific measurements): no record was
  committed after the last bundle. It re-runs from the start on the next
  baseline.
- **Two exploratory candidates on R09's,** A1's build streams ported
  (`f487833e`) and a window least-rotation kernel (`3e4fe770`), local
  commits that were never pushed. Neither had a timed result. They are
  not part of R10.

## The diagnostic before this declaration (disclosed)

**When and how.** 2026-10-16, about 01:55 UTC. One process per arm on
`M1`'s first target at three sizes, under `taskset -c 2` on an otherwise
idle machine, not through `isolated_bench run`. It is a diagnostic, not a
measurement. The binaries are the arms above.

**What it found.** Every recovered logarithm was the same in both arms.
The cold times are the prediction table above. What it set is the
prediction; the rows, the rule and the A/A requirement were not chosen
from it.

## Commands

```
# the base
git worktree add v3-wt 995ea2071cc7453877a503d30eb561ae82cddab9
git show ba345bbcf:Cargo.lock > v3-wt/Cargo.lock
cargo build --release --bin ic          # in v3-wt; SHA-256 00cff0b0…

# the candidate, from main at 92246530c
cargo build --release --locked --bin ic # SHA-256 d8daac9a…

# the round
for step in manifest pin aa compare holdout extend callgrind; do
  icprog run r10 $step --base <v3's ic> --cand <main's ic> \
    --isolate <isolated_bench> \
    --base-commit 995ea2071cc7453877a503d30eb561ae82cddab9 \
    --cand-commit 92246530c4f333724374e70897330bc38fa159a1
done
icprog analyse r10 > research/ic_tool_program/rounds/R10-main-head/analysis.json
```
