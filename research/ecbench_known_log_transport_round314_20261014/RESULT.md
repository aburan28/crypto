# Round 314: exact-log transport does not reach rho parity

Date: 2026-10-14

**The post-holdout hash schedule improves the sequential transported-log
candidate by 21.3--36.8%, but its complete cold counted cost remains
2.121--3.267 times the matched strong Pollard-rho reference.  It also does
not transfer to the fixed P-256 Dickson factor base.**  This is a reproducible
negative result, not a parity claim.

The experiment follows the frozen [preregistration](PREREGISTRATION.md), its
pre-holdout [K rule](FROZEN_RULE.json), and the sole permitted
[post-holdout iteration](ITERATION_PROTOCOL.md).  The candidate is measured
through `ecbench` from cold factor-base construction to exact scalar replay.
All 960 final-session records verify; all 800 measured records replayed
identically in the saved audits.  The four operation-count A/A comparisons
are exactly 1.  Every execution earned L0, so wall time is descriptive and is
not used below.

## Boundary table

The single unit is `S = cold ecbench.gae / sqrt(r)`.  Every row is a lower
bound because native field arithmetic, hashing, allocation and modular
combination are counted but unpriced.  Ratios are paired ratios of total
counted work over 40 measured one-target runs: eight public targets and five
rounds per target.

| subgroup, in increasing `r` | variant | class | mean S | candidate / rho | correctness |
|:--|:--|:--|--:|--:|:--|
| `icv1-f2m37-tm534059-32aad96b`, log2 r=27.781 | strong signed-Frobenius rho | reference | 0.281047 | 1.000 | 40/40 verified |
| same | pair claw | generic control | 0.406623 | 1.447 | 40/40 verified |
| same | sequential transported logs | baseline | 0.757512 | 2.695 | 40/40 verified |
| same | hash transported logs | engineering | 0.596163 | **2.121** [1.908, 2.409] | 40/40 verified; A/A=1 |
| `icv1-f2m43-tm998717-e2e742b0`, log2 r=32.113 | strong signed-Frobenius rho | reference | 0.170048 | 1.000 | 40/40 verified |
| same | pair claw | generic control | 0.335075 | 1.970 | 40/40 verified |
| same | sequential transported logs | baseline | 0.575686 | 3.385 | 40/40 verified |
| same | hash transported logs | engineering | 0.440667 | **2.591** [2.190, 3.059] | 40/40 verified; A/A=1 |
| `icv1-f2m47-t22705043-f4e44623`, log2 r=36.636 | strong signed-Frobenius rho | reference | 0.120642 | 1.000 | 40/40 verified |
| same | pair claw | generic control | 0.327921 | 2.718 | 40/40 verified |
| same | sequential transported logs | baseline | 0.623242 | 5.166 | 40/40 verified |
| same | hash transported logs | engineering | 0.394109 | **3.267** [2.578, 4.301] | 40/40 verified; A/A=1 |
| `icv1-f2m41-tm2308219-7f48b14a`, log2 r=39.000 | strong signed-Frobenius rho | reference | 0.163145 | 1.000 | 40/40 verified |
| same | pair claw | generic control | 0.330307 | 2.025 | 40/40 verified |
| same | sequential transported logs | baseline | 0.710814 | 4.357 | 40/40 verified |
| same | hash transported logs | engineering | 0.515723 | **3.161** [2.361, 4.389] | 40/40 verified; A/A=1 |

The hash schedule's paired cost relative to the sequential candidate is
0.787, 0.765, 0.632 and 0.726 in increasing subgroup order.  This is a real
constant-factor improvement.  It is not an index-calculus advance: the
candidate chooses scalar-labelled points and is therefore a generic
time--memory method.  The already registered pair claw is cheaper on every
cell.

## Scaling and dominant phases

The native analyzer fits `log2(metric) = intercept + exponent*log2(r)` over
the four frozen sizes.  These are finite-size fits, not asymptotic proofs.

| metric | fitted exponent | 95% t interval, df=2 | R2 |
|:--|--:|--:|--:|
| hash candidate cold counted GAE | 0.473794 | [0.335270, 0.612318] | 0.990850 |
| strong rho cold counted GAE | 0.417818 | [0.235452, 0.600184] | 0.979834 |
| hash pair-table entries | 0.491786 | [0.467495, 0.516076] | 0.999737 |
| hash target attempts | 0.255981 | [0.072722, 0.439240] | 0.947535 |
| peak process RSS KiB | 0.073087 | [0.017278, 0.128897] | 0.940741 |

The candidate work-exponent point estimate passes the preregistered 0.55
gate, but its interval includes values above 0.55 and the measured parity
gate fails.  The pair table consumes 60.8--69.0% of cold GAE; target PDP
consumes 19.8--36.8%.  The exact materialized table widths are 4,441,
18,894, 86,774 and 206,465 entries.  Peak process RSS is 7,316, 8,084,
10,488 and 12,996 KiB.  The algorithm materializes no on-disk table and
reports zero algorithmic disk traffic; session serialization is outside the
solve interval.

The growing table is the dominant obstruction even after removing the
sequential schedule's additive correlations.  The measured fit is consistent
with a square-root time--memory method, not evidence for a sub-square-root
ECDLP algorithm.

## Correctness and replay

Each reported decomposition is replayed as a full group sum.  Each
transported column log is checked by a charged scalar multiplication, and
each recovered discrete logarithm is checked by `[d]G = Q`.  Across the four
final sessions there are:

- 960/960 verified executions, including warmups;
- 800/800 measured executions reproduced bit-for-bit by `--replay-all`;
- 160 exact counted A/A pairs;
- zero wrong answers, exhausted runs, crashes or timeouts;
- zero observed false positives or false negatives on the checked instances.

The audit receipts are
`86b49fa611f08157f18b014c8b0d6855fccd0f89c659b7d62972af5e97dcd172`
(n37),
`0cc5b71fcd7eadac625b7d84cfff73f5e348105b19425c777e9dfdc9a1299c80`
(n43),
`f44247caba6ac1513c6d005a00033196d1ab37471616a2b913147890c5d69a16`
(n47), and
`58f93c2243202d87eb3e58b2bd66f488102fe20ca5877d499beb4da1aea30bc8`
(n41).  They are same-environment deterministic replays; CI supplies the
separate Linux replay before merge.

## P-256 transfer and degree decision

The motivating object remains
`icv1-fp256-t89188191154553853111372247798585809583-f188c491` with
`FB1h2f8621cda105`.  Arbitrary points `[s]G` have known logs because their
construction supplies `s`; that does not place them in this fixed Dickson
base.  Conversely, producing `s` for an existing FB1 point is the ECDLP being
solved.

The repository's exact
[P-256 transport certificate](../icv1-fp256-t89188191154553853111372247798585809583-f188c491_fb_log_transport_round21_20261006/RESULT.md)
has SHA-256
`4aaf5fde72257ed1f262d5674ab2fe81c7b4f90698fa3767615adfed6524e757`.
It exhausts 8,413,312 small-multiple images and 5,986,344 completed
Dickson-block images, replays every candidate, finds zero transport
relations, and leaves all 131,458 logarithm dimensions independent.
Therefore the exact P-256 membership/log-transport gate fails and no
full-depth unplanted P-256 attempt is authorized.

This method bypasses polynomial solving, so it establishes no new Groebner
degree of regularity.  For the fixed P-256 base, the independent certificate
still gives structured residual degree 4 and leaves the unsplit S17 degree of
regularity unknown.  A low local degree does not remove the 131,458-dimensional
global logarithm problem.

## Promotion decision

| gate | result |
|:--|:--|
| every record verifies and every measured record replays | pass |
| exact counted A/A | pass |
| degree-37 hash candidate / rho <= 1.30 | **fail: 2.121** |
| four-size work exponent <= 0.55 | pass by point estimate: 0.473794 |
| cold setup retained and storage reported | pass |
| exact transport into the fixed P-256 factor base | **fail: rank 0, K=131,458** |
| all conjunctive gates | **fail** |

Round 314 therefore closes as **engineering improvement; negative parity and
P-256 transfer result**.  Its one preregistered post-holdout iteration has
been consumed.  Another algorithmic branch needs a new preregistration and
fresh hidden targets; retuning against these sessions would not be valid
evidence.

## Reproduction

The deterministic analysis is
[`ANALYSIS.json`](ANALYSIS.json), SHA-256
`e887ca37106e24b34fc089f9fe13e8f2bd048780118bbf5495065d423b6e16c2`.

```bash
cargo test --release --lib transported_log
cargo build --release --bin ecbench --example known_log_transport_round314_analyze

target/release/ecbench verify \
  --dir research/ecbench_known_log_transport_round314_20261014/sessions/iteration-n37 \
  --replay-all --exit-code
# Repeat for iteration-n43, iteration-n47 and iteration-n41.

mkdir /tmp/round314-analysis-replay
target/release/examples/known_log_transport_round314_analyze \
  /tmp/round314-analysis-replay/ANALYSIS.json
cmp research/ecbench_known_log_transport_round314_20261014/ANALYSIS.json \
  /tmp/round314-analysis-replay/ANALYSIS.json
```
