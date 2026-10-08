# The first two challenge verdicts: toy prime and toy Koblitz, epoch 1

**Result.** Two standing challenges from `docs/bounds/challenges/` were run
once each, at epoch 1, and judged with every measured run replayed. On the
prime toy domain, `bsgs.negation` against the reference walk `rho.negation`
reads **`trade`**: `0.8250 [0.725, 0.945]` of the incumbent's operations with
`15.79 [13.38, 19.06]` times its table. On the Koblitz toy domain,
`rho.signed_frobenius` against the strong single-target reference
`rho.signed_frobenius_strong` reads **`advances` at the constant level**:
`0.2991 [0.235, 0.402]` of the incumbent's operations and `0.560 [0.483,
0.656]` of its table. Each verdict wrote one bound record. Both are what the
seed frontier already showed as unpaired entries, now as paired verdicts on
frozen workloads; neither moves an exponent, and neither says anything above
`r = 2^24`.

Nothing here is index calculus and nothing here is a speedup over the
generic floor. The BSGS row is the time–memory trade `docs/ecbench/README.md`
§7 describes; the Koblitz row prices the strong walk's fixed set-up at toy
sizes, where it is 71 % of that walk's total.

| | |
|---|---|
| protocol | [`docs/bounds/README.md`](../../docs/bounds/README.md) §§5–6; the challenges are the frozen protocol, the epoch fixes the targets |
| challenges | [`prime-toy-reference.json`](../../docs/bounds/challenges/prime-toy-reference.json) `ECCH1h904625668686` (incumbent `rho.negation`, bound `ECBND1h7b69b9787056`); [`koblitz-toy-reference.json`](../../docs/bounds/challenges/koblitz-toy-reference.json) `ECCH1h76d71da2a89c` (incumbent `rho.signed_frobenius_strong` at registry defaults, bound `ECBND1h1573ea580e89`) |
| epoch | 1 for both; nobody had run either challenge before |
| specs | [`specs/prime-toy-reference-e1.json`](specs/prime-toy-reference-e1.json) `ECS1h6ae84963ddfc`, [`specs/koblitz-toy-reference-e1.json`](specs/koblitz-toy-reference-e1.json) `ECS1h97880b56eeb4`, each written by `ecbench challenge spec` and checked against the challenge by the verdict (`spec_matches: true`) |
| sessions | [`sessions/prime-toy-reference-e1`](sessions/prime-toy-reference-e1) `ECBS1h3a5708ad901b`; [`sessions/koblitz-toy-reference-e1`](sessions/koblitz-toy-reference-e1) `ECBS1h5eaf9eef9f25`; each `complete`, 384 records (96 warm-up + 288 measured), 384 verified |
| binary | the pinned release `ecbench`, SHA-256 `58c0eea59b37c190b86f7729fb5f613ee79d9a99a2e08b6df670f0379a519016`, run from outside a git worktree, so the sessions record `git_commit: null`; the hash is the identity |
| host | Intel Xeon @ 2.10 GHz, 4 logical CPUs, 1 NUMA node, virtual machine (KVM), Linux 6.18 x86-64; the sessions' own `host.json` carries class `ECBENV2h393e9c141834`, which the verdicts cite; [`HOST.txt`](HOST.txt) is the capsule printed afterwards from the shell |
| isolation | `--cpus none` on a container without a CPU reservation: **every measured run is L0**. The challenge asks for L2, and that gates wall time only; no axis, fit or acceptance rule reads a clock (README §8 rule 4), so the verdicts are unaffected and no wall-clock figure is reported below |
| audits | run by the verdicts themselves with `--replay-all`: [`audit-prime-toy-reference-e1.json`](audit-prime-toy-reference-e1.json) `ok`, 288 of 288 replays `identical`; [`audit-koblitz-toy-reference-e1.json`](audit-koblitz-toy-reference-e1.json) `ok`, 288 of 288 replays `identical`. An independent `ecbench verify --replay-all --exit-code` on each directory agrees |
| verdicts | [`verdict-prime-toy-reference-e1.json`](verdict-prime-toy-reference-e1.json) `ECVD1h0b0f5e0d7f67`; [`verdict-koblitz-toy-reference-e1.json`](verdict-koblitz-toy-reference-e1.json) `ECVD1h53a3dfe36118`; both re-derive byte for byte (`docs/bounds/refit.sh`) |
| new bounds | [`challenge-prime-toy-reference-e1-bsgs-negation.json`](../../docs/bounds/records/challenge-prime-toy-reference-e1-bsgs-negation.json) `ECBND1he0aea8671f10`, `improves_on` `ECBND1h7b69b9787056`; [`challenge-koblitz-toy-reference-e1-rho-signed-frobenius.json`](../../docs/bounds/records/challenge-koblitz-toy-reference-e1-rho-signed-frobenius.json) `ECBND1hbc8602528f48`, `improves_on` `ECBND1h1573ea580e89` |
| tables | [`table-prime-toy-reference-e1.md`](table-prime-toy-reference-e1.md), [`table-koblitz-toy-reference-e1.md`](table-koblitz-toy-reference-e1.md), written by `ecbench table` |
| comparisons | one `ecbench compare --a incumbent --b candidate --save` per session, under each session's `comparisons/` |
| logs | [`run-prime-toy-reference-e1.log`](run-prime-toy-reference-e1.log), [`run-koblitz-toy-reference-e1.log`](run-koblitz-toy-reference-e1.log) (`--quiet`, one line each) |

## The question

Does the candidate do what the incumbent does with fewer operations, less
memory, or both, on the challenge's frozen curves and on targets nobody saw
before the epoch was named? A `compare` inside a session of one's own design
is a measurement; only a verdict against a sealed challenge moves the
frontier. The two reference challenges were the natural first ones to run:
each asks the question every later candidate on its domain will be asked.

Each spec interleaves three arms on the same 32 one-target workloads (four
curves, eight planted targets each, three measured rounds after one warm-up):
`incumbent`, `candidate`, and `incumbent-aa`, an A/A control of the incumbent
with shared seeds. The control's paired ratio is `1.0000 [1.000, 1.000]` in
both sessions, as a deterministic method under shared seeds must give.

## Prime, toy: `bsgs.negation` against `rho.negation`

Curves `icv1-fp18-tm337-d28d5e09`, `icv1-fp20-t727-cd198a38`,
`icv1-fp22-tm1385-475dcb5f`, `icv1-fp24-t1059-6df599df` (`r = 2^17.7` to
`2^24.0`), floor `√(π/4) = 0.886` in `S`. Mean `S` over 24 verified runs per
cell, 95 % two-stage bootstrap intervals, `ref` the ratio to `incumbent` on
the same curve.

| arm | fp18 (`2^17.7`) | fp20 (`2^19.1`) | fp22 (`2^21.2`) | fp24 (`2^24.0`) |
|---|---:|---:|---:|---:|
| `incumbent` `rho.negation` | S 1.305 [1.139, 1.478], floor×1.473 | S 1.289 [1.048, 1.560], floor×1.455 | S 1.247 [1.006, 1.503], floor×1.408 | S 1.143 [0.844, 1.484], floor×1.290 |
| `candidate` `bsgs.negation` | S 1.111 [0.969, 1.257], floor×1.254, **ref 0.851** | S 1.055 [0.900, 1.219], floor×1.191, **ref 0.819** | S 0.987 [0.801, 1.181], floor×1.113, **ref 0.791** | S 0.960 [0.765, 1.164], floor×1.083, **ref 0.840** |
| `incumbent-aa` control | ref 1.000 | ref 1.000 | ref 1.000 | ref 1.000 |

Verdict `ECVD1h0b0f5e0d7f67`, statement verbatim:

> trade (better on ops; worse on memory): bsgs.negation / rho.negation = 0.8250 [0.725, 0.945] in ecbench.gae over 96 pairs on 4 size(s) (memory worse, uncharged indistinguishable (reported, not deciding)); bounded: unpriced work on one or both arms; α 0.460 against 0.444; A/A control 1.0000 [1.000, 1.000]

| axis | incumbent | candidate | Σ candidate / Σ incumbent | 95 % | reads |
|---|---:|---:|---:|---|---|
| `ops` (× floor) | 1.406 | 1.160 | 0.8250 | [0.725, 0.945] | better, decides |
| `memory` (entries / √r) | 0.032 | 0.502 | 15.79 | [13.38, 19.06] | worse, decides |
| `uncharged` (/ √r) | 1.054 | 1.014 | 0.962 | [0.822, 1.130] | indistinguishable, reported only |

Per curve the paired ratio is `0.851 [0.704, 1.055]`, `0.819 [0.622, 1.057]`,
`0.791 [0.559, 1.078]`, `0.840 [0.610, 1.140]`: every point estimate below
one, no single curve's interval excluding it; the pooled 96 pairs do. The
stage rows say where it lands: the candidate spends 52 % of its total in
`setup` (the baby-step table, `4.58 ×` the incumbent's set-up) and 48 % in
`search` (`0.447 ×` the incumbent's). Fits: candidate `α 0.460 [0.377,
0.502]`, incumbent `0.444 [0.338, 0.530]`, `exponent_moved: false`. The
incumbent measured `1.406 × [1.270, 1.531]` the floor here against its
recorded `1.466 × [1.278, 1.653]`: inside, so the recorded bound and this
session agree.

New bound `ECBND1he0aea8671f10`: `bsgs.negation`, `1.160 ×` the floor
`[1.068, 1.248]`, `0.502 √r` entries, `α 0.460 [0.377, 0.502]`, four sizes,
96 verified runs, `bounded` (inserts and lookups counted, not charged). It
names `ECBND1h7b69b9787056` under `improves_on` and joins the prime toy
frontier beside the two unpaired `bsgs.negation` entries.

## Koblitz, toy: `rho.signed_frobenius` against `rho.signed_frobenius_strong`

Curves `icv1-f2m29-tm40309-30c52b96` (`r = 2^15.4`),
`icv1-f2m17-tm101-00378d4e` (`2^16.0`), `icv1-f2m19-t797-b6cf2467`
(`2^17.0`), `icv1-f2m23-t5197-69e76b73` (`2^21.0`); floor `√(π/4n)`, so
`0.165`, `0.215`, `0.203`, `0.185` in `S`. The incumbent is the strong
lockstep walk at its registry defaults (`dp_bits=4`, `lanes=32`,
`step_cap_factor=2000`); the candidate is the lean tuned walk
(`cap_multiple=64`).

| arm | f2m29 (`2^15.4`) | f2m17 (`2^16.0`) | f2m19 (`2^17.0`) | f2m23 (`2^21.0`) |
|---|---:|---:|---:|---:|
| `incumbent` `rho.signed_frobenius_strong` | S 4.538 [4.395, 4.682], floor×27.58 | S 4.066 [3.940, 4.194], floor×18.92 | S 3.082 [2.964, 3.200], floor×15.16 | S 1.088 [1.007, 1.173], floor×5.89 |
| `candidate` `rho.signed_frobenius` | S 1.173 [0.956, 1.680], floor×7.13, **ref 0.259** | S 1.447 [0.851, 2.641], floor×6.73, **ref 0.356** | S 0.760 [0.663, 0.917], floor×3.74, **ref 0.247** | S 0.440 [0.319, 0.619], floor×2.38, **ref 0.405** |
| `incumbent-aa` control | ref 1.000 | ref 1.000 | ref 1.000 | ref 1.000 |

Verdict `ECVD1h53a3dfe36118`, statement verbatim:

> advances (better on ops, memory) at the constant level: rho.signed_frobenius / rho.signed_frobenius_strong = 0.2991 [0.235, 0.402] in ecbench.gae over 96 pairs on 4 size(s) (memory better, uncharged better (reported, not deciding)); bounded: unpriced work on one or both arms; α 0.229 against 0.126; A/A control 1.0000 [1.000, 1.000]

| axis | incumbent | candidate | Σ candidate / Σ incumbent | 95 % | reads |
|---|---:|---:|---:|---|---|
| `ops` (× floor) | 16.886 | 4.996 | 0.2991 | [0.235, 0.402] | better, decides |
| `memory` (entries / √r) | 0.030 | 0.017 | 0.560 | [0.483, 0.656] | better, decides |
| `uncharged` (/ √r) | 1.106 | 0.483 | 0.437 | [0.250, 0.742] | better, reported only |

Per curve: `0.259 [0.209, 0.372]`, `0.356 [0.209, 0.662]`, `0.247 [0.214,
0.298]`, `0.405 [0.299, 0.550]`; every interval excludes one. The stage rows
locate the move in `setup`: the incumbent spends 71 % of its total there (32
lane scalar multiplications, a fixed charge of roughly 750 to 1 000
group-addition equivalents) and the candidate 36 %, paired ratio `0.163`; the
`search` phase reads `0.664`. Fits: candidate `α 0.229 [-0.013, 0.491]`,
incumbent `0.126 [0.097, 0.266]`, `exponent_moved: false`. Both exponents are
the signature of a fixed cost at `2^15` to `2^21`, not a law (README §4), and
the level named is `constant`. The incumbent measured `16.886 × [9.156,
25.275]` the floor against its recorded `16.794 × [9.107, 25.310]`: inside.

New bound `ECBND1hbc8602528f48`: `rho.signed_frobenius`, `4.996 ×` the floor
`[2.985, 7.374]`, `0.017 √r` entries, `α 0.229 [-0.013, 0.491]`, four sizes,
96 verified runs, `bounded` (canonicalisations counted, not charged). It names
`ECBND1h1573ea580e89` under `improves_on`. On the rebuilt page the Koblitz
toy frontier holds three `rho.signed_frobenius` entries and nothing else; the
strong walk's toy bound is dominated by all three.

## Class

**Accounting** (AGENTS.md §3), both. No algorithm changed; two registered
methods were paired on frozen workloads for the first time and the frontier
records what the pairing showed. The prime `trade` is a BSGS table bought with
a `0.50 √r` table against the walk's `0.03 √r`, the textbook time–memory
trade. The Koblitz `advances` is
the strong walk's set-up priced at sizes where `√r` is a few hundred: an
engineering fact about the reference at toy size, and the reason the IC
measurement rules keep `rho.signed_frobenius_strong` as the reference by rule
rather than by this page. A verdict does not change that rule.

## What this does and does not say

- **Toy tier.** Fields of 18 to 24 bits and degrees 17 to 29. Nothing above
  `r = 2^24` is measured, and `toy` never speaks for `medium`: on the medium
  Koblitz domain (`research/ecbench_pair_claw_20261003`, degrees 41 to 61) the
  strong walk reads `1.089 × [0.927, 1.261]` the floor and the lean walk
  `1.188 × [0.954, 1.464]` over three degrees; no medium verdict exists.
- **One epoch each.** A second epoch draws fresh targets from the same nonce;
  verdicts are never averaged across epochs. The challenge ids and the epoch
  are the identity of these results.
- **Constant level only.** Neither candidate's `α` interval is disjoint from
  its incumbent's over the four sizes; no exponent moved.
- **No wall-clock figure.** All 288 measured pairs in each session are L0, so
  `compare` reports the wall ratio as descriptive only (medians `0.77` and
  `1.20`, the lean Koblitz walk slower in time while cheaper in operations).
  Operation counts are the result; they survive hardware and reproduce bit for
  bit (576 of 576 replays identical across the two sessions).
- **Bounded totals.** Both verdicts carry `bounded: true`: table inserts,
  lookups, canonicalisations and partition hashes are counted and not charged
  on one or both arms, reported under `uncharged`, and in neither verdict does
  the unpriced work shift against the candidate (`uncharged_shift: false`).

## Reproduce

```sh
cargo build --release --bin ecbench
E=./target/release/ecbench; B=research/ecbench_bounds_challenges_20261006
$E challenge spec --challenge docs/bounds/challenges/prime-toy-reference.json --candidate '{"id":"bsgs.negation"}' --epoch 1 --out /tmp/spec-prime.json
cmp /tmp/spec-prime.json $B/specs/prime-toy-reference-e1.json
$E verify --dir $B/sessions/prime-toy-reference-e1 --replay-all --exit-code
$E challenge verdict --challenge docs/bounds/challenges/prime-toy-reference.json --dir $B/sessions/prime-toy-reference-e1 --epoch 1 --replay-all --bounds docs/bounds/records --root . --out /tmp/verdict-prime.json --exit-code
cmp /tmp/verdict-prime.json $B/verdict-prime-toy-reference-e1.json
```

The same four lines with `koblitz-toy-reference`, `'{"id":"rho.signed_frobenius"}'`
and `koblitz-toy-reference-e1` reproduce the second verdict.
`docs/bounds/refit.sh` runs both verdicts and re-derives every bound record.
A new epoch runs `ecbench run --spec` on a fresh `challenge spec --epoch N`
into a new directory under `sessions/`; these two directories are evidence and
are never edited.
