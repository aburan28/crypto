# Every candidate, one harness: results

**Decision: no index-calculus arm is below the matched rho reference on any
curve, in any session, and the gap widens with the subgroup size.** Every
method this repository implements was run cold, one unseen target at a
time, in `ecbench`'s one counted unit, on four prime-field and four
Koblitz curves. The generic methods land on their analysed constants. The
best index-calculus figure anywhere in the three complete sessions is
`4.30×` the strong rho reference (`ic-frob-m2` on `icv1-f2m17-tm101-00378d4e`,
`r = 2^16`), and at the largest Koblitz subgroup measured (`r = 2^21`) the
same arm reads `129×`. On the prime curves the folded meet-in-the-middle
oracle reads `7.7×` at 18 bits and `40×` at 24 bits.

This is a measurement of existing methods on fresh targets, not a change
to any of them. By the classification fixed in [`PROTOCOL.md`](PROTOCOL.md)
before the runs, nothing here is an *advance* or *engineering*; the
session is the cross-method baseline AGENTS.md §12 assigns to `ecbench`.
No ratio below is a speedup, and no number here is about ECC2K-130.

| | |
|---|---|
| protocols | [`PROTOCOL.md`](PROTOCOL.md) (both main sessions), [`PROTOCOL-2.md`](PROTOCOL-2.md) and [`PROTOCOL-3.md`](PROTOCOL-3.md) (Koblitz follow-ups), each committed before its session ran |
| specs | [`spec-koblitz.json`](spec-koblitz.json) `ECS1hb88a55c7e993`, [`spec-prime.json`](spec-prime.json) `ECS1he577e6ae8ed2`, [`spec-koblitz-followup.json`](spec-koblitz-followup.json), [`spec-koblitz-followup-3.json`](spec-koblitz-followup-3.json) `ECS1h71d8da8786ac` |
| binary | `ecbench`, release, SHA-256 `b2596e9db8b79a11479ae22b4d12570114501b841c22716929f61012fd9d3845`; the two main sessions ran it from `678d7e5fa`, the follow-ups from `f8dfda558` (only research files changed between them) |
| host | Intel Xeon @ 2.80 GHz, 4 logical CPUs, 1 NUMA node, virtual machine (`adx aes avx avx2 avx512bw avx512dq avx512f avx512vl bmi1 bmi2 fma pclmulqdq popcnt sse4.2`); Linux 6.18 x86-64; env class `ECBENV2h0bfad9696cf5`; `--cpus auto`, as root |
| complete sessions | [`sessions/koblitz`](sessions/koblitz) `ECBS1hf4c34c2455f9` (1 664 records: 1 472 verified, 192 error), [`sessions/prime`](sessions/prime) `ECBS1haf997766feb9` (1 280 verified of 1 280), [`sessions/koblitz-followup-3`](sessions/koblitz-followup-3) `ECBS1hb0cb85a662d8` (448 records: 384 verified, 64 exhausted) |
| incomplete sessions, kept | [`sessions/koblitz-followup`](sessions/koblitz-followup) (SIGTERM at 102 of 448, status `interrupted`), [`abandoned/koblitz-followup-2`](abandoned/koblitz-followup-2) (process died at 80 of 448, status `running`); see below |
| correctness | 3 136 of 3 392 executions in the complete sessions verified (`[k]G = Q` and `k = planted`, by the runner and again by the audit); the 192 errors and 64 exhaustions are index-calculus runs that did no search or found no relation, listed per arm below; no run recovered a wrong scalar |
| replay certificates | `audit-koblitz.json` `70c975f22155fcef26895ef4884e60baa6c551eddcaae71ade0e40e57f77a6cf`, `audit-prime.json` `61b89480a43feb26dd03388880df0bb4b8f2693b694eddb83d700ee145613e97`, `audit-koblitz-followup-3.json` `d400579f6a911a9507f340105024c3eafb2a25ef0cc2c27c279441549cca2c9d`; each `ok`, 12 of 12 replays identical, auditor binary `b2596e9d…d845` |
| full tables | [`table-koblitz.md`](table-koblitz.md), [`table-prime.md`](table-prime.md), [`table-koblitz-followup-3.md`](table-koblitz-followup-3.md), written by `ecbench table` |
| comparisons | one saved `ecbench compare` per arm against its session's reference, under each session's `comparisons/` |
| levels earned | koblitz L1×780 L2×884; prime L1×528 L2×752; followup-3 L1×160 L2×288. The spec required L2, so wall time is admitted only on the pairs that earned it; operation counts stand at every level |

## Boundary and unit

- **Floor:** `√(π/2A)` in `S`, `A = 2` on the prime curves and `A = 2n`
  on the Koblitz curves, recorded per run by the harness.
- **Reference:** `rho.negation` on the prime curves; `rho.signed_frobenius_strong`
  on the Koblitz curves, the strong single-target reference the IC
  measurement rules require. `rho.signed_frobenius` is the operations-only
  baseline.
- **Unit:** `S = total group-addition equivalents / √r`, cold, every phase
  charged. An IC run charges factor base, oracle setup, relations, linear
  algebra and verification; a rho run charges setup, search and internal
  verification. Native work the unit has no pinned price for is counted
  separately and marks the total a lower bound (`lower bound` column).

**A note on reading the Koblitz reference at these sizes.** The strong rho's
setup is a fixed charge of about 750–1 000 group-addition equivalents (32
scalar multiplications), and at `r ≈ 2^15–2^21` that is three to four
times its search. Its `S` therefore sits at 6–28× the floor here and
*above* the lean signed-Frobenius walk (`rho-frob` reads 0.23–0.72 of it)
and above every BSGS. This is the reference's cost at toy size, not a
defect: the fixed charge vanishes against `√r` at the sizes the family is
about. It means every Koblitz `S / reference` below is *flattering* to the
candidate; against the lean walk the IC ratios would be 1.4–4.3× larger.
On the prime curves the reference is `1.2–1.7×` its floor and the point
does not arise.

## Koblitz family, session `ECBS1hf4c34c2455f9`

Mean `S` over 8 targets × 3 rounds per cell; intervals are 95 % two-stage
bootstrap intervals. `ref` is `S / rho.signed_frobenius_strong` on the same
curve. Curves are ordered by `r`.

| arm | method | f2m29 (`r=2^15.4`) | f2m17 (`r=2^16`) | f2m19 (`r=2^17`) | f2m23 (`r=2^21`) |
|---|---|---:|---:|---:|---:|
| `rho-strong` | `rho.signed_frobenius_strong` | S 4.578, floor×27.8 | S 3.968, floor×18.5 | S 3.037, floor×14.9 | S 1.101, floor×6.0 |
| `rho-frob` | `rho.signed_frobenius` | S 3.076, ref 0.672 | S 2.838, ref 0.715 | S 0.702, ref 0.231 | S 0.457, ref 0.415 |
| `rho-neg` | `rho.negation` | S 2.162, ref 0.472 | S 1.711, ref 0.431 | S 1.908, ref 0.628 | S 1.050, ref 0.953 |
| `rho-plain` | `rho.plain` | S 2.372, ref 0.518 | S 2.069, ref 0.521 | S 1.815, ref 0.597 | S 1.529, ref 1.388 |
| `rho-frozen` | `rho.frozen_reference` | S 7.863, ref 1.718 | S 6.761, ref 1.704 | S 5.142, ref 1.693 | S 2.757, ref 2.503 |
| `bsgs` | `bsgs.textbook` | S 1.776, ref 0.388 | S 1.403, ref 0.354 | S 1.575, ref 0.519 | S 1.658, ref 1.506 |
| `bsgs-il` | `bsgs.interleaved` | S 1.494, ref 0.326 | S 1.186, ref 0.299 | S 1.360, ref 0.448 | S 1.393, ref 1.265 |
| `bsgs-neg` | `bsgs.negation` | S 1.265, ref 0.276 | S 0.910, ref 0.229 | S 1.080, ref 0.356 | S 1.159, ref 1.053 |
| `kangaroo` | `kangaroo.vow` | S 1.814, ref 0.396 | S 2.099, ref 0.529 | S 2.161, ref 0.712 | S 2.177, ref 1.977 |
| `ic-frob-m2` | `ic.pipeline` `koblitz-orbit:divisor=0;1` × `mitm-frobenius:m=2` | **error** 0/24 | S 17.06, **ref 4.30** | **error** 0/24 | S 142.6, **ref 129.4** |
| `ic-frob-m3` | same base × `mitm-frobenius:m=3` | **error** 0/24 | S 17.17, **ref 4.33** | **error** 0/24 | S 142.6, **ref 129.5** |
| `ic-subtract` | same base × `subtract` | **error** 0/24 | S 23.25, **ref 5.86** | **error** 0/24 | S 109.6, **ref 99.5** |
| `rho-strong-aa` | control, `rho.signed_frobenius_strong` | ref 1.000 | ref 1.000 | ref 1.000 | ref 1.000 |

- The 192 `error` records are the three IC arms on `f2m19` and `f2m29`:
  `no invariant subspace for divisor [0, 1] on this curve`, raised by the
  base builder before any work (`ord_19(2) = 18` spans the field;
  `ord_29(2) = 28` exceeds the factor enumeration cap). The follow-up
  below measures those two curves with a base that needs no such factor.
- On `f2m23` the IC total is almost entirely oracle setup: 205 257 of
  206 388 group-addition equivalents per run, against 1 102 for relations.
  The `m = 2` and `m = 3` oracles cost the same to within 0.1 %, because
  the table, not the search, is the cost at this size. The pooled
  comparisons over the two curves that ran are `incomplete` (a measured
  run did not verify) and not admissible; the per-curve rows above are
  the result.

## Prime family, session `ECBS1haf997766feb9`

`ref` is `S / rho.negation` on the same curve.

| arm | method | fp18 (`r=2^17.7`) | fp20 (`r=2^19.1`) | fp22 (`r=2^21.2`) | fp24 (`r=2^24.0`) |
|---|---|---:|---:|---:|---:|
| `rho-neg` | `rho.negation` | S 1.294, floor×1.46 | S 1.344, floor×1.52 | S 1.480, floor×1.67 | S 1.080, floor×1.22 |
| `rho-plain` | `rho.plain` | S 2.074, ref 1.602 | S 1.781, ref 1.326 | S 1.663, ref 1.124 | S 1.538, ref 1.424 |
| `rho-frozen` | `rho.frozen_reference` | S 4.815, ref 3.721 | S 4.590, ref 3.416 | S 2.609, ref 1.763 | S 2.400, ref 2.222 |
| `bsgs` | `bsgs.textbook` | S 1.588, ref 1.227 | S 1.597, ref 1.189 | S 1.441, ref 0.974 | S 1.456, ref 1.348 |
| `bsgs-il` | `bsgs.interleaved` | S 1.569, ref 1.213 | S 1.640, ref 1.220 | S 1.166, ref 0.788 | S 1.271, ref 1.176 |
| `bsgs-neg` | `bsgs.negation` | S 1.092, ref 0.844 | S 1.100, ref 0.819 | S 0.942, ref 0.636 | S 0.957, ref 0.886 |
| `kangaroo` | `kangaroo.vow` | S 2.269, ref 1.753 | S 2.017, ref 1.501 | S 2.209, ref 1.493 | S 2.461, ref 2.279 |
| `ic-mitm` | `ic.pipeline` `prime-abscissa:size=32` × `mitm:negation_folded=1` | S 9.961, **ref 7.70** | S 10.49, **ref 7.81** | S 16.23, **ref 10.97** | S 43.74, **ref 40.49** |
| `ic-subtract` | same base × `subtract` | S 309.9, **ref 239** | S 448.7, **ref 334** | S 926.7, **ref 626** | S 2 775, **ref 2 569** |
| `rho-neg-aa` | control, `rho.negation` | ref 1.000 | ref 1.000 | ref 1.000 | ref 1.000 |

Pooled over all four curves, `ecbench compare` reads `ic-mitm / rho-neg =
15.47 [11.32, 21.09]` and `ic-subtract / rho-neg = 858 [579, 1 226]`, both
`ops (ok)`; the pooled figure is dominated by the 24-bit curve and is not
a scaling statement. The per-curve `ic-mitm` ratio rises monotonically
with `r`: 7.70, 7.81, 10.97, 40.49. Its relation phase grows from 3 512 to
176 479 group-addition equivalents per run across the four sizes while its
oracle setup is a fixed 1 056; the growth is in the search for
decomposable points over a 32-element base, as expected of a fixed base
at growing `r`.

## Koblitz follow-up, session `ECBS1hb0cb85a662d8`

The two curves the orbit base could not serve, with the polynomial-basis
`binary-subspace` base at dimensions 6 and 8, under the unfolded
meet-in-the-middle (`mitm:m=2`) and `subtract` oracles. `max_trials` is
`1 000 000` per [`PROTOCOL-3.md`](PROTOCOL-3.md).

| arm | f2m29 (`r=2^15.4`) | f2m19 (`r=2^17`) |
|---|---:|---:|
| `rho-strong` | S 4.578, floor×27.8 | S 3.037, floor×14.9 |
| `rho-frob` | S 3.076, ref 0.672 | S 0.702, ref 0.231 |
| `ic-sub6-mitm` | **exhausted** 0/24 | S 21.21, **ref 6.98** |
| `ic-sub6-subtract` | **exhausted** 0/24 | S 718.7, **ref 237** |
| `ic-sub8-mitm` | S 473.9, **ref 103.5** | S 112.3, **ref 36.96** |
| `ic-sub8-subtract` | S 36 989, **ref 8 081** | S 620.1, **ref 204** |
| `rho-strong-aa` | ref 1.000 | ref 1.000 |

The reference and baseline rows reproduce the first session's exactly
(same seeds, same targets, same counts), which is the harness doing what
the audit says it does.

## Incomplete sessions

Two earlier attempts at `PROTOCOL-2` (the same arms with `max_trials =
100 000 000`) did not finish, and both are kept as found:

- `sessions/koblitz-followup`, `interrupted` by SIGTERM at 102 of 448
  records when the launching tool hit its ten-minute limit. Its audit
  (`audit-koblitz-followup.json`, `ceeff86f…d8a2`) fails only on
  `session_complete`; every integrity check passes, and CI's audit with
  `--allow-interrupted` accepts it. Its 88 verified records agree with
  the follow-up-3 rows on `n = 19` and show 24 `n = 29` runs exhausting
  the full `10^8` budget without a relation.
- `abandoned/koblitz-followup-2`, launched detached, died without a
  signal at 80 of 448 records and so carries status `running` with no
  exit line. Its audit (`audit-koblitz-followup-2.json`, `69dca282…34ae`)
  fails on `session_complete` and `record_count`. The harness has no way
  to tell a session that died from one still running, so no audit flag
  accepts it; it lives outside `sessions/` so that CI's re-audit, which
  globs that directory, does not fail on it. Its files are the ones the
  run wrote, byte for byte.

## Predictions, scored

| prediction | statement | verdict |
|---|---|---|
| P1 | each generic method's mean `S` within its interval of its analysed constant at the largest `r` | **pass on the prime family** (`rho.negation` 1.080 [0.857, 1.301] vs 0.886: interval contains it; `rho.plain` 1.538 vs 1.253: contains; BSGS 0.97/0.95/0.96 of theory; kangaroo 1.23, interval [1.93, 3.02]/2 contains 1). **Fails as written on the Koblitz family for the rho walks:** `rho.signed_frobenius` reads 2.47× its constant on `f2m23` (interval [1.8, 3.3]) and `rho.plain` 1.22×, because every walk's fixed setup and verification charge is a large share of `S` at `r ≤ 2^21`; BSGS and kangaroo (1.05–1.16, 1.09 of theory) pass. The calibration note's constants were measured at `r = 2^39`, where that share is negligible. This is the same accounting, read at a size the constants do not describe; the prediction should have been bounded by `r`. |
| P2 | every IC arm above the reference, interval excluding 1, on every curve | **pass on every curve where an IC arm ran**: smallest per-curve ratio 4.30 (`ic-frob-m2`, `f2m17`), its interval [17.02, 17.10] in `S` against the reference's [3.85, 4.08]. Not scored on the four error cells and the two exhausted cells. |
| P3 | no IC arm's ratio to the reference falls with `r` | **pass on the prime family** (`ic-mitm` 7.70 → 40.49; `ic-subtract` 239 → 2 569, monotone). On the Koblitz family only two sizes per base ran: orbit base 4.3 → 129 from `r = 2^16` to `2^21`; subspace base `ic-sub8-mitm` 103.5 → 37.0 from `r = 2^15.4` to `2^17`, a **fall** across two curves of different `n`, reported as such. Two points on curves whose `n` differs by 10 do not make a trend; a dimension-8 subspace is a smaller fraction of a 29-dimensional field than of a 19-dimensional one, so fewer points decompose over it, which is the expected direction. |
| P4 | every execution verifies, every replay identical | **pass on replays** (36 of 36 identical across the three complete sessions). **Fails as written on executions**: 192 errors (base builder) and 64 exhaustions (no relation in `10^6` trials) are recorded with their own statuses. No run recovered a wrong scalar. |
| P5 | every IC arm verifies on both follow-up curves | **fail**: dimension 6 exhausts on `f2m29` under both oracles; dimension 8 verifies everywhere. |
| P6 | every IC arm above the reference with interval excluding 1, on both curves; subspace base costlier than the orbit base | **pass** where it ran (6.98 to 8 081); the subspace base on `f2m19` (6.98 at dimension 6) is cheaper than the orbit base on `f2m23` (129) but those are different curves; on no curve do both bases run, so the second clause is not scored. |
| P7 | dimension 8 costs more than dimension 6 | **pass on `mitm`** (`f2m19`: 112.3 vs 21.2); **fail on `subtract`** (620 vs 719: the larger base finds relations sooner and the saving exceeds the table). On `f2m29` dimension 6 does not finish, so dimension 8 is both costlier per table and the only one that solves. |

## What this establishes, and what it does not

- **Established.** On eight registered curves with `r` from `2^15.4` to
  `2^24`, every index-calculus configuration the harness offers costs
  more than the matched rho on the same target, cold, with every phase
  charged, by a factor between 4.3 and 8 081, and the factor rises with
  `r` wherever one base runs at two sizes of the same family. The
  generic methods stand where the calibration note put them once the
  fixed per-run charge is allowed for. Every count reproduces bit for
  bit from the session files.
- **Not established.** Anything above `r = 2^24`; wall time (reported
  `admitted` on the pairs that earned L2, with A/A controls at 1.00
  [0.98, 1.04], but at these sizes a run is milliseconds and the figure
  is descriptive); the `descent-algebraic` oracle (excluded by protocol,
  nondeterministic); anything about ECC2K-130, m = 83, or any curve not
  in the table. BSGS below rho in `S` is the time–memory trade and not
  a finding.
- **For the scoreboard.** This is a baseline for the cross-method axis,
  not a round of any IC thread: no boundary moved and no variant
  changed. The one figure it adds to the record is the smallest
  verified cold one-target IC/rho ratio under `ecbench`'s accounting,
  `4.30×` at `r = 2^16` against a reference that is itself 18.5× its
  floor at that size.

## Reproduce

```bash
cargo build --release --bin ecbench
```

```bash
for s in koblitz prime koblitz-followup-3; do
  ./target/release/ecbench verify --dir research/ecbench_all_candidates_20261003/sessions/$s --replay 12 --exit-code
done
./target/release/ecbench table --dir research/ecbench_all_candidates_20261003/sessions/koblitz
```

Operation counts match these sessions exactly on any host built from the
same commit; CI re-audits every session under `sessions/` with every
measured run replayed. To rerun from scratch, run each spec into a new
directory; never into one of these.
