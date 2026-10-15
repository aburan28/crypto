# Phase 0 pilot: findings (2026-10-15)

**What ran.** Nine `ecbench` sessions, one or two public *train* targets per
arm (target seed `20261015100`), counts only (`--cpus none`, L0), on this
repository's cloud container (4 vCPU, 15 GiB). Binary SHA-256 in
[`ecbench.sha256`](ecbench.sha256); method inventory in
[`methods.txt`](methods.txt). Specs `spec-*.json`; sealed sessions under
[`sessions/`](sessions/); every record laid out in [`TABLE.md`](TABLE.md)
(by `summarise.py`, recorded values only). **These are feasibility and
pricing probes, not measurements: one or two targets, no intervals, and
nothing here is a comparison result.**

| session | curve | arms | outcome |
|:--|:--|--:|:--|
| `pK13` | Koblitz `E_0`, n=13 (r ≈ 2^11.0) | 7 × 2 targets | 12 verified, 2 exhausted (negfold, full-rank linear algebra) |
| `fK19` | Koblitz n=19 (r ≈ 2^17.0) | 7 × 1 | 7 verified, 19 s |
| `fK23` | Koblitz n=23 (r ≈ 2^21.0) | 7 × 1 | 7 verified, 6 min (F4 arm dominant) |
| `pP20` | `prime_search` 20-bit | 5 × 2 | 10 verified |
| `fP24`, `fP28` | PKM ladder curves | 6 × 1 | all verified, < 5 s |
| `fP32`, `fP36`, `fP40` | `find_a3_curve` rungs | 5 × 1 | rho/BSGS/kangaroo verified; both IC arms exhausted at 2·10⁶ trials with a 64-point base |

## Findings that change the protocol

1. **Almost nothing is a complete cost, including the references.** Every
   rho, BSGS and kangaroo record is a lower bound (canonicalisations,
   partition hashes, table inserts and queries unpriced). The IC arms are
   complete only where the repository has a *pinned calibration*
   (`calibration_pinned_units`: ns per double, lookup, row op, AS solve,
   word XOR, Frobenius, canonicalisation): negfold and F4 are complete at
   n=13 and n=23 and lower bounds at n=19, which has no calibration. So
   "complete IC vs lower-bound rho" comparisons in earlier studies mixed a
   wall-time-calibrated conversion on one side with nothing on the other.
   → **Amendment A3**: report every arm with (a) charged `S` and (b) a
   sensitivity band pricing each unpriced native counter at 0, 0.1 and 1
   GAE; rank only where the conclusion holds across the band.
2. **Some arms do not emit the counters they leave unpriced.**
   `ic-compact` (12 counters) and both prime IC arms (sqrt solves,
   lookups, row ops, Legendre symbols, inversions) report the names but
   not the counts (`n/e` in `TABLE.md`), so no band can be computed.
   → **Implementation gap** in `src/cryptanalysis/ecbench/` (outside this
   study's reserved paths): emit the counts. Until fixed these arms are
   reported as lower bounds only.
3. **Calibration coverage is per curve.** Pinned conversions exist for some
   sizes only (n=13, 23 here; not 19). → Either calibrate every ladder size
   on one fixed host before Phase 1 or report the uncalibrated sizes in
   the band only.
4. **Full-rank linear algebra breaks negfold.** With
   `incremental-gauss-full-rank-checked`, negfold exhausted 200,000 trials
   at n=13 on both targets (`S` 4,555); with target-pin linear algebra it
   verifies. → **Amendment A2** (applied): target-pin linear algebra for
   every IC arm, the scalar verified outside the solver.
5. **`prime_search` cannot build the prime ladder** (O(p) point count;
   planning a 28-bit curve had not finished after several minutes and was
   stopped). → **Amendment A1** (applied):
   `examples/ic_frontier_prime_curves.rs` emits the deterministic `a = −3`
   P-256-class rungs of `prime_fast::find_a3_curve` (Hasse-interval BSGS
   order) at 20–40 bits in < 1 s ([`../prime_curves.jsonl`](../prime_curves.jsonl)).
6. **IC parameters must scale with the curve.** A fixed 64-point prime base
   exhausts from 32 bits on. For two-summand IC the optimum base is
   N ≈ (2r)^{1/3} (relations need ≈ 2r/N trials, the pair table N²), which
   predicts cost ∝ r^{2/3} against rho's r^{1/2}. → **Amendment A4**: size
   rules, not tuned constants — N = ⌈(2r)^{1/3}⌉ for prime m=2 arms;
   Koblitz dimension/columns chosen per size on the train panel (Phase 0
   step 3, still to run) with the rule recorded.
7. **The inherited-F4 arm is too slow to ladder.** One execution took about
   6 minutes at n=23 (charged `S` 85,071, complete); 32 executions per size
   would take hours and n ≥ 29 is out of reach. → F4 runs at n ≤ 23 only,
   as a stage reference.
8. **BSGS has the lowest charged `S` on several curves** (it trades √r
   memory for time). → **Amendment A5**: the frontier is a time–memory
   Pareto set; retained memory is reported for every arm and BSGS is
   never ranked against a constant-memory method on time alone.

## Not yet done (Phase 0)

- Step 3, parameter selection on the train panel for each Koblitz size
  (n ∈ {13, 17, 19, 23, 29, 31, 37, 41}) and the prime size rule check.
- Calibrations for the uncalibrated ladder sizes (finding 3).
- The counter-emission fix (finding 2), which needs a scope expansion into
  `src/cryptanalysis/ecbench/`.
- Phase 2 (wall time, second-host replay) needs a RunPod pod; this
  container has no `RUNPOD_API_KEY`.
