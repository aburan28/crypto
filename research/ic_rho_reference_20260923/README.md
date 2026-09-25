# The matched rho reference: calibration, evaluation, re-pricing

Evidence for §18 of
[`research/notes/index-calculus/RESEARCH_IC_BOUNDARY_LEDGER.md`](../notes/index-calculus/RESEARCH_IC_BOUNDARY_LEDGER.md).
The accounting contract
([`ec_index_calculus_contract.json`](../index_calculus_baseline_20260914/ec_index_calculus_contract.json),
`comparison_contract.rho`) prices a whole attack against "the eligible
automorphism-aware rho reference". Through §17 the prime and random binary
rows of the boundary ledger, and every `ic bench` row, divided by a walk that
used no automorphism. This directory holds the evidence for the matched walk
that replaces it.

It also holds the re-pricing of every frozen report whose ratio to rho the
scoreboard or the ledger quotes. The boundaries, the calibration protocol,
six targets and the abandon conditions were written into §18.1–§18.2 and
committed before anything here ran.

## The walks

| name | `A` | what | code |
|:--|--:|:--|:--|
| `frozen-plain` | 1 | the walk every row was priced against through §17: 16 jumps, two scalar multiplications per walk start, distinguished points every `√r/8` | `ic_boundary::rho_reference` (unchanged) |
| `tuned-plain` | 1 | the tuned walk on points: stride starts (one addition each), distinguished points every `r^{1/4}`, a jump table sized to `r` by joint double-and-add | `rho_reference_walk(…, RhoWalk::plain())` |
| `negation` | 2 | the same walk on `{P, −P}`: canonical representatives, the Wiener–Zuccherato look-ahead, short-cycle detection and a deterministic doubling escape — **the matched reference** on prime and random binary curves | `rho_reference_negation` |
| `signed-frobenius` | 2n | the repository's signed-Frobenius walk, priced as the Koblitz regime prices it | `ic_boundary::signed_frobenius_rho` |

Every walk charges every group operation it performs: the table, the starts,
the walk with its rejected look-ahead additions and escape doublings, and the
verification `[d]G = Q`. The tuned walks' `counters` say where each operation
went, and the tests check that the counters add up to the total.

## Layout

| path | what | ledger |
|:--|:--|:--|
| `calibration/prime.json`, `char2.json` | sixteen curves no evaluation uses (seed `0xCA11B`), `J ∈ {4, 8, 16}` for both tuned walks, 128 runs each: the jump-count rule | §18.3 |
| `evaluation/ladder.json` | sixteen fresh curves (seed `0xE7A1`, every prime curve generated), the three walks paired, 128 runs each | §18.4–§18.5 |
| `reprice/*.json` | the Round-5 ladder, headline and holdout, and the seven §17.11 whole-method reports, each re-priced; each first replays the frozen walk on the recorded seeds | §18.6 |
| `diagnostic/koblitz.json` | seven Koblitz curves: the frozen, negation and signed-Frobenius walks paired (not a declared target) | §18.8 |
| `analyse.py`, `analysis.json` | grades the six targets and produces every table §18 and the scoreboard quote | §18.4 |
| `run.sh` | every command, in order | §18.10 |

## Results, in brief

- **All six declared targets met.** The exponent target is met on the
  pooled fit only, at its lower edge (`α = 0.452` against `0.45`), because
  set-up is still visible at these sizes.
- **Every run verified** (21,096), and the frozen walk replayed exactly on
  all 616 recorded runs of the ten re-priced reports.
- **The matched walk sits `1.02`–`1.20×` its own floor from `2^{20}`.** Its
  whole `S` is `0.98`–`1.20` there. The paired `√2` is `1.362`, 95 %
  interval `[1.296, 1.430]`.
- **Every re-priced `vs rho` rises `2.1`–`6.5×`**, mostly because of the
  frozen walk's set-up. The scoreboard's best prime and binary rows go from
  `3.55×` and `16.0×` to `13.5×` and `38.2×`, and §17.11's pair table from
  `3.5×`, `29×` and `86×` to `21.0×`, `133×` and `257×`. The frozen ladder's
  one row below its rho (`0.84×` at 12 bits) reads `5.5×`.
- **Left open here, settled in §19** (`research/ic_rho_koblitz_20260923/`):
  the Koblitz reference charges nothing for its per-step canonicalisation
  (§18.8).  Priced, a canonical step costs `2.74–2.92` batched additions;
  but the collection thread's figures had a larger error the other way,
  32 targets against one, and against batch rho they read `5.6–8.6×`
  higher (`1.17×` → `6.51×`).

## Reproducing

```
cargo build --release --bin ic
sh research/ic_rho_reference_20260923/run.sh
python3 research/ic_rho_reference_20260923/analyse.py
```

`ic` never overwrites a report, so move the frozen files aside first. Every
report names the binary's blake3 and the commit it ran from, and every
re-pricing report names its source file's blake3.
