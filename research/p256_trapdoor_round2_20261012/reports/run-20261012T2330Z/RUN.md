# P-256 trapdoor search, round 3 — run report

- **Run id.** run-20261012T2330Z (UTC start 2026-10-12T23:30Z)
- **Code.** `src/cryptanalysis/p256_trapdoor_round3.rs` at this PR's commit.
- **Command.** `cargo run --release --example p256_trapdoor_round3 -- 1048576`
- **Inputs.** 2^20 deterministic splitmix64-derived 160-bit seeds (counters 0..2^20); X9.62 A.3.3.1 derivation with P-256's `p` and `a = −3`; roots pooled (`±b`), 1,048,170 b-samples from 1,048,576 seeds (the QR filter rejects ~0.02%... recorded: 2^20 − 406 seeds rejected).
- **Runs requested / completed.** 1 / 1, exit 0, no outlier flagged.

## Measured (copied from round3_records.json)

| Probe | Published | Honest (n = 1,048,170) | corrected p |
|---|---|---|---|
| QR(b) | true | 0.500 | 1.000 |
| cubic residuosity of b | false | 0.332 | 1.000 |
| Hamming weight of b | 134 | 128.00 ± 8.00 | 1.000 |
| CF max quotient of b/p | 1275 (11 bits) | ≥ published: 0.5092; ≤: 0.9565; max seen 30 bits | 1.000 |
| Hamming weight of r | 125 | 127.50 ± 7.99 | 1.000 |
| QR(r) | true | 0.500 | 1.000 |

Seed-preimage phrase probe: 1000 curated phrases, no preimage of
`c49d360886e704936a6678e1139d26b7819f7e90` (expected; bounty unclaimed).

## Graph

`graph.svg` plots the recorded published values against the recorded
honest summary statistics. Recorded values only.

## Observations

- No outlier: the published `b` is at the honest median on every
  probed statistic (see note §5 "Round 3").
- Interpretation belongs to the note; this report records measured
  values only.
