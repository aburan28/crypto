# Trait detection on isogeny walks, queued

Status: **tooling; no ECDLP cost measured**

Date: 2026-10-04

Result class: not an index-calculus measurement.  No IC/rho ratio and no
scoreboard row.  The curves are those of
[`research/isogeny_walk_20261004`](../isogeny_walk_20261004/README.md), and
their walks were rebuilt byte for byte (routes hashes below).

## What changed

The walker now runs trait detection on every curve it finds, in two
scopes (`src/cryptanalysis/isogeny_walk/traits.rs`).

- **Class audits.** These are the repository's audits that depend only on
  `p` and `#E`, which every curve in an `F_p`-isogeny class shares:
  - `ecc_safety`: ten checks, from order size and smoothness through MOV,
    anomalous, Weil descent and the Cl(O) orbit;
  - the structural report: `|t² − 4p|` and the twist order trial-factored,
    `#E(F_{p^k})`, twist leakage and the phase-10 parity;
  - the Petit–Kosters–Messeng signals.

  They run once on the root.  The two that read a model and generator are
  re-run on walked curves, and the record states whether every verdict
  matched the root's.
- **Curve detectors.** Each reads one curve's model and generator and
  writes one `trait_status` entry: `non_singular`, `generator_valid`,
  `a_minus_3_model`, `qr_prefix_64` and `coefficient_bits` (bit lengths of
  `min(v, p − v)` for `a` and `b`).  A new detector added to
  `default_detectors()` runs everywhere.
- **Queueing.** `isogeny_walk plan` writes one `taskq.task-spec/v1` per
  shard.  Each job rebuilds the deterministic walk at a pinned commit and
  runs `isogeny_walk traits --shard i --of N` into `$TASKQ_OUTPUT_DIR`.
  `isogeny_walk collect` merges the shards and refuses:
  - a missing or repeated shard;
  - a curve covered twice or not at all;
  - shards built from different walks;
  - edited shard output.

  The specs pass taskq's own `normalize_spec`; this was checked ad hoc, not
  in CI.

## Runs

Each run was split into 4 or 8 shards, run locally exactly as a worker
would run them, and then collected.  No taskq worker was available on this
host (no Redis), so the queue path itself was exercised only up to spec
validation.

| run | curves | shards | routes SHA-256 (matches #1336) | invariance sample | all `ecc_safety` pass |
|:--|--:|--:|:--|:--|:--|
| [`p256-ell31-radius2`](runs/p256-ell31-radius2/) | 84 | 4 | `9241deccb957…` yes | 8 of 8 identical | yes |
| [`p224-ell61-radius2`](runs/p224-ell61-radius2/) | 93 | 4 | `1c84fe54c11b…` yes | 8 of 8 identical | yes |
| [`p256-ell61-20k`](runs/p256-ell61-20k/) | 20,000 | 8 | `05f7df41c391…` yes | 8 of 8 identical | yes |
| [`p224-ell61-20k`](runs/p224-ell61-20k/) | 20,000 | 8 | `e88d7c763a46…` yes | 8 of 8 identical | yes |

Each directory holds `collect.json` (per-trait distributions, the merged
`traits.jsonl` hash) and `class_audits.json`.  The radius-2 runs also hold
`walk.json`.  The 20k `traits.jsonl` files (about 12 MB each) are not
committed; `collect.json` pins their SHA-256, and the commands below
regenerate them.

### Class audits

| | P-256 | P-224 |
|:--|:--|:--|
| `ecc_safety` | all 10 pass | all 10 pass |
| PKM overall score (1 = resistant) | 0.775 | 0.647 |
| PKM special-prime score, Solinas weight | 0.50, 5 | 0.20, 3 (near `2^k`) |
| `\|t² − 4p\|` bits, residue after trial division | 258, 255-bit composite | 226, 182-bit composite |
| twist residue | 241-bit prime | 212-bit composite |
| max twist leak bits | 16 | 13 |
| phase-10 parity blocked | yes | yes |

P-224's low special-prime score is the PKM criterion's known flag on its
modulus `2²²⁴ − 2⁹⁶ + 1` (`pkm_criterion.rs` cites it).  It depends on `p`
alone and is therefore the same on every curve in the class.  The PKM
module measures vulnerability *signals*, not an attack.

### Curve detectors

| | P-256 20k | P-224 20k |
|:--|--:|--:|
| `non_singular` | 20,000 / 20,000 | 20,000 / 20,000 |
| `generator_valid` | 20,000 / 20,000 | 20,000 / 20,000 |
| `a_minus_3_model` | 10,195 (51.0%) | 4,937 (24.7%) |
| smallest `b`, bits of `min(b, p − b)` | 241 | 208 |
| expected smallest `b` over 20,000 uniform draws | ≈ 240.7 | ≈ 208.7 |

In the 93-curve P-224 run, walk node 19 is an `a = −3` model with a
208-bit `b`.  That is a `2⁻¹⁵` event per curve, or about 0.3% somewhere
among 93 curves.  It is one observation of a property of the recorded
model (`icwalk-canon/v1` picks the smaller of `b` and `−b`), not of the
curve's group.  A small `b` is not a known weakness.  In the 20k run the
same curve is still the smallest, at the expected minimum of 20,000 draws
(≈ 208.7 bits).

## Reproduce

```bash
cargo build --release --bin isogeny_walk
B=./target/release/isogeny_walk
$B walk --curve p224 --max-ell 61 --max-curves 20000 --no-class-audits --out W
for i in $(seq 0 7); do $B traits --curve p224 --dir W --shard $i --of 8 --class-audit-sample 8 --out S$i; done
$B collect --out T S0 S1 S2 S3 S4 S5 S6 S7
# queued instead: $B plan --curve p224 --max-ell 61 --max-curves 20000 --commit <sha> --shards 8 --out specs
```

## Not established

- Any DLP weakness.  Each class audit gives one verdict for the whole
  class.
- Any per-curve IC difference.  No detector here measures factor-base
  yield or solving degree at 224 or 256 bits.  The screening statistic for
  PR #1330's pilot must still be frozen before reading walk output.
- A live taskq run.  The specs are validated against taskq's normaliser
  only.
