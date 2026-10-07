# ECC2K-130 DP ingest: Rust port protocol

## Hypothesis

The pure logic of `ecc2k130/aws/dp_ingest.py` — record decoding, key
classification, envelope checks, coverage counting, and the weight-32 cutoff
verdict — can move to a standalone Rust crate under `ecc2k130/dp-ingest/`
without changing any root-crate dependency, so the IC autolab frozen
`Cargo.lock` is untouched. The Python program stays the live daemon until a
later phase proves Postgres/S3 parity.

## Boundary

| kind | quantity |
|:--|:--|
| **Reference** | `ecc2k130/aws/dp_ingest.py` on `main`, held by `test_dp_ingest.py` |
| **Floor** | not a performance claim; this phase is a port |

This phase makes no end-to-end speed or rate claim. `S` and speedup stay unset.

## Frozen inputs

- Record layout: 32-byte little-endian `(seed, k0, k1, k2)`.
- Magics: `ECC2KDP2` (refuse), `ECC2KDT3` + `(3, 32)` header (strip after hash).
- Campaign weight constants: `CAMPAIGN_DP_WEIGHT = 32`,
  `CAMPAIGN_ITER_PER_DP_LOG2 = 28.41`, `DP_RATIO_MIN_RECORDS = 50000`,
  `DP_RATIO_TOLERANCE_LOG2 = 1.0`, `FIELD_BITS = 131`.
- Key shapes: legacy `dp/slot-N/<epoch>-<offset>.bin` and orbit
  `dp/slot-N/<32-hex>-<offset>-<64-hex>.bin`.
- Offline tests ported from `test_dp_ingest.py` classes `KeyShapes`,
  `Decoding`, `CutoffVerdict`, and the envelope/table-v3 body checks that need
  no network.

## Success

1. `cargo test` and `cargo clippy --all-targets -- -D warnings` pass in
   `ecc2k130/dp-ingest/` on the pinned toolchain 1.98.
2. Every ported offline test agrees with the Python values it was taken from
   (cutoff tests use the same seed numbers as `CutoffVerdict`).
3. Root `Cargo.toml` / root `Cargo.lock` / the frozen IC lock are unchanged.
4. `dp_ingest.py` is still the live ingest; this PR does not flip systemd,
   Lambda or `ingest.sh`.

## Stop

- A root dependency would be needed to make the pure tests pass.
- Binary COPY or status JSON formatting cannot be made to agree with Python
  in a later phase (declared here so Phase B knows the abandon condition).

## Cost accounting

Phase A: offline unit tests only. No fleet GPU time. No RDS.

## Phases

| phase | delivers | lives where |
|:--|:--|:--|
| **A (this PR)** | pure library + offline tests + CLI stub | `ecc2k130/dp-ingest/` |
| **B** | Postgres/S3 adapters + parity harness against Python on fixtures | same crate; ephemeral Postgres |
| **C** | cutover of systemd / Lambda / `ingest.sh`; Python kept for rollback | deploy scripts |

## Classification

**engineering** — the algorithm is unchanged; the language of the offline
half moves. No advance against a floor.
