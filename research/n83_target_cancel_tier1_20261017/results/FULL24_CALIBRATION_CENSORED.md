# N83 cancellation: 24-point cloud panel and calibration gate

Cooperative two-worker target cancellation reduced the observed target-online
geometric mean from 219.617063 ms with join to 59.758863 ms across the 24
sealed public points. All 24 paired cancellation intervals were shorter than
their join intervals; the geometric paired join/cancellation ratio was
3.675054. These are descriptive measurements from cloud run `39048976945`,
attempt 1. The preregistered 24-point Tier 1 selector remains
`BOUNDED_UNKNOWN`: strata 0 and 1 failed the per-pair 10% before/after
calibration gate. The clean strata 2 and 3 remain admitted six-point results.

| Stratum | Point indices | Cloud status | Join geometric online ms | Cancellation geometric online ms | Paired join/cancellation | Faster cancellation points | Same-point rho/cancellation | Failed calibration pairs |
| ---: | --- | --- | ---: | ---: | ---: | ---: | ---: | --- |
| 0 | 0–5 | `BOUNDED_UNKNOWN` | 170.601142 | 67.362880 | 2.532569 | 6/6 | 40.085954 | point 1: 14.637%; point 3: 10.338% |
| 1 | 6–11 | `BOUNDED_UNKNOWN` | 289.882754 | 34.390063 | 8.429259 | 6/6 | 59.371109 | point 7: 10.270% and 10.336% |
| 2 | 12–17 | `STRATUM_ADMITTED` | 125.404958 | 39.286846 | 3.192034 | 6/6 | 64.131893 | none |
| 3 | 18–23 | `STRATUM_ADMITTED` | 375.098477 | 140.122705 | 2.676929 | 6/6 | 18.944648 | none |
| All, descriptive | 0–23 | `BOUNDED_UNKNOWN` | 219.617063 | 59.758863 | 3.675054 | 24/24 | 41.236519 | four pairs on two hosts |

Each cloud analysis records `PASS_INDEPENDENT_REPLAY`, `PASS` topology,
resource, event-order, and telemetry checks. Strata 0 and 1 record
`NOISY_OR_INCOMPLETE`, `NOT_ADMITTED`; the others record
`PASS_10_PERCENT_PER_PAIR`, `PASS_TIER1`. Their per-point timings,
calibration readings, exact resource counts, and cold costs are preserved in
the versioned analysis JSON. Cloud `PASS_INDEPENDENT_REPLAY` is the provider
verifier's result; a separate local replay of stratum 2 reproduced its
analysis byte for byte. The remaining local archive replays are tracked
separately and are not silently assumed by this summary.

The join candidate is
`IC1N83Ckb1a1fb164672PDP4rootRCguidedLAgaussTDdirectISO0h33c9c3b2f933`;
the cancellation candidate is
`IC1N83Ckb1a1fb164672PDP4rootRCguidedLAgaussTDdirectISO0h58ec337f71f7`.
Both use the same K992/root11 factor base with 164,672 usable points and
992 folded columns. Each online interval includes target query,
decomposition, relation check, descent, and producer recovery check after
target-independent setup. Same-point rho is independently solved for each
candidate's target. The paired ratio above uses the geometric mean of
per-point ratios, with candidate IDs used to align each point's arms.

| Stratum | Provider artifact ID | Raw ZIP SHA-256 | Versioned receipts |
| ---: | ---: | --- | --- |
| 0 | 12108389377 | `8428d634647950ae1533307ddaccd5049b75d5db7ef4151c7055102509519c93` | `stratum0_cloud_analysis.json`, `stratum0_provider_receipt.json` |
| 1 | 12109415285 | `d028fe586e27ab89d0c9159feef20dafbe927a6cc89d03b2c9bebf8d1afc788c` | `stratum1_cloud_analysis.json`, `stratum1_provider_receipt.json` |
| 2 | 12108835537 | `b1d0248f993b8fd6deefaec5cfa024d65db9ba1d46c3d1135860402f87115181` | `stratum2_analysis.json`, `stratum2_provider_receipt.json`, `stratum2_local_replay_receipt.json` |
| 3 | 12108209419 | `081f576f5e0fb030749e700529dc30c407ee854674882b5b1af538a78667a75c` | `stratum3_cloud_analysis.json`, `stratum3_provider_receipt.json` |

The next confirmation uses a new frozen 24-point domain and independently
verified same-point pairs on four Tier 1 instances. Calibration repetitions
and the 10% threshold must be fixed before launch. The full selector can
choose a winner only when all four host strata pass their gates.
