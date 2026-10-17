# V2 instrumentation and replay freeze

Preregistered in `AMENDMENT_V2.md` at commit `2546920f04`. The V1 `runs/` and `SUMMARY.json` are retained. V2 uses the same public `INPUTS.jsonl`, modulus choices, selected H receipts, target, cap, and seed schedule. The V2 native executable adds phase counters and writes only to `runs_v2/`.

| Artifact | SHA-256 |
|:--|:--|
| `INPUTS.jsonl` | `ef5cd2ac1d8a2cb9aa3b7edc6f5b92408fafc00e7253b97dddc0d3555998ddc5` |
| `SELECTION_61.json` | `911d9a45e554812a26e554965f117fa4397bc6e2bb9317b2dc8ba9c65d28c9c4` |
| `SELECTION_271.json` | `958e910e8ad22092ffbef9a97b9a4636382a261273db0bd4e971b243dab210bc` |
| `examples/torsion_pairing_transfer.rs` V2 | `382ea60f8746c5145ebffe6d25f6809a520b09974f9f5e294da9d96b96d82de7` |
| release executable V2 | `9daceafea0dce4811299ff87a6f4ed3818dc1b49ed5410dffb2133fa56af804d` |
| `PREPARE_V2_61.json` | `fb2c885b509d051e97dd4ac99bf9f7ad72b4e110a3939c38c38afc00e3933ef7` |
| `PREPARE_V2_271.json` | `f0a125f5938ceb10420a5e1aebc6eb68d4a8dbc270b5140c9bc8efc713724f2d` |

Built with `cargo build --release --example torsion_pairing_transfer --offline`, rustc 1.93.1 on Darwin 25.6.0 arm64. Both V2 `prepare` commands exited zero. A structural equality check compared the V1 and V2 tower coefficient, Rabin certificate, complete selected-point receipt, selected H, α, β, full extension-operation counts, and base-field-operation counts; it returned true in both cells. This screening did not run the multiplicative DLP solver. The V2 `prepare` receipts additionally expose candidate-search, Miller, final-exponentiation and pairing-check charges.

After this freeze commit, execute exactly `sh research/isogeny_torsion_fibers_20261012/pairing_transfer_20261015/RUN_V2.sh`. The script checks all five source/input/selection/binary hashes before writing a new run. Its six pairings, six independent scalar replays, six bilinearity checks and twelve rho controls use the same schedule as V1. Do not replace a failed V2 row with a later silent rerun.
