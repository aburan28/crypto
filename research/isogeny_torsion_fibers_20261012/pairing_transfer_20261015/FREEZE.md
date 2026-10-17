# Pairing-transfer input and point-selection freeze

Frozen before any paired multiplicative DLP solve. The protocol was committed at `d573c9c95ee4f79e2d56ce534926715a889aefc6`. The two solver rows in `INPUTS.jsonl` contain public curve coordinates and no answer scalar; their source is `../extension_field_20261014/FIXTURES.jsonl` at SHA-256 `1a63cc22a22e0fe27d44359325d502038d26d3ac616374fc5558c041939218b7`.

| Artifact | SHA-256 |
|:--|:--|
| `INPUTS.jsonl` | `ef5cd2ac1d8a2cb9aa3b7edc6f5b92408fafc00e7253b97dddc0d3555998ddc5` |
| `examples/torsion_pairing_transfer.rs` | `a9ffa68da164d7d3fb72f48665f935971c2a340f59d70b6d9f5a76d0ee7c4ff1` |
| `SELECTION_61.json` | `911d9a45e554812a26e554965f117fa4397bc6e2bb9317b2dc8ba9c65d28c9c4` |
| `SELECTION_271.json` | `958e910e8ad22092ffbef9a97b9a4636382a261273db0bd4e971b243dab210bc` |
| built release example | `ba39e9fc96fd8c82c912c4dd57f90423a3a867ac0cae0fe4215214a7805d2a9a` |

Build: `cargo build --release --example torsion_pairing_transfer --offline` with rustc 1.93.1 (Homebrew, aarch64-apple-darwin) on Darwin 25.6.0, arm64 T6041. Preparation commands were:

```sh
target/release/examples/torsion_pairing_transfer prepare research/isogeny_torsion_fibers_20261012/pairing_transfer_20261015/INPUTS.jsonl 0 > research/isogeny_torsion_fibers_20261012/pairing_transfer_20261015/SELECTION_61.json 2> research/isogeny_torsion_fibers_20261012/pairing_transfer_20261015/SELECTION_61.stderr
target/release/examples/torsion_pairing_transfer prepare research/isogeny_torsion_fibers_20261012/pairing_transfer_20261015/INPUTS.jsonl 1 > research/isogeny_torsion_fibers_20261012/pairing_transfer_20261015/SELECTION_271.json 2> research/isogeny_torsion_fibers_20261012/pairing_transfer_20261015/SELECTION_271.stderr
```

The selections are deterministic and carry the Rabin irreducibility checks, full curve order, first accepted `r`-torsion point, and both pairing values. Both preparation commands exited zero; both stderr files are empty. The earlier `PREPARE_DIAGNOSTIC_*.json` records came from a pre-freeze scratch build and are retained as diagnostics only. They have no exact source snapshot and cannot support timing claims or replace the frozen selections.

The committed input, implementation, and selection receipts are the only admissible inputs to the timed pairing runs. The timed runs redo construction and selection and charge them. The fixture's answer scalar is available only for an external post-solve comparison, never to the native solver.
