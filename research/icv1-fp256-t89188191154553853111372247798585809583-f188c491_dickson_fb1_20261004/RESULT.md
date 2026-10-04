# Native ICV1/FB1 reproduction result

Date run: 2026-10-04

Protocol: [PROTOCOL.md](PROTOCOL.md)

## Verdict

**PASS.** The native Rust replay reproduced every frozen identity, count, and
digest for the depth-18 `dickson-torus` factor base on
`icv1-fp256-t89188191154553853111372247798585809583-f188c491` (P-256).
The result is a correctness and storage result.  It is not an ECDLP solve,
an F4 measurement, or a speed/exponent claim.

| quantity | frozen | native result |
|---|---:|---:|
| depth | 18 | 18 |
| root exponent | `0x30000000000000000000` | `0x30000000000000000000` |
| terminal trace | 0 | 0 |
| negation-folded columns | 131,239 | 131,239 |
| materialised signed points | 262,478 | 262,478 |
| legacy enumeration SHA-256 | `e52a7b606604641efc99d5a0ef933bb258457e4b65f3c3112be88d792d9e217e` | same |
| sorted wide point-key SHA-256 | `8fa207da35a4426d3992a1b8e9b08ccbfc8795208f450cf76b44bd47d05bd915` | same |
| factor-base identity | `FB1h2ea06bef7f7a` | `FB1h2ea06bef7f7a` |
| full FB1 SHA-256 | `2ea06bef7f7aa68ce83ad0889116a621885ed03937b5de5481ace9e1dc3626ad` | same |
| native verification | required | pass |
| database rows | 262,478 | 262,478 |

The registered identities carried by the dump are
`EC1P256Cp256h0523b774e066` and
`urn:ec-record:1:sha256:0523b774e0666cf37ff4dc9cf452eae024192a7f89bdad10ca226d6e64d4dd70`.

## Relation-count diagnostic

For 17 distinct columns with independent signs, the exact candidate domain is

```text
374163014528376139389697548555821950542279161732309414538616264465591578132480
```

Dividing by the P-256 subgroup order gives Poisson mean
`3.2313348613016615` and heuristic probability
`1-exp(-lambda) = 0.9604952697587432`.  The saturated Boolean selector ideal
has the structural upper bound `degree of regularity <= 18` (`m+1`).  This is
not an observed Gröbner-basis degree: no summation-polynomial/F4 solve was run,
so the measured degree of regularity remains unknown.

## Native command

```bash
cargo build --release --bin ecbench
./target/release/ecbench fb-wide \
  --curve icv1-fp256-t89188191154553853111372247798585809583-f188c491 \
  --factor-base dickson-torus:depth=18 \
  --out /tmp/icv1-fp256-t89188191154553853111372247798585809583-f188c491--FB1.factor-base.json \
  --relation-length 17 \
  --verify
./target/release/ecbench db sql \
  /tmp/icv1-fp256-t89188191154553853111372247798585809583-f188c491--FB1.factor-base.json \
  | sqlite3 -bail /tmp/p256-wide-fb.db
```

The full JSON dump is deterministic derived data and is deliberately not
committed.  Its reproducibility manifest is:

| item | value |
|---|---|
| schema | `ecbench.factor_base_dump/v1-wide` |
| bytes | 66,020,135 |
| SHA-256 | `2c35ed7db8b5c774101c1667dbe7b643f7f6c212799e63c48fc0bc45ca762279` |
| source SHA-256 | `c8bf22bb46bd7de23f32faa587328f053bc790a69b4b7762e641ca99ce23f13a` |
| release binary SHA-256 | `fbca2da0ed4ffafc0482bcf75437ca4176234c40e08333633a5efac96de8f3e8` |

The command above is its extraction/regeneration procedure from the committed
source.  The checked-in depth-9 vector pins the enumeration and identity
encoding, and the CLI integration test constructs, verifies, parses, and
converts a small wide dump to SQL.

## Verification receipt

- `cargo test p256_dickson_factor_base --lib`: 3 passed;
- `cargo test wide_factor_base_dump_loads --lib`: 1 passed;
- `cargo test --test ecbench a_wide_p256_factor_base_dump_round_trips_through_sql`: passed;
- full `--verify` rebuild: passed;
- SQLite ingestion: one factor base and 262,478 point rows;
- host: AMD EPYC 9V74, Linux 6.18.44 x86_64;
- toolchain: `rustc 1.99.0 (b940084d7 2026-09-28)`.

No timing is reported.  This run changes no scoreboard or leaderboard result
because it creates no `ic.pipeline` session, matched rho reference, or
whole-pipeline measurement.
