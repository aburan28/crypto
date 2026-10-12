# Native ICV1/FB1 reproduction result

Date run: 2026-10-04

Protocols: [PROTOCOL.md](PROTOCOL.md) and
[PROTOCOL_STANDALONE.md](PROTOCOL_STANDALONE.md)

## Verdict

**PASS.** The native Rust replay reproduced every frozen identity, count, and
digest for the depth-18 `dickson-torus` factor base on
`icv1-fp256-t89188191154553853111372247798585809583-f188c491` (P-256).
The result is a correctness and storage result.  It is not an ECDLP solve,
an F4 measurement, or a speed/exponent claim.

The first unmerged prototype added an `ecbench` subcommand.  CI rejected that
delivery surface because historical n37 evidence seals the `ecbench` source,
not because a target mismatched.  The preregistered standalone replay then
produced the same 66,020,135-byte JSON and left the frozen source files
identical to the PR base.

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
cargo build --release --bin p256_factor_base
./target/release/p256_factor_base \
  --curve icv1-fp256-t89188191154553853111372247798585809583-f188c491 \
  --factor-base dickson-torus:depth=18 \
  --out /tmp/icv1-fp256-t89188191154553853111372247798585809583-f188c491--FB1.factor-base.json \
  --sql-out /tmp/icv1-fp256-t89188191154553853111372247798585809583-f188c491--FB1.factor-base.sql \
  --relation-length 17 \
  --verify
sqlite3 -bail /tmp/p256-wide-fb.db \
  < /tmp/icv1-fp256-t89188191154553853111372247798585809583-f188c491--FB1.factor-base.sql
```

The full JSON dump is deterministic derived data and is deliberately not
committed.  Its reproducibility manifest is:

| item | value |
|---|---|
| schema | `ecbench.factor_base_dump/v1-wide` |
| bytes | 66,020,135 |
| SHA-256 | `2c35ed7db8b5c774101c1667dbe7b643f7f6c212799e63c48fc0bc45ca762279` |
| library source SHA-256 | `c44d0d95e1d6ccb68c6bbc456a52091d98081b8704fe19e15affadd9a03f6f40` |
| CLI source SHA-256 | `1cda9ac323a697c74a0c5776374c36b81d8b469b79093122fa4b9a1620228e76` |
| release binary SHA-256 | `c7d7fb9a1db48ba3b0ef021dd33ad3894394bbf51c6d8f6ecb141bae5722408e` |

The generated SQL loader is 77,212,627 bytes with SHA-256
`f84c90287838ecd70e4c85c2aa0f15c9d39cc831dafb9fe06574064b5ca841b3`.

The command above is its extraction/regeneration procedure from the committed
source.  The checked-in depth-9 vector pins the enumeration and identity
encoding, and the CLI integration test constructs, verifies, parses, and
converts a small wide dump to SQL.

## Verification receipt

- `cargo +1.98 test --lib p256_dickson_factor_base`: 4 passed;
- `cargo +1.98 test --test p256_factor_base`: passed;
- `cargo +1.98 clippy --all-targets -- -D warnings`: passed;
- full `--verify` rebuild: passed;
- SQLite ingestion: one factor base and 262,478 point rows;
- frozen `ecbench` source diff against the PR base: empty;
- host: AMD EPYC 9V74, Linux 6.18.44 x86_64;
- toolchain: `rustc 1.98.1 (48a229cea 2026-09-01)`.

No timing is reported.  This run changes no scoreboard or leaderboard result
because it creates no `ic.pipeline` session, matched rho reference, or
whole-pipeline measurement.
