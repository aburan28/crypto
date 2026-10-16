# Exact S3 root-table Bloom filter in the online Koblitz example

`KIC_ROOT_BLOOM=1` adds a three-hash Bloom prefilter to the exact canonical
root table in [`koblitz_orbit_dlp_fast_online.rs`](koblitz_orbit_dlp_fast_online.rs).
The table inserts every exact key into the filter. A negative filter answer
skips a table probe; a positive answer still goes through the exact table,
which retains the existing first-witness and shift rules. In
`koblitz_orbit_dlp_fast_online.rs` the option is off by default; in
`koblitz_orbit_dlp_fast.rs` it is on by default since 2026-10-15
(`KIC_ROOT_BLOOM=0` disables it there). The summary records
`root_bloom_enabled` and `root_bloom_bytes`.

This is a transfer of the opt-in root filter in
[crypto PR #1659](https://github.com/aburan28/crypto/pull/1659) to the
separate tracked online-timing example. The branch starts from crypto main
`4b1cc9324513166e9219cf166be08ca0725f5f9b`. Its method differs from
the conjugate-x A1 filter in
[cryptanalysis PR #561](https://github.com/aburan28/cryptanalysis/pull/561):
this filter screens already canonicalized S3 root keys, while A1 screens
x-coordinate orbit candidates before orbit-key construction.

The release example test command passed one table test:

```sh
CARGO_TARGET_DIR=/private/tmp/crypto-s3-bloom-target-20261010 \
  cargo test --release --offline --example koblitz_orbit_dlp_fast_online
```

The test inserts 4,096 deterministic keys into both filtered and unfiltered
tables, checks every inserted key and its first value, then compares 100,000
other lookups. The release example binary was also run twice with
`construct:23:1:8`, scalar fixture `17`, rank seed `20261010`, and
`KIC_RANK_THREADS=1`, changing only `KIC_ROOT_BLOOM=0/1`. The normalized
target records matched with SHA-256
`cd7045c0b714950b1c01ff31befc179d9a4f9b11d6837924b982f7255d0f8545`;
the normalized summaries matched with SHA-256
`d785a8839eff917e4453ab3663ab1faacc3fd3735a676c9c5b671193b8238db1`.
Normalization excludes timing fields, monotonic clock timestamps, peak RSS,
and the filter fields. Both runs had rank 8, 71 probes, the same ordered
point indices, x codes and intermediates, and verified scalar 17. The
enabled filter occupied 4,096 bytes.

This fixture is an integration and correctness check. The current N23 wall
measurements are ordinary-host diagnostics. A disjoint one-target panel on
an isolated CPU host, with exact preparation memory and online phase costs,
is the next measurement gate. The repository-wide release library test was
attempted on the same main lineage for PR #1659 and has failures in
unchanged library code; a focused rerun reproduced the N71 construction
failure at `src/cryptanalysis/koblitz_fast_arith.rs:1207`.
