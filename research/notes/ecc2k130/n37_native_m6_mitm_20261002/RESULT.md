# Frozen n37 descendant-native six-sum: full rank, incomplete target coverage

**Decision: `DIRECT_NATIVE_K42_M6_N37_L1024_BATCH_NO_GO`.** The fixed 42-column
descendant-native base reaches independently verified rank 42, but direct
six-summand decomposition recovers only **917/1024** frozen public target
logs. The other **107** are complete misses of the exact 3+3 oracle, not
timeouts or unverified answers. This rejects the *bare direct* fixed-base
batch route on this Q block. A residual shift, changed base, higher arity,
or other search policy must be separately frozen and charged. The result
does not reject those variants or establish a full-method IC/rho cost ratio.

The [preregistered protocol](PROTOCOL.md), committed as `3805bcd0`, pinned
the base and inputs before the measured run. Source was commit
`b675df8560f1b98cf1d85299291fdc3c10e11e84` on PR #1254, built with
Rust 1.93.1 on macOS arm64. The fixed base is the 42 leaf points in
[`NATIVE42.json`](../n37_native_basis_bridge_20261002/NATIVE42.json), derived
from the degree-73 descendant without consulting any Q. The producer read
only public target coordinates; the separate published fixture scalars were
read *afterward* by the independent replay.

| Arm | Formal sum capacity `V/r` | Frozen Q hits / 1024 | Descriptive 95% Wilson interval under uniform-Q sampling | Verified target logs | Decision |
|:--|--:|--:|--:|--:|:--|
| Signed native ≤5 | 0.161008 | 146 (14.26%) | 12.25–16.53% | diagnostic only | Consistent with the 16.10084% uniform-population support ceiling; finite samples can exceed that ceiling. |
| Signed native ≤6 | 2.288827 | 917 (89.55%) | 87.53–91.28% | 917 | Rank succeeds, direct 1024-target batch fails on 107 complete misses. |

The complete half table has 105,995 raw signed multisets of length at most
three and 102,391 distinct group sums. Independent source-curve reachability
recomputed the same set sizes `[1,85,3613,102391]` for arities 0–3. The
relation phase took 49 deterministic target-blind probes: 42 supported,
all 42 independent, and seven complete misses. A separate dense modular
solve reproduced all 42 base logs. Every reported relation and all 917
reported target logs were replayed using the *source-curve pullbacks* and
checked by full-point scalar multiplication. All 917 recovered values also
equal the published fixture scalars. No scalar is reported for the 107 misses.

[`RAW.json.gz`](RAW.json.gz) is the 82,629-byte gzip archive, SHA-256
`6a0e8c9243aae9360048fc0918d55463ea09b980caeae6d62669855bad043346`.
It expands to the 1,585,747-byte producer record with SHA-256
`2cc0275cfe5b2fc0bef4f84bd3dce9105b7d9576ed26340c8a5b404a16cd9d7e`.
It retains every relation probe, rank step, witness, target decision,
recovered scalar, source coordinate, phase timer, and logical operation
count. [`REPLAY.json`](REPLAY.json), SHA-256
`34d3a15ead7e547f2275835352d7cdb6edba4eeeb538201097c77f37c55ce095`,
is the independent receipt. Its oracle constructs at-most-three reachable
sets by iterative source-curve group-law closure, checks complete 2+3 and
3+3 membership for every Q and relation probe, then resolves the relation
matrix independently. It does not call the producer's half-table or lookup.

The verified input SHA-256 values are:

| Input | SHA-256 |
|:--|:--|
| Frozen native base | `bd4bd8af982bcc65234ae1fae1ad52eb0ed41a1c203bf4d2556fd8e87f788a4c` |
| n37/L1024 campaign specification | `da958a3f1117dd1b88703fed2055c5c9f64255516e1e920f35d193bf53c6a88d` |
| Public block-00 points | `187ec04fe50326bbb2f17dadf76056841f04af37a8abb8ff7f502fcd531711ad` |
| Degree-73 archive | `eb4773d556886672e8b80735fc486c16b73765fd5069ff81133100bfc3eafa90` |
| Published scalar fixture, replay only | `124c13b3daba0477b471466f95c4a9e6f9ce9fe01d5ea67b6582469cfb9d5697` |

The measured build used the repository-root `Cargo.lock` with SHA-256
`b28d3c2d81146a40d00df85f45c6460d9e2bc5c26875307a06cfd15efff99365`
(29,115 bytes). Because the root lockfile is ignored by this repository,
the exact bytes are preserved as [`Cargo.lock`](Cargo.lock) in this evidence
directory. Copy it to the repository root before using `--locked`.

This was **one diagnostic process**, not a matched performance campaign.
[`TIMING.txt`](TIMING.txt) records 4.07 s process wall, 2.66 s user and
0.05 s system for the combined five-sum control **and** six-sum candidate.
The instrumented six-sum query phase took 2.929 s and the five-sum control
0.354 s; these are descriptive single-run stage times, not a cold candidate
speedup. The same process also charged input validation and degree-73
transport (37.3 ms), table construction (77.7 ms), relation collection
(123.0 ms), dense solve/base verification (0.23 ms), and target-log
verification (11.8 ms). Internal timing before output serialization/write
was 3.550 s. The process timer includes output. The replay process cost
14.42 s wall and is a verification receipt, not attack work.

The producer counted 105,910 logical point additions building the table,
1,507,265 in relation lookups, 31,296,573 in six-sum target lookups,
3,699,712 in the separate five-sum control, 4,059 explicit scalar
multiplications across setup/relations/verification, and witness additions
in the raw record. Excluding the five-sum control, the instrumented oracle
and external witness replays requested 32,920,936 logical additions. Scalar
multiplication internals, field-basis work and isogeny polynomial work are
included in wall time but not converted to group-addition equivalents.
**Matched same-host signed-Frobenius batch rho, complete native cost `S`,
and speedup remain unset.** The earlier source-policy cold CPU/rho ratio
3.031 is not this arm's comparator.

Reproduce from this PR head, whose measured implementation was commit
`b675df8560f1b98cf1d85299291fdc3c10e11e84`:

```sh
cp research/notes/ecc2k130/n37_native_m6_mitm_20261002/Cargo.lock Cargo.lock
cargo test --locked --lib native_signed_mitm
cargo build --release --locked --example n37_native_m6_mitm --example n37_native_m6_replay
gzip -dc research/notes/ecc2k130/n37_native_m6_mitm_20261002/RAW.json.gz > /tmp/n37-native-frozen-raw.json
/usr/bin/time -p target/release/examples/n37_native_m6_mitm /tmp/n37-native-raw.json
jq -S 'del(.phase_wall_ms)' /tmp/n37-native-raw.json > /tmp/n37-native-counts.json
jq -S 'del(.phase_wall_ms)' /tmp/n37-native-frozen-raw.json > /tmp/n37-native-frozen-counts.json
cmp /tmp/n37-native-counts.json /tmp/n37-native-frozen-counts.json
target/release/examples/n37_native_m6_replay /tmp/n37-native-raw.json /tmp/n37-native-replay.json
jq -S 'del(.raw_result_sha256)' /tmp/n37-native-replay.json > /tmp/n37-native-replay-normalized.json
jq -S 'del(.raw_result_sha256)' research/notes/ecc2k130/n37_native_m6_mitm_20261002/REPLAY.json > /tmp/n37-native-frozen-replay.json
cmp /tmp/n37-native-replay-normalized.json /tmp/n37-native-frozen-replay.json
```

The normalized comparisons omit only wall-time fields and the receipt's
hash of the timed raw file. The exact archived files and hashes above are
the replay target for this one run. A second local execution from the
committed source passed both normalized byte-for-byte comparisons.

**Next evidence-ranked direction.** First preregister a bounded,
target-blind residual schedule on a *new disjoint Q block*: for an initially
unsupported Q, query `[a]G+Q` at fixed offsets, charge every failed oracle
call, and derive `log Q` only from a verified witness and known `a`.
That directly addresses the measured 107/1024 support gap. Then compare its
cold complete n37 single/batch CPU with same-Q strong batch rho and the
source, transported and pullback policies; carry only surviving policies to
n41/n53. A larger K, arity seven, or a symbolic solver needs a separately
charged admission gate because increasing formal capacity alone has not
produced complete target coverage here. This n37 result supplies no n131
transfer or Certicom attack claim.

Class: **scoped support no-go / accounting**. End-to-end crossover: **unset**.
