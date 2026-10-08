# Frozen n37 descendant-native residual search: all b01 logs recovered

**Decision: `N37_NATIVE_K42_M6_BOUNDED_RESIDUAL_BATCH_CORRECTNESS_PASS`.**
On the preregistered, orbit-disjoint b01 public-Q block, the fixed 42-column
descendant-native signed six-sum oracle recovered and independently verified
**1024/1024** logarithms. The direct query supported 939 targets; the 85
direct misses were all resolved by the frozen target-independent shifts.
Seventy-nine needed one shift and six needed two. No target used shifts 3–16.
This is a correctness and finite-block support result. It does **not**
establish a method speedup, a single-target IC/rho crossover, transfer to
n41/n53/n131, or an ECC2K-130 attack.

The [protocol](PROTOCOL.md) was committed as `b0bc1899` before this Q block
was queried. The measured native producer and independent replay were both
built from source commit `880c0eefe45a68b0dfd9ead8fb15930e484096ae`,
with the exact `Cargo.lock` preserved in the prerequisite
[n37 native m6 evidence](../n37_native_m6_mitm_20261002/Cargo.lock)
(SHA-256 `b28d3c2d81146a40d00df85f45c6460d9e2bc5c26875307a06cfd15efff99365`).
The solver read only the b01 public point file. The separate published
fixture scalars were read afterward by the replay, which checked all 1024.

| Frozen arm / Q block | Public Q | Direct support | Shifted recovery | Verified logs | Unresolved | Decision |
|:--|--:|--:|--:|--:|--:|:--|
| Prior b00 bare direct 3+3, PR #1254 | 1024 | 917 | not attempted | 917 | 107 | Bare direct batch no-go. |
| b01 direct query under this protocol | 1024 | 939 | not yet attempted | 939 | 85 | Direct alone still incomplete. |
| b01 direct plus at most 16 frozen shifts | 1024 | 939 | 85 | 1024 | 0 | Bounded residual batch passes this block. |

The b00 and b01 rows use different, frozen Q blocks. Their wall times and
per-target support counts are not paired performance observations. The b01
run rebuilt the table and relation system cold: 105,995 raw signed half
entries collapsed to 102,391 distinct half sums; 49 target-blind relation
probes yielded rank 42, with seven complete misses. The shift list was
generated once from frozen SplitMix64 seed `0x6e33375f72657331`; all 16
`[a_j]G_leaf` points were precomputed and charged, even though only the
first two shifts were needed. The observed target attempt histogram was
939 with one query, 79 with two, and six with three, for **1,115** total
oracle queries. Every answer passed full-point verification on both source
and leaf curves.

[`RAW.json.gz`](RAW.json.gz) preserves every input digest, selected shift,
relation, rank step, target attempt, hit, complete miss, witness, scalar,
operation count, and phase timer. Its deterministic gzip archive is 105,650
bytes, SHA-256
`e153c41a49833eeac8d75a3ca177c8328723679055cfb5744d27c99a84f8be1d`;
it expands to 1,798,008 bytes, SHA-256
`0ef64b4f965cc316ae2abd59fbdcb949655321bf9a584297dec0aa1bfca7263c`.
The independent [`REPLAY.json`](REPLAY.json) is 1,188 bytes, SHA-256
`fd06ba7b4a1cdeacbe2d2b5563064f6adba80ceb828664b7c60083ed0c00bd9e`.
The replay built source-curve reachable sets by iterative group-law closure,
checked complete 3+3 membership for every direct and shifted query, verified
every reported source-pullback witness, solved base logs independently,
checked every recovered Q by source scalar multiplication, and then compared
to the published fixture. Its reachable set sizes through three factors
were `[1,85,3613,102391]`; all support decisions and aggregate counts
matched the producer.

| Verified input | SHA-256 |
|:--|:--|
| Frozen native base | `bd4bd8af982bcc65234ae1fae1ad52eb0ed41a1c203bf4d2556fd8e87f788a4c` |
| n37/L1024 frozen specification | `da958a3f1117dd1b88703fed2055c5c9f64255516e1e920f35d193bf53c6a88d` |
| Public b01 points | `78553fdff5ae66521d3a1052978962258e78d48c92285df6e1df962034b43717` |
| Degree-73 archive | `eb4773d556886672e8b80735fc486c16b73765fd5069ff81133100bfc3eafa90` |
| b01 published scalar fixture, replay only | `3453995c0423e4911ad4a6afa7cfe50ec1d6b0e6c67a771bda590f1d6ed680df` |

This was **one diagnostic process** on macOS 26.6 arm64 with Rust 1.93.1.
[`TIMING.txt`](TIMING.txt) records 2.37 s process wall, 1.71 s user,
and 0.01 s system. Instrumented target queries took 1.750 s, but that
stage time is not an attack speedup. The cold process included 26.9 ms of
setup/transport, 58.8 ms of table construction, 79.3 ms of relation
collection, 0.15 ms of dense solve/base verification, and 1.762 s of target
queries and verification, plus output serialization/write. The separate
[`REPLAY_TIMING.txt`](REPLAY_TIMING.txt) records 9.54 s wall for independent
verification; it is not attack work.

The producer counted 105,910 table-build additions, 1,507,265 relation
query additions, 29,461,763 direct query additions, 2,559,946 residual
query additions, 91 explicit shift additions, 6,227 external coefficient
replay additions, and 6,227 oracle witness additions across relation,
direct, and residual queries. Together these are **33,647,429 logical
addition requests** under the producer's counting convention; 4,289
explicit scalar multiplications are separately recorded. Scalar
multiplication internals and field-basis/isogeny arithmetic are included in
wall time but not converted to group-addition equivalents. The matched
same-host rho cost, common-unit `S`, and speedup remain **unset**.

Reproduce from this PR head (measured implementation at `880c0eef`):

```sh
cp research/notes/ecc2k130/n37_native_m6_mitm_20261002/Cargo.lock Cargo.lock
cargo check --locked --offline --example n37_native_m6_residual --example n37_native_m6_residual_replay
cargo build --release --locked --offline --example n37_native_m6_residual --example n37_native_m6_residual_replay
gzip -dc research/notes/ecc2k130/n37_native_residual_20261002/RAW.json.gz > /tmp/n37-residual-frozen-raw.json
/usr/bin/time -p target/release/examples/n37_native_m6_residual /tmp/n37-residual-raw.json
cmp <(jq -S 'del(.phase_wall_ms)' /tmp/n37-residual-raw.json) <(jq -S 'del(.phase_wall_ms)' /tmp/n37-residual-frozen-raw.json)
target/release/examples/n37_native_m6_residual_replay /tmp/n37-residual-raw.json /tmp/n37-residual-replay.json
cmp <(jq -S 'del(.raw_result_sha256)' /tmp/n37-residual-replay.json) <(jq -S 'del(.raw_result_sha256)' research/notes/ecc2k130/n37_native_residual_20261002/REPLAY.json)
```

The normalized comparisons omit only the producer's phase wall times and
the receipt's hash of that timed raw file. A second local producer run from
the committed source passed the normalized byte-for-byte comparison; its
independent replay passed the corresponding normalized receipt comparison.

**Next gate.** Measure this complete native b01 batch against a cold,
same-host, same-Q signed-Frobenius batched rho run with a frozen resource
envelope, A/A noise control, and interleaved repeats. Even a favorable
batch result would remain separate from the repository's required
one-target same-Q online IC/rho comparison. Then compare equal-useful-size
source, transported, descendant-native, and pullback policies, and carry
only surviving, fully charged approaches to n41/n53. No n131 transfer can
be inferred from this result alone.
