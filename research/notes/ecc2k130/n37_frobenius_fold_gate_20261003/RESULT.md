# n37 source-support Frobenius-folded pair-table gate

Decision: **`FOLD_EXACTNESS_AND_SETUP_ADMITTED_N37_SOURCE`** at the frozen
support and point-only regression sample. The [protocol](PROTOCOL.md) was
committed as `1354075ad247e4652ef7c16b08cd10c58bad2819` and its
[PR #1281](https://github.com/aburan28/crypto/pull/1281) opened before
measurement. Both fixed source supports built a signed-Frobenius-folded
table with **66,822 counted group-addition requests**, as preregistered,
versus **4,831,386** for the complete nondecreasing pair table: **72.3023×
fewer setup additions**. The native producer and a separate general-law
replay agreed on all **8,192** m≤2/m≤3 hit/miss decisions over the two
2,048-target arms, the unique table sizes, every query count and every
reported witness and relation row. The independent Linux CI replay is
pending; this local result is not yet a merged claim.

| Frozen source support | Folded unique sums | Table additions | m≤2 hits / 2,048 | m≤3 hits / 2,048 | m≤3 query additions | Minimum retained table bytes |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Original | 64,467 | 66,822 | 43 | 2,048 | 99,917 | 1,491,456 |
| Pullback | 64,530 | 66,822 | 51 | 2,048 | 102,767 | 1,492,464 |

Each table also performed 66,780 nonidentity canonicalisations, 111,888
Frobenius maps and 111,888 index lookups while building its orbit map.
Those native operations are **counted but unpriced**. The m≤3 query path
performed another 99,917 or 102,767 canonicalisations, respectively, and
101,965 or 104,815 lookups including direct-base probes. The retained-byte
number is a strict lower bound for the explicit orbit-index array and the
logical map entries; hash-table capacity, allocator overhead and the frozen
support are excluded. The source's 64,467 folded entries compare with
4,773,667 distinct full-point table sums in the earlier complete-table
gate; the pullback has 64,530 versus 4,778,329. These are representation
counts, not measured RAM reduction factors.

The separate replay rebuilt the 66,822 pairs on each support using the
general binary-curve group law and keyed each nonidentity sum by the least
of its 37 Frobenius-conjugate abscissae. This is independent of the
producer's normal-basis rotation key and `FrobeniusPairTable` implementation.
It recomputed each m≤2 and m≤3 query, checked counters, compared decisions
with the independently replayed complete-table reference, and verified
every reported witness's full-point sum and frozen `(column, coefficient)`
row. The m≤2 controls include 2,005 original and 1,997 pullback misses;
the m≤3 sample is saturated. Agreement on these historical targets is a
correctness regression, not a new target-population yield estimate.

Two local producer runs had identical non-timing canonical JSON, SHA-256
`3cfb28e21dfcb32a252b32bb96619a4e72bb0a5997a4128e019fc0afde494e6b`.
Raw run A is [`RESULT_A.json.gz`](RESULT_A.json.gz), 184,684 bytes,
gzip SHA-256 `b1ba2093452ae74effd152ad6f15af2ce3957828820213918b7c896c42175c5c`,
expanded SHA-256 `f1bcf4fcea33a90b6b2d9d3becf2eab96c16484a9114a06d03137c7f5ba50d75`.
Raw run B is [`RESULT_B.json.gz`](RESULT_B.json.gz), 184,687 bytes,
gzip SHA-256 `7ff327b23e4d84c703e10386049047a84f6856b04c74e656b0f19e1a73d19b50`,
expanded SHA-256 `251aee18844a30765f02049c87c858bfd0549fb28941ee12ac04e3814040e83e`.
The [A](REPLAY_A.json) and [B](REPLAY_B.json) replay receipts each report
PASS. [Evidence hashes](EVIDENCE.json), [host facts](HOST.json) and
[mutation failures](MUTATIONS.json) accompany the runs. Changing a witness
index, hit flag, setup count or support coordinate caused replay to fail.

On this unisolated macOS arm64 L0 host, table construction took 30.1–37.8 ms
across the two runs and supports. Those wall observations describe this
implementation and host only. The support arrays were previously selected;
the cold source-base scan, leaf-seed selection, isogeny transport, exact
pullbacks, relation collection, rank, target-log recovery, table memory
overhead, native canonicalisation and field work are not fully charged.
**Complete `S`, matched-rho ratio and end-to-end speedup remain null.** In
particular, this gate does not establish a single-target or batched
ECC2K-130 advantage, a degree-263 descendant yield, or n131 transfer.

The next admitted experiment is a fixed-input `ic.pipeline` plug-in with a
fresh point-only Q panel, full-rank verified recovery and same-Q
signed-Frobenius rho. Its frozen-base input and every unpriced native count
must be explicit. A later cold native constructor must charge source-base
selection, isogeny work and the alternative descendant-native scalar
action before comparing policies or declaring a crossover. The large
setup-addition reduction makes batch amortisation worth testing, but the
earlier complete-table failure and this scoped success do not settle full
IC/rho economics.

Reproduce both raw runs and the independent replay on a checkout with the
frozen lockfile:

```sh
cp research/notes/ecc2k130/n37_native_m6_mitm_20261002/Cargo.lock Cargo.lock
cargo build --release --locked --example n37_frobenius_fold_gate --example n37_frobenius_fold_replay
gzip -cd research/notes/ecc2k130/n37_frobenius_fold_gate_20261003/RESULT_A.json.gz > /tmp/n37-fold-A.json
target/release/examples/n37_frobenius_fold_replay /tmp/n37-fold-A.json /tmp/n37-fold-A-replay.json
cmp /tmp/n37-fold-A-replay.json research/notes/ecc2k130/n37_frobenius_fold_gate_20261003/REPLAY_A.json
```

The [Linux workflow](../../../../.github/workflows/n37-frobenius-fold-gate.yml)
also regenerates the non-timing evidence and replays both preserved runs.
