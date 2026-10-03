# n37 four-policy natural-target PDP and rank gate

Decision: **`M3_SATURATED_NONDISCRIMINATING`** on the two frozen point-only
blocks. The [protocol](PROTOCOL.md) was committed as
`10043b865919e206291f0780ec987cf4e514cec0` and its PR opened before the
new measurement. Both seed selections gave exact at-most-three-summand
witnesses for **all 2,048** public targets. Both reached full 42-column rank
and recovered and verified all 2,048 logarithms. Therefore these n37 blocks
cannot distinguish three-summand yield between the fixed equal-useful-size
bases. This is a completed exact-PDP and rank result, **not** a cold IC/rho
speed comparison or an ECC2K-130 attack-speed claim.

Every policy uses the same 42 log columns and 1,554 signed support classes
(3,108 physical points) admitted by the earlier
[support gate](../n37_four_policy_support_20261003/RESULT.md). Each complete
table enumerated the identity, 3,108 singletons and 4,831,386 nondecreasing
pairs. A proved two- or three-summand miss required the complete table and
all 3,109 residual positions; no timeout or unknown was counted as a miss.
Original↔transported and descendant-native↔pullback had identical canonical
witness indices, hit flags, rows, rank trajectories and recovered scalars,
as required by the bijective degree-73 transport on this prime subgroup.

| Fixed policy | Distinct table sums | m≤2 hits / 2,048 | m≤3 hits / 2,048 | Probes to rank 42 | Dependent rows | Verified logs / 2,048 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Original source | 4,773,667 | 43 | 2,048 | 67 | 25 | 2,048 |
| Transported leaf | 4,773,667 | 43 | 2,048 | 67 | 25 | 2,048 |
| Descendant-native leaf | 4,778,329 | 51 | 2,048 | 98 | 56 | 2,048 |
| Pullback source | 4,778,329 | 51 | 2,048 | 98 | 56 | 2,048 |

The two-summand control found **43 source-only** and **51 native-only**
targets, with no overlap on these blocks. Its exact conditional two-sided
McNemar p-value is
`9319056435053552377244630248 / 2^94 ≈ 0.470492`;
the 8/2,048 absolute hit difference fails the frozen 21/2,048 and `p<0.01`
selection-lead rule. The uniform-population two-summand formal-count ceiling
is about 2.09646%; the observed 43 and 51 hits in a finite fixed sample do
not violate that population bound. All three-summand hit flags agree, so the
corresponding exact p-value is one. The formal three-multiset count is
`5,008,536,820 ≈ 21.72r`; that counting capacity does not prove complete
group coverage, and this experiment tested only the 2,048 frozen targets.

The source-selected base reached rank 42 after 67 fixed target-blind
SplitMix64 probes and the descendant-native base after 98. All relation
witnesses and modular rows verified. The difference is a **single-stream
rank diagnostic**, not a population estimate or a cost advantage. A fresh
preregistered rank-stream panel would be needed to decide whether it is
stable. Every source and leaf target logarithm from each arm was checked
by scalar multiplication in the producer. The independent replay rebuilt
complete hash tables for the two source supports, independently reduced
every rank row, verified each source witness with the general binary-curve
group law, checked every reported source logarithm with the general scalar
law, and checked the two paired leaf arms through the full-point isogeny.

The source-matched final raw manifest is
[`RESULT_CI_FIXED.json.gz`](RESULT_CI_FIXED.json.gz), 222,574 bytes with
SHA-256 `c84066a521137cfa61a2f490ed18a3a957b79da6e27f587ad1579b350cb9764a`.
It expands to 5,104,424 bytes with SHA-256
`2c8f2cc0767e177d21b3a0c62fbd50a40992150d1dd3bee7c07f9e629c647ea1`.
The independent [`REPLAY_CI_FIXED.json`](REPLAY_CI_FIXED.json) has SHA-256
`7b79fd455401c36ffaa4b115593f6c820af7e43de58b886d3ba6ce71e3d9f2cd`.
CI exposed unused shared helper items under strict Clippy. The source-only
revision scoped a dead-code allowance to these standalone examples and removed
one needless borrow; it changed no search or decision logic. The prior PASS
[`RESULT.json.gz`](RESULT.json.gz) and
[`RESULT_INITIAL.json.gz`](RESULT_INITIAL.json.gz) are retained with their
own replay receipts. All three canonical non-timing JSON projections have
the same SHA-256,
`5f8f00f7e9ca5ddbed87d88e0a1d04a55b11436bf52e52d6cae6ba798dd07b00`.
The initial producer source digest was not retained, so the source-matched
post-lint run is decisive. All input, output, code, mutation and lockfile hashes are in
[`EVIDENCE.json`](EVIDENCE.json).

Replay failed as required when a scratch result changed one witness index,
even after changing its paired isogenous row too; when a paired rank row or
hit flag changed; and when the reported decision changed. Temporarily
changing the first b03 public point caused the pinned input SHA-256 check to
fail, and the original point file was restored with its frozen hash. The
failure receipts are retained beside the PASS receipts.

Each policy's table construction used 4,831,386 full-point additions. The
source-selected table retained at least 77,401,664 bytes in its point/code
arrays; the post-lint producer's peak RSS before JSON output was 135,331,840
bytes. On this [L0 macOS host](HOST.json), its four-arm diagnostic took
15.32 descriptive seconds, including the intentionally expensive two-summand control.
The source-selected m≤3 queries made 99,917 table lookups across 2,048
targets; native-selected queries made 102,767. These are stage counts and
descriptive wall times. The archived support excludes a fresh source-base
scan, leaf-seed selection, scalar orbit closure, exact pullbacks, and matched
rho; this gate did rebuild the degree-73 map during input validation.
**Complete `S`, rho ratio and speedup remain
null.** A single target's online cost and a batch's amortized setup cannot
be inferred from this mixed control run.

The next yield experiment should move to a larger Koblitz subgroup (n41,
then n53) with equal useful support and newly frozen targets; repeating seed
selection on b03/b04 would contaminate the comparison. Separately, the
native `[lambda]` closure must be priced in fresh cold ecbench sessions,
since its 1,512 scalar actions may outweigh any relation benefit. The
degree-263 endomorphism-order change remains a hypothesis about action and
solver cost; this n37 saturation does not establish a degree-263 yield or
an n131 transfer.

Reproduce the decisive replay from a checkout with the frozen inputs:

```sh
cp research/notes/ecc2k130/n37_native_m6_mitm_20261002/Cargo.lock Cargo.lock
cargo build --release --locked --example n37_four_policy_pdp --example n37_four_policy_pdp_replay
gzip -cd research/notes/ecc2k130/n37_four_policy_pdp_20261003/RESULT_CI_FIXED.json.gz \
  > /tmp/n37-four-policy-pdp-result.json
printf '%s  %s\n' \
  '2c8f2cc0767e177d21b3a0c62fbd50a40992150d1dd3bee7c07f9e629c647ea1' \
  /tmp/n37-four-policy-pdp-result.json | shasum -a 256 -c
target/release/examples/n37_four_policy_pdp_replay \
  /tmp/n37-four-policy-pdp-result.json /tmp/n37-four-policy-pdp-replay-fresh.json
cmp /tmp/n37-four-policy-pdp-replay-fresh.json \
  research/notes/ecc2k130/n37_four_policy_pdp_20261003/REPLAY_CI_FIXED.json
```

The [CI workflow](../../../../.github/workflows/n37-four-policy-pdp.yml)
performs this on Linux. To rerun the producer, give it a new output path;
it refuses overwrites. It reads no fixture scalar.
