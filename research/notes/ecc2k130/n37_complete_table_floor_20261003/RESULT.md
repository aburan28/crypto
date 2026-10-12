# n37 complete three-summand table: cold setup floor

Decision: **`REJECT_COMPLETE_TABLE_COLD_N37` in charged group-addition
equivalents for the frozen 42-column, complete pair-table oracle.** The
preregistered [protocol](PROTOCOL.md) was committed as
`d40917c016643dee76ed926db5432de75dc7f154` and [PR #1280](https://github.com/aburan28/crypto/pull/1280)
opened before these runs. The independent Rust counter replayed every
nondecreasing pair of each policy through `CountedGroup`. Two runs produced
byte-identical receipts. The source, transported, descendant-native and
pullback policies each required **4,831,386** pair additions before any
relation or target solve, `S_floor = 318.155256` on the n37 prime subgroup.

Eight distinct new public `Q` were solved by the native strong signed-
Frobenius rho reference in two measured rounds. All 16 measured runs, and
all eight warmups, recovered a scalar independently checked by `[d]G = Q`.
`ecbench verify --replay-all` reproduced all 16 measured records exactly.
The largest rho charge was **5,652.494 GAE** (`S = 0.372227`), so even that
run used **854.735 times fewer charged GAE** than constructing one frozen
complete pair table. This exceeds the preregistered 100-fold threshold on
every run. No table setup, relation collection, rank solve, isogeny or
descent cost was added to the IC floor; those omissions make it lower.

| n37 arm or boundary | One-target charged GAE | `S = GAE/√r` | Scope |
| --- | ---: | ---: | --- |
| Generic rho floor, signed Frobenius `A=74` | about 2,212.5 | 0.145695 | Analytic average boundary, not a run |
| Strong rho, 16 public measured runs | mean 3,970.431; max 5,652.494 | mean 0.261460; max 0.372227 | Verified whole rho solve charge; unpriced native work remains |
| Original complete-table setup | at least 4,831,386 | at least 318.155256 | Mandatory pair enumeration only |
| Transported complete-table setup | at least 4,831,386 | at least 318.155256 | Same count, leaf arithmetic |
| Descendant-native complete-table setup | at least 4,831,386 | at least 318.155256 | Same count, before scalar orbit closure |
| Pullback complete-table setup | at least 4,831,386 | at least 318.155256 | Same count, before exact inverses |

The per-run table-floor/rho-GAE ratios range from **854.735 to 2,373.569**.
Every row, target coordinate, method ID, run ID, rho cost, isolation grade,
unpriced counter and ratio is retained in [DECISION.json](DECISION.json).
The counter's ordered pair-sum digests are in the two identical
[COUNT_A.json](COUNT_A.json) and [COUNT_B.json](COUNT_B.json). The four
digests differ because each policy has its own point set or curve
representation; each has the same counted pair cost. The frozen support
manifest was checked at SHA-256
`8eac2ae4b8d8f0fc452b7cd7c8edbd3183558e55454abbd67349665bd64610d5`.

The sealed [rho session](sessions/strong-rho-e1/) has ID
`ECBS1h509a4a823bd8`, spec ID `ECS1h999d25a394f1`, and records SHA-256
`4678b123907d33af3eb4964841347ff6a317db7f5383a78f278bc8cab4070932`.
It was recorded with binary SHA-256
`b9fe48bdc8ec4b02e8c59832131b2bf87873c9abe1afc254cf0f0203262cf033`
on a macOS aarch64 host (environment class `ECBENV2hd82681268e96`).
The [audit receipt](AUDIT.json) has SHA-256
`346345851722800f53d3448b85b8a0cebf7d9e3adf6efc1b1813e87ca1e4a273`
and zero problems. The two count receipts each have SHA-256
`f895700386ef5280fba175699c0f257a0702f377a779bdc00b4ad6b39c67a3c8`;
the [decision receipt](DECISION.json) has SHA-256
`1521ba98667186096cf3846f504a46edab4bd8e0db3752e4984e7f06d5b03598`.
This supersedes the original receipt SHA-256
`448be4d380abfd8387f60e80fa0edf1183fe3b07bf46a64da3352542bf7e992a`
after commit `f11fe93ee75a8901873108025eaf79875bd215a9` enabled
`serde_json/float_roundtrip`: one `rho_s` field now preserves the exact value
`0.37222665730402227` already present in the sealed source record rather than
the prior one-ULP-lower parse. No count, charge, ratio at reported precision,
or decision changed.
The producer and analyzer source SHA-256 digests are respectively
`9021b123eaccfd20272d29842a946fbbf8e0a880cbb4d3cad8ff931ff548ba7c`
and `163165c09619035d90c884b12c47a832784cf9eec6d915303d285b3a5023918b`.
Their measured-host release binaries have SHA-256 digests
`40898817d2bab1b0de61ef46ee6b866a28ef6515ac673a5acdc77e1c0cf7b052`
and `be2c616559e11ce67f7fd62242144d7b13f31a30437c0dbed5128f215e442cc5`,
respectively.
The PR's Linux replay recomputes all four pair-sum streams and replays all
16 measured rho runs on another host; the committed receipt is not silently
treated as an independent-host certificate before that check passes.

This is **not an end-to-end IC/rho ratio, wall-time result, or general
index-calculus no-go**. `ecbench` charges strong rho's scalar
multiplications at `1.5 log₂(r)` additions each, and records
canonicalizations, hashes and table work as unpriced. The table counter
charges its exact add calls but does not price field-level work, memory,
sorting or other IC stages. All session records earned L0 on macOS; no
wall-clock figure is admitted. A batched implementation can amortize the
table and requires its own matched batched-rho experiment. Against the
observed mean independent single-target rho charge, the table alone would
need at least about 1,217 targets merely to amortize below that rho
charge per target, before IC rank or descent and before any rho batch reuse.
That is a prioritization bound, not a measured batch crossover.

For the primary one-target route, stop investing in the **eager complete
pair table** at this fixed support size. The next evidence-ranked direction
is a target-demanded or algebraic three-summand oracle whose cold setup
does not enumerate `Θ(P²)` pairs, with complete relation/rank and descent
costs charged in native `ecbench` on fresh public targets. Continue the
equal-useful-support n41/n53 tests with fresh maps and held-out targets;
the degree-73 n37 kernel cannot simply be reused in those fields. The
degree-263 endomorphism-order change remains a separate n131 transfer
hypothesis, with no measured PDP or attack-speed benefit established here.
The newly merged [n31 cofactor-balance predictor](../../koblitz-isogeny/yield-spread-results-20261003.md)
does not separate these n37 supports: every point in all four fixed bases
was verified in the order-`r` subgroup, so its cofactor component `[r]F`
is zero for every factor. Its useful follow-on is where candidate bases
include different cofactor cosets, not this equal-subgroup-support panel.

Reproduce from a checkout with the frozen support files and lockfile:

```sh
cp research/notes/ecc2k130/n37_native_m6_mitm_20261002/Cargo.lock Cargo.lock
cargo build --release --locked --example n37_complete_table_floor --example n37_complete_table_floor_analyze --bin ecbench
note=research/notes/ecc2k130/n37_complete_table_floor_20261003
target/release/examples/n37_complete_table_floor /tmp/n37-floor-replay.json
cmp /tmp/n37-floor-replay.json "$note/COUNT_A.json"
target/release/ecbench verify --dir "$note/sessions/strong-rho-e1" --replay-all --exit-code --out /tmp/n37-rho-replay.json
target/release/examples/n37_complete_table_floor_analyze "$note/COUNT_A.json" "$note/COUNT_B.json" "$note/sessions/strong-rho-e1" "$note/AUDIT.json" /tmp/n37-decision-replay.json
cmp /tmp/n37-decision-replay.json "$note/DECISION.json"
```

The count and analyzer refuse to overwrite an output path; choose unused
paths for a repeat. To rerun the rho producer rather than replay the
committed session, use a new `--out` directory with the frozen
[`rho-public-spec.json`](rho-public-spec.json); never edit the sealed
session.
