# P-256 low-width Dickson-union screen, round 295: result

## Verdict

The screen constructs an actual complete Dickson factor base at Round 294's
theoretical minimum fixed-arity width.

`FB1hc72514a2a8d3` is the exact union of a 136-column depth-8 fibre and a
28-column depth-6 fibre.  Its 164 columns are distinct, all 328 signed points
replay, and the canonical `ecbench.factor_base_dump/v1-wide` artifact is
byte-identical on independent rebuild.  Fixed arities 109 and 110 cover the
P-256 subgroup.  The exact favourable collision boundary is **13.920747 times
rho**, improving the smallest previous complete covering Dickson fibre's
17.734068 boundary by **1.273931 times**.

This reaches the best width permitted by the exact coefficient-domain count,
but it does not reach rho.  The two component chains retain local generator
degree two; the degree of regularity of a global 109/110-term union system is
unset.  No relation solver was evaluated, zero relations are reported, the
relation false-positive/false-negative gate is unset, and no unplanted P-256
relation was attempted.  The factor base is retained; promotion fails.

## One boundary, one unit

The unit is favourable P-256 group-addition equivalents divided by
`sqrt(n)`, with rho fixed at `S=1.3`.  Every candidate receives a perfect free
decomposition oracle, success probability one, one operation per collision
sample, and no setup, solving, replay, matrix, recovery, or representation
cost.  These are generic-collection lower boundaries, not end-to-end costs.

| variant | complete native base | columns `B=K` | first covering arity | exact collision / rho | balanced two-list / rho | result |
|:--|:--:|--:|--:|--:|--:|:--|
| Pollard rho | yes | - | - | **1.000000** | **1.000000** | reference |
| Round-294 theoretical minimum | no inventory | 164 | 109 | 13.920747 | 19.701921 | counting floor |
| prior depth-9 `FB1h91cc12460e5a` | yes | 266 | 59 | 17.734068 | 25.091548 | superseded native minimum |
| **Round-295 `FB1hc72514a2a8d3`** | **yes** | **164** | **109** | **13.920747** | **19.701921** | **factor-base advance; promotion rejected** |
| registered depth-18 `FB1h2f8621cda105` | yes | 131,458 | 17 | 394.425280 | 557.802111 | comparison |

The ratio improvement from the prior native minimum is
`17.73406835956266 / 13.920747397073491 = 1.273930763465389`.  No other
complete independent-log base can improve this boundary by reducing `B`:
Round 294 proves that every `B<164` fixed-arity domain is smaller than the
P-256 subgroup.

## Frozen screen and selected FB1

The native runner built and exactly rebuilt 256 hash-derived depth-8 fibres
and 256 hash-derived depth-6 fibres, then evaluated all 65,536 pairs.  All 512
components passed; none was rejected.  Pair widths range from 123 through
196.  Exactly 23,905 pairs cover in at least one fixed arity, 2,275 have width
164, and every pair has zero duplicate columns.

The frozen ranking selects:

| component | index | FB1 | columns | root exponent | terminal |
|:--|--:|:--|--:|:--|:--|
| depth 8 | 0 | `FB1h061e868987d2` | 136 | `0x3d31fded66eb449bdba975` | `0xd280a56fdc5f82bb38bb4b9f0646a63f916aec63586aaac640427ea9e380d0e8` |
| depth 6 | 9 | `FB1h0eb207d95d5e` | 28 | `0x3e2334694ac173022f97679` | `0x42ab26df6d850f09c48b015c21c81fc22c4d0f567a7619fb169a1b6a13bbacd1` |

Their canonical union identity is:

```text
FB1hc72514a2a8d3
FB1 SHA-256      c72514a2a8d3641918a319e46060d7ffbb981376a7be53319b5d5dc99a9c15d5
point SHA-256    f7b8491caa0848163832e34f5febbd4c9c1a18332a99c29f14d45772cf28e899
columns / points 164 / 328
```

The checked-in factor-base artifact follows the repository's wide dump schema,
uses family `dickson-torus-union`, and carries both ordered component FB1s,
depths, indices, root exponents, terminals, and the standard
`prime-affine-x-plus-one-shift-sign/be33/v1` point-key encoding.

## Exactness, degree, and gates

The screen performs 1,028 component builds including verification and winner
replay, and checks 82,356 signed point rows.  The independent complete run is
byte-identical for both result JSON and selected factor-base JSON.  Exact
binomial recurrence covers all 73 observed widths with 46,862 counted
big-integer operations and reproduces `sum 2^m C(B,m)=3^B` at each width.

| gate | status | evidence |
|:--|:--:|:--|
| dependency and component identities | pass | Round-294 hash checked; 512/512 components rebuilt |
| exact width 164 realized | pass | 2,275 frozen pairs; selected union 136+28 with zero overlap |
| native boundary below depth-9 control | pass | 13.920747 vs 17.734068 times rho |
| boundary at or below rho | **fail** | favourable lower boundary is 13.920747 times rho |
| structured residual degree at most 5 | **fail / unset** | component generator degree 2 is not global solving degree |
| zero relation false positives and false negatives | **fail / unset** | no relation selector was run |
| usable relation below `2^103` | **fail / unset** | no complete global solver |
| complete collection below `2^120` | **fail / unset** | independent-log lower boundary already misses rho |
| projected materialized storage below `2^50` | **fail / unset** | 83,981-byte FB1 dump does not price a 109/110-term solver |
| promoted | **no** | zero relations; no complete DLP |

The exact factor-base-width route is now closed at its optimum.  A successor
cannot gain another factor by selecting fewer independent columns.  It must
instead produce multiple independent rows per non-generic global solve, or
prove exact logarithm transport that lowers `K`; Rounds 21--29 remain the
control against disguising scalar-orbit rho as such transport.

## Resources, artifacts, and reproduction

The canonical isolated run used reserved CPU 4 on the AMD EPYC 9V74 host.  It
completed successfully and uncontended in 7.067164 s wall, 7.040940 s user and
0.024070 s system, peaking at 35,268 KiB RSS.  Both successful isolated runs
are retained.

- canonical result JSON: 431,043 bytes, SHA-256
  `9550f2bbaceb9e297fb35c12e79480581ca80bd58b03a8493e86ea7c1eda91a5`;
- canonical selected FB1 wide dump: 83,981 bytes, SHA-256
  `d27516ca40a612ecf3ebabfa8ae04776084c7948110e20f221da438e1e64d8f8`;
- two-run isolation JSONL: 4,255 bytes, SHA-256
  `b99b222ca682491eae9dfed7aa611f8aab44f7bcd3352600f22624ebc2da1d1b`;
- semantic evidence SHA-256:
  `a1e6b5bae14680d23d69e0af3b2cb6a857b362bb91a60e27cec22a4a137fe760`.

```bash
cargo test --release --bin p256_dickson_union_screen
cargo clippy --release --bin p256_dickson_union_screen -- -D warnings
cargo build --release --bin p256_dickson_union_screen --bin isolated_bench
target/release/isolated_bench run --wait --cpus 4 \
  --out research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_dickson_union_round295_20261007/isolation.jsonl \
  --label p256-dickson-union-round295-canonical -- \
  target/release/p256_dickson_union_screen \
  --round294 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_variable_arity_round294_20261007/variable-arity-result.json \
  --out research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_dickson_union_round295_20261007/dickson-union-result.json \
  --factor-base-out research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_dickson_union_round295_20261007/selected-factor-base.json
```
