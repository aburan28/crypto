# P-256 decorrelated hash-jump factor-base screen, round 293: result

## Verdict

Round 293 reaches the registered selector-stage boundary but does **not**
achieve index-calculus parity with Pollard rho.

The best candidate, `LDHJR4093h94827a2592b7`, exactly covers all 48,889,856
states emitted by its complete `2^17` sign orbit.  Its ideal-local ratio is
`0.966922` rho and its exact-delta-frequency-corrected selector lower bound is
`0.997995` rho.  This is the first exact full-coverage row below one in this
P-256 selector series.

The margin is only `0.200466%`, and the construction obtains it by replacing
the two-delta algebraic band with 4,093 hash-derived jump deltas.  The resulting
131,458 x-coordinates are an explicit random set whose exact membership
polynomial has degree at least 131,458; no same-family structured residual
degree at most five exists.  The monotone path itself advances by the constant
translation `+G`, the generic transition already classified by Round 29.
Relation solving, collection, rank, sparse linear algebra and scalar recovery
remain unpriced.  The candidate therefore fails promotion and no unplanted
full-depth relation was attempted.

## One boundary, one unit

The table uses charged P-256 addition equivalents per exact distinct selector
state, normalized to rho.  It is a construction-stage lower bound, not a cold
one-target DLP ratio.

| variant | rare delta classes | emitted states | exact distinct | distinct fraction | ideal-local / rho | delta-corrected / rho | algebraic reading | result |
|:--|--:|--:|--:|--:|--:|--:|:--|:--|
| Pollard rho | generic | complete method | complete method | 1.000000 | **1.000000** | **1.000000** | complete-method boundary | reference |
| Round 189 measured two-delta | 1 | 29,753,344 | 13,096,682 | 0.440175 | 2.200455 | 2.320918 | arithmetic-band collisions | rejected |
| Round 189, capacity-226 optimistic ceiling | 1 | 29,753,344 | at most 19,714,708 | at most 0.662605 | at least 1.461785 | at least 1.541810 | frozen-capacity certificate | closed |
| **two-delta global capacity sweep** | 1 | best relaxation at C=63: 8,388,608 | at most 8,388,608 | at most 1.000000 | at least 0.979408 | **at least 1.033025** | all optimistic capacities 0–306 | **all tuples closed** |
| `LDHJR1021h062c8464da75` | 1,021 | 200,540,160 | 114,352,128 | 0.570221 | 1.692270 | 1.705516 | incomplete decorrelation | rejected |
| `LDHJR2047hbc24b14a54ad` | 2,047 | 100,663,296 | 76,021,760 | 0.755208 | 1.278578 | 1.298802 | incomplete decorrelation | rejected |
| **`LDHJR4093h94827a2592b7`** | **4,093** | **48,889,856** | **48,889,856** | **1.000000** | **0.966922** | **0.997995** | explicit high-degree set; generic translation | **selector pass, promotion fail** |

The measured selector ratio moves from 2.320918 to 0.997995, but this is
classified as **relabelling**, not an end-to-end advance: the missing work has
been moved into an unstructured membership constraint and a generic
translation process.  Cold one-target `S` and a complete speed ratio remain
unset.

## Tuple-independent closure of the two-delta family

For the frozen Round 189 base with `R=7021`, every coefficient satisfies

```text
R a_j = R alpha + rho_j (mod n),  0 <= rho_j < B.
```

For a selected column of common-run capacity `c_j`, the anchor-free scaled
sign-boundary weight is `v_j=2 rho_j+R c_j`.  Exact replay of all 131,458
coefficients finds zero scaling or endpoint failures and confines the observed
weights to the registered band of width `B+R`.

Within Hamming layer `k`, complement symmetry confines the scaled starts to
width `min(k,17-k)(B+R)`.  Splitting that band by residue modulo `R` gives the
preregistered optimistic ceiling

```text
min(C(17,k)(C+1), R(C + ceil(min(k,17-k)(B+R)/R) + 2)).
```

At Round 189's measured capacity 226, summing all 18 layers while deliberately
ignoring cross-layer overlap gives at most 19,714,708 distinct states, or
66.260478% of emissions.  That fixed-capacity relaxation costs at least
1.541810 rho after the frozen two-delta correction.

The initial receipt applied that number too broadly to all tuple capacities.
The corrected native certificate therefore sweeps every integer capacity from
zero through the absolute 17-slot maximum 306, including unattainable values
to remain optimistic.  Its global minimum occurs at capacity 63: even granting
perfect coverage of all 8,388,608 emissions, the ideal-local ratio is 0.979408
but the delta-corrected floor is **1.033025 rho**.  Consequently no tuple
selector, exhaustive or heuristic, can repair this two-delta base.

The coefficient digest matches Round 189, the cycle closes, and all 131,458
endpoint identities replay.  Certificate digest:
`06f9759a48f6adc3f16457e15706c3518dd8cbd9c98f8dfcb769b13779847ac5`.

## Registered hash-jump construction

All three candidates keep mechanical long runs of `+G`.  Their first `R-1`
rare deltas are deterministic SHA-256 reductions and the final delta exactly
closes the P-256-order coefficient cycle.  Every candidate was accepted on
its first jump and anchor attempt.  All 394,374 coefficients are distinct and
nonidentity, unique up to sign, and all 394,374 edges plus all 394,374
common-run endpoints replay with zero failures.

| candidate | coefficient SHA-256 | rare-delta SHA-256 | tuple draw / cutoff / capacity | tuple SHA-256 |
|:--|:--|:--|:--|:--|
| `LDHJR1021h062c8464da75` | `062c8464da75a330b5aa3ae6a0412ba77599e05f1a8f0c4afee9ae1471f87f13` | `c1585ac72b913d8542f142c6657cfadde85a265f2173ddd0e7d9a453e95d08ed` | 867 / 1,488 / 1,529 | `6f1a373364f1b0df353c8535a6b5b60bcb23a4f5d00be0ee435f11c3a637c1b9` |
| `LDHJR2047hbc24b14a54ad` | `bc24b14a54ada8e9dba3ea3885b6f1504d152bc79c2cad392e68561e0e54a3dd` | `1ce5e1d3c28a10075e1651a74383c76641e172c2fc0e93380c4d55bca260c105` | 366 / 742 / 767 | `f4dd025304a69e4e66f1e5a4a4acd5b1bf28384c3523fa017d758c807168ec4f` |
| `LDHJR4093h94827a2592b7` | `94827a2592b735b1e216fa3a3cab1300bfdf65fc2726644ed375668fafa2bb36` | `a57b419ef79f0c2584d36c9f47417180eae4ab004fe037ca6f7b1db842104cfd` | 10 / 372 / 372 | `87a9f88d9abca32977f9275ad6b9df744c85dfa5c5845fb23a10989ad3647ff7` |

The corrected boundary uses the exact oriented-delta inventory

```text
p_delta = ((B-R)^2 + R) / B^2,
hidden_ratio = 0.964336477130181 / sqrt(p_delta).
```

The respective hidden ratios are 0.971885, 0.979590 and 0.995326.  No rare
delta duplicates the common delta or another registered rare delta.

## Exactness controls and deterministic replay

The two complete toys reproduce 66 distinct states from 80 emissions at
`p=257` and 592 from 672 at `p=65537`.  For every hash-jump candidate, native
direct enumeration at sign depths 8, 10 and 12 uses the same selected columns
with each slot capacity truncated to at most two.  All nine native comparisons
match exact interval union.  Across all eleven controls there are zero false
positives and zero false negatives.

An independent second native execution produced the same two-delta
certificate, candidate identities, tuples, exact state counts, boundary
ratios, support digests, gates and semantic evidence digest.  After removing
timing-only fields, both canonical cores have SHA-256
`5a88d413cc58aa591cd96dff0309c66f2c7e4b5d77d52efdc79ecc7d4400a612`.

## Degree and end-to-end obstruction

Unique-up-to-sign P-256 points have distinct short-Weierstrass x-coordinates.
Thus every accepted hash-jump candidate has 131,458 distinct x-values.  In the
absence of a separately verified defining variety, the exact univariate
membership polynomial for that explicit set has degree at least 131,458.
Hashing the rare jumps supplies statistical-looking subset weights, not a
low-degree relation system.

The best selector therefore fails the requirements that motivated the factor
base:

| promotion gate | status | obstruction |
|:--|:--|:--|
| exact selector coverage and corrected ratio below rho | passed | 1.000000 coverage, 0.997995 corrected lower bound |
| zero FP/FN and exact replay | passed | eleven complete controls; all coefficient edges/endpoints |
| peak materialized storage below `2^50` | passed | 46,649,344-byte process high-water mark; 22,166,512-byte largest logical candidate |
| same-family structured residual degree at most 5 | **failed** | degree unset; explicit membership degree at least 131,458 |
| cost per usable relation below `2^103` | **failed / unset** | no relation oracle or usable-row probability |
| 138,031 independent rows below `2^120` | **failed / unset** | collection, duplicates, rank and sparse solve absent |
| demonstrably non-generic complete path | **failed** | the in-segment transition is constant group translation |
| complete one-target DLP below rho | **failed / unset** | no relation solve, matrix solve or scalar recovery |
| promoted | **no** | selector success alone is insufficient |

The reasonable next target is now sharply defined: retain roughly four
thousand independent-looking boundary classes **and** exhibit a verified
low-degree global defining structure whose relation solver, collection and
linear algebra fit inside the 0.2005% selector margin.  The hash-jump family
does the first by destroying the second.

## Resources, artifacts and reproduction

The corrected isolated execution took 1.145222 s wall, 1.114530 s user and
0.029397 s system on reserved CPU 4 of the AMD EPYC 9V74 host.  It was
uncontended and peaked at 45,412 KiB RSS.  The JSONL also retains the earlier
capacity-226 successful run (0.966049 s wall, uncontended) whose result is
preserved as `hash-jump-result-capacity226-v1.json`.  A still earlier zero-work
launch failed before process creation because the standalone binary had not yet been built;
the isolation wrapper emitted no JSONL entry for a spawn failure.  That event
is preserved separately in `launch-failure.json` and is not counted as an
experiment execution.

The corrected canonical result is 681,841 bytes with SHA-256
`02626b9147e024fd222ad6505688c5113a955a0e68b73afa1e7600cda04c6c97`.
The preserved 681,492-byte capacity-226 receipt has SHA-256
`a7bb2a796275ca50de7792cdc0bc0a47a3222a97a741d4b260d591f75de04b4b`.
The two-successful-run isolation JSONL is 5,061 bytes with SHA-256
`82fd495cc455f806bf8e304fda6b8b17d564421c986a2781a4f8276f33bb2b85`.
The 683-byte launch-failure receipt has SHA-256
`7cb9d7dd87b09a6061255587490efb16d0b5f041f89065ec98bd58ed9d9e005e`.
Semantic evidence SHA-256:
`222fac00fbc81d4cbe49b6be449f30e40a4a6d330097f90e6c6dfc81b1f965ab`.

```bash
cargo test --release --bin p256_hash_jump_screen
cargo clippy --release --bin p256_hash_jump_screen -- -D warnings
cargo build --release --bin p256_hash_jump_screen --bin isolated_bench
target/release/isolated_bench run --wait --cpus 4 \
  --out research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_hash_jump_round293_20261007/isolation.jsonl \
  --label p256-hash-jump-round293-global-capacity -- \
  target/release/p256_hash_jump_screen \
  --round189 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_executed_screen_round2_290_20261006/rounds/round189-r07021.json \
  --executed-manifest research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_executed_screen_round2_290_20261006/manifest.json \
  --round29 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_memoryless_closure_round29_20261006/memoryless-closure-result.json \
  --round30 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_hidden_state_round30_20261006/hidden-state-result.json \
  --round31 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_low_delta_round31_20261006/low-delta-result.json \
  --out research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_hash_jump_round293_20261007/hash-jump-result.json
```
