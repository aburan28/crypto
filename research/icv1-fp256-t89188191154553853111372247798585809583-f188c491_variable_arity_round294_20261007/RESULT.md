# P-256 Dickson variable-arity closure, round 294: result

## Verdict

Changing the relation arity does **not** bring an independent-log Dickson
factor base to Pollard-rho parity.

The exact signed distinct-column domain first covers the P-256 subgroup across
all arities at `B=162`, but no single fixed arity covers there.  The first
fixed-arity cover is `B=164`, at arities 109 and 110.  Even after granting a
perfect free decomposition oracle, success probability one, one operation per
collision sample, no setup, no verification and no linear algebra, the exact
negation-folded `K=164` collision boundary is **13.920747 times rho**.  The
balanced two-list boundary is 19.701921 times rho.

The smallest executed Dickson base with a covering fixed arity is the
terminal-zero depth-9 base `FB1h91cc12460e5a`: 266 independent columns, first
covering arity 59, and an exact favourable boundary of **17.734068 times
rho**.  The registered `FB1h2f8621cda105` first covers at arity 17 and has a
394.425280-times-rho boundary before the distinct-column correction, solver,
replay, rank allowance, sparse linear algebra or recovery.

This is a generic-collection closure for variable distinct signed arity, not
an impossibility theorem for a future non-generic algebraic solver.  No row
passes the boundary gate, structured variable-arity regularity remains unset,
zero relations are reported, and no unplanted full-depth P-256 relation was
attempted.

## One boundary, one unit

The unit is P-256 group-addition equivalents divided by `sqrt(n)`, with the
repository rho reference fixed at `S=1.3`.  Candidate ratios grant exact
coverage and omit every implementation cost.  They are lower boundaries, not
end-to-end measurements.

| variant | columns `B=K` | fixed arity | covers one fixed arity | exact K-th collision / rho | balanced two-list / rho | result |
|:--|--:|--:|:--:|--:|--:|:--|
| Pollard rho | - | - | yes | **1.000000** | **1.000000** | reference |
| all-arity counting threshold | 162 | max at 108 | no | 13.835474 | 19.581419 | `3^B >= n`, but no fixed arity covers |
| **global fixed-arity threshold** | **164** | **109 or 110** | **yes** | **13.920747** | **19.701921** | global optimistic minimum; rejected |
| Dickson depth 8, `FB1h3bafdc8978dd` | 135 | max at 90 | no | 12.628053 | 17.875308 | cannot cover P-256 |
| **Dickson depth 9, `FB1h91cc12460e5a`** | **266** | first at **59** | **yes** | **17.734068** | **25.091548** | best executed covering base; rejected |
| Dickson depth 10, `FB1h15aff3e341f4` | 494 | first at 45 | yes | 24.172704 | 34.194017 | rejected |
| Dickson depth 11, `FB1h39288302415d` | 1,029 | first at 36 | yes | 34.892057 | 49.350815 | rejected |
| Dickson depth 12, `FB1ha2314e64a1f9` | 2,035 | first at 31 | yes | 49.071256 | 69.401499 | rejected |
| registered `FB1h2f8621cda105` | 131,458 | first at 17 | yes | 394.425280 | 557.802111 | rejected |

Depth 7 has only 67 columns and also cannot cover the subgroup.  The exact
runner checks all widths 1--512, all arities at those widths, and each fixed
arity 2--256 against widths through 131,458.  It performs 2,127,402 counted
big-integer recurrence operations.  The ordered boundary-table digest is
recorded in the canonical JSON.

## Exact coefficient-domain boundary

For `B` independent logarithm columns and fixed distinct signed arity `m`,

```text
D(B,m) = 2^m * C(B,m).
```

The exact recurrence reproduces

```text
sum(m=0..B) D(B,m) = 3^B
```

at every checked width.  `3^161 < n <= 3^162`, so 162 columns are necessary
even when every arity is pooled.  At `B=162`, however, the largest single
fixed domain is the arity-108 domain and remains below `n`.  At `B=164`, the
arity-109 and arity-110 domains are equal:

```text
116383775127750444427048016195109203089116454475801022361283226351238392053760
```

which is just above the registered P-256 subgroup order.  This is deliberately
optimistic: it treats every coefficient pattern as a distinct useful group
element.

The exact negation-folded K-th cross-colour collision law imported from Round
24 is

```text
T_K = sqrt(2*n) * Gamma(K+1/2) / Gamma(K).
```

Here `K=B` because the geometric bases retain independent logarithms.  Setting
the disjointness probability to one further favours the candidate.  Thus the
13.920747 ratio at `B=164` is below the cost of the applicable collector, not
an estimate of an implementation.

## Accounting correction preserved

The preregistered H2 formula copied Round 21's superseded shortcut
`sqrt((pi/2)*K*n)`.  The first successful receipt therefore reported an
overly favourable 12.346348-times-rho value.  Dashboard reconciliation caught
the mismatch with Round 24 before publication.  The canonical v2 artifact
hash-checks Round 24, uses the exact gamma recurrence for every gate, and
retains the shortcut only in fields named `preregistered_superseded_*`.

Both v1 receipts are preserved.  The correction raises the global floor from
12.346348 to 13.920747 and cannot change a failure into a pass.

## Native Dickson inventory and exactness

The runner builds and independently rebuilds the terminal-zero Dickson bases
at depths 7--12.  It separately rebuilds the registered depth-18 nonzero coset
and exactly reproduces:

```text
FB1h2f8621cda105
factor-base SHA-256 2f8621cda105a51f8703dabe38dec87f3715ab855770acd7b2af0a14906f9e42
point-set SHA-256   70ab21cbda76205ff43ab6bf85c68607be9f6b89df6399e6fa98e3489b888cf1
columns / points    131458 / 262916
```

Fourteen native builds verify 270,968 signed points with zero failures.  The
canonical result replay is byte-identical, including its semantic evidence
digest `37c5392892a7ea29b7187c438322c1c6ef976cf8493f9f707cfc545834492d98`.

The imported local Dickson solving-degree maxima remain 3, 3 and 4, and the
split image gates remain degree 2.  Those measurements reach only small local
systems and split image trees through 16 summands.  They do not establish a
degree-at-most-five global system at the depth-9 covering arity 59 or the
global threshold arities 109/110.  Variable-arity unsplit degree therefore
remains null.

## Promotion gates

| gate | status | evidence |
|:--|:--:|:--|
| dependency hashes and exact integer identities | pass | five dependencies; all `3^B` identities exact |
| native factor-base verification | pass | 14 builds; 270,968 signed points; zero failures |
| any covering generic boundary at or below rho | **fail** | global minimum 13.920747 times rho |
| structured residual degree at most 5 | **fail / unset** | local imported degrees do not cover arity 59 or 109/110 |
| usable relation below `2^103` | **fail / unset** | no complete decomposition solver |
| complete collection below `2^120` | **fail / unset** | generic lower boundary already misses rho |
| projected materialized storage below `2^50` | **fail / unset** | no promoted complete solver/storage plan |
| promoted | **no** | zero relations; no complete DLP |

The correct next target is no longer another arity or Dickson depth.  It must
be a non-generic global algebraic solver that amortizes at least a 13.920747
factor at the theoretical `B=164` boundary, or 17.734068 on the smallest
executed covering Dickson base, while retaining exact membership, degree at
most five and the complete cost gates.  Exact log transport is the other
logical escape, but Rounds 21--29 show that the available P-256 transports
either leave `K=B` or reduce to scalar-orbit rho.

## Resources, artifacts, and reproduction

The corrected isolated execution used reserved CPU 4 on the AMD EPYC 9V74
host.  It completed successfully and uncontended in 21.114896 s wall,
21.042157 s user and 0.072106 s system, peaking at 203,024 KiB RSS.  The
largest deterministic logical point inventory is 66,088,215 bytes and the
algorithm writes zero disk bytes apart from its result artifact.

- canonical v2 JSON: 116,909 bytes, SHA-256
  `9716fc1ba7d24a023627e8701094ca870b2121a84f7a0148784f70565640815a`;
- v1 deterministic shortcut receipt: 95,512 bytes, SHA-256
  `a36b1ba88b567e9300b18a740c252aa8404a3fd5294f3676da5dc6ce4e2ea21a`;
- v1 process-telemetry receipt: 95,549 bytes, SHA-256
  `006aa9e11f20f749818764b7b332be22c0a6c1282115b508996cfc3e1e974988`;
- three-run isolation JSONL: 7,502 bytes, SHA-256
  `9a32f34edcfe3b9cbd68d81f50bb95283e0c6c3c38c634dc2f111de86819023d`.

```bash
cargo test --release --bin p256_variable_arity_closure
cargo clippy --release --bin p256_variable_arity_closure -- -D warnings
cargo build --release --bin p256_variable_arity_closure --bin isolated_bench
target/release/isolated_bench run --wait --cpus 4 \
  --out research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_variable_arity_round294_20261007/isolation.jsonl \
  --label p256-variable-arity-round294-exact-kth -- \
  target/release/p256_variable_arity_closure \
  --round21 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_fb_log_transport_round21_20261006/transport-result.json \
  --round24 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_negation_cross_colour_round24_20261006/negation-result.json \
  --round25 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_scalar_orbit_round25_20261006/scalar-orbit-result.json \
  --round29 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_memoryless_closure_round29_20261006/memoryless-closure-result.json \
  --round293 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_hash_jump_round293_20261007/hash-jump-result.json \
  --out research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_variable_arity_round294_20261007/variable-arity-result.json
```
