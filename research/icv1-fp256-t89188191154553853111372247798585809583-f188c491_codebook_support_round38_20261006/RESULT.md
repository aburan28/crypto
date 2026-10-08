# P-256 correlated pair-codebook support, round 38: result

Date run: 2026-10-06

No registered correlated codebook reaches rho after jointly charging lost
support and cold pair construction.  At 256 MiB, the best candidate keeps four
of eight pair slots hot and constructs four pairs.  Its deliberately
overcounted state support exceeds the P-256 group, but its arithmetic floor is
still `1.015685` rho.

Making all eight slots hot avoids construction but collapses the support upper
bound to `2^202.266` states.  The resulting target-retry lower bound makes that
row `1.496e16` times rho.  An all-hot codebook needs at least 441,075,261 signed
entries, or 28,228,816,704 bytes (`2^34.716`), merely for its optimistic bound
to touch parity.  That threshold counts ordered repeats, overlapping columns,
inconsistent signs, starts below cutoff and maximum-length paths as valid, so
it is necessary and emphatically not sufficient.

Probabilistic mixtures cannot improve the selected rows: their total
cost/coverage ratio is a coverage-weighted average of the pure-row ratios.
No correlated-start implementation, benchmark or full-depth unplanted
relation was attempted.

## Frozen boundary

Round 37 is imported by exact SHA-256
`8b86bbc1488e2c3f9d6bc8a6c4fc9cf0c8c21a1a440adcc8fecba5569c7ed687`.
The audit preserves:

```text
factor-base columns                       131,458
pair slots                                      8
signed singleton choices                  262,916
maximum path states                           307
complete signed-pair entries       34,562,148,612
Montgomery-affine entry                      64 bytes
corrected base ratio           0.998569034150286 rho
local-oracle ratio             0.964336477130181 rho
mean segment capacity          224.361249307511
```

For codebook size `K` and `h` hot slots, every row grants the support upper
bound

```text
binom(8,h) * K^h * 34,562,148,612^(8-h) * 262,916 * 307.
```

It then charges exactly `8-h` cold additions and no access, generation,
factor-base-read or validation cost.  The reported ratios are therefore lower
bounds under assumptions more favorable than an implementable selector.

## Selected frontier

| codebook | hot / cold slots | represented-state upper bits | target coverage upper | ratio before retry | complete / rho lower bound |
|:--|--:|--:|--:|--:|--:|
| 2 MiB | 2 / 6 | 271.124535 | 1.000000 | 1.024243 | **1.024243** |
| 64 MiB | 3 / 5 | 267.116061 | 1.000000 | 1.019964 | **1.019964** |
| **256 MiB** | **4 / 4** | **260.429516** | **1.000000** | **1.015685** | **1.015685** |
| necessary all-hot threshold | 8 / 0 | 255.997934 | 0.998569 | 0.998569 | **0.999999987** |

The necessary threshold contains 441,075,261 entries, 1.276180% of the
complete table.  Its 28.23-GB representation is below `2^50` bytes but is not
cache resident.  The tiny numerical margin below one comes from the integer
ceiling at the frozen base ratio; no actual coverage, matching or runtime is
established.

The 256-MiB sweep exposes the tradeoff directly:

| hot slots | cold additions | support upper bits | coverage upper | complete / rho lower bound |
|--:|--:|--:|--:|--:|
| 3 | 5 | 273.116061 | 1.000000 | 1.019964 |
| **4** | **4** | **260.429516** | **1.000000** | **1.015685** |
| 5 | 3 | 247.099114 | 0.002092 | 483.459 |
| 6 | 2 | 233.090640 | 1.269e-7 | 7.934e6 |
| 7 | 1 | 218.274811 | 4.401e-12 | 2.278e11 |
| 8 | 0 | 202.266337 | 6.677e-17 | 1.496e16 |

One additional hot slot moves work out of pair construction only by removing
too much start entropy.  The minimum is not hidden between rows: for branch
weights `q_h`, the mixed ratio is
`sum(q_h*C_h)/sum(q_h*epsilon_h)`, a weighted average of
`C_h/epsilon_h` with nonnegative weights `q_h*epsilon_h`.

## Exact controls

```text
sweep rows                              27
support monotonicity failures            0
construction-cost monotonicity failures  0
selection failures                       0
mixture dominance identity            true
toy cells                                6
toy formula failures                     0
toy upper-bound failures                 0
invalid toy sequences credited      61,334
threshold arithmetic failures            0
false positives                           0
false negatives                           0
```

The two tractable signed-pair universes enumerate every ordered sequence and
signed singleton.  The formula matches all sequence populations exactly and
deliberately credits 61,334 invalid overlapping-column sequences.

Deterministic digests:

```text
ordered sweep  1bffeb9aeb7db1468c43ce951d7cb5bd4d70a04f89e4bf26f84a59770bd73bd6
threshold      dba4dd6364eae4b1a170eb0a4bdd13353ae008f4910d815b0e597fa14ac31233
semantic       b454893599248b0e9bd795fab341a759e95396bf6e4862a93c4c618d481f034e
```

An independent release replay reproduced the complete artifact byte for byte.

## Requirement status

| requirement | status | evidence or gap |
|:--|:--|:--|
| exact support/cost sweep | verified | 27 arbitrary-precision rows; zero arithmetic, selection or toy failures |
| probabilistic mixtures | bounded | weighted-average identity proves no mixture beats the best pure row under the registered model |
| registered codebook at or below rho | failed | best 256-MiB lower bound is `1.015685` rho |
| complete measured time at or below rho | failed | analytical screen fails before implementation |
| storage below `2^50` | verified for necessary threshold | 28.23 GB, but cache residency and sufficiency fail |
| proved P-256 usable-relation probability | not established | represented-state count is an overcount, not achieved coverage |
| structured residual degree at most 5 | not established | no algebraic residual solve in this round |
| collection below `2^120` and per row below `2^103` | not established | no usable full-depth relation or collector |
| demonstrably non-generic end-to-end method | not established | no implemented correlated selector |

## Decision

Reject correlated pair codebooks at 2, 64 and 256 MiB.  Do not implement or
benchmark them: each registered depth fails under a support bound that already
credits invalid branches.  A further iteration must reuse pair work without
restricting independent start entropy—for example, exact batched pair-sum
generation amortized across otherwise independent starts—and must charge the
batch construction and storage explicitly.

## Reproduction

```bash
cargo test --bin p256_codebook_support
cargo run --release --bin p256_codebook_support -- \
  --round37 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_pair_entropy_round37_20261006/pair-entropy-result.json \
  --out research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_codebook_support_round38_20261006/codebook-support-result.json
```

Canonical result artifact: `codebook-support-result.json`, 23,352 bytes,
SHA-256
`4297328eee07822331c68bddd65d8b0155ba1c72c6e0dda044dc97d42aba5cce`.
