# Residual-column closure allowance for retained N83 pilot pairs

The [saved calculation](closure.json) covers all 2,232 admissible
pair/shape addresses in the pinned [archive-backed support screen](../pilot-pair-support-screen-20261014/RESULT.md).
It gives an exact **partial-row planning threshold**, not a measured relation
rate or a factor-base choice. No N83 solver ran in this calculation.

For a base with `K` signed-Frobenius columns and a disjoint pool with `L`
residual signed-Frobenius columns, let `t` valid partial rows be represented
over the prime subgroup order. Their `t × L` residual coefficient matrix has
rank at most `L`; its left nullspace therefore has dimension at least `t-L`.
At `t = K+L`, there are at least `K` independent **combinations of input
rows that cancel all residual columns**. This is only a residual-closure
allowance. The projected small-base rows can be dependent or zero, a bounded
eliminator can stop, and none of the input rows is assumed to have been found.

The archived support count for each shape is
`M = C(166K+u-1,u) C(166L+v-1,v)` for `u` small and `v` residual summands.
Under the stated uniform nonidentity-query marginal, the expected number of
point-multiset decompositions in `A` queries is at most `A M/(r-1)`, with
`r = 2417851639230796216685689`. Thus
`ceil((K+L)(r-1)/(2M))` queries are necessary before Markov's upper bound
even permits a one-half chance of obtaining **the sufficient `K+L` input-row
threshold**. This is deliberately *not* a necessary query count for `K`
useful full rows: graph algebra can close earlier, and the count says nothing
about solver access or base-row independence. The worker's fixed first-trial
coefficients are outside the uniform-query premise.

For the frozen K600 sequential-small-base, 3+2 Gray/hash contrast:

| Pool | Residual columns `L` | Partial-row threshold `K+L` | Conditional half-threshold query floor |
| --- | ---: | ---: | ---: |
| Gray-prefix, seed `2026100801` | 97 | 697 | 39,468 |
| Hash, seed `2026100801` | 600 | 1,200 | 1,777 |

The hash pool needs a larger residual graph before the unconditional nullity
allowance applies, while its larger exact decomposition support lowers this
specific conditional query floor. Neither number predicts natural partial
yield or chooses the faster total-runtime configuration. The current paired
gate still needs measured solver access, valid natural partials, independently
checked rank, sparse linear algebra, a fresh individual logarithm and matched
rho under complete phase accounting.

The program pins the parent support receipt by SHA-256, recomputes every
point-multiset count, checks the tight integer boundary of all 2,232 floors,
and hashes the full sorted table in the compact output. Its focused controls
check the rank-nullity allowance by independent modular elimination on small
matrices, pin both selected-pair values, replay the saved JSON byte-for-byte,
and reject a modified parent. Reproduce with:

```sh
python3 research/koblitz_n83_factor_base_sweep_20261008/verification/pilot-pair-residual-closure-20261015/closure.py
python3 -m unittest discover -s research/koblitz_n83_factor_base_sweep_20261008/verification/pilot-pair-residual-closure-20261015 -p 'test_*.py' -v
```

This source-only screen adds zero solver-stage charge to the 23,100-second
local N83 ledger. The separately approved Gray/hash paired trial remains
unlaunched while the host fails its frozen quiet admission.
