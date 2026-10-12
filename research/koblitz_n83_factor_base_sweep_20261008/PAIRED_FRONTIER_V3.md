# N83 paired-base frontier: exact addressability, not a runtime result

The frozen v2 size frontier has 57,024,000 single-base/solver tuples but no
small-base/envelope or residual-suffix coordinate. This v3 design retains all
v2 ordinals unchanged and appends a separate finite grid for nested public-x
base recipes. The [generator](paired_frontier_v3.py) and [saved design](paired-frontier-v3.json)
bind the v2 file by SHA-256
`7526b6218dfe0b074368901f3869ddf4a978cb0bacffe8654056056700f8a972`.
The K=600 import receipt and prior K=1,182 search audit are also hash-pinned;
their limited evidence cannot promote a tuple to complete-runtime status.
The subgroup orders used by conditional support calculations are bound to
the original support receipt, SHA-256
`38abcfd36691b9565f97b6d86d4e89ce71ce0b18675586c65b40f945320b8c3f`.

Within each common curve arm, public-x policy, seed and closure, the 11 sizes
`32, 64, 128, 256, 600, 900, 1182, 2048, 4096, 8192, 16627` yield 55
ordered small/envelope choices. Across two arms, three policies, three seeds
and three closures, that is **2,970 paired specifications**. Each crosses
residual summands `0, 1, 2` and **51,840 solver settings**, including the
original v2 solver axes plus the explicit `ordered_first_pair` symmetry level.
The paired suffix contributes **461,894,400** addressable tuples; with the
unaltered v2 prefix, v3 has **518,918,400** ordinal cases. A case is decoded
in constant space, without writing one row per tuple.

This is exhaustive only over that declared finite Cartesian product. It does
not pair different policies, seeds, closures, or curve arms: those sets have
no established ordered suffix, so treating them as the current double-large-
prime domain would be mathematically wrong. It also omits other field bases,
source-point generators, continuous tuning values, solver versions and new
algorithms. The categorical symmetry axis also excludes simultaneous
compositions of two or more listed symmetry policies. A new level,
symmetry-composition axis or cross-stream pairing requires another versioned
design. The original WDSat, Gray and FES backend levels are crossed and
assigned explicit unavailable-adapter dispositions; addressability is not
execution. The zero-residual cases are retained and marked as redundant on
the pair axis; the new first-pair symmetry level itself was absent from v2.
Counted dispositions sum exactly to the paired tuple
population under one ordered rule per tuple.

Only **eight** paired tuples currently reach a bounded mixed-orbit partial
worker: primary `a=0`, `public_x_hash`, seed `2026100801`, signed-Frobenius
closure, small K=600 or K=1,182, envelope K=2,048, five or six summands,
one or two residuals, SAT-CDCL, zero split bits, exact first-pair ordering,
binary enumeration, one thread and full-width modular arithmetic. This is a
source-availability classification. The worker produces at most a bounded
partial-search observation; no tuple has a complete cold index-calculus
runtime. After this design was generated, K=600/K=2,048 gained one
[audited natural public-target trial](verification/k600-mixed-orbit-first-2400-20261012/RESULT.md):
it reached 5,000 conflicts with zero models or partials and
`UNKNOWN_search_cap`. K=1,182/K=2,048 also has one capped natural public-target
trial recorded as `UNKNOWN_search_cap`. These later receipts do not change the
saved design's source-availability counts or rank the bases.

For signed-Frobenius pairs only, a case also reports the conditional count
`C(B+s-1,s) C(L+r-1,r)` of unordered point multisets, where `B=166K_small`,
`L=166(K_envelope-K_small)`, `s=m-r`, and `r` is the residual-summand count.
Dividing that count by the nonidentity subgroup size gives a uniform-target
Markov hit ceiling, clipped at one. The K=600 and K=1,182 three-small/two-
residual cases both clip to one, which supplies no yield prediction. Cases
without signed-Frobenius closure leave this support field null rather than
inheriting an invalid `166K` point count. All point counts are conditional
on successful distinct-orbit construction and an exact nested suffix.

Two exact primary examples are ordinal `156849446` for K=600/K=2,048 and
`158249126` for K=1,182/K=2,048, both with five summands and two residuals.
Their conditional multiset counts are respectively
`4757348390396427029827200` and `13009394065076279872900344`. They can be
inspected and the saved design checked with:

```sh
python3 research/koblitz_n83_factor_base_sweep_20261008/paired_frontier_v3.py --case 156849446
python3 research/koblitz_n83_factor_base_sweep_20261008/paired_frontier_v3.py --check
python3 -m unittest discover -s research/koblitz_n83_factor_base_sweep_20261008 -p test_paired_frontier_v3.py -v
```

The seven focused tests verify the frozen v2 ordinal prefix, every pair's
parent-base indices, paired-case boundary decoding, exact Cartesian and
disposition totals, the eight
source-available recipes, conditional support counts, null support for an
incompatible closure and rejection of changed parent evidence. The saved
design's SHA-256 is
`f7b860c9a652154541c4e6cf7d52ad13bd5a7436a41d38855d206cde36d4b0ac`.
No factor base was constructed or uploaded by this design extension, and the
generator made no SAT, rank, linear-algebra, individual-log, rho or
total-runtime measurement. The separately authorized K=600 trial used the
stored K=600 and K=2,048 objects and charged the cumulative local ledger to
17,400 seconds. Its one capped solver call still leaves useful relation yield
and the selected complete-runtime base unresolved.
