# Incremental pullback improvements, always paired with the benchmark

This round tests two separate engineering candidates against the frozen
coefficient-pullback implementation at commit
`b8bea555e04968f68edf170f50f19bfab38a3c03`:

1. Normalize the coefficient equation and cancel redundant terms before
   constructing the restricted binary-linear solve.
2. Cache coordinates during Gaussian elimination while retaining every
   counted field addition. This changes representation overhead, not the
   mathematical operation vector.

[RESULTS.md](RESULTS.md) contains the frozen comparison and verdict.
[comparison.json](comparison.json) is the machine-readable audit and aggregate
data. The prior experiment and its results remain intact. These are binary
a=0 Nagao experiments; no prime-field or F4 implementation is included.

## Frozen contract and benchmark gate

[contract.json](contract.json) was written before measuring either candidate.
[targets.json](targets.json) fixes all 12 prior benchmark targets plus 12 new
holdouts. Both corpora have two uniform and two known-decomposable targets at
each of (n,d)=(11,5),(23,7),(29,8). The fresh targets exclude the historical
solver_02 through solver_09 and prior pullback targets. All targets are frozen
before candidate validation or timing. They become development data after
this round; subsequent tuning needs new holdouts.

The stage benchmark runs all five variants on each target in first-relation
and complete-enumeration modes: **240 cells**. It reruns the unmodified
pullback reference, both incremental variants, symmetric S4, and chained S3.
The Semaev controls use the existing CryptoMiniSat adapters. All solvers have
the same cooperative three-second all-phase budget. Timeouts, partial
enumerations and budget overruns are retained. There is one repetition, so
wall times remain descriptive without a runtime-speedup confidence claim.

The original WDSat suite uses a=1 and another input protocol. Following the
parent benchmark contract and AGENTS.md §8, this is the equivalent matched
a=0 suite; it leaves the WDSat corpus unchanged. `compare.py` checks that
completed baseline replays reproduce the prior exact solution sets and field
counters. A replay censored by the time limit is reported separately.

The declared diagnostic gate is at least 10% fewer multiplications for the
normalized solver on **each** complete, equal-output eleven-bit enumeration
corpus. Both corpora must contain all four completed target pairs. Cached and
uncached normalized solvers must have identical field counts on those pairs.
This gate is not the repository's calibrated full-cost gate.

The conditioned branch count remains O(2^(2d)). Elimination uses O(d²) field
additions and O(d) field multiplications per branch; arithmetic itself depends
on n. For B signed factor points and uniform affine targets, the same support
ceiling holds: p <= min(1,8*C(B/2,3)/(#E-1)). Domain exclusions only decrease p.
No exponent is fitted and no calibrated full-DLP floor is derived here.
S, cost/rho and cost/floor stay null until a common operation conversion prices
all field, scalar and control costs. The classification is **engineering
candidate, pending normalized-cost evidence**, not an advance.

## Cancellation that removes branch multiplications

Use the same target (r,s), curve y²+xy=x³+1, support V, chosen h2=b²+b+r,
nonzero b, and conditioned root z as the previous pullback. Here z != 0,r,h2.
Let t=r+z, u=h2+z, and I_u=image(w -> w²+u*w on V). The previous solver
searches v in I_u using

```
C2*v²+C1*v = rhs
C2 = t*gamma²
C1 = b*z*gamma
rhs = K+t*delta²+b*z*delta
gamma = t/(b*r)
delta = (t*c0+K)/(b*r)
c0 = r*h2+z*u
K = H_(a=0)(z)
```

Since t != 0, divide the equation by t. The linear coefficient becomes z/r,
which is independent of b and h2. The relation b*r*delta=t*c0+K gives

```
gamma²*v²+(z/r)*v = delta²+b*delta+c0
```

The residual cubic has h1=r*h2 and h0=b*c+r²*h2 at a=0, where
c=r²+b*(r+s). Therefore t*c0+K=b*c+r*z*u. Define target constants

```
D = (r+s)/r
J = D²+D+r
eta = z*u/b
delta = r+b*D+eta
```

Expanding the right-hand side, the terms b*eta=z*u cancel the z*u in c0,
and h2=b²+b+r cancels the remaining b*r terms. The final restricted equation is

```
gamma²*v²+(z/r)*v = eta²+b²*J
gamma = (1+z/r)/b
a = gamma*v+delta
```

`NormalizedSearch` hoists D,J once per target, z/r and 1+z/r once per root,
and r+b*D and b²*J once per b. All are counted as actual field operations.
This removes repeated cubic evaluation and redundant products from each
branch. It retains the old setup, support-image basis, coefficient deduplication
and exact point recovery. The exceptional r=0 chart uses the previous solver.
The transformation changes neither the allowed support nor relation yield.

`CachedSearch` then stores the coordinate vector of each elimination column.
When it subtracts a column, it still calls `f.add` for both the field value
and its preimage. An XOR updates the cached coordinates of that same value;
it does not replace or hide counted field arithmetic. Pivot coordinates no
longer require repeated field-to-coordinate conversion. This is a separate
ablation: the operation vectors must match the normalized solver exactly.

## Correctness and complete scalar recovery

Before any stage measurement, validation enumerates all 43 affine five-bit
targets and compares the full (a,b,z) candidate sets across all three solvers.
Both candidates then completely enumerate the projected point decompositions
and compare with the independent group pair oracle. The coordinate-cache
variant must match the normalized field counts, including the exceptional
target chart. Larger measured solutions are checked against exact group
ground truth; completed enumerations must match it in full.

The cold ECDLP benchmark executes the **unmodified** previous `e2e.py` function
in a separate module namespace. The driver substitutes only its decomposition
callable and restores the binding after each run. Target generation, cofactor
projection, factor columns, relation collection, scalar matrix, individual
descent and final verification are identical code. There is no shared warm
precomputation. The run is sequential; the isolated module binding is not a
thread-safe concurrent adapter.

There are 48 cold runs: two toy groups (orders 44 and 508, prime subgroups 11
and 127), six seeds, and four variants including matched rho. Seeds 91,137,211
rerun the frozen targets; 307,401,503 are fresh. Each run verifies the scalar
against the planted test value and [k]G=Q. The seed is used to create Q; its
secret is never supplied to the decomposition solver. The harness projects
entire relations by cofactor 4 before scalar elimination, and handles equal,
opposite and identity projected factor points correctly.

Fresh seeds do not guarantee distinct targets in these tiny groups. The
five-bit fresh-seed panel has two distinct Q values, one overlapping the
frozen panel; the nine-bit panel has two distinct Q values, both unseen.
All seed runs remain in the evidence. The comparison explicitly reports
these collisions rather than treating six seeds as six distinct targets.

The post-run comparison audits **every** IC oracle call, including empty
answers, against the independent pair oracle. It independently verifies every
scalar, checks prior cold-baseline field/scalar counter replays, and reports
whether the paired variants follow identical attempt paths. Different valid
first relations are permitted, but their downstream costs remain in the cold
totals. Missing completions block an unqualified end-to-end claim.

Field add/mul/square and scalar modular add/mul/inversion-call vectors are
separate. Setup includes exhaustive toy curve-order checking; it is charged
but is not a large-group scaling model. Python control, coordinate conversion,
bookkeeping and the cache's bit operations appear only in cold wall time.
No calibrated total-operation speedup follows from these partial vectors.

## Rerun

Use Python 3.12 with `pycryptosat==5.14.7`. From the repository root:

```sh
python research/nagao_relations/pullback_incremental/run.py --targets research/nagao_relations/pullback_incremental/targets.json --output /tmp/pullback-increment-new
python research/nagao_relations/pullback_incremental/compare.py --input /tmp/pullback-increment-new --output /tmp/pullback-increment-comparison-new.json
```

Both commands refuse to overwrite evidence. The first command includes
validation, all 240 solver cells, and all 48 cold DLP runs. `--validate-only`
performs just the exhaustive validation. `--freeze /tmp/new-targets.json`
recreates this contract's target file; changing the holdout seed requires a
new versioned contract before further tuning.

An initial runner failed before any solver measurement because the inherited
module `optimized` shadowed the new file. Its traceback, source and provenance
remain in `results/` and `campaign.stderr`. The runner now loads the candidate
by an explicit file path and unique module name. The successful campaign is
`results_v2/`; candidate arithmetic and the predeclared contract were unchanged.
