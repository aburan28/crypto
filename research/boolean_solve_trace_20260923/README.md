# Complete generic Boolean-solving traces

The previous goal turn established conditional gains in a specialized linear-
tail kernel, with small-case guard failures and no whole-solver measurement.
That was progress, not completion of the broader goal. This study measures
complete generated Boolean solves with identical search rules around the
different kernels. It also retains a search-only alternative, since reducing
node count can still increase total cost.

## Scope and solver contract

This is a standalone mathematical driver, **not the repository inherited-F4
solver**. It accepts only generated public Boolean systems. No curve, imported
target, discrete-log scalar, private key, relation collection or rho workflow
is present. A whole-solve result here remains specific to this driver and
fixture distribution; it is not a production/index-calculus claim.

The solver specializes squarefree polynomials exactly, cancelling duplicate
terms by parity. It checks zero/constant equations, reduces explicit affine
equations, applies forced assignments and solves wholly affine residuals
directly. Kernel-backed arms then compute the same complete degree-3 affine
tail and apply its forced assignments. Otherwise they branch on the most
frequent variable, breaking ties by lowest index and trying zero first.

Both Boolean branches are exhaustive. Every propagated equation is either an
exact specialization or a row-space consequence of valid products. UNSAT is
returned only after a sound contradiction or both complete branches fail.
Node and kernel-resource limits return **UNKNOWN**, and are never reinterpreted
as UNSAT. SAT assignments are evaluated independently on the original equations.

The flat, sparse-bucket and hybrid kernel arms must produce the same result,
model, logical search counters and stable trace digest. The trace digest is a
reproducibility aid, not a collision-resistant certificate. Search-only follows
the same non-kernel rules and may traverse a different tree. Its completed
outcome independently checks the status of kernel-backed runs. An unverified
UNSAT or any censored sample blocks an unqualified performance promotion.

## Fixed workloads and costs

`protocol.json` freezes four sizes (12/16/20/24), three families, discovery seeds
17/937, two new holdout seeds and four balanced arm-position repetitions. Each
system has n+2 sparse quadratic equations. Planted witnesses are public fixture
construction data and are never given to the solver. Cross-planted systems
include an affine consequence hidden between quadratic generators. Unplanted
systems have independent constants and no assumed status.

The discovery-only `resource_probe_01` established a feasible execution budget
before the holdout protocol. Its source and executable hashes, command, output
and runtime are retained. It is not promotion evidence. The complete grid has
48 cells, a deterministic 200,000-node cap, explicit kernel caps and bounded
worker/campaign runtimes. No fixture is screened or dropped for being slow.

For each alternative against flat, every family at n=16/20/24 must have a 95%
paired-bootstrap lower bound above 2.0. The hybrid also faces the pointwise
fastest of search, flat and bucket. Every arm on the complete grid must finish
with verified results before unqualified promotion. Completion costs remain
null for censored cells; their observed work is retained. Intervals describe
the fixed seeds/repetitions, not a population or independent reproduction.

Cold totals include solver creation, system copying, specialization, cheap
propagation, recursive search, kernel work, solver/context destruction and
independent result validation. Diagnostic-record formatting is outside timing.
Fixture and search-reference preparation are common supplied-input work outside
arm timing and inside process receipts. Kernel time and calls by active-variable
count expose whether the faster kernel matters to complete solves.

`KERNEL_SOURCE.json` binds the kernel derivation. The generic generator cap is
raised to 36 for these square/overdetermined systems. Solve-lifetime caches have
explicit limits. Four tests compare all backends against Boolean enumeration
on small systems, check models, trace equality, specialization and censored
limits. No production solver code is changed.

```sh
python3 research/boolean_solve_trace_20260923/run.py --out research/boolean_solve_trace_20260923/run_01
python3 -m unittest discover -s research/boolean_solve_trace_20260923 -p 'test_*.py'
```

Every output directory is new. Frozen sources/protocols, executable hashes, raw
fixtures/samples, process receipts, test output and a manifest are retained.
Evidence replays use temporary copies. The previous kernel studies remain
unchanged, and the full thread goal stays open unless the broader requirements
are actually verified.

## Selective degree-2 successor

The first complete grid finishes every solve with verified results. Its hybrid
kernel does not improve total solving, and search-only often wins despite many
more nodes. `run_01` remains immutable.

`word_protocol.json` fixes a new policy before fresh holdouts: invoke degree-2
inference only when at most ten variables remain active. The space then has at
most `1 + 10 + C(10,2) = 56` monomials and fits in one u64. The word implementation
generates every admitted product, cancels by XOR, eliminates the quadratic
columns and canonically reduces the affine tail. It includes all required
multiples of linear generators; it is not merely cancellation of input rows.

Two arms use exactly this policy: `small_flat` runs the existing general flat
kernel at degree 2, and `word_tail` uses compressed active-variable coordinates
and the one-word implementation. They must produce identical complete search
traces, counters and models. Their policy differs from always-degree-3 inference,
so equality is required within policy groups, not across them.

The successor retains all four previous methods, uses six rotated repetitions
and fresh holdouts 20261010/1310719. The new full-solve gate requires every
n=16/20/24 family comparison to have a 95% lower bound above 2.0 against the
fastest prior whole-solve method, and above 1.05 against the same-policy flat
control. All 18 comparisons and every completion/verification gate must pass.
Six Rust tests now cover the added solver modes, degree-2 kernel equality on
embedded variable sets, and complete-search trace equality for the new pair.

```sh
python3 research/boolean_solve_trace_20260923/run.py --protocol research/boolean_solve_trace_20260923/word_protocol.json --out research/boolean_solve_trace_20260923/run_02
```

## Ordered specialization successor

The second run improves the same-policy flat implementation by about 1.13–1.39x,
but does not beat the strongest prior solver. `merge_protocol.json` therefore
targets specialization in search-only, without changing its mathematical search.
It uses fresh holdouts 20261011/1441793 and seven balanced arm-position repetitions.

Setting variables to zero only removes terms and preserves monomial order.
When exactly one variable is set to one, terms containing it and terms not
containing it form two ordered sublists. Removing the common variable from the
first sublist preserves relative degrees and the symmetric differences used
by DegRevLex. A linear symmetric-difference merge therefore gives the exact
canonical specialization. The merge writes backward into the original buffer
so unread terms are not overwritten, then compacts the surviving suffix.
Multiple true assignments retain the sorting fallback. Equation order and
deduplication are unchanged.

The search and merge_search arms must agree on every result/model, logical
counter and trace digest. `specialized_terms` remains a logical input-term
count, not a calibrated instruction count. The new gate requires a 2x lower
confidence bound against the fastest of all six prior complete-solve methods
in every n=16/20/24 family cell. All earlier gates and failed runs remain.
Eight Rust tests now include 6,912 polynomial/partial-assignment checks and
larger complete-trace comparisons for the specialization change.

```sh
python3 research/boolean_solve_trace_20260923/run.py --protocol research/boolean_solve_trace_20260923/merge_protocol.json --out research/boolean_solve_trace_20260923/run_03
```
