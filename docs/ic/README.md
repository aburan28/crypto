# ic research framework

The standalone ic executable inspects elliptic-curve parameters and runs
bounded, reproducible index-calculus experiments on internally generated
known-answer or public hash-derived Koblitz targets. Imported parameter-file
points are used only for mathematical validation.

**Agent scoreboard:** per-stage records and next targets to beat live in
[`BOUNDARY_TARGETS.md`](./BOUNDARY_TARGETS.md) and
[`boundary_targets.json`](./boundary_targets.json) (binary, Koblitz, prime;
`schema_version` 2). Beat claims must include the ledger's **measurement
schema** fields — including **FFD / degree of regularity** on algebraic
`decomposition` frontiers — or they fail closed.

**Autolab runner:** agents push those beats with the local control plane at
[`research/sat_factor_base_review_20260908/autolab/`](../../research/sat_factor_base_review_20260908/autolab/)
(`boundary_autolab.py`). It pins the ledger, fail-closed validates measurement
reports, and launches the public-synthetic `koblitz_rank_fixture` /
`koblitz_rho_fixture` producers for the priority Koblitz `vs_rho` rungs.

## Commands

    cargo build --release --bin ic
    ./target/release/ic --help
    ./target/release/ic list

A bare curve name is an inspection shortcut:

    ./target/release/ic ecc2k-130
    ./target/release/ic inspect --curve ecc2k-95 --json
    ./target/release/ic inspect --curve secp256k1
    ./target/release/ic inspect --file docs/ic/prime-example.json

Named profiles are ECC2K-130, ECC2K-95, sect163k1, and secp256k1. Their data
comes from the repository's ECC2K client definitions and existing curve
constructors. ECC2K-130 has an abstract normal-basis profile: the inspector
checks its Koblitz group-order recurrence but explicitly reports that the
coordinate representation and challenge points have not been imported.
It does not substitute polynomial-basis coordinates.

## Generated instances and larger synthetic curves

    ./target/release/ic
    ./target/release/ic --solver sat
    ./target/release/ic run --degree 11 --curve-a 1 --known-log 53 --solver enumerate
    ./target/release/ic run --degree 9 --random-target --seed 42 --json
    ./target/release/ic run --degree 13 --curve-a 0 --solver pair-table --summands 3

The default run uses K_0 over GF(2^9) and known logarithm 53. The solver
constructs the target from the declared known answer, then checks both the
recovered value and the group identity. No public point file can be passed
to the run command.

The staged `workflow` command also accepts a scalar-blind public target:

```json
{"targets":[{"public_hash_seed":29}]}
```

This form hashes the curve identity, seed, and counter to an abscissa,
chooses a lift from the next digest bit, and applies the public cofactor.
It constructs and records no target scalar. Both individual descent and
the signed-Frobenius rho baseline accept a recovered value only after
checking `[d]G = Q`. The older `known_log` and `random_seed` forms retain
their expected-scalar check as an additional known-answer control.

The following knobs are recorded in each run report:

- degree: the field degree n; odd values 3 through 63 for the Koblitz
  family (a curve is usable only when its largest prime factor exceeds
  the cofactor, which above 23 holds for degrees 29 (curve-a 1), 31, 37,
  and 39 (curve-a 0)), or k times an odd number for a subfield curve;
- curve-a: 0 or 1, with b fixed to 1, for a Koblitz curve; the
  coordinate of a in the subfield basis otherwise;
- subfield: the degree k of the subfield the curve is defined over,
  q = 2^k (default 1, the Koblitz family; up to 8), see below;
- curve-b: the coordinate of b in the subfield basis (default 1);
- known-log: a positive scalar smaller than the selected subgroup order;
- random-target: draw a known-answer scalar reproducibly from seed;
- seed: seed for the generated fixture and relation sampler;
- factor-index: candidate in the legacy degree-ord_n(2) factor family;
- factor-base: a recipe file written by `ic search`, replacing factor-index;
- summands: factor-base points per relation, 2 (default), 3, or 4;
- max-trials: 1 through 1000000;
- solver: groebner, sat, enumerate, or pair-table;
- batch: targets decomposed per parallel batch (0 = CPU count);
- control: legacy accounting, see below.

Solvers answer the same question, "is this target a sum of `summands`
factor-base points", and every returned decomposition is re-added in the
group before it becomes a relation:

- enumerate: ordered-tuple search, `|F|^(m-1)` group operations per target;
- pair-table: meet in the middle over a table of all `|F|(|F|+1)/2` pair
  sums built once per run — one lookup per target for two summands,
  `|F|` for three, `|F|²` for four (16 bytes per table entry);
- groebner: the Weil-restricted Semaev system reduced by matrix-F4;
- sat: the same system, CDCL with native parity rows.

Not every degree/coefficient combination has a usable subgroup. A valid
curve does not guarantee successful collection or an invertible relation
system. Factor-base materialization is capped at 4096 abscissae
(dimension 12 for the legacy family). Fixture generation and inspection do
not require a usable factor base.

The run display follows factor-base construction, an optional pair table,
relation collection, linear algebra, and verification. By default relation
columns are signed Frobenius orbits merged by their cofactor projection, the
relation matrix is kept in reduced echelon form over Z/rZ as relations
arrive, and the run stops the moment the scalar is pinned — dependent
relations are counted, never padded. An inconsistent relation or a pinned
scalar that fails `[d]G = Q` invalidates the run outright. `--control`
restores the earlier accounting for matched comparisons: one column per
Frobenius orbit, no projection merge, a fixed surplus of relations, and a
single solve at the end. Incomplete results exit unsuccessfully and remain
incomplete in JSON reports.

### Subfield curves beyond the Koblitz family

    ./target/release/ic run --degree 14 --subfield 2 --curve-a 0 --curve-b 2 --solver pair-table
    ./target/release/ic search --degree 22 --subfield 2 --curve-a 1 --curve-b 3 --family divisor
    ./target/release/ic logs --degree 14 --subfield 2 --curve-a 0 --curve-b 2 --database logs14.json

The Galbraith–Granger–Merz–Petit construction needs only that the curve
be defined over a subfield: with `--subfield k` the synthetic curve is
`y² + xy = x³ + a x² + b` with `a, b ∈ GF(2^k) ⊂ GF(2^n)`, `n = k · e`
and `e` odd, and the `2^k`-power Frobenius `π` plays the role squaring
plays on a Koblitz curve. `a` and `b` are named by their coordinates in
an `F_2`-basis of the subfield (the kernel of `X^{2^k} + X`), so
`--subfield 1 --curve-a a --curve-b 1` is exactly `K_a`. Point counting
goes through `#E(GF(2^k))`, found by enumeration, and the trace
recurrence `s_i = t·s_{i−1} − q·s_{i−2}`; `λ` is the root of
`λ² − tλ + q` with `π(G) = [λ]G`.

The invariant factor bases are the kernels of `q`-linearised
polynomials `Σ c_i X^{q^i}` with `c_i ∈ GF(q)`, classified by the
irreducible factors of `x^e − 1` over `GF(q)` (Cantor–Zassenhaus over
`GF(q)`, carried out inside `GF(2^n)`); a factor of degree `d` gives a
subspace of `2^{kd}` abscissae whose points fall into orbits of length
dividing `e`. Recipe indices refer to that factor list, which for `k = 1`
is the familiar `F_2` list in the same order, so every Koblitz recipe,
document and report is unchanged. Documents record `subfield` and
`curve_b` (omitted when 1) and are bound to them. Two things differ
from the Koblitz case in practice: an even `n` is allowed (the
Artin–Schreier solve for even degree is a linear solve, not the
half-trace), and when `x^e − 1` splits into binomials `x^d − c` over
`GF(q)` the invariant subspaces are multiplicative cosets whose
inverses land in the reciprocal factor's subspace, so a curve with
`Tr(a) = 1` and `b = 1` has no points over them at all — the search
scores such bases at zero and a run over one reports a base with no
usable columns, so pick `b` (or `a`) accordingly, as the examples above
do.

## Searching for a factor base

    ./target/release/ic search --degree 15 --curve-a 1 --spec-out fb15.json
    ./target/release/ic run --degree 15 --curve-a 1 --factor-base fb15.json --solver pair-table
    ./target/release/ic search --degree 31 --curve-a 0 --summands 3 --family divisor --max-dimension 11

`search` scores factor bases by the number the pipeline actually pays for:
the expected trials to collect a determining system,
`(columns + 1 + extra) / coverage`, where coverage is the fraction of
subgroup targets that decompose into `summands` base points. Coverage is
measured exactly on one shared target set — the whole subgroup when
`r − 1 ≤ --exhaustive-cap` (default 4096), otherwise `--targets` seeded
samples — by enumerating every witness of every target through the pair
table. Candidates come from three families:

- factor: the legacy single irreducible factors (`--factor-index`);
- divisor: every product of irreducible factors of `x^n − 1` whose degree
  lies in `[--min-dimension, --max-dimension]`, the complete list of
  Frobenius-stable linear subspaces;
- union: Frobenius closures of random seed spaces of dimension
  `--union-min-seed` to `--union-max-seed`, `--union-samples` per dimension
  plus the standard basis.

Each candidate is also tried after 2-torsion saturation (when the cofactor
is even; `--no-saturate` skips it) and after greedy orbit pruning
(`--no-prune` skips it): signed orbits are dropped while doing so lowers the
expected trial count, which is an exact recount over the witness list
rather than a re-search. The pruned base stays Frobenius- and
negation-closed, so every relation identity survives.

The best `--validate-top` candidates are then validated by real child runs
on `--holdout` fresh known-answer fixtures with `--solver` (default
pair-table); the selected candidate is the fastest one that verified every
holdout. `--spec-out` saves its recipe, a small JSON document bound to the
degree and coefficient it was found on, which `run --factor-base` replays.
Status is `complete` only with a validated winner; `--validate-top 0`
reports `unvalidated` with the census ranking alone.

The census is exact on its target set and ignores per-trial oracle cost;
the report carries `enumeration_ops_per_trial`, `pair_table_lookups_per_trial`
and `sat_variables` as cost proxies, and the validation runs measure the
wall time that combines both. A selected candidate is the best validated
observation on these fixtures, not a global optimum.

Why this matters: `ic run --degree 15 --curve-a 1` with the legacy family
collects no relation at all in 20 000 trials — the 31-point base never
reaches the order-211 subgroup with two summands — while the search finds,
scores and validates bases covering all 210 targets in a few seconds.

## Factor-base logarithms and individual-logarithm descent

    ./target/release/ic logs --degree 9 --curve-a 0 --solver pair-table --database logs9.json
    ./target/release/ic solve --degree 9 --curve-a 0 --logs logs9.json --known-log 53
    ./target/release/ic logs --degree 31 --curve-a 0 --summands 3 --factor-base fb31.json --database logs31.json
    ./target/release/ic solve --degree 31 --curve-a 0 --summands 3 --logs logs31.json --random-target --seed 7

Like a number-field-sieve pipeline, `ic` separates the two costs the
plain `run` conflates. `run` bakes the target into every relation
(`R = [a]G + [b]Q`) and rebuilds the whole relation matrix per target.
Instead:

- `ic logs` **precomputes**, once per curve and factor base, the
  discrete logarithm of every relation column — the factor-base
  logarithm database. It draws `R = [a]G` probes (no target), decomposes
  each over the factor base, rewrites it as `Σ_o c_o x_o ≡ h·a (mod r)`
  over the projected columns, and once the rows determine every column
  reads the whole logarithm vector off in one solve. Every column log is
  certified by `[x_o]G == R_o` before the database is written; a
  database that fails that check is never emitted.
- `ic solve` **descends** a target with a single relation, reusing the
  database. It draws `R = [a]G + [b]Q` until one decomposes, giving
  `h·a + h·b·d ≡ Σ_o c_o x_o (mod r)`; with the column logs known, the
  scalar `d = log_G Q` falls out of one modular inverse, and the
  recovered `d` is re-checked as `[d]G == Q` before it is returned. On
  load the whole database is re-verified against the reconstructed curve,
  so a tampered or mismatched database is rejected, not trusted.

The database is a JSON document bound to its degree, coefficient,
subgroup order and factor-base spec; `solve` rejects it on any other
curve. Only the projected column representation is used (the `ic`
default): its columns are canonical cofactor projections `R_o ∈ ⟨G⟩`, so
each column logarithm is a genuine, self-certifying discrete log.

The precomputation reaches whatever the factor base and summand count
support. On `K_0/2^31` over a search-selected dimension-11 base
(`--summands 3`, 35 columns) it completes in about 15 s; each subsequent
target then needs one or two relations. The `logs` trial budget
(`--max-trials`, default 200000) bounds the search for a full-rank
relation set; a base whose coverage cannot determine every column
reports `incomplete` rather than emitting an unverified database.

### Linear algebra: relation filtering and block Wiedemann

    ./target/release/ic logs --degree 31 --curve-a 0 --summands 3 --factor-base fb31.json --database logs31.json --linear-algebra sparse --block-size 4
    ./target/release/ic logs --degree 31 --curve-a 0 --summands 3 --factor-base fb31.json --database logs31.json --linear-algebra dense

Each relation has at most `m` nonzero entries, so the relation matrix is
sparse in exactly the way a number-field-sieve matrix is. By default
(`--linear-algebra sparse`) `ic logs` solves it the way CADO-NFS does:

1. **Filtering** — duplicate rows are dropped; a column occurring in only
   one row (a *singleton*) is removed with that row and recovered later
   by back-substitution; surplus rows beyond a small excess are removed,
   choosing rows whose removal cascades through the weight-2 columns
   (the clique rule); light columns are merged away by structured
   Gaussian elimination under a fill-in bound. Every elimination is
   recorded, so the eliminated logarithms are reconstructed exactly from
   the core solution (back-substitution, then propagation through the
   original rows, then a small dense residual if anything is left).
2. **Block Wiedemann** — the reduced core is made square by folding its
   excess rows into random earlier rows, homogenised to `M (x, 1)ᵀ = 0`,
   and a kernel vector is read off the Krylov sequence `X Mⁱ Y` of
   `block_size × block_size` blocks through a matrix Berlekamp–Massey
   step (a shifted minimal approximant basis). Only sparse
   matrix-times-block products touch the matrix, in parallel over rows.

The sparse path never attempts a solve before every column occurs in
some row, and the solution is checked against every relation before
the group certification `[x_o]G == R_o` runs. `--linear-algebra dense`
keeps the reference behaviour: full big-integer elimination after every
new relation. Both paths certify the same database (the `ic` tests
compare them); the report's `linear_algebra` object records the mode,
the attempts, the time, and for the sparse path the filtering counts
and the Wiedemann run (`core_dimension`, `sequence_length`, products).

## Running the pipeline as a resumable workflow

    ./target/release/ic workflow --params wf.json --dir runs/k0n31
    ./target/release/ic workflow --params wf.json --dir runs/k0n31 --stop-after logs
    ./target/release/ic workflow --params wf.json --dir runs/k0n31          # resumes

Like a number-field-sieve run, `ic workflow` executes the pipeline as
stages whose outputs live on disk, so a run can be stopped, inspected
and resumed without redoing finished work:

1. **select** — the factor base, either an explicit recipe or the
   best-by-census candidate of the factor-base search; written as
   `factor_base.json`.
2. **collect** — relations, in work units (see below); each unit is
   written as `relations/unit-NNNNN.json`.
3. **logs** — the units are merged, every relation re-verified in the
   group and deduplicated, and the factor-base logarithm database solved
   (`logs.json`), every column certified by `[x]G == R`. If the
   relations do not yet determine every column, further units are
   collected up to `collection.max_units`.
4. **solve** — each target descended with one relation reusing the
   database; `solutions.json` is rewritten after every target, so an
   interrupted run resumes at the first unsolved one. The pair table is
   built once per process and shared by collection and descent.

`state.json` records a BLAKE3 digest of the parameter file and each
stage's status. A rerun in the same directory reloads existing
artifacts, re-verifies them against the reconstructed curve (a stale or
tampered artifact is an error, never trusted), and continues from the
first incomplete stage; a parameter file whose digest differs is refused
so one directory never mixes two experiments. Artifacts are written
atomically. `--stop-after select|collect|logs|solve` ends the run early.

### Distributed relation collection

    ./target/release/ic workflow --params wf.json --dir runs/k0n31 --collect-units 0-3    # worker A
    ./target/release/ic workflow --params wf.json --dir runs/k0n31 --collect-units 4-7    # worker B
    ./target/release/ic workflow --params wf.json --dir runs/k0n31                        # merge, solve, descend

Relation collection is embarrassingly parallel, and the workflow splits
it the way a sieve is split into `q`-ranges. The probe scalar of trial
`t` depends only on the parameter seed and `t`, so the probe sequence
is one fixed, reproducible sequence; a **work unit** `k` is the slice
`[k · unit_trials, (k + 1) · unit_trials)` of it. Any process with the
parameter file (and `factor_base.json`, when the base came from the
search) can run a set of units with `--collect-units` — it writes
`relations/unit-NNNNN.json` for each and stops — and the files from
several workers or machines are simply placed in the run directory. The
driver without `--collect-units` collects whatever units of the first
`collection.units` are still missing itself, then merges everything
present. Inside a unit the trials run in parallel over the cores.

A unit costs what its trials cost and nothing more: the collector — whose
point index map is keyed by big integers and takes about as long to build
as a short unit takes to run — is built once for the whole stage rather
than once per unit, so splitting the same work into eight units instead
of one no longer adds to the bill.

Two setup costs alongside it were the same shape. Selecting a subgroup
base rebuilt the whole base after every batch of eight abscissae, which
is quadratic in the abscissae; it now skips the rebuilds that could only
have come back short, since an abscissa carries at most two points. And
the projected signed-orbit map — every base point multiplied by the
cofactor, then each one's whole Frobenius orbit walked — ran in
big-integer arithmetic; it now runs in single words where the field fits,
which it does for every degree this pipeline reaches. At degree 53 on a
15264-point base, selecting goes from 8.60 s to 0.94 s and the orbit map
from 8.16 s to 0.048 s, and the base and its columns are unchanged. The
pair table's 6.7 s is real `|F|²/2` work and is untouched; it is now the
dominant fixed cost, and unlike the scanning it does not grow with `r`.

A relation file carries only the probe scalar and the factor-base point
indices of each relation, bound to the parameter digest, curve, factor
base and summand count. On merge every relation is re-verified in the
group (`[a]G == Σ P_i`, exactly `m` indices, all in range) and exact
duplicates are dropped, so a corrupt or forged file cannot poison the
database — a rejected relation is counted, never used — and a file from
another run or base is ignored, not merged. Partition invariance is
tested: the union of any set of units equals the relations of a
single-process run over the same range, and `ic logs` itself now draws
its probes from the same sequence in parallel batches.

Parameters (`collection`, all optional):

    "collection":{"unit_trials":4096,"units":4,"max_units":64}

`unit_trials` probes per unit; `units` the number the driver collects
before the first solve; `max_units` the most it may collect when the
relations do not yet determine every column. The report's collect stage
lists the units present, run and ignored; the logs stage reports the
relations loaded, rejected and deduplicated and which units were used.

A parameter file (schema_version 1):

    {"schema_version":1,"name":"k0n31","curve":{"degree":31,"curve_a":0},
     "summands":3,"solver":"pair_table","seed":1,"max_trials":200000,
     "linear_algebra":{"mode":"sparse",
                       "sparse":{"wiedemann":{"block_m":4,"block_n":4},
                                 "filter":{"target_excess":32,"merge_max_weight":8}}},
     "collection":{"unit_trials":4096,"units":4,"max_units":64},
     "factor_base":{"mode":"search","family":"divisor","min_dimension":5,
                    "max_dimension":11,"targets":256,"saturate":false},
     "targets":[{"known_log":"654009"},{"random_seed":7},{"random_seed":8}]}

`factor_base.mode` is `spec` (with a recipe as written by `ic search`)
or `search` (the census search's knobs; the best candidate is taken
without child validation). Each target is a synthetic known-answer
instance: `known_log` names the scalar, `random_seed` draws one
reproducibly. Solver is `pair_table` (default), `enumerate`, `groebner`
or `sat`. `linear_algebra` is optional: `mode` is `sparse` (default;
filtering + block Wiedemann, with every knob of the sparse solver under
`sparse`) or `dense`. `collection` sizes the work units (above). The
report lists every stage with whether it ran or was reused — the collect
stage with its units, the logs stage with its verification counts and
linear-algebra statistics — and every solution with its expected and
recovered scalar.

### The ρ baseline (`vs_rho`)

    "baseline":{"rho":true,"rho_seed":5931241263826122831,"rho_max_iterations":268435456}

With `baseline.rho` the driver runs the signed-Frobenius Pollard ρ
(`koblitz_signed_frobenius_rho_with_progress`, which already carries the
Koblitz automorphism discount: it walks the `A = 2n` classes) on the
same known-answer targets, in the same process, after the descent, and
writes `baseline.json` plus a `vs_rho` block in the report with the
three timing classes the boundary ledger distinguishes:

- **charged** — per-target descent wall (one decomposition and a lookup
  against the reused database) against the ρ walk on the same target;
- **amortised** — select + collect + logs + pair-table wall divided over
  the targets, plus the descent;
- **whole process** — the precompute counted once against the ρ total.

Each comes with the ratio `ρ / IC` and a boolean verdict; a verdict is
`true` only when every target was solved *and* ρ-verified. The ledger's
`vs_rho` row asks for exactly these fields (`timing_class`,
`automorphism_discount`, both costs, `claim_boundary`); the report
carries them, but promoting a row still goes through the ledger's own
autolab and independent-replay process.

The baseline is a real opponent, not a formality, so read it first:

- It walks the `A = 2n` signed Frobenius classes, by the
  distinguished-point method — a direct-mapped cache of recent points
  for the short cycles, a sparse table of stored points for the
  long-range collision.
- Fruitless cycles (a cycle whose jumps cancel, which the negation map
  makes common) are escaped by doubling the cycle's own smallest state,
  so the walk stays a deterministic map. `charges.fruitless_cycles`
  counts them.
- Walks are stepped in batches sharing one field inversion
  (`parallel_walks`, default 32, scaled down on instances whose whole
  walk is shorter than the setup would cost).
- Its step count tracks `√(πr/2) / √(2n)`; a run far above that is a
  broken baseline, and a `vs_rho` verdict built on one means nothing.
  `rho.verified` must equal the target count — **a ρ that fails to
  recover its logarithms makes every ratio in the block meaningless**,
  which is exactly how an earlier revision of this baseline produced a
  spurious charged crossover at `n = 41`.

### Choosing a factor base

Five families are recipes (`ic search --family`): `factor`, `divisor`
and `union` are linear — an invariant subspace, a divisor of `x^e − 1`,
or a union of Frobenius translates of a seed span — and `subgroup` is
not. The linear families exist because the algebraic oracles need them:
Semaev's polynomials and the Weil descent are written over a subspace.
The pair-table oracle is a meet-in-the-middle search and needs no
structure at all, and for it the structure is a cost, not a feature.

`subgroup` draws abscissae pseudo-randomly and keeps `[h]P` for the
cofactor `h`, so every base point lies in the prime-order subgroup the
targets live in. Sums of such points cannot leave that subgroup, and the
measured decomposition rate lands on the `|F|³/(3!·r)` a random base
would give — at degree 41, one target in 28 against one in 273 for the
subspace union of the same size, with the same column count. Selection
uses no logarithm: `[h]P` is in the subgroup whatever `P` is.

Use it with the pair-table solver. `groebner` and `sat` need a linear
domain and will refuse.

**How big?** With a subgroup base the work per target falls as
`r/|F|²` — four times cheaper per doubling — while the pair table grows
as `|F|²`, so the answer depends on how many targets share the database.
At degree 41: 5248 points costs 3.0 s of precompute and 16.2 ms a
target, 10496 points costs 9.0 s and 7.2 ms, and the two cross at about
665 targets. `docs/ic/runs/koblitz-base-size-20260912.json` has the
sweep and both end-to-end runs. Past about 10500 points at that degree
the decomposition rate saturates and further growth only makes each
trial dearer.

`PairSumTable::build` keeps a base inside 4 GiB, and past that budget it
changes representation rather than refusing. The full table stores each
pair as `(packed sum, i, j)`, sixteen bytes; the **compact** one stores a
bucketed hash of the sum in a `u32` — never a false negative, a false
positive about one time in `2²⁸` — and recovers the summands of a hit by
one `|F|`-long scan, since `target − P_i` is a base point exactly when
`i` is a summand. That scan asks the group rather than the table, so a
false positive costs an empty scan and never a wrong answer. Hits are
rare, so it is paid about once per relation rather than once per probe.

At a fixed budget `B` the base is `|F| = √(2B / bytes per pair)`, and the
descent needs `2r/|F|²` probes, so the width of a stored pair is
proportional to the descent's cost. Four and a half bytes instead of
sixteen is a base of 42302 points instead of 23169 at 4 GiB, and 3.3
times fewer probes. Below the budget nothing changes: a base that fits
with its summands keeps them. Below even the compact size the build
still refuses with a number instead of an allocation.

`docs/ic/params/k0n61-subgroup-wide.json` is the largest rung this family
offers: `K_0/F_{2^61}`, a 48-bit subgroup, `r = 162 888 033 982 417`, on
a 36112-point compact base. It solves 32 of 32 with a 53.1 ms descent
against ρ's 3.248 s — charged 61.2, amortised 1.29, and all three
verdicts true at 32 targets.
`docs/ic/runs/koblitz-degree61-20260912.json` records it, and says what
it does to the earlier reach projection: that projection put the charged
ratio at 5.1 at 48 bits and wanted 315 targets, because it was measured
on constants that have since moved three times.

`docs/ic/params/k0n53-subgroup-wide.json` is the degree-53 rung on a
36464-point compact base: the descent falls from 50.4 ms a target to
16.3 ms and the charged ρ/IC ratio rises from 25.0 to 75.2, while the
precompute rises from 16.5 s to 87.8 s. That is a trade, and
`docs/ic/runs/koblitz-compact-pair-table-20260912.json` prices it — about
2200 targets before the wider base is the cheaper one.

**`descent_summands`** lets the descent ask for a different number of
summands than collection, which shares only the base and its pair table.
Collection wants few probes (each costs two scalar multiplications) so
it takes three; the descent walks its probes by `+G` in blocks sharing
one inversion, which makes a probe cheaper than the lookup after it, so
two wins. At degree 41 that is 14.35 ms a target against 8.63 on a
5248-point base, and 7.20 against 4.65 on a 10496-point one.

### Ledger rungs, ready to run

`docs/ic/params/k0n{31,37,39,41}.json` are the four Koblitz rungs of the
boundary ledger as parameter files — 32 known-answer targets each, the
ρ baseline on, collection units sized to the base:

    ./target/release/ic workflow --params docs/ic/params/k0n41.json --dir /tmp/n41

`docs/ic/params/k0n53-subgroup.json` is the largest subgroup this family
offers (44 bits, 38× the degree-41 rung). It solves 32 of 32 with a
49 ms descent against ρ's 1.216 s, and
`docs/ic/runs/koblitz-degree53-and-reach-20260912.json` carries it
together with the projection of where the advantage runs out: the
charged class survives to about a 52-bit subgroup, and the whole-process
class needs a batch of targets that grows with `r` (16 at 40 bits, 96 at
45, 834 at 50). Relation collection is what stops it, not the descent.

`docs/ic/params/k0n{31,37,39,41}-subgroup.json` are the same four rungs
with subgroup bases; `docs/ic/runs/koblitz-subgroup-bases-20260912.json`
records them, and the charged ρ/IC ratio there crosses 1 at degrees 37,
39 and 41 (8.3 at 41). Read that file's `what_this_is_not` before
quoting it — in particular, its degree-41 whole-process verdict is a
bulk statement about 32 targets, not a single-instance one.

`docs/ic/runs/koblitz-scaling-20260911.json` records two consecutive
series of all four, with per-stage timings, filter and Wiedemann
statistics, and both ρ and IC verification counts. No rung crosses: the
charged ρ/IC ratio is below 1 at every one. Read its
`what_this_is_not` before quoting any number from it.

## Random fixtures and custom parameters

    ./target/release/ic generate --degree 11 --curve-a 1 --seed 42 --out fixture.json
    ./target/release/ic inspect --file fixture.json --json

Generation randomizes the known scalar and target, not every field or
curve coefficient. The degree/coefficient select the curve; the existing
constructor chooses its field polynomial and subgroup generator.

A parameter document has schema_version 1 and these fields:

- name
- field
- a and b
- subgroup_order and cofactor
- optional generator and point coordinate pairs
- optional fixture containing known_log and seed

All integers except degree, polynomial_terms, schema_version, and seed
are strings, in decimal or hexadecimal with a 0x prefix. Unknown fields
are rejected. The input file limit is 1 MiB.

Field forms:

    {"kind":"prime","modulus":"97"}
    {"kind":"binary","degree":9,"polynomial_terms":[9,1,0]}
    {"kind":"binary_abstract","degree":131,"representation":"normal-basis metadata"}

Polynomial terms include the leading degree and constant term; they
must be unique. The inspector verifies irreducibility using polynomial
GCD and Frobenius identities. Abstract binary profiles accept only
Koblitz base-field coefficients a in {0,1}, b=1, and no coordinates.

See prime-example.json and binary-example.json for complete examples.

A factor-base recipe has schema_version 1, degree, curve_a, and a spec:

    {"schema_version":1,"degree":15,"curve_a":1,"spec":{"kind":"divisor","indices":[0,2]}}

Spec kinds are factor, divisor, frobenius_union, two_torsion_saturated, and
pruned (a parent spec plus the canonical abscissa of every retained signed
orbit). A recipe is rejected on any other curve.

## Meaning of validation

The inspector reports every check separately:

- strict numeric and canonical coordinate ranges;
- nonsingularity;
- exact irreducibility for a supplied binary polynomial;
- fixed-base Miller-Rabin screens for prime fields and subgroup orders;
- supplied points on the curve and annihilation by the declared order;
- supplied known-answer identities;
- exact Koblitz cardinality using an integer trace recurrence;
- the necessary Hasse bound for a declared cardinality.

Miller-Rabin passes are not primality proofs. An annihilation check does
not independently prove exact point order. Hasse's bound is necessary,
not a cardinality certificate. General prime/binary point counting,
missing points, and unsupported coordinate representations are reported
as not_checked. The aggregate checks_passed status means none of the
performed checks failed, not that all possible properties were proved.

Binary inspection supports degrees 2 through 571. Prime-field inspection
supports moduli at most 512 bits. Curve names are descriptive; they do
not change these validation rules.

## Comparing factor bases

    ./target/release/ic compare --degree 7 --curve-a 1 --samples 3 --holdout 2 --seed 42 --out comparison.json

Comparison uses the degree-ord_n(2) irreducible-factor family already
implemented by the materialized builder. It is not a search over every
possible factor base — `search` is. It launches one child process at a
time, each using only generated known-answer inputs. All candidates
receive the same training fixtures and solver settings.

A candidate must complete and verify every training run before it is
eligible. The training winner has the smallest observed median process
wall time; it must then verify separate generated holdout fixtures.
Failed and timed-out attempts remain in the report. If no candidate
passes, selected_factor_index is null and the result is inconclusive.

Comparison limits are 16 candidates, 1..8 training samples, 1..8 holdout
samples, and a per-child watchdog of 1..120 seconds (default 30).
Standalone run limits are algorithmic; the comparison watchdog applies
to child runs launched by compare.

Selection time includes process startup, curve and factor-base
construction, relation collection, matrix solving, verification, and
report emission. Child completion is polled every 10 ms; very small
timing differences should not be interpreted as meaningful. The report
also includes the total comparison time, so selection and unsuccessful
candidate costs remain visible. A selected candidate is a bounded
observation, not a global optimum, scaling claim, or challenge result.

## Reports and resource accounting

Add --json for machine-readable output or --out PATH to save the JSON
report while retaining the human progress display. Output files are
created exclusively: an existing path is never overwritten.

Run, comparison and search reports include normalized parameters,
known-answer fixtures, seeds, solver settings, status, counts, stages,
elapsed time, platform, package version, and a BLAKE3 fingerprint of the
executable. Run reports also carry the factor-base recipe and its size,
per-stage timing, and the relation accounting (independent, dependent and
inconsistent relations). Inspection reports retain normalized inputs and
the input-file hash. Generated parameter documents remain directly
importable.

On macOS and Linux, CPU time and peak resident memory are sampled with
getrusage before report emission. Peak RSS is normalized to bytes.
Each comparison or validation run is a fresh process; its resource
counters therefore do not inherit earlier candidates' high-water marks.
Parent counters cover the parent only. Unsupported measurements are null,
never zero. A watchdog-killed child has no complete resource report.

Exit status is zero for completed synthetic operations and inspections
whose performed checks pass. Invalid inputs, incomplete experiments,
inconclusive comparisons and searches, and output-file errors are
unsuccessful. Clap usage errors use its standard nonzero exit status.

## Verification

    cargo test --release --test ic_framework --test ic_progress
    cargo test --release --lib koblitz_

Tests cover named profiles, custom prime curves, generated-fixture
round trips, reproducibility, malformed and ambiguous parameters,
resource reporting, a degree-11 synthetic run under both accountings,
every solver with three summands, comparison eligibility, separate holdout
inputs, a search that succeeds where the legacy family yields nothing, the
binding of recipes to their curve, and exclusive artifact creation. The
library tests cross-check the pair table, the exact census, orbit pruning
and the incremental relation solver against exhaustive search and the
dense modular solver.
