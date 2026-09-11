# ic research framework

The standalone ic executable inspects elliptic-curve parameters and runs
bounded, reproducible index-calculus experiments on internally generated
known-answer Koblitz instances. Imported points are used only for mathematical
validation.

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

The following knobs are recorded in each run report:

- degree: odd values 3 through 41; a curve is usable only when its largest
  prime factor exceeds the cofactor, which above 23 holds for degrees 29
  (curve-a 1), 31, 37, and 39 (curve-a 0);
- curve-a: 0 or 1, with b fixed to 1;
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
