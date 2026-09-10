# ic research framework

The standalone ic executable inspects elliptic-curve parameters and runs
bounded, reproducible index-calculus experiments on internally generated
known-answer Koblitz instances. Imported points are used only for mathematical
validation.

**Agent scoreboard:** per-stage records and next targets to beat live in
[`BOUNDARY_TARGETS.md`](./BOUNDARY_TARGETS.md) and
[`boundary_targets.json`](./boundary_targets.json) (binary, Koblitz, prime).

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

The default run uses K_0 over GF(2^9) and known logarithm 53. The solver
constructs the target from the declared known answer, then checks both the
recovered value and the group identity. No public point file can be passed
to the run command.

The following knobs are recorded in each run report:

- degree: odd values 3 through 23, within the existing library cap of 24;
- curve-a: 0 or 1, with b fixed to 1;
- known-log: a positive scalar smaller than the selected subgroup order;
- random-target: draw a known-answer scalar reproducibly from seed;
- seed: seed for the generated fixture and relation sampler;
- factor-index: candidate in the existing materialized factor-base family;
- max-trials: 1 through 20000;
- solver: groebner, sat, or enumerate.

Not every degree/coefficient combination has a usable subgroup. A valid
curve does not guarantee successful collection or an invertible relation
system. Factor-base dimension is capped at 12 before materialization.
Fixture generation and inspection do not require a usable factor base.

The run display follows factor-base construction, relation collection,
linear algebra, and verification. Every synthetic run requires the
relation matrix and uses fixed-surplus collection. Incomplete results
exit unsuccessfully and remain incomplete in JSON reports.

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
possible factor base. It launches one child process at a time, each using
only generated known-answer inputs. All candidates receive the same
training fixtures and solver settings.

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

Run and comparison reports include normalized parameters, known-answer
fixtures, seeds, solver settings, status, counts, stages, elapsed time,
platform, package version, and a BLAKE3 fingerprint of the executable.
Inspection reports retain normalized inputs and the input-file hash.
Generated parameter documents remain directly importable.

On macOS and Linux, CPU time and peak resident memory are sampled with
getrusage before report emission. Peak RSS is normalized to bytes.
Each comparison run is a fresh process; its resource counters therefore
do not inherit earlier candidates' high-water marks. Parent counters
cover the parent only. Unsupported measurements are null, never zero.
A watchdog-killed child has no complete resource report.

Exit status is zero for completed synthetic operations and inspections
whose performed checks pass. Invalid inputs, incomplete experiments,
inconclusive comparisons, and output-file errors are unsuccessful.
Clap usage errors use its standard nonzero exit status.

## Verification

    cargo test --release --test ic_framework --test ic_progress

Tests cover named profiles, custom prime curves, generated-fixture
round trips, reproducibility, malformed and ambiguous parameters,
resource reporting, a degree-11 synthetic run, comparison eligibility,
separate holdout inputs, and exclusive artifact creation.
