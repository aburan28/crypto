# Native ISO-1 census, audit, and continuation worker

The absolute-Frobenius quotient needs exactly `(p^4+3*p^2+8)/12`
point counts for `p > 3`. Its orbit weights recover all `2*p^4+2*p^2`
normalized weak representatives and both trace signs. The standalone
[theorem](../THEOREM.tex) proves the quotient and the norm-one conductor
obstruction. The retained native arithmetic kernels are extracted from
commit `62becf9572fe74cbe8b3d8cebee3bf8a240708a1`; see
[source provenance](SOURCE_PROVENANCE.txt). This focused Cargo package
builds those kernels independently of other library modules.
The [field receipt](../census_field_receipt.txt) records explicit moduli
and Frobenius multipliers. At p = 53 the tower is
`u^2=2`, `theta^3=1+u`, and `theta^(p^2)=(26+24*u)*theta`.

From the repository root:

```sh
cargo build --release --locked --manifest-path research/iso1_weak_classes_20261007/census_runtime/Cargo.toml --bins
cargo test --release --locked --manifest-path research/iso1_weak_classes_20261007/census_runtime/Cargo.toml --bins
RAYON_NUM_THREADS=4 research/iso1_weak_classes_20261007/census_runtime/target/release/iso1_class_census 53 p53.csv --derive-twists --absolute-orbit-quotient
ISO1_P=53 gp -q -f research/iso1_weak_classes_20261007/gp_large_prime_control.gp > p53_gp.csv 2> p53_gp.stderr
test ! -s p53_gp.stderr
research/iso1_weak_classes_20261007/census_runtime/target/release/iso1_census_audit p53.csv samples p53_gp.csv
```

The default trace counter validates a Hasse-interval baby-step/giant-step
result on two random points. Its complete trace rows remain probabilistic
computational labels. Independent GP controls test positive membership;
weighted totals, twist symmetry, Hasse bounds, discriminants, conductors,
orbit call counts, and the proved congruence detect structural failures.
An audit failure preserves the attempted result and stops the queue.

For independent exact counting at `p <= 13`, add `--exact-squares`.
This builds the field-square table once and counts every field element
for every selected orbit. The completed p = 7 run and complete GP orbit
census agree on every weighted trace row:

```sh
RAYON_NUM_THREADS=4 research/iso1_weak_classes_20261007/census_runtime/target/release/iso1_class_census 7 p7.csv --derive-twists --absolute-orbit-quotient --exact-squares
ISO1_P=7 gp -q -f research/iso1_weak_classes_20261007/gp_full_absolute_orbit.gp > p7_gp.csv 2> p7_gp.stderr
test ! -s p7_gp.stderr
research/iso1_weak_classes_20261007/census_runtime/target/release/iso1_census_audit p7.csv orbits p7_gp.csv
```

The queue takes an output directory and an ordered prime list:

```sh
research/iso1_weak_classes_20261007/census_runtime/target/release/iso1_census_queue /private/tmp/iso1-larger-primes-next 59 61 67
```

It runs sequential primes with four census workers, runs 5,000 independent
GP controls per completed prime, validates each CSV, and saves command,
source, executable, environment, exit-status, and checksum receipts. It
compresses successful CSVs with zstd, tests the archive, and checks that
the decompressed SHA-256 equals the raw CSV before removing its own raw
copy. `queue_status.tsv` records starts, completions, and failures.

An already running p = 53 census can be followed using
`--after-pid PID CENSUS.csv CONTROL.csv` before the remaining primes.
The worker waits for that process to exit, then requires the completed
CSV and its independent controls to pass the audit before proceeding.
This mode needs permission to inspect the specified process with `ps`.
Do not overwrite its frozen census or audit executables while it runs.

Before each prime, the storage guard requires one GiB plus twice
`140*(p^3+1)` bytes, a conservative raw-CSV estimate. A refusal is
recorded as `STOPPED_STORAGE_PREFLIGHT`; later primes remain unfinished.
The [continuation protocol](../CONTINUATION_PROTOCOL.md) preserves the
requested complete prime range. Wall seconds in receipts account for
resources on a contended host; the call-count theorem supplies the
algorithmic reduction.

## Resumable runs (2026-10-16)

Set `ISO1_CHECKPOINT=FILE` on `iso1_class_census` or on
`iso1_lambda_census` (with `ISO1_OMEGA_ONLY=1`) to make a run resumable. Each
finished job or block appends one line to `FILE`. A rerun with the same
arguments and `FILE` reloads those lines, skips their work and writes the
same output. Each job seeds its own RNG, so results do not depend on
scheduling. A torn final line from an interrupted write is ignored. Verified
by killing runs part-way and resuming:

- A p = 37 certified census, killed after 1,876 of 2,738 jobs, resumed to a
  CSV byte-identical to `p37_twist_derived.csv`. Its point-count and F_p
  multiplication totals were identical to the uninterrupted run's.
- A p = 13 Ω census, killed after 20 blocks, resumed to a CSV byte-identical
  to `../residual_model_20261016/data/om13.csv`.

The certified queue (`../certified_20261016/queue.sh`) uses this from p = 71
on, because the cloud container that runs it restarts every few tens of
minutes.
