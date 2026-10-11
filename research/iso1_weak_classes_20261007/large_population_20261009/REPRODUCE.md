# Replaying the cubic population and quadratic scope controls

Build the focused native package, keeping the frozen run directories intact:

```sh
study=research/iso1_weak_classes_20261007
cargo build --release --locked --manifest-path "$study/report_check_runtime/Cargo.toml" --bins
```

The native checker validates source freezes, raw process hashes, exact model
identities, summaries, and family-specific labels. Run it on each archive:

```sh
pop="$study/large_population_20261009"
for run in validation_p7_run1 pilot_run1 positive_controls_run1 population_run1; do
  "$study/report_check_runtime/target/release/iso1_population_check" "$pop" "$pop/$run"
done
```

Replay the saved coordinates and every retained route with compiled PARI/GP:

```sh
ISO1_POPULATION_STUDY="$pop" gp -q -f -s 256M "$pop/replay_models.gp"
ISO1_POPULATION_STUDY="$pop" gp -q -f -s 256M "$pop/quadratic_replay.gp"
gp -q -f -s 256M "$pop/structure_exact.gp"
gp -q -f -s 256M "$pop/torsion_exact.gp"
```

The first replay reconstructs 1,326 records and 72 isogeny edges and emits
75 constructive Hilbert–90 coordinate maps for retained positive cubic
models. The second independently enumerates all affine x values, plus two
infinity points, for each of the three frozen quadratic-family controls.
Their exact traces are −38, −10, and 610.

The local population driver uses the recorded Homebrew GP executable path.
On this host, a fresh independently recorded run has the following command;
choose a new output directory rather than reusing a frozen one:

```sh
"$study/report_check_runtime/target/release/iso1_population_run" "$pop" /private/tmp/iso1-population-fresh population
```

The native driver records its binary, GP binary, protocol and source hashes;
it runs four GP workers and enforces the registered caps. `small`, `pilot`,
and `controls` are separate phases, with their own new output directories.
Linux CI replays fixed field moduli and records with its own compiled GP;
it does not substitute a Linux execution for the archived macOS run.

The class labels and admission rates refer to the cubic norm-one family.
The quadratic scope controls are a separate constructed cohort. See
[the prior-work comparison](PRIOR_WORK.md) for the corrected interpretation
of the full Joux–Vitse family.
