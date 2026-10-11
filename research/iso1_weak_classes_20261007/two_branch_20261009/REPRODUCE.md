# Reproducing the two-branch study

Run from the repository root. Arithmetic runs use the installed compiled
PARI/GP and native Rust tools; the selected Cargo workspace imports the existing
SHA and identity-certificate modules without the sparse checkout's missing
root library. The GP process runner uses `/opt/homebrew/bin/gp` on this host.
For another host, invoke the GP scripts with its installed `gp` as below.

```sh
cargo test --release --locked --manifest-path research/iso1_weak_classes_20261007/report_check_runtime/Cargo.toml --bins
cargo test --release --locked --manifest-path research/iso1_weak_classes_20261007/report_check_runtime/Cargo.toml --lib
cargo build --release --locked --manifest-path research/iso1_weak_classes_20261007/report_check_runtime/Cargo.toml --bins
```

The following independent native replay checks raw stdout/stderr hashes,
archived arithmetic sources and loaded-launcher source bindings, both census
branches, all attempted construction statuses, CM/control agreement, the
512-source counts, and the historical cubic row counts at p7, p13 and p37.
Its generated summary and CSVs are derived views of the raw receipts.

```sh
research/iso1_weak_classes_20261007/report_check_runtime/target/release/iso1_two_branch_check research/iso1_weak_classes_20261007/two_branch_20261009
research/iso1_weak_classes_20261007/report_check_runtime/target/release/iso1_two_branch_algebra research/iso1_weak_classes_20261007/two_branch_20261009/identity_certificates.json
research/iso1_weak_classes_20261007/report_check_runtime/target/release/iso1_two_branch_exact research/iso1_weak_classes_20261007/two_branch_20261009/p7_run1/stdout.txt
```

The complete p7 arithmetic replay needs 48,118,441 affine-x evaluations.
The two-torus GP enumeration emits exact field moduli, coefficient vectors,
orbit sizes, point counts, all ordinary class labels and a completion marker.
Requiring empty stderr is necessary: GP can exit with code zero after a script
error. The native wrapper requires both empty stderr and the completion marker.

```sh
ISO1_P=7 ISO1_BRUTE=1 gp -q -f -s 256M research/iso1_weak_classes_20261007/two_branch_20261009/torus_census.gp
ISO1_P=37 ISO1_BRUTE=0 gp -q -f -s 256M research/iso1_weak_classes_20261007/two_branch_20261009/torus_census.gp
gp -q -f -s 256M research/iso1_weak_classes_20261007/two_branch_20261009/algebra.gp
gp -q -f -s 256M research/iso1_weak_classes_20261007/two_branch_20261009/cm_controls.gp
gp -q -f -s 256M research/iso1_weak_classes_20261007/two_branch_20261009/cm_microfactor.gp
gp -q -f -s 256M research/iso1_weak_classes_20261007/two_branch_20261009/cm_degree6.gp
```

The source fixtures are encoded explicitly in `p7_sources.gp`,
`large_sources.gp`, and `general_sources.gp`, independently of generation
seeds or the GP version's default field polynomial. They retain the original
source moduli, roots and orders from the prior population panel. Existing
run directories must not be overwritten; use a fresh output path for a new
execution. Every native process receipt preserves a cap or failure.

```sh
ISO1_INPUTS=research/iso1_weak_classes_20261007/two_branch_20261009/p7_sources.gp gp -q -f -s 256M research/iso1_weak_classes_20261007/two_branch_20261009/construct.gp
ISO1_INPUTS=research/iso1_weak_classes_20261007/two_branch_20261009/replay_inputs.gp gp -q -f -s 256M research/iso1_weak_classes_20261007/two_branch_20261009/replay.gp
ISO1_INPUTS=research/iso1_weak_classes_20261007/two_branch_20261009/general_replay_inputs.gp gp -q -f -s 256M research/iso1_weak_classes_20261007/two_branch_20261009/replay_general.gp
ISO1_INPUTS=research/iso1_weak_classes_20261007/two_branch_20261009/general_controls.gp ISO1_GENERAL_CONTROL_INPUTS=research/iso1_weak_classes_20261007/two_branch_20261009/general_control_inputs.gp gp -q -f -s 256M research/iso1_weak_classes_20261007/two_branch_20261009/construct_general.gp
```

The construction variants are frozen separately. The adaptive protocols
document all selected fixtures and changed edge orders. In the prioritized
follow-ups the inherited raw `edges7` column is an aggregate of odd degrees
other than 3 and 5; exact route degrees remain in every kernel row. In the
initial degree-2/3 runner `tested` is the queue cursor, one past the tested
count on a component closure. The main report uses `visited` and exact edge
counts. The `odd_run2` diagnostic retains its earlier aggregate counter; the
corrected `odd_run3` provides the degree-specific validation result.

Seven identity records and native mutation controls accompany the mathematical
proofs. Record `run_v1_snapshot.txt` is the actual launcher source for the
first panel. Fifteen original source-observation mismatches are explicitly
corrected in `provenance_correction.json`; original receipts are preserved.
Later runs archive their matching driver source automatically. The complete
evidence seal checks all retained files, excluding itself and the transient
TeX log.

```sh
cargo run --release --locked --manifest-path research/iso1_weak_classes_20261007/report_check_runtime/Cargo.toml --bin iso1_two_branch_seal -- research/iso1_weak_classes_20261007/two_branch_20261009 --check
```

The editable report is `REPORT.tex`. Its tables and SVGs are generated from
`summary.json` by `iso1_two_branch_present`; `render_svg.swift` produces the
complete SVG viewport with native WebKit. Quick Look clipped wide previews on
this host, so its initial thumbnails were replaced by the full render.

```sh
tectonic --keep-logs --outdir research/iso1_weak_classes_20261007/two_branch_20261009 research/iso1_weak_classes_20261007/two_branch_20261009/REPORT.tex
```
