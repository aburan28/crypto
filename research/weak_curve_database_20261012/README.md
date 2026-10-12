# Weak-curve evidence database

This release makes the existing weak-curve research queryable at three distinct
levels: exact elliptic models, ordinary isogeny classes, and mechanisms that
can transfer a specified subgroup problem. The SQLite database contains
**339,740 class observations**, **814 exact synthetic models with 814 linked
property records**, **21 family or condition entries with 21 proof obligations**, **12 construction
recipes**, **10 theorem entries**, **4 calibrated or refuted heuristics**, and
**8 hypothesis entries**. The class rows split by mechanism and evidence
grade: 87,486 complete exact two-branch ordinary rows for
p=7,11,13,17,19,29,37, and 252,254 earlier cubic norm-one diagnostic rows
at p=41,43,47. The latter include nonordinary Hasse candidates and use a
randomized point counter; zero labels are diagnostic. These sets have
different branch scopes and are never pooled as a universal prevalence
estimate.

The model inventory includes 537 records over 192-, 204-, 216-, 228-,
240- and 252-bit fields. Of these, 384 are unclassified by the present
genus-three model index and 153 are outside that index's domain;
selected-subgroup advantage is unmeasured for all 537. At 17 bits,
24 exact models carry a positive cubic or quadratic covering-model label.
The model count therefore measures indexed identities and controls;
covering support and selected-subgroup outcomes are separate database fields.

[Open the database](weak_curves.sqlite), [read the schema](schema.sql), or
[inspect the class graph](class_counts.svg). The graph is generated from
[CLASS_SUMMARY.tsv](CLASS_SUMMARY.tsv), with the two scopes labeled separately.
[CONSTRUCTION.md](CONSTRUCTION.md) gives a complete candidate construction
path and a recorded p7 model; [construction_route.svg](construction_route.svg)
shows the witness gates. [THEOREMS_AND_HEURISTICS.md](THEOREMS_AND_HEURISTICS.md)
separates proof, necessary screen, statistical model and falsified rule.
[PROOF_MAP.md](PROOF_MAP.md) links Lean declarations to their exact scope.
[cubic_identity.certificate.json](cubic_identity.certificate.json) is the
native certificate `IDC1h3561e7fee541a297` for the denominator-cleared
algebraic identity; the Lean theorem states the stronger group implication.
[SOURCE_AUDIT.md](SOURCE_AUDIT.md) records the published-class import backlog.

## Query contract

The primary key for a model is its exact ICV1 identity. A class row is keyed
by (characteristic, trace, mechanism_scope); its field and trace domain are in
the associated input receipt. An isogenous representative shares a trace but
needs its own model and route witness. A class-support label does not imply
that every model in the class has the covering. A curve model with a cover
does not by itself establish a cost advantage for any selected subgroup.

The main tables are family, curve_model, curve_property, class_observation,
theorem, heuristic, conjecture, mechanism, family_proof_obligation,
source_ref, and input_file. Each
imported row retains its source file ID; class rows also retain the source
line. Every source file is SHA-256 locked in [INPUT_LOCK.json](INPUT_LOCK.json).
Generated database, visuals, report, certificate and Lean source hashes are
in [ARTIFACTS.sha256](ARTIFACTS.sha256).
The build checks exact trace domains, input hashes, model/trace/order
arithmetic, unique identities, p7 oracle agreement and model/property
one-to-one correspondence. The check command reimports into a temporary
database and compares a canonical digest of every table row.

Example queries from this directory:

    sqlite3 -readonly weak_curves.sqlite "SELECT characteristic,mechanism_scope,rows,positive_ordinary_rows,zero_ordinary_rows FROM class_summary ORDER BY characteristic"
    sqlite3 -readonly weak_curves.sqlite "SELECT icv1,trace,curve_order,subgroup_advantage FROM curve_model WHERE cubic_model_label='yes' ORDER BY icv1 LIMIT 10"
    sqlite3 -readonly weak_curves.sqlite "SELECT id,name,validation_status FROM mechanism ORDER BY id"
    sqlite3 -readonly weak_curves.sqlite "SELECT id,lean_status,open_bridge FROM theorem ORDER BY id"
    sqlite3 -readonly weak_curves.sqlite "SELECT family_id,proof_state,next_bridge FROM family_proof_obligation ORDER BY family_id"

## Build and verification

Run from this directory:

    cargo test --release --locked
    cargo run --release --locked -- build
    cargo run --release --locked -- check
    cargo run --release --locked -- summary
    cd lean
    lake build WeakCurveDatabase
    lake env lean AxiomAudit.lean
    lake env leanchecker --fresh --verbose WeakCurveDatabase
    cd ../..
    cargo run --release --example weak_curve_cubic_identity -- check research/weak_curve_database_20261012/cubic_identity.certificate.json

The database is built from fixed public research files in this repository.
The Rust importer accepts only build, check, and summary; it accepts no
external curve, key, target, or input path. The theorem project pins Lean
4.34.1 and Mathlib commit d13f23b723b8a846827a245b89c10fc7d3f11612.
Local validation results and retained development failures are in
[VALIDATION.md](VALIDATION.md).

## Scope of the catalogue goal

The database presently indexes the exact local corpus and the published
mechanism registry, with attribution. The 2013 Iijima–Momose–Chao
classification reports additional covering cases; its tables need
normalization against the 21 entries before those cases become records.
A claim about being the largest database needs a comparable census of other
public catalogues under the same definition of a distinct model, class and
family. The next import priority is the published Type I/II case table,
then odd-prime polarized covers and complete p71/p73/p79 two-branch class
oracles. The [protocol](PROTOCOL.md) fixes evidence and sampling gates.
