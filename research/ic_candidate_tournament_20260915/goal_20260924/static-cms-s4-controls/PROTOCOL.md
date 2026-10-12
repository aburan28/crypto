# Source-receipted static CryptoMiniSat on six disclosed n17 controls

The [exact-yield audit](../generic-exact-yield-audit/RESULTS.md) found four
n17a1 ordinary queries with exact three-point decompositions that the F5
worker missed under its registered bound (trials 4, 67, 71, 75). Freeze all
four as known-positive correctness controls and trials 0 and 3 as
exact-negative controls in the [panel](panel.json). These six points were
disclosed before this pilot. Their outcome-selected allocation **does not**
estimate natural relation yield, and they are distinct from the four n19a0
points in the earlier [loader-failed pilot](../cms-s4-controls/RESULTS.md).
No fresh target or sealed confirmation set is opened.

The question is whether the repository's *wide symmetrised S4* XOR-DIMACS
export, solved by an independently source-receipted static CryptoMiniSat
5.14.7 build, can return a valid group decomposition for these four known
positives. This differs from the worker's chained-S3 F5 system, so it is a
solver/encoding feasibility probe, not a paired timing comparison.

Before any exporter or SAT instance, verify the panel hash, the parent exact
result and worker archive, all six scalar-to-point calculations and exact
labels, the pinned exporter source, and the Phase-B terminal bundle that
contains the completed static CryptoMiniSat build receipt. Require that the
receipt binds the executable hash and the declared source and dependency
commits. Copy that executable into the final output directory and run a
metered `--version` **from that copied path** before any instance dispatch.
The copy and its source path must have the registered digest; no `@rpath`
Homebrew dynamic library may be required. Retain the copied executable,
receipt, verifier result, preflight process streams, and a host record.

Build `koblitz_pdp_export` from the exact admitted Rust/dependency source.
Run it with `--export-only`, the six public points and a public nonce. It may
not receive the exact labels. Check field, curve, factor-base geometry,
S4 representation, unused internal solvers, target and all three exported
source-system files before starting SAT. Run the copied static CryptoMiniSat
once per valid export, in registered order, with one thread, random seed 1,
one returned model, at most one million conflicts, and a 120-second watchdog.
The exporter has a separate 60-second watchdog. Preserve all failed attempts.

For SAT, parse the model and verify the full XOR-DIMACS formula. Decode the
first 18 source bits through the recorded basis, enumerate rational point
lifts in the verified geometric base, and independently group-check the
three-summand sum against the supplied point. A formula model without a
point witness is `SOURCE_MODEL_NONLIFTING`; a witness on an exact-negative
control or UNSAT on an exact-positive control is a contradiction requiring
investigation. A timeout, solver error, or partial model proves neither SAT
nor UNSAT. Scalar replay is for fixture/evidence verification only, outside
the stage process timings.

Stop after the six registered attempts. The no-retry rule also applies if a
preflight or exporter fails. Archive the exact processes and publish the
outcome without choosing further points in this panel. A positive result
can motivate a separately registered integration into full natural relation
collection, rank/linear algebra and one-target recovery. This pilot itself
establishes no complete IC solve, online time or speedup.
