# Native n17 prepared public-point boundary, version 1

This source change follows ordinary-controller PR #1349. It closes a specific
interface gap: the consumed control adapter only accepts its old preparation
seal and disclosed target, and its independent verifier runs outside the
producer interval. Those controls and their archived sources remain closed.

The new `PreparedN17Target::from_ordinary_math` interface accepts the target-free
mathematical input emitted by native ordinary preparation. It reconstructs the
exact educational n17 curve, 63 geometric points, 62 usable cofactor images,
29 signed/Frobenius columns and log table. Every column is replayed with both
the producer and the existing independent shift-and-add arithmetic. It accepts
no target, imported historical state seal or known target scalar. Table replay
alone cannot admit preparation provenance, natural yield or a complete solver.
The outer controller must bind the original source-frozen preparation audit.

The F5 source entry accepts one supplied affine public point. It builds the
target-independent solver context before entering the online interval. The
fixed source uses three summands, MatrixF5 degree 3, node budget 8192, at most
eight queries, no direct collisions and the observed-query API. Range, curve
and subgroup checks occur inside the online interval before solver invocation.
Every failed query remains in the report. Independent witness replay and
independent final scalar replay occur before the interval closes. Five IC
online phases retain their exclusive clock costs; rho remains absent. Missing
phases on failed runs stay null rather than becoming inferred zeros.

## Source controls, before any scientific invocation

Hypothesis: retained mathematical data can cross the new preparation/target
boundary without loosening an old controller, and a corrupt witness, scalar,
base, log, target, missing failure or unexplained early stop cannot become a
verified result. The reference is the retained disclosed n17 preparation and
three-query target transcript. Its Python provenance and previously measured
costs remain historical; the controls neither replay the solver nor reuse the
old registration. Synthetic phase markers in a transcript control establish
clock closure only and are not F5 or relation-yield measurements.

Success requires all constructor, corruption, failure preservation, public
point, bounded plan and in-interval verification controls to pass. Stop on any
failed assertion or compile error, retain the failure, fix the source and rerun
only these source controls. No target generation or new solver search is part
of this protocol. Use Rust through the repository's busy wrapper; ordinary
host test durations are not isolation receipts or comparative measurements.

## Scientific gates still pending

The source API itself grants no dispatch authority. A separate versioned
one-use execution controller and published external seal must bind the actual
prepared state, independently admitted preparation and every executed source
(including `src/bin/icprog/oracle.rs` and `json.rs` imported unchanged), native
build/dependencies, resource limits, target exposure and all original receipts.
Source/library callback reports cannot supply that custody. The accepted
external CryptoMiniSat target frontend still requires separate integration;
the built-in Rust SAT engine is not a replacement.

Only after both exact pipeline preparations and new complete solver adapters
pass their source-bound audits may a new frozen native `ecbench` protocol
introduce fresh paired one-target workloads, the incumbent and strong same-point
rho. Canonical EC1/IC1/workload/run manifests, exposure tracking, independently
verified complete intervals, memory and operation accounting remain required.
Controlled CPU speedup also requires the host isolation/noise receipt. This
module emits null identities and speedup, and false runtime, fresh, headline,
promotion and full-goal flags. It does not admit new natural yield, performance
or a global winner, and never reopens an old confirmation set.
