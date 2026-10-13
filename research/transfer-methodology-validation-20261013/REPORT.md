---
title: "ECDLP transfer methodology: known-case and net-cost validation"
date: "2026-10-13"
author: "aburan28/crypto research record"
geometry: margin=0.6in
fontsize: 10pt
colorlinks: true
linkcolor: blue
urlcolor: blue
---

# Outcome

The fail-closed transfer methodology and native validator classify six known
public cases consistently with their frozen source evidence. All six
assessments pass the version-three machine contract, and ten semantic controls
pass, including rejection of unnamed solvers, stage-only promotions, unpriced
work, inconsistent phase totals, unsupported measurements, and incorrect net
ratios. Two fixtures support a valid subgroup transfer, two remain unknown at
the transfer gate, and two are refuted. None supports a modeled or measured
computational advantage.

The extension now defines what “a weaker destination” means operationally. A
candidate must name the destination solver, show that it reaches the
transported order-`r` component, and compare three arms in one unit: the
strongest direct source solver, a generic solver on the destination image, and
the specialized destination pipeline with all nine transport/solve/recovery
phases charged. A destination-stage gain is recorded but cannot promote the
net-transfer verdict.

The replay also found and repaired a pre-existing build defect in
`tools/ecdlp-surfaces`: the checked-in binary module was not exported from the
crate root, so the package's own validation script could not compile the n=83
checker. After adding the missing module export, the existing package reported
22 tests passed, zero failed, Clippy clean, and byte-for-byte replay of its
committed examples.

![The diagnostic first establishes the transfer and named solver, then admits
only a three-arm, cold end-to-end net comparison.](validation-flow.pdf){width=100%}

# Requested scope and evidence

| Requirement | Status | Evidence |
|---|---|---|
| Put the methodology in `aburan28/crypto` | Verified complete | `docs/transfer-assessment/METHODOLOGY.md` |
| Make it machine-checkable | Verified complete | `tools/transfer-assessment`, v3 schema, template, ten semantic controls |
| Run it on known transfers | Verified complete for six bounded public cases | `assessments/` and `validation-summary.json` |
| Include covers, surfaces, higher genus, and lattices | Verified complete as diagnostic coverage | quartic image, product surfaces, genus-3 cover, torsion-free lattice control |
| Define when to use the repository skills | Verified complete | `docs/transfer-assessment/MEASUREMENT.md` and the `ecdlp-surfaces` skill handoff |
| Measure whether a destination solver is actually weaker | Protocol and machine gates complete; no new solver run performed | named-solver record, three arms, nine phases, net ratio, replay/interval/scaling gates |
| Establish a new faster ECDLP reduction | Not established and not implied | all advantage verdicts are `unknown` or `not_applicable` |
| Validate arbitrary covers/lattices/surfaces | Open beyond the six cases | each new candidate still needs its own identities, map, subgroup, solver, and cost evidence |

The validation profile is an assessment consistency check. It does not make an
official state transition, independently prove every cited artifact, or turn
template completion into evidence.

# Pinned basis

Repository revision:
`3d4a891bd07000cf63935f5c876a22ff5e193315`.

The cases cite these exact paths at that revision:

- [`tools/ecdlp-surfaces/examples/transfer_certificate.json`](https://github.com/aburan28/crypto/blob/3d4a891bd07000cf63935f5c876a22ff5e193315/tools/ecdlp-surfaces/examples/transfer_certificate.json);
- [`tools/ecdlp-surfaces/examples/frobenius_result.json`](https://github.com/aburan28/crypto/blob/3d4a891bd07000cf63935f5c876a22ff5e193315/tools/ecdlp-surfaces/examples/frobenius_result.json);
- [`tools/ecdlp-surfaces/examples/binary_n83_a0.json`](https://github.com/aburan28/crypto/blob/3d4a891bd07000cf63935f5c876a22ff5e193315/tools/ecdlp-surfaces/examples/binary_n83_a0.json);
- [`tools/ecdlp-surfaces/README.md`](https://github.com/aburan28/crypto/blob/3d4a891bd07000cf63935f5c876a22ff5e193315/tools/ecdlp-surfaces/README.md); and
- [`research/hyperelliptic_cover_infrastructure_20261006/README.md`](https://github.com/aburan28/crypto/blob/3d4a891bd07000cf63935f5c876a22ff5e193315/research/hyperelliptic_cover_infrastructure_20261006/README.md).

The lattice case is a known-false type control from the methodology: an
ordinary additive lattice `Z^m` is torsion-free, so a homomorphism from the
finite cyclic group `Z/rZ` has zero image. It does not reject finite lattice
quotients, tori, endomorphism-kernel lattices, or relation lattices; those must
be typed and audited separately.

# Exact classifications

| Case | Strongest established stage | Valid transfer | Algorithmically relevant | Modeled / measured advantage |
|---|---|:---:|:---:|:---:|
| Small-prime quartic image | explicit coordinate isomorphism, subgroup and recovery checks | supported | supported for generic group algorithms | unknown / unknown |
| n=83 inclusion `E_0 -> E_0 x E_1` | explicit factor inclusion and projection recovery | supported | supported for generic group algorithms | unknown / unknown |
| Compatible degree-four Frobenius polynomial | necessary isogeny-class filter only | unknown | unknown | unknown / unknown |
| n=83 genus-3 cover | curve construction and `C -> E` map; exact n=83 Jacobian pullback unexecuted | unknown | unknown | unknown / unknown |
| n=83 proposed `E_0 -> E_1 x E_1` | order and Frobenius obstruction | refuted | refuted | not applicable |
| proposed `Z/rZ -> Z^m` | torsion obstruction | refuted | refuted | not applicable |

This is the intended discrimination:

1. A faithful representation may be a valid transfer without an advantage.
2. Frobenius compatibility is necessary evidence, not a homomorphism.
3. A curve map `C -> E` is not the required executable subgroup map
   `E -> J(C)`; the n=83 pullback remains unknown until replayed.
4. Group-order, Frobenius, and torsion contradictions reject impossible paths
   before solver work.
5. Generic rho is applicable to both valid image subgroups, but that says only
   that the DLP survives transport; it does not make either destination weaker.

# Net-advantage decision rule

The measurement record contains three complete arms:
`source_reference`, `destination_generic_control`, and
`destination_candidate`. Its primary scope is one cold target. Candidate cost
is the sum of route search, construction, forward map, solver setup,
relation/walk, linear algebra, individual-log/descent, source recovery, and
verification. The only net ratio eligible for promotion is
`source_reference.total / destination_candidate.total`.

The destination generic control answers whether any observed change is merely
representation overhead. The specialized arm answers whether index calculus,
Weil descent, or another named method exploits the destination. Both remain
secondary to the net ratio. Amortization, online-only intervals, and individual
phase ratios are recorded separately and cannot promote the cold verdict.

The skill handoff is deliberate: `transfer` and `ecdlp-surfaces` find and
certify routes; `curve-algebra` establishes relevant endomorphism/Frobenius
actions; `ecbench-measure` runs the matched cross-method comparison; and
`ecbench-bounds` tests an admitted measurement against the established
frontier. `report-evidence` and `research-visuals` govern every run and result.

# Native validator and controls

`tools/transfer-assessment` parses every assessment, binds each supported or
refuted state and every cost arm to a declared source id, and enforces
cumulative verdict gates:

- `valid_transfer=supported` requires supported map correctness, subgroup
  preservation, and scalar recovery;
- `algorithmically_relevant=supported` also requires a named solver with
  supported applicability, image access, search/relation compatibility, and
  individual-log/recovery gates;
- modeled advantage additionally requires the cold one-target scope, all three
  arms, a common unit, nine-phase totals, no unpriced work, a matching ratio
  greater than one, correctness, resources, and bound evidence;
- measured advantage additionally requires measured arms, five verified pairs,
  four sizes, a replay receipt, and a 95% interval entirely above one;
- novelty cannot be promoted without the earlier cumulative verdicts; and
- the diagnostic profile refuses an official status change.

Ten Rust tests cover a valid unknown record, a supported transfer, missing
transfer gates, unbound evidence, missing solver details, a complete modeled
control, stage-only and unpriced rejection, inconsistent phase/ratio
arithmetic, a complete measured control, and missing measured-evidence gates.
Result: 10 passed, 0 failed.

The six study assessments then produced
`transfer-assessment-validation-summary-v1` with
`assessments_checked=6`, `all_valid=true`, and no validation errors.

# Reproduction

Environment:

```text
Linux 6.18.44 x86_64 GNU/Linux
rustc 1.99.0 (b940084d7 2026-09-28)
cargo 1.99.0 (5f94df478 2026-08-27)
Graphviz 2.42.4
Pandoc 3.1.11.1
```

Commands:

```sh
cargo test --release --locked \
  --manifest-path tools/transfer-assessment/Cargo.toml

cargo clippy --release --locked \
  --manifest-path tools/transfer-assessment/Cargo.toml -- -D warnings

npx --yes ajv-cli@5.0.0 validate --spec=draft2020 \
  -s docs/transfer-assessment/assessment.schema.json \
  -d docs/transfer-assessment/assessment-template.json \
  -d 'research/transfer-methodology-validation-20261013/assessments/*.json'

cargo run --release --locked \
  --manifest-path tools/transfer-assessment/Cargo.toml -- \
  --output research/transfer-methodology-validation-20261013/validation-summary.json \
  research/transfer-methodology-validation-20261013/assessments/*.json

tools/ecdlp-surfaces/check.sh

cargo test --release --lib --quiet
```

The first surface-package replay at the pinned base failed during compilation:

```text
error[E0432]: unresolved import `ecdlp_surfaces::binary`
```

After `pub mod binary;` was restored in `tools/ecdlp-surfaces/src/lib.rs`, the
same command passed 9 binary tests and 13 math tests, Clippy, and committed
example regeneration. The repair changes module reachability only; the fixture
mathematics and frozen JSON were not edited.

The final repository-wide run on the solver-aware extension, with all
compile-time fixtures materialized and the repository's local-socket tests
permitted, reported 4,058 passed, zero failed, and 114 ignored in 145.78
seconds. The first sparse/sandboxed attempt from the original methodology run
reported 4,043 passed and 15 failed: six failures were missing sparse-checkout
fixtures, and nine were local-socket `Operation not permitted` errors. After
materializing the exact two remaining fixture directories, the six affected
tests passed; all 37 `pollard_collab` tests passed with local sockets enabled;
then the clean full run above passed. These were access/setup failures, not
mathematical evidence and not changes to the tested implementation.

# Graph scope

The boundary ledger, index-calculus scoreboard, progress timeline, curve
registry, cover catalog, and existing surface diagrams were checked for
applicability. None changes: this PR adds no new curve identity, map edge,
admitted performance measurement, solver advantage, or official claim. The
study-local diagram records only diagnostic classifications.

# Limits and next discriminating check

This round validates methodology behavior on bounded known cases, not general
completeness. It contains no performance benchmark, so no `benchmarks/` entry
or PNG performance chart is required. No new or substantially revised
mathematical theorem is asserted, so there is no new Lean theorem to map.

The next higher-genus workstream has three ordered gates:

1. Execute the native n=83 `OrdinaryBinaryCover` pullbacks for the exact pinned
   `P` and `Q`; verify image order `r`, degree-three pull/push composition,
   exceptional fibers, and source scalar recovery.
2. Select and implement a destination solver adapter that can reach that exact
   Jacobian image. Before n=83 performance work, validate factor-base,
   relation, linear-algebra, individual-log, norm/projection, and false-relation
   controls on bounded public examples.
3. Run the three-arm cold comparison over at least four faithful toy sizes,
   preserve each phase and failure, then apply the relevant n=83 family gate.
   A stage crossover is not a method crossover; an n=83 result is not an
   automatic n=131 result.

Completing gate one would support the transfer, not an advantage. Completing
gate two would support solver applicability, not net speed. Only gate three
can populate a measured advantage record.

Conductor scope coordination could not be performed: `conductor` was absent
from this execution environment and the required command exited 127. Work was
kept to `docs/transfer-assessment/`, `tools/transfer-assessment/`, this study
directory, the `ecdlp-surfaces` skill handoff, and the one-line surface module
repair disclosed above.
