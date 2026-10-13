---
title: "ECDLP transfer methodology: known-case validation"
date: "2026-10-13"
author: "aburan28/crypto research record"
geometry: margin=0.75in
fontsize: 10pt
colorlinks: true
linkcolor: blue
urlcolor: blue
---

# Outcome

The fail-closed transfer methodology and native validator classify six known
public cases consistently with their frozen source evidence. All six
assessments pass the machine contract, and four deliberately invalid semantic
controls pass by being rejected. Two fixtures support a valid subgroup
transfer, two remain unknown at the transfer gate, and two are refuted. None
supports a modeled or measured computational advantage.

The replay also found and repaired a pre-existing build defect in
`tools/ecdlp-surfaces`: the checked-in binary module was not exported from the
crate root, so the package's own validation script could not compile the n=83
checker. After adding the missing module export, the existing package reported
22 tests passed, zero failed, Clippy clean, and byte-for-byte replay of its
committed examples.

![The diagnostic accepts complete transfers, preserves unknowns after merely
necessary filters or unexecuted subgroup maps, and rejects algebraically
impossible destinations.](validation-flow.pdf){width=100%}

# Requested scope and evidence

| Requirement | Status | Evidence |
|---|---|---|
| Put the methodology in `aburan28/crypto` | Verified complete | `docs/transfer-assessment/METHODOLOGY.md` |
| Make it machine-checkable | Verified complete | `tools/transfer-assessment`, schema, four invalid controls |
| Run it on known transfers | Verified complete for six bounded public cases | `assessments/` and `validation-summary.json` |
| Include covers, surfaces, higher genus, and lattices | Verified complete as diagnostic coverage | quartic image, product surfaces, genus-3 cover, torsion-free lattice control |
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

# Native validator and controls

`tools/transfer-assessment` parses every assessment, binds each supported or
refuted state to a declared source id, and enforces cumulative verdict gates:

- `valid_transfer=supported` requires supported map correctness, subgroup
  preservation, and scalar recovery;
- `algorithmically_relevant=supported` also requires a named applicable solver;
- modeled advantage additionally requires complete accounting;
- measured advantage additionally requires scaling support and measured data;
- novelty cannot be promoted without the earlier cumulative verdicts; and
- the diagnostic profile refuses an official status change.

Four Rust tests exercised one valid unknown record, one fully supported
transfer, promotion without the three transfer gates, and unbound evidence plus
an unsupported measured-advantage promotion. Result: 4 passed, 0 failed.

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

The final repository-wide run, with all compile-time fixtures materialized and
the repository's local-socket tests permitted, reported 4,058 passed, zero
failed, and 114 ignored in 217.16 seconds. The first sparse/sandboxed attempt
reported 4,043 passed and 15 failed: six failures were missing sparse-checkout
fixtures, and nine were local-socket `Operation not permitted` errors. After
materializing the exact two remaining fixture directories, the six affected
tests passed; all 37 `pollard_collab` tests passed with local sockets enabled;
then the clean full run above passed. These were access/setup failures, not
mathematical evidence and not changes to the tested implementation.

# Limits and next discriminating check

This round validates methodology behavior on bounded known cases, not general
completeness. It contains no performance benchmark, so no `benchmarks/` entry
or PNG performance chart is required. No new or substantially revised
mathematical theorem is asserted, so there is no new Lean theorem to map.

The next check for the higher-genus route is concrete: execute the native n=83
`OrdinaryBinaryCover` pullbacks for the exact pinned `P` and `Q`, verify image
order `r`, the degree-three pull/push composition including exceptional fibers,
and source scalar recovery. A passing transfer still would not establish an
advantage until a destination solver and complete matched cost comparison are
attached.

Conductor scope coordination could not be performed: `conductor` was absent
from this execution environment and the required command exited 1 with no
output. Work was kept to `docs/transfer-assessment/`,
`tools/transfer-assessment/`, this study directory, and the one-line surface
module repair disclosed above.

# Graph scope

The boundary ledger, index-calculus scoreboard, progress timeline, curve
registry, cover catalog, and existing surface diagrams were checked for
applicability. None changes: this PR adds no new curve identity, map edge,
admitted performance measurement, solver advantage, or official claim. The
study-local diagram records only diagnostic classifications.
