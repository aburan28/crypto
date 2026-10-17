# ECDLP audit control plane

This directory holds the machine-readable planning layer for the proposed
[comprehensive elliptic-curve private-key recovery audit](../../research/notes/ecdlp-general/RESEARCH_COMPREHENSIVE_ECDLP_AUDIT_PLAN.md).
The precise scope is authorized research into unauthorized recovery routes;
most protocol and implementation incidents do not solve the original ECDLP and
must remain separately classified.

## Current artifacts

| Artifact | Role | Authority |
|:--|:--|:--|
| [`ECDLP_AUDIT_PLAN.mock.json`](ECDLP_AUDIT_PLAN.mock.json) | Exact eight-workstream plan, dependencies, acceptance gates, taxonomy skeleton, and open obligations | **non-authoritative mock**; cannot promote claims |
| [`ecdlp-audit-plan.schema.json`](ecdlp-audit-plan.schema.json) | Draft 2020-12 structural contract | validates the mock only |
| [`ecdlp-audit`](../../src/bin/ecdlp_audit.rs) | `verify-plan` checks planning identity, references, hashes, status, gates, confinement, and claim boundaries; `verify-fixed-f8-order` independently replays one exact point-count/zeta certificate | the plan subcommand cannot promote claims; the order subcommand certifies only the stated fixed-fixture cardinality |
| [`program-dependencies.svg`](program-dependencies.svg) | Dependency and promotion-boundary diagram | visual twin of the plan |
| [`execution-boundary.svg`](execution-boundary.svg) | Current cover/Jacobian evidence boundary | checkpoint visual, not a completeness claim |
| [`execution checkpoint`](../../research/notes/ecdlp-general/RESEARCH_ECDLP_AUDIT_EXECUTION_CHECKPOINT_20261014.md) | Locally reproduced first implementation and audit checkpoint | all 20 scientific gates remain open |
| [`fixed-F8 order certificate`](../curves/fixed-f8-genus3-order-certificate.json) | Implementation-separated extension-count and zeta replay for the exact genus-three fixture | certifies `#Jac(H)(GF(8))=486` only; no arithmetic or norm-map oracle |

Draft workstream contracts are separated by responsibility:

| Directory | Workstream slice | Boundary |
|:--|:--|:--|
| [`ledger/`](ledger/) | WS3 identity audit | proposed codec and open obligations only |
| [`evidence/`](evidence/) | WS4 E0–E7 profiles and transition fixtures | no authoritative transition engine |
| [`routes/`](routes/) | WS5 route-registry and mock matrix contracts | candidate inventory; universal denominator unfrozen |
| [`corpus/`](corpus/) | WS6 protocol, schema, pointer seed, and negative fixtures | no systematic search executed |
| [`accounting/`](accounting/) | WS7 cross-route receipts and structural controls | invented fixture values only |
| [`governance/`](governance/) | WS8 policy, roles, and workflow drills | proposed; no protected approval state |

Run the planning checks from a full matching repository checkout:

```bash
cargo run --release --bin ecdlp-audit -- verify-plan \
  --input docs/ecdlp-audit/ECDLP_AUDIT_PLAN.mock.json \
  --repo-root .
```

A locally built executable can invoke the same command directly as
`target/release/ecdlp-audit verify-plan ...`.

Replay the separately scoped fixed-fixture certificate with:

```bash
cargo run --release --bin ecdlp-audit -- verify-fixed-f8-order \
  --input docs/curves/fixed-f8-genus3-order-certificate.json \
  --repo-root .
```

This command does not import production curve, Jacobian, or transfer
arithmetic. It certifies only `#Jac(H)(GF(8))=486` for the exact recorded
polynomial-basis fixture; it does not certify a cover map, group presentation,
norm, catalog row, solver, or security conclusion.

The release workflow is configured to build `ecdlp-audit`, start it with
`--help` on each native Rust target, and stage this directory in the Rust
archives as non-authoritative examples. That configuration is not evidence of
a successful CI run or published package; check the workflow run and release
assets for the relevant commit.

The archive stages the order certificate and supports
`verify-fixed-f8-order`. The staged audit directory is not a self-contained
`verify-plan` checkout. The mock
binds each historical source baseline to an exact Git
`(commit, repository path, raw-blob SHA-256)` tuple. A pinned source is replayed
from that immutable Git object; it deliberately does not assert that a mutable
worktree copy still equals the historical baseline. Every deliverable with an
artifact separately binds its current repository-relative path and exact raw
file bytes by SHA-256. Use the archive copy for inspection only. Run semantic
verification from the matching full repository with `--repo-root`, and do not
interpret missing archive paths as a scientific or completeness result.

The mock intentionally has:

- `authoritative=false` and `promotion_allowed=false`;
- no authoritative subjects or assessments;
- an unfrozen universal denominator;
- all 20 execution gates open and explicit; and
- five top-level taxonomy skeletons plus a 43-family candidate registry; the
  frozen, cited route population remains a WS5/WS6 deliverable.

It must not be renamed to `master-ledger.json` or treated as that ledger. WS3
requires reviewed identity rules, a frozen population manifest, an additive
generator, an independent semantic verifier, migration receipts, and
governance approval before an authoritative ledger exists.

## Status semantics

Claim status is one of `supported`, `refuted`, `unknown`, or
`not_applicable`. Evidence maturity E0–E7 is separate. An explicit `unknown`
counts as a classified cell but not as scientific support, safety, or
completion. A failed bounded search supports only its exact bounded result.
Workstream requirement status is derived from its deliverable and required-gate
states, and every required gate must be referenced by at least one deliverable.
An artifact digest establishes byte identity only; it is not scientific or
semantic validation of the artifact's claims.

Executable work defaults to passive public evidence and public
synthetic/known-answer fixtures. No live-target execution, secret material, or
production exploitation is authorized by these files.
