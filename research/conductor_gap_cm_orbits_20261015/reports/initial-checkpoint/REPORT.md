# Native conductor-gap controls and larger explicit floors

15 October 2026. This report is cumulative; sealed run directories and the
historical evidence archive are immutable. The 100-pair solve campaign is
running. Proposed or pending experiments are not counted as completed solves.

## Completed observations

All 8,000 construction-only trials finished: 100 gaps, 20 shuffled repetitions,
three scientific arms and a duplicate scalar A/A arm. The backend compared
entire crater point/column/coefficient vectors, in addition to recording
supplementary fingerprints. All 2,000 paired map/scalar comparisons agreed.
Each trial constructed 8,192 points from four seed columns and 2,048 levels.

| Arm | Runs | Mean thread CPU per construction |
|---|---:|---:|
| Crater, degree-5 CM map | 2,000 | 0.005104398 s |
| Same crater, scalar action | 2,000 | 0.040143001 s |
| Floor, scalar action | 2,000 | 0.040216660 s |
| Crater scalar A/A | 2,000 | 0.040501120 s |

These are descriptive L0 observations under concurrent work, not an isolated
wall speedup. Native operation counters retain unpriced modular exponentiation,
hash, allocation and elimination work. The construction-stage difference does
not establish a difference in total ECDLP hardness.

Raw data: [orbits.jsonl](results/orbits.jsonl), native-derived
[summary](results/orbit-summary.json), [conditions](results/CONDITIONS.md),
[pre-run hashes](results/orbit-pre-run.sha256), and [native source](native/folded.cpp).

The new square-root Vélu codomain implementation matched both coefficients of
all 100 previously constructed direct codomains. The reference check retained
every row in [sqrt-velu-reference.jsonl](prepared/sqrt-velu-reference.jsonl).
A separate independent PARI run proved primality, checked both endpoint
orders using Hasse-unique prime-order witnesses, and checked the prime-kernel
witness for all 100 pairs: [certificates](results/pair-certificates.jsonl).

The independent checker also certified 15 newly constructed floor models:
five D = −19 fields through 91 bits and ten D = −43 controls. These include
prime conductor gap 17,592,186,219,863 over
p = 1547425079971670554692473983. Exact coefficients, subgroup generators,
orders and discriminants are in [floor-certificates-v2.jsonl](results/floor-certificates-v2.jsonl).
Classification uses the ordinary CM order lattice and the class-number-one
maximal-order j invariant; this mathematical bridge is not claimed as a Lean
formalization. No explicit full isogeny point map or dual is certified here.

The 95/96-bit constructor rows remain pending. The first process was stopped
after the 91-bit result when combined memory reached 15,982,964,736 bytes of
a 17,179,869,184-byte cgroup limit. [Stop receipt](results/larger-construction-stop.json).
Its first output schema marks an absent direct reference as true; those values
are not validation receipts. The corrected version emits null and an explicit
`reference_supplied` field. All new endpoint claims use independent PARI checks.

## Native end-to-end harness

`ic.cm_folded_native_v1` runs the existing compiled native relation backend
through ecbench, with an executable-content hash required before launch.
Rust independently checks every factor-base scalar, every retained relation,
matrix rank, and recovered public point. No validation scalar is an input.
All reusable setup is included in supplementary cold cost. The five exclusive
online phases include independent replay and sum to the complete online clock.
Unpriced work stays explicit; zero charged GAE is not a complete zero cost.

The first small control failed all six jobs before solving because measured
children lost the supplied public fixture. That session is retained in
`sessions/control`. `CurveSpec::for_child` fixes measured execution, profiling
and replay. The next control verified all six jobs, including exact map/scalar
query-count agreement in both rounds: `sessions/control-fixture-fix`.

The primary campaign requests 300 IC attempts and 200 independent same-point
rho solves across the original 200 public regression fixtures. Its eight
shards run four co-scheduled workers at L0, with no shared target tables.
Source, executable and spec hashes were saved before launch in
[campaign-pre-run.sha256](results/campaign-pre-run.sha256). Do not interpret
a running shard, a censored target or an incomplete cost vector as a win.

## Requirement ledger

| Requested work | Current evidence/status |
|---|---|
| A: 100-pair repeated three-arm orbit comparison | Complete at L0; exact records agree; duplicate scalar controls included. Isolated cold/warm timing remains open. |
| B: 200-target folded campaign plus controls | Running, with 300 IC and 200 rho attempts requested. No final solve-rate or speed ratio yet. |
| C: fresh targets, repeated query seeds, coordinate/floor/family replication | Not yet executed. Ten second-family endpoints are constructed and independently certified. |
| D: separate structural and generic effects | Same-crater orbit equality established. Full solve comparison and replication pending. No prime-field F4/SAT complexity claim. |
| E: native fast codomain construction and larger gaps | All 100 direct references match; 15 new floor models independently certified through 91 bits. 95/96-bit rows and the original non-rational-kernel 96-bit floor remain open. |
| PR, reports, required checks and canonical indexes | In progress. Passing controls and opening a PR do not complete the research ledger. |

## Proof and provenance

The three evaluated polynomial identity certificates are
`IDC1hd41ecc586e59b223` (CM norm numerator),
`IDC1hc6bfc34b0ebbfa79` (paired second moment), and
`IDC1hd387bbbb1058c248` (paired third moment). Each has 64 issuer points and
64 additional domain-separated checker points. See [certificates](certificates/).

Lean 4.34.1 checks conditional orbit-equality and norm-bound lemmas with no
axioms. [Proof map](lean/README.md) lists the unformalized CM bridges explicitly.

The original six-page PDF and source evidence remain in
[historical-evidence-summary.zip](inputs/historical-evidence-summary.zip).
The initial source revision was e15e29d. Conductor coordination was attempted
but returned `not logged in`; no shared task lease is claimed. Work is on an
isolated feature branch. No remote scientific registry or paper is modified.
