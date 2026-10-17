# Complete ordinary isogeny classes, 2026-10-16

This is an enumeration and correctness protocol. It makes no performance or
cryptanalytic improvement claim. `inputs/` retains the exact five user cases
(the filenames preserve the historical A–E labels); each native plan emits the
ICV1 identity of the exact input model. No new curve is substituted for a case.

## Requested scope and acceptance evidence

| Requirement | Acceptance evidence | Initial state |
| --- | --- | --- |
| Every F_p-isomorphism class in each of the five fixed-trace classes | Certified root order, exact class-number sum, distinct fixed-trace vertices, connected kernel-certified paths | Only historical case B is fully enumerated |
| Every stratum | Conductor for every vertex and exact observed/expected population equality | Pending native full-class replay |
| Every isogeny | A finite prime-degree set must be named; every rational cyclic kernel at every vertex in that set, including loops and parallel kernels | Arbitrary degrees are an infinite family; finite scope is recorded separately from vertex completeness |
| Reusable CLI | Curve JSON input; plan/build/resume/status/verify/solve; `crypto isogeny class` and `isogeny class` adapters | In implementation |
| Solve ECDLP on the graph | Correct target equation for each requested target and algorithm at every vertex | Full workload remains pending; single-target BSGS is a separate CLI capability |
| Rho/BSGS/PH, direct EC factor bases and applicable finite-field IC | Frozen matched ecbench inputs, operation counters, end-to-end timing, relation yields and verified logarithms; min/max/mean over explicitly identified populations | Historical B samples exist outside this repository; no new qualified comparison in this PR |
| Durable dataset access | Content hashes, bytes, manifests and store receipts; historical Drive links retained | Historical B Drive archive exists; new large outputs require a working durable store |

## Construction contract

Ordinary prime fields `3 < p < 2^61`, `D_K < -4`, `|D_K| <= 10^7`;
prime degrees at most 31, different from p. Unsupported cases fail explicitly.
The source order must be certified by the existing native order/twist backend,
or exact enumeration on small fields. Merely killing sample points is insufficient.
The conductor of Z[pi] is **not** asserted to be the source endomorphism conductor.
For nonzero ordinary trace and D_K<-4, j identifies the F_p-isomorphism class
within that trace class. Exceptional-j and supersingular cases need another backend.

Enumerate every rational kernel with the existing native division-polynomial
backend (rational cubic roots for degree 2). Preserve normalized quotient kernels
and an explicit F_p isomorphism to the stored target model. Check the exact
rational-map curve identity. Store each complete vertex expansion atomically in
SQLite with full synchronous durability; interruption rolls back only that expansion.
Do not load a whole-class visited set into RAM. Degree and input configuration
is frozen inside the database. Resume requires the frozen source-code hash.

An exhausted frontier certifies only the selected-degree connected component.
Complete class coverage additionally requires the exact expected class population.
Automatic degrees contain conductor primes and a conservative Minkowski generating
set for the maximal-order class group. An explicit degree set can omit required
connections and remains a component even when its frontier is empty.

Replay regenerates the complete outgoing kernel set, verifies every rational map,
every target isomorphism, root reachability and database integrity, then labels
strata by distance from the floor for each conductor prime and compares every
population with the ring-class formula. This is a separate replay pass using shared
arithmetic, **not an independently implemented cryptographic verifier**.

The database is a lossless finite-degree multigraph: kernel plus source model
reconstructs its normalized rational map; target scaling identifies the chosen
target model. Degrees, twists, loops and distinct kernels are not collapsed.

## Controls and stop conditions

Before attempting large generation, run the native unit/integration controls:
pure 2-volcano, pure 3-volcano, mixed 2/3 class, repeated/parallel kernels,
checkpoint/resume equivalence, missing-degree component, exact all-model census
at small p, corrupt-map rejection, and supplied-order/domain rejection.
Add single-target BSGS planted-log and outside-subgroup controls.

Freeze all five plans and compare class/stratum counts to the historical request.
Preflight storage records **24 bytes per vertex for a three-u64 a,b,j layout**;
SQLite indexes, edges, journals, free space and DLP evidence are additional.
This is a layout baseline, not an information-theoretic bound; compressed or
generative representations may be smaller. Do not launch a full run whose
chosen representation exceeds available storage.
Retain bounded prefixes explicitly as partial; they do not discharge full generation.
Failures, unsupported inputs and resource limits remain evidence, never zero costs.

No timing admission, speed comparison, IC relation-yield experiment, or new
mathematical theorem is attempted in this enumeration PR. Those require their
existing ecbench/isolation/proof protocols and separate execution evidence.
The class-number and volcano results used here are established results, also
documented in `src/cryptanalysis/isogeny_degree_search/class_number.rs` and
`docs/endomorphism-rules.md`; this PR checks their implementation on the controls.

## Environment and workflow gaps

Base: `082d5507859860f8a1ac1f20318c49e964d2af02`. Workspace: 4 CPU quota,
16 GiB memory, 32 GiB filesystem, separate 8.8 GiB temporary mount.
The mandated Conductor invocation returned exit 127 (`conductor: command not
found`); `.conductor/WORKFLOW.md` is absent from the base tree. No reservation or
task is invented. Scope: the new native package, two CLI adapters and lockfiles,
this study, documentation, and the named problem record.
