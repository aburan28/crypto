# Frozen follow-up plan: native-XOR n83 rotated-m10 gate

## Motivation and separation from prior outcomes

PR #1343's complete-per-edge ordinary-CNF circuit and PR #1346's smaller
inductive-validity ordinary-CNF circuit both returned `UNKNOWN` at their
separately frozen 120-second CaDiCaL caps.  The latter removed about 30% of the
variables and clauses and 42% of peak RSS, but still supplied no public
relation.

This follow-up changes the solving representation.  It keeps the frozen
inductive DAG and emits each Boolean XOR gate as one CryptoMiniSat extended-
DIMACS parity row instead of four Tseitin clauses.  AND gates retain their
three exact CNF clauses.  The target, factor slots, complete group-law
branches, validity policy, and primary inputs do not change.

The experiment can establish an exact export reduction or one verified PDP
relation.  It cannot establish a relation yield, full ECDLP, or a crossover
against rho.

## Development observations disclosed before the public run

During implementation, the deterministic planted instance was emitted and
validated against its known complete model.  It retained 693,872 variables,
replaced 1,980,768 XOR Tseitin clauses by 495,192 native parity rows, and used
587,049 remaining CNF clauses.  The resulting file was 22,859,716 bytes.

One ten-second CryptoMiniSat development smoke on that planted instance ended
`INDETERMINATE`.  Its log reported that the connected 465,344-row residual
parity matrix exceeded the default 100,000-row/column limits and that zero
Gauss-Jordan matrices were enabled.  This is not hidden or retuned.  The
frozen gate therefore tests compact native parity storage, parity propagation,
and CryptoMiniSat preprocessing/search under its disclosed defaults; it does
not claim that a giant dense Gaussian matrix is active.  No public-target XOR
instance was emitted or solved during development.

## Correctness gates before the public attempt

1. Run all seven native module tests, including exhaustive type-II-ONB n=3
   field laws, every valid n=3 point-addition case without intermediate
   validity checks, exact n83 rho-scalar replay, and both planted n83 variants.
2. On a toy DAG containing an XOR and an AND, require the streaming extended-
   DIMACS verifier to accept the true full assignment and reject independently
   corrupted XOR and AND outputs.
3. Run the frozen CryptoMiniSat binary on committed three-variable parity
   controls.  The exact `z=a xor b` positive control must exit 10/SAT and its
   one-bit-corrupted twin must exit 20/UNSAT.
4. Emit the deterministic planted n83 XOR certificate.  Stream-check every
   parity row and CNF clause against the complete model, then replay every DAG
   gate, all ten curve points, and the numerical group sum in a fresh process.
5. Require the public inductive DAG prefix to remain
   `e8715d805517fc35cb70f6b96c5d387a28579c0800e259f69e642c91e5636f88`.

Any failed control or digest drift is a hard stop.

## Frozen public instance and solver

Use the exact `EC1N83Ckb1h876c2921cb64` public Q, type-II optimal normal
basis, ten slots `[9,9,9,8,8,8,8,8,8,8]`, inductive-validity circuit, and
8,000,000-node cap from PR #1346.  Use the official Linux amd64 release asset
for CryptoMiniSat 5.16.0, one thread, random seed 1, one requested solution,
default preprocessing/Gaussian settings, and `--maxtime 120`.

Run exactly one public-Q solver attempt.  Preserve the exporter receipt,
extended DIMACS digest, full stdout/stderr, exit code, elapsed time, solver
statistics, and peak memory.  If SAT, extract the complete model, stream-check
it against every emitted constraint, rebuild the DAG in a fresh process, and
independently verify all factors and their group sum.

The outer cold-process cap is 900 seconds, the file cap is 2,000,000,000
bytes, and the address-space/RSS cap is 12 GiB.  Do not change solver options,
representation, target, or timeout after observing the public result.

## Decision rule

- `SAT_VERIFIED_RELATION` requires both extended-DIMACS validation and full
  numerical replay.  It permits a separately frozen disjoint-target yield
  panel; it is not a DLP or rho crossover.
- UNSAT rejects this exact tuple domain for Q, subject to the correctness
  controls; it is not a general no-go result for high-arity descent.
- INDETERMINATE, timeout, cap, crash, or replay failure is a bounded negative
  for this native-XOR/CryptoMiniSat route.

No faster-than-rho claim is allowed without relation collection, rank, factor
logarithms, held-out target recovery, and fully charged cold cost below the
matched strong rho baseline of 201,733,439,488 walk iterations.

Degree 51 is composite (`51=3*17`) and has proper subfields.  It remains a
separate control case and cannot establish the claimed prime-extension n83
result.
