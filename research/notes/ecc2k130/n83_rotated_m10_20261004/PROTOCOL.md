# Frozen protocol: exact-public-target n83 rotated-m10 complete-point gate

## Question and claim boundary

Test one concrete high-arity relation mechanism on the exact public n83
Koblitz target already solved by the strong signed-Frobenius rho campaign.
The mechanism chooses ten finite points on

`E0: y^2 + x*y = x^3 + 1 / GF(2^83)`

whose x coordinates lie in Frobenius-rotated normal-basis subspaces of
dimensions `[9,9,9,8,8,8,8,8,8,8]`, then constrains their complete group-law
sum to the target.  The 83 selector coordinates are disjoint and exhaustive.

This is a **PDP representation and bounded natural-target solver gate**.  It
is not a relation-yield, factor-log, discrete-log, asymptotic, or faster-than-
rho result.  A SAT witness only admits a later full-pipeline experiment.

## Exact instance and independent baseline

- Cryptanalysis curve ID: `EC1N83Ckb1h876c2921cb64`.
- Crypto registry curve ID: `icv1-f2m83-tm6151469093347-debefd74`.
- Basis: the type-II optimal normal basis with
  `gamma_i = zeta^i + zeta^-i`, `zeta^167 = 1`.
- Curve order: `9671406556923184866742756 = 4*r`.
- Prime subgroup order: `r = 2417851639230796216685689`.
- Public generator and target are the coordinates frozen in `INPUT.json`.
- The independent rho replay recovered
  `467066815623456506232910`, verified `[k]G=Q`, and charged
  `201733439488 = 2^37.553659...` walk iterations across all contributing
  workers.  Its recorded work is the comparison baseline; this gate does not
  manufacture a wall-time value absent from that receipt.

The native `check` command must first verify that G and Q are on E0, `[r]G=O`,
and the independent rho scalar reproduces Q.  Any failure is a hard stop.

## Rotated domain

Let `c_0,...,c_82` be the Frobenius order of normal coordinates, beginning at
`gamma_1`, so squaring maps `c_j` to `c_(j+1 mod 83)`.  Slot `i` contains
coordinates `c_(10*j+i)` for `0 <= j < d_i`, with

`(d_0,...,d_9) = (9,9,9,8,8,8,8,8,8,8)`.

The concatenated positions must equal every integer from 0 through 82 exactly
once.  Each factor has these selector bits as x and 83 free y bits, is finite,
and must satisfy the full curve equation.  This domain has `2^83` raw x tuples;
that is only a capacity observation, not a target-support probability.

## Complete relation

Nine affine additions join the ten factors from left to right.  Every edge
contains all copy, inverse, doubling, and generic branches, canonical infinity
encoding, curve validity for both inputs and the result, and an existential
83-bit slope.  The final result is wired to the fixed finite target Q.  The
Boolean DAG has hash-consed constants, primary inputs, binary XOR, and binary
AND only.  Tseitin export uses four clauses per XOR, three per AND, units for
false, true, and the chain output, and no solver-specific XOR extension.

The expected primary-input count is

`83 + 10*83 + 8*(2*83+1) + 9*83 = 2996`.

Any accounting drift is a hard stop.  The DAG node cap is 8,000,000; the
ordinary DIMACS cap is 2,000,000,000 bytes.  Cap failures remain results.

## Frozen stages

Run stages in this order on one clean source commit:

1. **Identity replay.** Run `n83_rotated_m10 check` and archive the JSON.
2. **Planted certificate.** Deterministically select the first liftable x from
   each slot under the source-defined scan, build their target, emit the exact
   n83 DAG/CNF and a full known Boolean model, then independently replay every
   gate, every factor equation, and the ten-point sum.  No solver result is
   required for this positive control.
3. **Natural public target.** Cold-build the same circuit with Q fixed, write
   ordinary DIMACS, and invoke the committed CaDiCaL 3.0.1 binary for exactly
   120 wall-clock seconds, one process and its default single solver thread.
   Preserve CNF and DAG digests, byte/variable/clause counts, build and write
   wall time, solver stdout/stderr, exit status, model if any, and peak RSS if
   the supervisor supplies it.
4. **Independent replay.** If SAT, rebuild the DAG from the frozen source,
   replay the complete Boolean model, decode all ten factors, check each point
   on E0, and recompute their group sum as Q.  A solver SAT line without this
   replay is failure.  Preserve UNKNOWN/time-limit and UNSAT honestly.

The whole cold stage has a 900-second external wall cap and 12-GiB address/
RSS cap.  CaDiCaL's internal 120-second cap is additionally mandatory.  Do not
retune preprocessing, clauses, slot order, target, timeout, or solver after
seeing the natural-target outcome.

## Decision rule

- `SAT_VERIFIED_RELATION`: admit a separate preregistered relation-yield and
  full-DLP experiment.  Make no speed claim.
- `UNSAT`: this exact domain does not represent the frozen target; stop this
  arm.
- timeout, cap, crash, missing/full-model failure, or replay failure: record a
  bounded negative/inconclusive gate and stop this arm.

To claim a full method beats rho, a later frozen experiment must collect
enough independent relations, solve factor logarithms, recover at least one
held-out exact target scalar, and report complete cold operation counts and
resources—including base construction, every failed PDP query, solver work,
linear algebra, and verification—against a matched strong rho arm.  Online
reuse may be reported separately but cannot replace cold accounting.

## Degree 51 disclosure

Extension degree 51 is composite (`3*17`) and has proper subfields; it is not
a prime-degree transfer case.  Any n51 subfield/Weil-descent route must be a
separate experiment and cannot be used as evidence for n83.

