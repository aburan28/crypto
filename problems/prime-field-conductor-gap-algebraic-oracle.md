# Prime-field conductor-gap algebraic decomposition oracle

- **Status:** open
- **Kind:** implementation gap
- **Affected work:** [large-gap IC measurements](../research/ecbench_conductor_gap_20261012/README.md), [complete instruction census](../research/ecbench_conductor_gap_ir_20261014/PROTOCOL.md)
- **Next action owner:** unassigned; native IC oracle workstream
- **Last updated:** 2026-10-14

## Problem

The native ecbench pipeline cannot yet compare algebraic F4/SAT IC on the
archived prime-field conductor-gap pairs because its algebraic decomposition
adapters accept binary curves and Boolean systems, while its prime adapter
does not implement that oracle.

The affected exact models are `icv1-fp29-t47-14392892` and
`icv1-fp29-t47-1f934956` at ell = 8209, and
`icv1-fp37-t277-91fd6e29` and `icv1-fp37-t277-090dc414` at ell = 131101.
Exact fields, subgroups, generators, targets and representation preimages
are retained in the linked studies. The implemented MITM comparison cannot
discharge an algebraic-method comparison or an intrinsic-hardness question.

## Evidence and uncertainty

- **Observed:** on base revision `275beb93f3bf33ade6c19f19e8b336027a00b297`,
  `src/cryptanalysis/ecbench/methods.rs::solve_ic_prime` selects `subtract`,
  `mitm` or `mitm-fiber`; an algebraic name is rejected. The binary adapter
  selects `descent-algebraic` and its solver plug-ins. The solver interface in
  `src/cryptanalysis/ic_framework/solvers.rs` accepts a Boolean system.
  This is source inspection, not a failed prime-field F4/SAT experiment.
- **Derived or proved:** no theorem of mathematical impossibility is established.
- **Suspected explanation:** no conjecture is needed to identify this adapter
  gap; a prime-field encoding and native solver connection remain missing.
- **Missing capability:** an executable prime-field decomposition oracle with
  explicit factor-base membership, witness verification and complete counters.
- **What this evidence does not establish:** absence of an algebraic advantage,
  security of either pair, or inability to implement a suitable solver.

## What would unblock it

- [ ] Add a native prime-field algebraic oracle to the same ecbench pipeline;
  expose every parameter that changes the decomposition search.
- [ ] Verify decompositions against direct group arithmetic and an exhaustive
  small-field/MITM control, retaining exceptional and no-solution cases.
- [ ] Record native encoding, solving and witness costs, including failures;
  do not convert partial solver counters into complete-cost claims.
- [ ] Freeze and run matched one-target F4/SAT and rho workloads on the exact
  ell = 8209 and 131101 pairs, with audit/replay evidence and resource limits.

These checks would establish the implemented comparison in its measured
scope. Representation controls, additional independent isogeny classes,
explicit target transport and cryptographic-scale conclusions remain separate.

## Next discriminating step

- **Action:** implement and validate one bounded two-summand prime-field
  decomposition encoding before attempting the large-gap benchmark.
- **Inputs and controls:** a fully enumerable small prime field and a fixed
  factor base; compare every returned witness to the native MITM/direct control.
- **Execution limits:** freeze field, base, queries, solver parameters and
  watchdog in a versioned protocol before executing; no campaign spending cap
  is introduced by this record.
- **Outcomes:** correct witness-complete bounded results permit a measured
  large-gap follow-up; budget exhaustion retains an implementation/resource
  observation; wrong witnesses block any inference until corrected.
- **Recheck trigger:** a native prime-field oracle passes the controls above.

## Resolution and history

Unresolved. The instruction census measures the existing MITM implementation;
it supplies no F4/SAT result and does not resolve this problem.
