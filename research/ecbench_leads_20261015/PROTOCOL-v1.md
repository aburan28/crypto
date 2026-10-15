# Follow-up to the ten EC index-calculus questions

Frozen before measurement, 2026-10-15. Parent evidence is
[the ten-question campaign](../ecbench_open_questions_20261014/REPORT.md).
The request is to pursue its leads using ecbench; the original ten mathematical
resolution criteria remain in force. This iteration addresses implementation
and discriminating experiments, not a promise to solve open mathematics.

## Leads and frozen tests

1. **Correct algebraic shape (Q2/Q6).** Add a versioned PDP3 adapter reporting
   the actual chained system: m*ell+(m-2)*n Boolean variables and (m-1)*n
   coordinate equations. Preserve the old adapter and all old sessions. Compare
   old/v2 inherited-F4 and F6-IC on the same public targets, degrees
   13/19/23/31, ell=3, m=3, degree cap4, node budget4096, eight trials,
   two targets and one measured round. Independently construct the systems
   and check metadata against their stored variables/equations. Identical
   recovery outcomes and group counts are required; metadata may differ.

2. **Width-generic ordinary IC (Q3/Q7/Q10).** Implement a new method identity,
   `ic.orbit_rank_v1`, with unknown factor-base logarithms, exact point witnesses,
   modular rank collection and target pinning. Its factor base scans raw x from1,
   projects by the cofactor and retains distinct signed Frobenius orbits. Compare
   exact full-point and signed-Frobenius quotient pair tables on an identical
   base. Each quotient witness retains seed columns, relative Frobenius phase
   and sign, then lifts to and checks the full point equation. This is an
   enumerated table certificate, not a low-degree algebraic correspondence.
   Known-log seed bases are excluded. All point additions/doublings and scalar
   multiplications use counted ledgers; native work remains explicit and unpriced.
   Stop at the explicit pair-visit and query caps; retain counts on exhaustion.

   Narrow controls: E0 degrees13/19/23/31, columns4, raw-x cap4096,
   pair-visit cap1000000, query cap2048, four public targets, one round;
   exact and orbit tables, plus the existing counted orbit pipeline and strong
   single-target rho/A-A controls. Held-out degree37 uses columns4/8,
   query cap2048 and the same pair cap, four targets, one round.

   Mandatory wide diagnostic: exact registered E0 degree83, modulus
   z^83+z^45+z^2+z+1, subgroup order2417851639230796216685689, cofactor4.
   Orbit tables with columns4/16, raw-x cap4096, query cap256 and pair cap1000000,
   two public targets, one round. Also run the previous unsupported counted
   pipeline once, preserving its error as a control. A completed bounded search
   is not the N83 recovery gate. No N131 transfer is admitted. A full N83 rho
   run is not scheduled: at its existing minimum nonzero step factor it exceeds
   this bounded local campaign; no ratio is emitted for this diagnostic.

3. **Prime factor-base balance (Q1/Q10).** Reuse the original prime models
   (16/18/20/22-bit constructor seed59297), but fresh public targets. Compare
   negation-folded two-summand MITM at base sizes16/32/64/128 and the existing
   size64 three-summand direct-subtraction control, with matched rho.negation
   and A/A. Eight independent single-target workloads, one round, query
   cap100000. This tests the observed setup/yield tradeoff; these dimensions
   cannot establish an asymptotic exponent.

Every spec uses target_seed=measurement_seed=20261015, public targets,
warmup0, alternating order, required isolationL2 and watchdog180 seconds per
execution. Local runs use --cpus none --allow-busy (L0). Counts are admissible
as diagnostics; wall time is not a performance result. No runs share setup
across targets. The same-point reference is strong rho on Koblitz and negation
rho on prime curves. The harness birthday reference sqrt(pi/(2A)) is not
claimed as an exact universal lower bound.

## Evidence, decisions and limits

Freeze source/spec hashes in a commit before measurement. Preserve immutable
sessions, errors, unsuccessful attempts and raw commands. Audit every verified
run; replay failed deterministic outcomes separately. Compare per curve and
per target, with missing/censored costs left missing. S is counted GAE/sqrt(r).
Unpriced native costs prevent full-cost ratios; a quotient of counted lower
bounds is not a bound on complete-cost ratios. No exponent fitting or speedup
claim is allowed. New mathematical claims require faithful proof and Lean
formalization where feasible; this iteration claims executable constructions
and bounded observations only.

Report all ten question statuses, including untouched original obligations.
Q4/Q5 nonlinear/quasi-subfield scalable families and Q8/Q9 paid isogeny/norm
transport retain the preceding evidence and remain open. Missing independent
runner/GPU access remains an access gap. Each session and the results receive
source-linked Markdown/PDF reports and SVG/150dpi PNG graphs beside the sealed
records. Update affected canonical views, retain the old measurements and
publish commits on the existing ready PR1853. Conductor command was attempted
and is unavailable (exit127); no reservation is claimed.
